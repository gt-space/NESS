"""
structural.py - liner stress and low-cycle-fatigue life from a solved RegenCircuit.

Takes the outputs of RegenCircuit.solve_circuit (wall temperatures, coolant and gas pressures, channel
geometry) and computes, at every station along the circuit:

  1. Thermal stress  - restrained through-wall temperature gradient (biaxial)
  2. Mechanical stress - pressure bending of the liner between ribs (fixed-fixed strip),
                         liner hoop membrane stress, and (optionally) outer-jacket hoop stress
  3. Combined elastic stress state at four liner locations and its von Mises equivalent
  4. Neuber-corrected local elastic-plastic stress/strain for the 0 -> operating -> 0 duty cycle
  5. Cycles to failure (Coffin-Manson-Basquin with Morrow mean stress correction)

Sign convention: tension +, compression -. P_c > P_gas pushes the liner toward the gas side.

Liner locations (rib = fixed edge of the strip, mid = middle of the channel span):

    hot_mid   : gas face, channel mid-span     (pressure bending tensile, thermal compressive)
    hot_rib   : gas face, over the rib edge    (pressure bending compressive, thermal compressive)
    cold_rib  : coolant face, rib edge         (pressure bending tensile, thermal tensile)
    cold_mid  : coolant face, channel mid-span (pressure bending compressive, thermal tensile)

Simplifications to be aware of (all documented in the report printed by .summary()):
  * Thermal stress is the fully-restrained *gradient* stress E*alpha*dT / (2(1-v)); the mean-temperature
    hoop restraint ("dog-house" ratcheting) and 2-D rib/land temperature variation are not modeled.
  * The liner hoop membrane term treats the liner as an independent thin ring (conservative: the ribs
    and jacket share that load in reality).
  * One duty cycle = unstressed ambient -> full-power steady state -> unstressed. No creep, no hold-time
    damage, no thermal-transient overshoot.
  * Material strength/fatigue data come from materials.py; if `structural_data_verified` is False on the
    material, the numbers there are placeholders.
"""

import numpy as np
from scipy.optimize import brentq

PA_TO_KSI = 1.0 / 6894757.293168
PA_TO_MPA = 1.0e-6

LOCATIONS = ("hot_mid", "hot_rib", "cold_rib", "cold_mid")
LOCATION_LABELS = {
    "hot_mid": "gas face, mid-span",
    "hot_rib": "gas face, rib edge",
    "cold_rib": "coolant face, rib edge",
    "cold_mid": "coolant face, mid-span",
}

LIFE_CAP = 1.0e12   # cycles reported when the strain amplitude is below the CMB curve everywhere


def von_mises(s1, s2, s3):
    """Von Mises equivalent of three principal stresses."""
    return float(np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2)))


def neuber_first_load(sigma_e, E, K_cyc, n_cyc):
    """
    Local stress / strain after the first load excursion from zero.

    Neuber:  sigma * eps = sigma_e^2 / E        (sigma_e = elastic-equivalent peak stress, >= 0)
    Curve:   eps = sigma / E + (sigma / K')^(1/n')       (Ramberg-Osgood, cyclic)
    """
    if sigma_e <= 0.0:
        return 0.0, 0.0
    target = sigma_e ** 2 / E

    def g(sig):
        return sig * (sig / E + (sig / K_cyc) ** (1.0 / n_cyc)) - target

    sig = brentq(g, 0.0, sigma_e, xtol=1.0, rtol=1e-10)   # sigma_local <= sigma_e always
    return sig, sig / E + (sig / K_cyc) ** (1.0 / n_cyc)


def neuber_range(dsigma_e, E, K_cyc, n_cyc):
    """
    Local stress range / strain range for a full reversal (Masing doubling of the cyclic curve).

    Neuber:  dsigma * deps = dsigma_e^2 / E
    Curve:   deps = dsigma / E + 2 (dsigma / (2 K'))^(1/n')
    """
    if dsigma_e <= 0.0:
        return 0.0, 0.0
    target = dsigma_e ** 2 / E

    def g(ds):
        return ds * (ds / E + 2.0 * (ds / (2.0 * K_cyc)) ** (1.0 / n_cyc)) - target

    ds = brentq(g, 0.0, dsigma_e, xtol=1.0, rtol=1e-10)
    return ds, ds / E + 2.0 * (ds / (2.0 * K_cyc)) ** (1.0 / n_cyc)


def cycles_to_failure(eps_a, sigma_m, E, sigma_f, b, eps_f_prime, c):
    """
    Coffin-Manson-Basquin with Morrow mean stress correction, solved for cycles N_f:

        eps_a = (sigma_f' - sigma_m)/E * (2N)^b + eps_f' * (2N)^c

    Returns N_f (cycles). Values are bounded to [0.5, LIFE_CAP].
    """
    def eps_of_2N(log10_2N):
        two_n = 10.0 ** log10_2N
        return ((sigma_f - sigma_m) / E) * two_n ** b + eps_f_prime * two_n ** c

    lo, hi = 0.0, np.log10(2.0 * LIFE_CAP)       # 2N from 1 reversal to 2*LIFE_CAP
    if eps_a >= eps_of_2N(lo):
        return 0.5
    if eps_a <= eps_of_2N(hi) or (sigma_f - sigma_m) <= 0.0:
        return LIFE_CAP if (sigma_f - sigma_m) > 0.0 else 0.5
    log10_2N = brentq(lambda x: eps_of_2N(x) - eps_a, lo, hi, xtol=1e-10)
    return min(10.0 ** log10_2N / 2.0, LIFE_CAP)


class LinerStructuralAnalysis:
    """
    Stress and low-cycle-fatigue analysis of the regen liner (hot wall) for a solved RegenCircuit.

    Parameters
    ----------
    regen : RegenCircuit
        A circuit on which solve_circuit() has already been run.
    Kt : float
        Local stress concentration factor applied to the elastic stress at the liner/rib fillet
        (1.0 = none). With Kt = 1 the Neuber step is a pure plasticity correction.
    life_safety_factor : float
        Factor on cycles: allowable cycles = N_f / life_safety_factor.
    required_cycles : int or None
        Mission cycle requirement; the summary reports pass/fail against allowable cycles.
    include_liner_hoop : bool
        Include the liner hoop membrane stress -(P_c - P_gas) r_m / t_w (conservative).
    t_jacket : float or None
        Outer jacket thickness [m]. If given, the jacket hoop stress P_c r / t_jacket is checked against
        yield at the bulk coolant temperature (it does not enter the liner life).
    """

    def __init__(
        self,
        regen,
        Kt=1.0,
        life_safety_factor=4.0,
        required_cycles=None,
        include_liner_hoop=True,
        t_jacket=None,
    ):
        if not hasattr(regen, "T_hw_arr"):
            raise RuntimeError("Run regen.solve_circuit(...) before the structural analysis.")
        self.regen = regen
        self.material = regen.material
        self.Kt = float(Kt)
        self.life_safety_factor = float(life_safety_factor)
        self.required_cycles = required_cycles
        self.include_liner_hoop = bool(include_liner_hoop)
        self.t_jacket = None if t_jacket is None else float(t_jacket)
        self.warnings = []

        o, i = regen.circuit_outlet, regen.circuit_inlet
        eng = regen.engine
        self.z = np.asarray(regen.circuit_z, dtype=float)                       # [in]
        self.t_w = float(regen.t_w)                                             # [m]
        self.C_w = np.asarray(regen.C_w_arr, dtype=float)                       # [m]
        self.C_h = np.asarray(regen.C_h_arr, dtype=float)                       # [m]
        self.r_gas = np.asarray(eng.Contour_r, dtype=float)[o:i + 1] * 0.0254   # [m] gas-side liner radius
        self.T_hw = np.asarray(regen.T_hw_arr, dtype=float)                     # [K]
        self.T_cw = np.asarray(regen.T_cw_arr, dtype=float)                     # [K]
        self.T_c = np.asarray(regen.T_c_arr, dtype=float)                       # [K] bulk coolant
        self.P_c = np.asarray(regen.P_c_in_arr, dtype=float)                    # [Pa] coolant, station inlet
        self.P_g = np.asarray(eng.P, dtype=float)[o:i + 1] * 6894.76            # [Pa] gas static

        if not getattr(self.material, "structural_data_verified", True):
            self.warnings.append(
                f"'{self.material.material_name}' strength/fatigue data in materials.py are PLACEHOLDERS - "
                "replace with test or vendor data before trusting any stress margin or cycle count."
            )

    # ------------------------------------------------------------------ material
    def _props(self, T):
        props = self.material.structural_props(T)
        missing = [k for k in ("E", "alpha", "v", "sigma_y", "sigma_uts", "eps_f", "n_cyc") if props.get(k) is None]
        fat = props.get("fatigue")
        if fat is not None:
            missing = [k for k in missing if k not in ("sigma_uts", "eps_f")]
        if missing:
            raise ValueError(
                f"Material '{self.material.material_name}' is missing {missing} needed for the structural "
                "analysis. Add them in materials.py (see Material.set_strength_tables / structural_props)."
            )
        if not props["in_table"]:
            self._T_outside_table = max(getattr(self, "_T_outside_table", 0.0), float(T))
        return props

    @staticmethod
    def _cmb_coefficients(props):
        """CMB constants: user-supplied `fatigue` dict, else Manson's Universal Slopes from UTS/E/eps_f."""
        fat = props.get("fatigue")
        if fat is not None:
            return fat["sigma_f"], fat["b"], fat["eps_f_prime"], fat["c"]
        sigma_f = 1.9 * props["sigma_uts"]               # sigma_f' = 1.9 UTS
        eps_f_prime = 0.76 * props["eps_f"] ** 0.6       # eps_f'   = 0.76 eps_f^0.6
        return sigma_f, -0.12, eps_f_prime, -0.6

    # ------------------------------------------------------------------ solve
    def solve(self):
        n = len(self.z)
        loc_arr = lambda: {loc: np.zeros(n) for loc in LOCATIONS}

        # Stress components [Pa]
        self.dP = self.P_c - self.P_g
        self.r_mean = self.r_gas + 0.5 * self.t_w
        self.sigma_bend_rib = 0.5 * self.dP * (self.C_w / self.t_w) ** 2     # fixed-fixed strip, edge moment
        self.sigma_bend_mid = 0.5 * self.sigma_bend_rib                      # mid-span moment is half the edge
        self.sigma_hoop_membrane = (
            -self.dP * self.r_mean / self.t_w if self.include_liner_hoop else np.zeros(n)
        )
        self.sigma_thermal_hot = np.zeros(n)     # gas-face thermal stress (compressive, negative)
        self.sigma_jacket_hoop = None

        self.sigma_theta = loc_arr()
        self.sigma_z = loc_arr()
        self.sigma_r = loc_arr()
        self.sigma_vm_elastic = loc_arr()       # signed elastic equivalent, before Kt
        self.sigma_local = loc_arr()            # Neuber local peak stress (signed)
        self.eps_amp = loc_arr()                # strain amplitude of the stabilized loop
        self.sigma_mean = loc_arr()
        self.sigma_yield = loc_arr()
        self.cycles = loc_arr()
        self.yield_ratio = loc_arr()            # elastic equivalent * Kt / yield

        for k in range(n):
            T_mean = 0.5 * (self.T_hw[k] + self.T_cw[k])
            pm = self._props(T_mean)
            dT = self.T_hw[k] - self.T_cw[k]
            sig_th = pm["E"] * pm["alpha"] * dT / (2.0 * (1.0 - pm["v"]))     # magnitude, biaxial
            self.sigma_thermal_hot[k] = -sig_th

            # location -> (face temperature, bending sign/factor in theta, thermal sign, radial stress)
            specs = {
                "hot_mid": (self.T_hw[k], +self.sigma_bend_mid[k], -sig_th, -self.P_g[k]),
                "hot_rib": (self.T_hw[k], -self.sigma_bend_rib[k], -sig_th, -self.P_g[k]),
                "cold_rib": (self.T_cw[k], +self.sigma_bend_rib[k], +sig_th, -self.P_c[k]),
                "cold_mid": (self.T_cw[k], -self.sigma_bend_mid[k], +sig_th, -self.P_c[k]),
            }
            for loc, (T_face, bend, th, s_r) in specs.items():
                p = self._props(T_face)
                s_theta = bend + self.sigma_hoop_membrane[k] + th
                s_z = p["v"] * bend + th                       # cylindrical bending: axial = v * tangential
                vm = von_mises(s_theta, s_z, s_r)
                sign = np.sign(s_theta + s_z + s_r) or np.sign(s_theta) or 1.0
                self.sigma_theta[loc][k] = s_theta
                self.sigma_z[loc][k] = s_z
                self.sigma_r[loc][k] = s_r
                self.sigma_vm_elastic[loc][k] = sign * vm

                # --- Neuber + Ramberg-Osgood (cyclic) ---
                sig_y_cyc = p["sigma_y"] * p["cyc_yield_ratio"]
                K_cyc = sig_y_cyc / 0.002 ** p["n_cyc"]
                sigma_e = self.Kt * vm
                sig_peak, _ = neuber_first_load(sigma_e, p["E"], K_cyc, p["n_cyc"])
                d_sig, d_eps = neuber_range(sigma_e, p["E"], K_cyc, p["n_cyc"])
                sig_peak *= sign
                sig_mean = sig_peak - sign * 0.5 * d_sig          # loop spans [sig_peak, sig_peak - sign*d_sig]
                eps_a = 0.5 * d_eps

                sigma_f, b, eps_f_p, c = self._cmb_coefficients(p)
                Nf = cycles_to_failure(eps_a, sig_mean, p["E"], sigma_f, b, eps_f_p, c)

                self.sigma_local[loc][k] = sig_peak
                self.eps_amp[loc][k] = eps_a
                self.sigma_mean[loc][k] = sig_mean
                self.sigma_yield[loc][k] = p["sigma_y"]
                self.cycles[loc][k] = Nf
                self.yield_ratio[loc][k] = sigma_e / p["sigma_y"]

        if getattr(self, "_T_outside_table", 0.0) > 0.0:
            self.warnings.append(
                f"Wall temperature up to {self._T_outside_table:.0f} K is outside the strength table of "
                f"'{self.material.material_name}'; the nearest table-end strengths were used (unreliable)."
            )

        # Outer jacket hoop stress (informational)
        if self.t_jacket is not None:
            r_jacket = self.r_gas + self.t_w + self.C_h + 0.5 * self.t_jacket
            self.sigma_jacket_hoop = self.P_c * r_jacket / self.t_jacket
            self.jacket_yield = np.array([self._props(T)["sigma_y"] for T in self.T_c])

        # Governing results per station and overall
        stack = np.vstack([self.cycles[loc] for loc in LOCATIONS])           # (4, n)
        self.cycles_min = stack.min(axis=0)
        self.governing_loc = np.array([LOCATIONS[j] for j in stack.argmin(axis=0)])
        self.allowable_cycles = self.cycles_min / self.life_safety_factor
        self.k_crit = int(np.argmin(self.cycles_min))
        return self

    # ------------------------------------------------------------------ reporting
    def summary(self):
        if not hasattr(self, "cycles_min"):
            self.solve()
        k = self.k_crit
        loc = self.governing_loc[k]
        print("\n#### ---- LINER STRESS & LIFE ---- ####")
        print(f"Material: {self.material.material_name}   Kt: {self.Kt}   Life safety factor: {self.life_safety_factor}")
        print(f"Max |pressure bending| at rib [ksi]: {np.max(np.abs(self.sigma_bend_rib)) * PA_TO_KSI:.2f}")
        print(f"Max liner hoop membrane  [ksi]: {np.max(np.abs(self.sigma_hoop_membrane)) * PA_TO_KSI:.2f}")
        print(f"Max thermal (gas face)   [ksi]: {np.max(np.abs(self.sigma_thermal_hot)) * PA_TO_KSI:.2f}")
        worst_vm = max(np.max(np.abs(self.sigma_vm_elastic[l])) for l in LOCATIONS)
        print(f"Max elastic von Mises    [ksi]: {worst_vm * PA_TO_KSI:.2f}")
        print(f"Max (Kt*elastic stress) / yield: {max(np.max(self.yield_ratio[l]) for l in LOCATIONS):.2f}")
        print(f"Critical station: z = {self.z[k]:.3f} in, {LOCATION_LABELS[loc]}")
        print(f"  T_hw = {self.T_hw[k]:.0f} K, T_cw = {self.T_cw[k]:.0f} K, dP = {self.dP[k] * PA_TO_KSI * 1000:.0f} psid")
        print(f"  Neuber local stress [ksi]: {self.sigma_local[loc][k] * PA_TO_KSI:.2f}, "
              f"mean [ksi]: {self.sigma_mean[loc][k] * PA_TO_KSI:.2f}, strain amplitude: {self.eps_amp[loc][k]:.5f}")
        print(f"  Cycles to failure: {self._fmt_cycles(self.cycles_min[k])}   "
              f"allowable (/{self.life_safety_factor:g}): {self._fmt_cycles(self.allowable_cycles[k])}")
        if self.required_cycles is not None:
            ok = self.allowable_cycles[k] >= self.required_cycles
            print(f"  Required cycles: {self.required_cycles}  ->  {'PASS' if ok else 'FAIL'}")
        if self.sigma_jacket_hoop is not None:
            ratio = self.sigma_jacket_hoop / self.jacket_yield
            print(f"Jacket hoop stress max [ksi]: {np.max(self.sigma_jacket_hoop) * PA_TO_KSI:.2f}  "
                  f"(max stress / yield = {np.max(ratio):.2f})")
        for w in self.warnings:
            print(f"WARNING: {w}")

    @staticmethod
    def _fmt_cycles(n):
        if n >= LIFE_CAP:
            return f">{LIFE_CAP:.0e}"
        if n <= 0.5:
            return "<1"
        return f"{n:,.0f}" if n < 1e6 else f"{n:.3e}"

    def plot(self, show=True):
        import matplotlib.pyplot as plt
        if not hasattr(self, "cycles_min"):
            self.solve()

        fig, axs = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
        z = self.z

        axs[0].plot(z, self.sigma_thermal_hot * PA_TO_KSI, label="Thermal, gas face")
        axs[0].plot(z, self.sigma_bend_rib * PA_TO_KSI, label="Pressure bending at rib")
        axs[0].plot(z, self.sigma_hoop_membrane * PA_TO_KSI, label="Liner hoop membrane")
        if self.sigma_jacket_hoop is not None:
            axs[0].plot(z, self.sigma_jacket_hoop * PA_TO_KSI, "--", label="Jacket hoop")
        axs[0].set_title("Stress components")
        axs[0].set_ylabel("Stress [ksi]")
        axs[0].legend()

        for loc in LOCATIONS:
            axs[1].plot(z, self.sigma_vm_elastic[loc] * PA_TO_KSI, label=LOCATION_LABELS[loc])
        axs[1].plot(z, self.sigma_yield["hot_mid"] * PA_TO_KSI, "k--", label="Yield at gas-face T")
        axs[1].set_title("Elastic von Mises equivalent (signed)")
        axs[1].set_ylabel("Stress [ksi]")
        axs[1].legend()

        for loc in LOCATIONS:
            axs[2].semilogy(z, np.maximum(self.cycles[loc], 0.5), label=LOCATION_LABELS[loc])
        axs[2].semilogy(z, np.maximum(self.allowable_cycles, 0.5), "k--", label="Allowable (with safety factor)")
        if self.required_cycles is not None:
            axs[2].axhline(self.required_cycles, color="r", linestyle=":", label="Required cycles")
        axs[2].set_title("Cycles to failure (Neuber + Coffin-Manson-Basquin)")
        axs[2].set_xlabel("Axial Position [in]")
        axs[2].set_ylabel("Cycles")
        axs[2].legend()

        fig.tight_layout()
        if show:
            plt.show()
        return fig

    def export_to_excel(self, filename="liner_stress_life.xlsx"):
        try:
            from openpyxl import Workbook
        except ImportError as exc:
            raise ImportError("openpyxl is required for Excel export. Install it with: pip install openpyxl") from exc
        if not hasattr(self, "cycles_min"):
            self.solve()

        wb = Workbook()
        ws = wb.active
        ws.title = "Liner Stress & Life"
        headers = [
            "Axial Position (in)", "Hot Wall T (K)", "Cold Wall T (K)", "dP coolant-gas (Pa)",
            "Channel Width (m)", "Thermal Stress Gas Face (Pa)", "Pressure Bending @Rib (Pa)",
            "Pressure Bending @Mid (Pa)", "Liner Hoop Membrane (Pa)",
        ]
        for loc in LOCATIONS:
            headers += [
                f"{loc} sigma_theta (Pa)", f"{loc} sigma_z (Pa)", f"{loc} elastic vM signed (Pa)",
                f"{loc} Neuber peak (Pa)", f"{loc} mean (Pa)", f"{loc} strain amp (-)", f"{loc} cycles (-)",
            ]
        headers += ["Min Cycles (-)", "Governing Location", "Allowable Cycles (-)"]
        ws.append(headers)
        for k in range(len(self.z)):
            row = [
                self.z[k], self.T_hw[k], self.T_cw[k], self.dP[k], self.C_w[k], self.sigma_thermal_hot[k],
                self.sigma_bend_rib[k], self.sigma_bend_mid[k], self.sigma_hoop_membrane[k],
            ]
            for loc in LOCATIONS:
                row += [
                    self.sigma_theta[loc][k], self.sigma_z[loc][k], self.sigma_vm_elastic[loc][k],
                    self.sigma_local[loc][k], self.sigma_mean[loc][k], self.eps_amp[loc][k], self.cycles[loc][k],
                ]
            row += [self.cycles_min[k], str(self.governing_loc[k]), self.allowable_cycles[k]]
            ws.append([float(v) if not isinstance(v, str) else v for v in row])
        wb.save(filename)
        print(f"Saved liner stress/life export: {filename}")
