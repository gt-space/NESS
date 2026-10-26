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

from redlines import draw_redline, span_text
from redline_explain import explain_structural, print_explanations

PA_TO_KSI = 1.0 / 6894757.293168
PA_TO_MPA = 1.0e-6

LOCATIONS = ("hot_mid", "hot_rib", "cold_rib", "cold_mid")
LOCATION_LABELS = {
    "hot_mid": "gas face, mid-span",
    "hot_rib": "gas face, rib edge",
    "cold_rib": "coolant face, rib edge",
    "cold_mid": "coolant face, mid-span",
}

FIRE_COUNTS = (1, 2, 3, 4, 5, 10, 20)   # firings reported by yielding_history()
MAX_FIRES = FIRE_COUNTS[-1]               # nothing is reported or plotted beyond this many fires
YIELD_STRAIN_THRESHOLD = 2.0e-4   # plastic strain counted as "measurable yielding" (10% of the 0.2% offset strain)

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
    mission_fires : int or None
        Hotfires the hardware must survive (verdict at the top of the summary). Defaults to required_cycles, else 10.
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
        mission_fires=None,
    ):
        if not hasattr(regen, "T_hw_arr"):
            raise RuntimeError("Run regen.solve_circuit(...) before the structural analysis.")
        self.regen = regen
        self.material = regen.material
        self.Kt = float(Kt)
        self.life_safety_factor = float(life_safety_factor)
        self.required_cycles = required_cycles
        # Number of hotfires the hardware must survive; defaults to required_cycles, else 10.
        self.mission_fires = int(round(mission_fires if mission_fires else (required_cycles if required_cycles else 10)))
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
        self.eps_p_first = loc_arr()            # plastic strain left by the first firing's load
        self.eps_p_cycle = loc_arr()            # plastic strain range of each later (stabilized) firing
        self.eps_f_loc = loc_arr()              # fracture ductility at the face temperature

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
                sig_peak, eps_first = neuber_first_load(sigma_e, p["E"], K_cyc, p["n_cyc"])
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
                self.eps_p_first[loc][k] = max(eps_first - abs(sig_peak) / p["E"], 0.0)
                self.eps_p_cycle[loc][k] = max(d_eps - d_sig / p["E"], 0.0)
                self.eps_f_loc[loc][k] = p["eps_f"] if p["eps_f"] is not None else np.nan

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

    # ------------------------------------------------------------------ redlines
    def _redline_masks(self):
        """Failure-limit crossings per station: {label: bool array}."""
        ratio = np.max(np.vstack([self.yield_ratio[loc] for loc in LOCATIONS]), axis=0)
        out = {
            "Elastic stress above yield (plastic cycling)": ratio > 1.0,
            "Cycles to failure < 1 (fails on first cycle)": self.cycles_min < 1.0,
        }
        if self.required_cycles is not None:
            out[f"Allowable cycles below required ({self.required_cycles:g})"] = self.allowable_cycles < self.required_cycles
        if self.sigma_jacket_hoop is not None:
            out["Jacket hoop stress above yield"] = self.sigma_jacket_hoop > self.jacket_yield
        return out

    def check_redlines(self):
        """List of violation messages (empty if every limit is respected)."""
        if not hasattr(self, "cycles_min"):
            self.solve()
        return [f"{label} at {span_text(self.z, m)}" for label, m in self._redline_masks().items() if m.any()]


    # ------------------------------------------------------------------ yielding vs number of fires
    def yielding_history(self, fires=FIRE_COUNTS):
        """
        How much permanent (plastic) deformation has built up after each number of fires.

        Per location and station: the first fire leaves plastic strain eps_p1 (Neuber local strain minus its
        elastic part). Every later fire adds a plastic strain range d_eps_p to the stabilized loop. The
        accumulated plastic strain (sum of |d eps_p|, forward + back) after n fires is

            P(n) = eps_p1 + (2 n - 1) * d_eps_p

        Returns a list of dicts, one per fire count, with the worst station/location, how many stations
        have yielded measurably (P >= YIELD_STRAIN_THRESHOLD), fraction of fracture ductility used and
        Miner fatigue damage n / N_f at the critical station.
        """
        if not hasattr(self, "cycles_min"):
            self.solve()
        n_st = len(self.z)
        rows = []
        for n in fires:
            P = {loc: self.eps_p_first[loc] + (2 * n - 1) * self.eps_p_cycle[loc] for loc in LOCATIONS}
            stack = np.vstack([P[loc] for loc in LOCATIONS])                 # (4, n_st)
            k_flat = int(np.argmax(stack))
            j, k = divmod(k_flat, n_st)
            loc = LOCATIONS[j]
            yielded = int(np.sum(np.any(stack >= YIELD_STRAIN_THRESHOLD, axis=0)))
            ductility = stack / np.vstack([self.eps_f_loc[l] for l in LOCATIONS])
            rows.append(dict(
                fires=n, peak_strain=float(stack[j, k]), z=float(self.z[k]), loc=loc, k=k,
                stations_yielded=yielded, stations=n_st,
                ductility_used=float(np.nanmax(ductility)),
                fatigue_used=float(n / self.cycles_min[self.k_crit]),
            ))
        return rows

    def print_yielding_history(self, fires=FIRE_COUNTS):
        rows = self.yielding_history(fires)
        k = self.k_crit
        loc = self.governing_loc[k]
        first = max(np.max(self.eps_p_first[l]) for l in LOCATIONS)
        per_fire = max(np.max(self.eps_p_cycle[l]) for l in LOCATIONS)
        print("\n---- YIELDING vs NUMBER OF FIRES ----")
        if rows[-1]["peak_strain"] < YIELD_STRAIN_THRESHOLD:
            print("No station reaches measurable yielding (plastic strain < "
                  f"{YIELD_STRAIN_THRESHOLD * 100:.3f}%) within {rows[-1]['fires']} fires; the liner stays essentially elastic.")
        else:
            if per_fire < 1e-6:
                story = (f"The first fire yields the liner (up to {first * 100:.3f}% plastic strain) and the residual stress then keeps "
                         "later fires elastic ('shakedown'): damage stops growing after fire 1.")
            else:
                story = (f"The first fire leaves up to {first * 100:.3f}% plastic strain, and every fire after that pushes the metal "
                         f"through a further +/-{per_fire * 100 / 2:.3f}% plastic swing, so deformation and fatigue damage keep "
                         "growing with each firing.")
            print(story)
        print(f"{'Fires':>6} | {'Peak plastic strain':>20} | {'Where':<30} | {'Stations yielded':>16} | {'Ductility used':>14} | {'Fatigue life used':>17}")
        for r in rows:
            where = f"z={r['z']:.2f} in, {LOCATION_LABELS[r['loc']]}"
            print(f"{r['fires']:>6} | {r['peak_strain'] * 100:>18.4f} % | {where:<30} | "
                  f"{r['stations_yielded']:>8} of {r['stations']:<5}| {r['ductility_used'] * 100:>12.1f} % | {r['fatigue_used'] * 100:>15.4g} %")
        print("Plastic strain = accumulated |d eps_p| from the Neuber/Ramberg-Osgood model; 'yielded' means >= "
              f"{YIELD_STRAIN_THRESHOLD * 100:.3f}%. Fatigue life used = fires / cycles-to-failure at the critical station "
              f"(z = {self.z[k]:.2f} in, {LOCATION_LABELS[loc]}). Ratcheting (net growth of the wall each fire) and creep are not modeled.")



    def ratchet_stations(self):
        """
        Two-bar (Bree-type) ratchet quantities at every station (see throat_ratcheting for the model):
        X = hoop pressure stress / hot-face yield, Y = thermal stress / yield, d_eps_r = net hoop strain
        gained per fire, regime in {"elastic", "shakedown", "ratcheting", "collapse"}, direction (-1 inward).
        """
        if hasattr(self, "_ratchet_cache"):
            return self._ratchet_cache
        n = len(self.z)
        keys = ("E", "alpha", "v", "sigma_y", "dT", "eps_T", "dP", "sigma_p", "X", "Y", "d_eps_r", "direction")
        out = {key: np.zeros(n) for key in keys}
        out["regime"] = np.empty(n, dtype=object)
        for k in range(n):
            p_hot = self._props(self.T_hw[k])
            p_mid = self._props(0.5 * (self.T_hw[k] + self.T_cw[k]))
            E, alpha, v = p_mid["E"], p_mid["alpha"], p_mid["v"]
            sigma_y = p_hot["sigma_y"]                       # hot-face yield: the controlling (weakest) layer
            dT = self.T_hw[k] - self.T_cw[k]
            eps_T = alpha * dT / (1.0 - v)
            dP = self.dP[k]
            sigma_p = abs(dP) * self.r_mean[k] / self.t_w     # primary hoop membrane stress, always evaluated here
            X = sigma_p / sigma_y
            Y = E * eps_T / (2.0 * sigma_y)                   # surface thermal stress / yield
            d_eps_r = max(0.0, eps_T - 4.0 * (sigma_y - sigma_p) / E)
            if X >= 1.0:
                regime = "collapse"
            elif X + Y <= 1.0:
                regime = "elastic"
            elif d_eps_r > 0.0:
                regime = "ratcheting"
            else:
                regime = "shakedown"
            for key, val in (("E", E), ("alpha", alpha), ("v", v), ("sigma_y", sigma_y), ("dT", dT), ("eps_T", eps_T),
                             ("dP", dP), ("sigma_p", sigma_p), ("X", X), ("Y", Y), ("d_eps_r", d_eps_r),
                             ("direction", -1.0 if dP > 0 else 1.0)):
                out[key][k] = val
            out["regime"][k] = regime
        self._ratchet_cache = out
        return out

    # ------------------------------------------------------------------ throat size vs fires (thermal ratcheting)
    def throat_ratcheting(self, fires=FIRE_COUNTS):
        """
        Approximate throat size after each number of fires due to thermal ratcheting.

        Thermal ratcheting: a steady primary stress (here the hoop pressure stress on the liner, dP*r/t,
        compressive when coolant pressure exceeds gas pressure) plus a stress that cycles every fire
        (the through-wall thermal gradient) can make the wall creep a little further in the same direction
        each fire, instead of just cycling back and forth.

        Estimate (classic two-bar / Bree model, elastic-perfectly-plastic):

            eps_T   = alpha * dT_wall / (1 - v)           differential thermal strain, hot face vs cold face
            X       = sigma_p / sigma_y                   primary (hoop) stress ratio
            d_eps_r = max(0, eps_T - 4 (sigma_y - sigma_p) / E)     net hoop strain gained EVERY fire

        The throat radius then changes by r(n) = r0 * exp(-s * n * d_eps_r), with s = +1 (shrinks) when the
        coolant pressure is above the gas pressure and s = -1 (grows) otherwise. Chamber pressure at a fixed
        mass flow scales as 1/At. This is a first-order bound for a smooth perfectly plastic wall; it ignores
        strain hardening, creep, erosion, and 3-D rib/crown bowing, so use it for trends, not for tolerances.

        Returns dict with regime, per-fire strain and one row per fire count (None if the throat is outside
        the regen circuit).
        """
        if not hasattr(self, "cycles_min"):
            self.solve()
        eng = self.regen.engine
        k = int(eng.throat_ind - self.regen.circuit_outlet)
        if k < 0 or k >= len(self.z):
            return None

        st = self.ratchet_stations()
        E, alpha, v = st["E"][k], st["alpha"][k], st["v"][k]
        sigma_y, dT, eps_T, dP = st["sigma_y"][k], st["dT"][k], st["eps_T"][k], st["dP"][k]
        sigma_p, X, Y, d_eps_r = st["sigma_p"][k], st["X"][k], st["Y"][k], st["d_eps_r"][k]
        regime = st["regime"][k]
        direction = st["direction"][k]

        r0 = float(self.r_gas[k])
        rows = []
        for n in fires:
            r_n = r0 * float(np.exp(direction * n * d_eps_r))
            d_area = (r_n / r0) ** 2 - 1.0
            rows.append(dict(
                fires=n, radius=r_n, d_radius=r_n - r0, d_area_pct=100.0 * d_area,
                d_Pc_pct=100.0 * ((r0 / r_n) ** 2 - 1.0),
            ))
        return dict(k=k, z=float(self.z[k]), regime=regime, X=float(X), Y=float(Y), eps_T=float(eps_T),
                    d_eps_r=float(d_eps_r), direction=direction, r0=r0, dT=float(dT),
                    sigma_p=float(sigma_p), sigma_y=float(sigma_y), rows=rows)

    def print_throat_ratcheting(self, fires=FIRE_COUNTS):
        res = self.throat_ratcheting(fires)
        print("\n---- THROAT SIZE vs NUMBER OF FIRES (thermal ratcheting) ----")
        if res is None:
            print("The throat is outside the regen circuit, so no throat ratcheting estimate is available.")
            return
        mm, um = 1000.0, 1.0e6
        print(f"At the throat (z = {res['z']:.2f} in): wall gradient {res['dT']:.0f} K, hoop pressure stress "
              f"{res['sigma_p'] * PA_TO_KSI:.1f} ksi vs yield {res['sigma_y'] * PA_TO_KSI:.1f} ksi (X = {res['X']:.2f}), "
              f"thermal stress / yield Y = {res['Y']:.2f}.")
        if res["regime"] == "elastic":
            print("Regime: ELASTIC. The throat wall never yields, so the throat does not change size from ratcheting.")
        elif res["regime"] == "collapse":
            print("Regime: PRIMARY STRESS ABOVE YIELD. The steady pressure load alone exceeds yield on the liner as an independent "
                  "thin ring, so the numbers below are a worst-case plastic-collapse trend, not a true ratcheting rate. In a real "
                  "engine the ribs and outer jacket carry part of that load, so this is very conservative. Check the pressure "
                  "load path (thicker liner / narrower channels / lower coolant pressure) before reading the table literally.")
        elif res["regime"] == "shakedown":
            print("Regime: SHAKEDOWN. The wall yields in the first firings, then settles into purely elastic cycling. "
                  "Ratcheting stops, so the throat size stays fixed after that.")
        else:
            way = "SHRINKS (closes)" if res["direction"] < 0 else "GROWS"
            print(f"Regime: RATCHETING. Each fire the wall yields a little further in the same direction instead of "
                  f"returning, gaining {res['d_eps_r'] * 100:.4f}% hoop strain per fire, so the throat {way} fire after fire. "
                  "The steady pressure squeezing the liner is what picks the direction; the hot/cold thermal cycle supplies the "
                  "energy to keep it moving.")
        print(f"{'Fires':>6} | {'Throat radius [mm]':>18} | {'Change [um]':>11} | {'Throat area change':>18} | {'Chamber pressure change':>23}")
        print(f"{0:>6} | {res['r0'] * mm:>18.4f} | {0.0:>11.2f} | {0.0:>16.3f} % | {0.0:>21.3f} %")
        for r in res["rows"]:
            print(f"{r['fires']:>6} | {r['radius'] * mm:>18.4f} | {r['d_radius'] * um:>11.2f} | "
                  f"{r['d_area_pct']:>16.3f} % | {r['d_Pc_pct']:>21.3f} %")
        print("Approximate: two-bar (Bree-type) ratchet strain for an elastic-perfectly-plastic wall. Ignores strain hardening, "
              "creep, erosion and 3-D bowing between ribs; chamber pressure change assumes constant propellant flow and c*.")



    # ------------------------------------------------------------------ mission verdict
    def mission_verdict(self):
        """
        Will the liner survive `mission_fires` hotfires? Four checks, each with a plain reason:
        fatigue life (with safety factor), fracture ductility, hot wall below melting, no throat ratcheting.
        Yielding every fire is reported as a note, not a failure, because low-cycle fatigue already accounts for it.
        """
        if not hasattr(self, "cycles_min"):
            self.solve()
        n = self.mission_fires
        row = self.yielding_history((n,))[0]
        rat = self.throat_ratcheting((n,))
        k = self.k_crit
        melt = getattr(self.material, "MeltingPoint", None)
        checks = []
        allow = float(self.allowable_cycles[k])
        checks.append(("Fatigue life", allow >= n,
                       f"allowable {self._fmt_cycles(allow)} cycles vs {n} needed (cycles to failure {self._fmt_cycles(self.cycles_min[k])}, "
                       f"{self.life_safety_factor:g}x safety factor)"))
        checks.append(("Fracture ductility", row["ductility_used"] < 1.0,
                       f"{row['ductility_used'] * 100:.1f}% of ductility used after {n} fires"))
        if melt is not None:
            margin = melt - float(np.max(self.T_hw))
            checks.append(("Hot wall vs melting", margin > 0,
                           f"peak {np.max(self.T_hw):.0f} K, {abs(margin):.0f} K {'below' if margin > 0 else 'ABOVE'} melting"))
        if rat is not None:
            last = rat["rows"][-1]
            checks.append(("Throat size", rat["regime"] != "ratcheting",
                           f"{rat['regime']}; radius {last['d_radius'] * 1e6:+.0f} um, area {last['d_area_pct']:+.2f}%, chamber pressure "
                           f"{last['d_Pc_pct']:+.2f}% after {n} fires"))
        notes = []
        ratio = max(float(np.max(self.yield_ratio[l])) for l in LOCATIONS)
        if ratio > 1.0:
            notes.append(f"liner yields every fire (worst stress {ratio:.2f}x yield, {row['peak_strain'] * 100:.2f}% plastic strain after {n} fires)")
        if rat is not None and rat["regime"] == "collapse":
            notes.append("throat hoop stress is above yield (plastic-collapse trend; very conservative model)")
        return dict(fires=n, survives=all(c[1] for c in checks), checks=checks, notes=notes)

    def print_mission_verdict(self):
        v = self.mission_verdict()
        print(f"\n--------- WILL IT SURVIVE {v['fires']} HOTFIRES?  {'YES' if v['survives'] else '** NO **'} ---------")
        for name, ok, why in v["checks"]:
            print(f"  [{'pass' if ok else 'FAIL'}] {name:<20}: {why}")
        for note in v["notes"]:
            print(f"  note: {note}")

    # ------------------------------------------------------------------ digestible datapoints
    def key_metrics(self):
        """The handful of numbers that summarize the liner's structural health, as an ordered dict."""
        if not hasattr(self, "cycles_min"):
            self.solve()
        k = self.k_crit
        loc = self.governing_loc[k]
        hist = {r["fires"]: r for r in self.yielding_history((1, 10, 20))}
        rat = self.throat_ratcheting((1, 10, 20))
        melt = getattr(self.material, "MeltingPoint", None)
        m = {
            "material": self.material.material_name,
            "mission_fires": self.mission_fires,
            "survives_mission": self.mission_verdict()["survives"],
            "critical_z_in": float(self.z[k]),
            "critical_location": LOCATION_LABELS[loc],
            "cycles_to_failure": float(self.cycles_min[k]),
            "allowable_cycles": float(self.allowable_cycles[k]),
            "required_cycles": self.required_cycles,
            "life_verdict": ("n/a" if self.required_cycles is None
                             else ("PASS" if self.allowable_cycles[k] >= self.required_cycles else "FAIL")),
            "max_stress_over_yield": float(max(np.max(self.yield_ratio[l]) for l in LOCATIONS)),
            "max_hot_wall_K": float(np.max(self.T_hw)),
            "melting_margin_K": None if melt is None else float(melt - np.max(self.T_hw)),
            "plastic_strain_pct_1_fire": hist[1]["peak_strain"] * 100,
            "plastic_strain_pct_10_fires": hist[10]["peak_strain"] * 100,
            "plastic_strain_pct_20_fires": hist[20]["peak_strain"] * 100,
            "stations_yielded_20_fires": f"{hist[20]['stations_yielded']} of {hist[20]['stations']}",
            "fatigue_life_used_pct_20_fires": hist[20]["fatigue_used"] * 100,
        }
        if rat is None:
            m["throat_regime"] = "outside circuit"
        else:
            last = rat["rows"][-1]
            m.update({
                "throat_regime": rat["regime"],
                "throat_strain_per_fire_pct": rat["d_eps_r"] * 100,
                "throat_radius_change_um_20_fires": last["d_radius"] * 1e6,
                "throat_area_change_pct_20_fires": last["d_area_pct"],
                "chamber_pressure_change_pct_20_fires": last["d_Pc_pct"],
            })
        return m

    def print_scorecard(self):
        """Compact at-a-glance block; the detailed tables follow it in summary()."""
        m = self.key_metrics()
        verdict = {"PASS": "PASS", "FAIL": "** FAIL **", "n/a": "(no requirement set)"}[m["life_verdict"]]
        self.print_mission_verdict()
        print("\n==================== STRUCTURAL SCORECARD ====================")
        print(f"  Material                 : {m['material']}")
        print(f"  Life vs requirement      : {verdict}   (allowable {self._fmt_cycles(m['allowable_cycles'])}"
              + (f" vs required {m['required_cycles']:g})" if m["required_cycles"] is not None else ")"))
        print(f"  Critical spot            : z = {m['critical_z_in']:.2f} in, {m['critical_location']}")
        print(f"  Cycles to failure        : {self._fmt_cycles(m['cycles_to_failure'])}")
        print(f"  Worst stress / yield     : {m['max_stress_over_yield']:.2f}x" + ("  (yielding every fire)" if m["max_stress_over_yield"] > 1 else "  (elastic)"))
        if m["melting_margin_K"] is not None:
            print(f"  Hot wall                 : peak {m['max_hot_wall_K']:.0f} K, {m['melting_margin_K']:.0f} K below melting")
        print(f"  Plastic strain (peak)    : {m['plastic_strain_pct_1_fire']:.3f}% after 1 fire, "
              f"{m['plastic_strain_pct_10_fires']:.3f}% after 10, {m['plastic_strain_pct_20_fires']:.3f}% after 20")
        print(f"  Stations yielded         : {m['stations_yielded_20_fires']} after 20 fires")
        print(f"  Fatigue life used        : {m['fatigue_life_used_pct_20_fires']:.4g}% after 20 fires")
        if m["throat_regime"] != "outside circuit":
            print(f"  Throat ratcheting        : {m['throat_regime'].upper()}, "
                  f"{m['throat_radius_change_um_20_fires']:+.1f} um radius / {m['throat_area_change_pct_20_fires']:+.2f}% area / "
                  f"{m['chamber_pressure_change_pct_20_fires']:+.2f}% chamber pressure after 20 fires")
        print("===============================================================")

    def plot_lifecycle(self, show=True):
        """Yielding and fatigue build-up vs number of fires, and where along the engine the yielding is."""
        import matplotlib.pyplot as plt
        if not hasattr(self, "cycles_min"):
            self.solve()
        grid = [int(v) for v in np.unique(np.round(np.logspace(0, np.log10(MAX_FIRES), 40)))]
        rows = self.yielding_history(grid)
        fires = np.array([r["fires"] for r in rows])
        marks = np.array(FIRE_COUNTS)
        mrows = self.yielding_history(FIRE_COUNTS)

        fig, axs = plt.subplots(2, 2, figsize=(13, 9))

        ax = axs[0, 0]
        ax.plot(fires, [r["peak_strain"] * 100 for r in rows], color="C3")
        ax.scatter(marks, [r["peak_strain"] * 100 for r in mrows], color="C3", zorder=3)
        ax.set_xscale("log")
        ax.set_title("Worst accumulated plastic strain vs fires")
        ax.set_xlabel("Number of fires")
        ax.set_ylabel("Plastic strain [%]")
        ax.grid(True, alpha=0.3)

        ax = axs[0, 1]
        ductility = np.array([r["ductility_used"] for r in rows]) * 100
        ax.plot(fires, ductility, color="C1", label="Ductility used")
        draw_redline(ax, fires, ductility, 100.0, "Ductility exhausted (100%)", draw_limit=True)
        ax.set_xscale("log")
        ax.set_title("Fracture ductility used vs fires")
        ax.set_xlabel("Number of fires")
        ax.set_ylabel("Percent of ductility [%]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="center left")

        ax = axs[1, 0]
        ax.plot(fires, [r["stations_yielded"] for r in rows], color="C0", label="Stations yielded")
        ax.set_xscale("log")
        ax.set_title("How much of the liner has yielded")
        ax.set_xlabel("Number of fires")
        ax.set_ylabel(f"Stations with plastic strain >= {YIELD_STRAIN_THRESHOLD * 100:.2f}%")
        ax.set_ylim(0, len(self.z) * 1.05)
        ax.grid(True, alpha=0.3)

        ax = axs[1, 1]
        for n, col in ((1, "C2"), (10, "C1"), (MAX_FIRES, "C3")):
            P = np.max(np.vstack([self.eps_p_first[l] + (2 * n - 1) * self.eps_p_cycle[l] for l in LOCATIONS]), axis=0)
            ax.plot(self.z, P * 100, color=col, label=f"After {n} fire{'s' if n > 1 else ''}")
        ax.axhline(YIELD_STRAIN_THRESHOLD * 100, color="k", linestyle=":", label="Measurable yielding")
        ax.set_yscale("log")
        ax.set_title("Where the yielding is (by axial position)")
        ax.set_xlabel("Axial Position [in]")
        ax.set_ylabel("Plastic strain [%]")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.suptitle(f"Lifecycle: {self.material.material_name} liner", fontsize=13, fontweight="bold")
        fig.tight_layout()
        if show:
            plt.show()
        return fig

    def plot_ratcheting(self, show=True):
        """Thermal ratcheting: Bree map, where it happens, and throat size / chamber pressure vs fires."""
        import matplotlib.pyplot as plt
        if not hasattr(self, "cycles_min"):
            self.solve()
        st = self.ratchet_stations()
        res = self.throat_ratcheting(tuple(range(0, MAX_FIRES + 1)))
        fig, axs = plt.subplots(2, 2, figsize=(13, 9))
        colors = {"elastic": "tab:green", "shakedown": "tab:olive", "ratcheting": "tab:red", "collapse": "k"}

        ax = axs[0, 0]
        x_max = min(max(1.2, float(np.max(st["X"])) * 1.05), 3.0)
        X = np.linspace(0, x_max, 300)
        ax.fill_between(X, 0, np.clip(1 - X, 0, None), color="tab:green", alpha=0.15, label="Elastic")
        ax.fill_between(X, np.clip(1 - X, 0, None), np.clip(4 * (1 - X), 0, None), where=X <= 1.0, color="tab:olive", alpha=0.15,
                        label="Shakedown")
        ax.fill_between(X, np.clip(4 * (1 - X), 0, None), 3.0, where=X <= 1.0, color="tab:red", alpha=0.12, label="Ratcheting")
        ax.fill_between(X, 0, 3.0, where=X >= 1.0, color="k", alpha=0.12, label="Primary stress > yield")
        for reg, col in colors.items():
            sel = st["regime"] == reg
            if sel.any():
                ax.scatter(st["X"][sel], st["Y"][sel], s=14, color=col, edgecolor="none")
        if res is not None:
            ax.scatter([res["X"]], [res["Y"]], s=140, marker="*", color="blue", edgecolor="white", zorder=5, label="Throat")
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, min(max(1.2, float(np.max(st["Y"])) * 1.1), 3.0))
        ax.set_title("Ratcheting map (every station is a dot)")
        ax.set_xlabel("X = hoop pressure stress / yield")
        ax.set_ylabel("Y = thermal stress / yield")
        ax.legend(loc="upper right", fontsize=8)

        ax = axs[0, 1]
        ax.plot(self.z, st["d_eps_r"] * 100, color="C3")
        ax.set_title("Net hoop strain gained per fire, by axial position")
        ax.set_xlabel("Axial Position [in]")
        ax.set_ylabel("Ratchet strain per fire [%]")
        ax.grid(True, alpha=0.3)
        if res is not None:
            ax.axvline(res["z"], color="b", linestyle="--", label="Throat")
            ax.legend()

        if res is None:
            for ax in (axs[1, 0], axs[1, 1]):
                ax.text(0.5, 0.5, "Throat outside the regen circuit", ha="center", va="center", transform=ax.transAxes)
        else:
            n = np.array([r["fires"] for r in res["rows"]])
            rad = np.array([r["radius"] for r in res["rows"]]) * 1000
            ax = axs[1, 0]
            ax.plot(n, rad, color="C0")
            ax.axhline(res["r0"] * 1000, color="k", linestyle=":", label="Original radius")
            if res["d_eps_r"] == 0.0:
                ax.text(0.5, 0.5, "No ratcheting at the throat:\nradius does not change", ha="center", va="center",
                        transform=ax.transAxes, fontsize=11, color="tab:green")
            ax.set_title(f"Throat radius vs fires  ({res['regime'].upper()})")
            ax.set_xlabel("Number of fires")
            ax.set_ylabel("Throat radius [mm]")
            ax.grid(True, alpha=0.3)
            ax.legend()
            ax = axs[1, 1]
            ax.plot(n, [r["d_area_pct"] for r in res["rows"]], color="C1", label="Throat area change [%]")
            ax.plot(n, [r["d_Pc_pct"] for r in res["rows"]], color="C3", label="Chamber pressure change [%]")
            ax.axhline(0.0, color="k", linewidth=0.8)
            ax.set_title("Effect on the engine")
            ax.set_xlabel("Number of fires")
            ax.set_ylabel("Change [%]")
            ax.grid(True, alpha=0.3)
            ax.legend()

        fig.suptitle("Thermal ratcheting (approximate, two-bar/Bree model)", fontsize=13, fontweight="bold")
        fig.tight_layout()
        if show:
            plt.show()
        return fig

    def export_csv(self, prefix="liner"):
        """Write plain CSV files (no extra packages): key metrics, yielding vs fires, throat vs fires, per-station."""
        import csv
        if not hasattr(self, "cycles_min"):
            self.solve()
        files = []

        def write(name, header, rows):
            path = f"{prefix}_{name}.csv"
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(header)
                w.writerows(rows)
            files.append(path)

        write("key_metrics", ["metric", "value"], list(self.key_metrics().items()))
        write("yielding_vs_fires",
              ["fires", "peak_plastic_strain_pct", "z_in", "location", "stations_yielded", "stations", "ductility_used_pct", "fatigue_life_used_pct"],
              [[r["fires"], r["peak_strain"] * 100, r["z"], LOCATION_LABELS[r["loc"]], r["stations_yielded"], r["stations"],
                r["ductility_used"] * 100, r["fatigue_used"] * 100] for r in self.yielding_history()])
        res = self.throat_ratcheting()
        if res is not None:
            write("throat_vs_fires", ["fires", "throat_radius_mm", "change_um", "throat_area_change_pct", "chamber_pressure_change_pct"],
                  [[0, res["r0"] * 1000, 0.0, 0.0, 0.0]] + [[r["fires"], r["radius"] * 1000, r["d_radius"] * 1e6, r["d_area_pct"], r["d_Pc_pct"]] for r in res["rows"]])
        st = self.ratchet_stations()
        write("per_station",
              ["z_in", "hot_wall_K", "cold_wall_K", "cycles_to_failure", "allowable_cycles", "governing_location",
               "stress_over_yield", "ratchet_X", "ratchet_Y", "ratchet_strain_per_fire_pct", "ratchet_regime"],
              [[self.z[k], self.T_hw[k], self.T_cw[k], self.cycles_min[k], self.allowable_cycles[k], self.governing_loc[k],
                max(self.yield_ratio[l][k] for l in LOCATIONS), st["X"][k], st["Y"][k], st["d_eps_r"][k] * 100, st["regime"][k]]
               for k in range(len(self.z))])
        print("Saved CSVs: " + ", ".join(files))
        return files

    # ------------------------------------------------------------------ reporting
    def summary(self):
        if not hasattr(self, "cycles_min"):
            self.solve()
        self.print_scorecard()
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
        violations = self.check_redlines()
        if violations:
            print("!!!! REDLINES EXCEEDED !!!!")
            for v in violations:
                print(f"  - {v}")
            print_explanations(explain_structural(self))
        else:
            print("No stress/life redlines exceeded.")
        self.print_yielding_history()
        self.print_throat_ratcheting()
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
        masks = self._redline_masks()
        if self.sigma_jacket_hoop is not None:
            draw_redline(axs[0], z, self.sigma_jacket_hoop * PA_TO_KSI, self.jacket_yield * PA_TO_KSI,
                         "Jacket yield", mask=masks["Jacket hoop stress above yield"], banner_slot=0)
        axs[0].set_title("Stress components")
        axs[0].set_ylabel("Stress [ksi]")
        axs[0].legend()

        for loc in LOCATIONS:
            axs[1].plot(z, self.sigma_vm_elastic[loc] * PA_TO_KSI, label=LOCATION_LABELS[loc])
        vm_max = np.max(np.vstack([np.abs(self.sigma_vm_elastic[loc]) for loc in LOCATIONS]), axis=0)
        draw_redline(axs[1], z, vm_max * self.Kt * PA_TO_KSI, self.sigma_yield["hot_mid"] * PA_TO_KSI,
                     "Yield at gas-face T", mask=masks["Elastic stress above yield (plastic cycling)"], banner_slot=0)
        axs[1].set_title("Elastic von Mises equivalent (signed)")
        axs[1].set_ylabel("Stress [ksi]")
        axs[1].legend()

        for loc in LOCATIONS:
            axs[2].semilogy(z, np.maximum(self.cycles[loc], 0.5), label=LOCATION_LABELS[loc])
        axs[2].semilogy(z, np.maximum(self.allowable_cycles, 0.5), "k--", label="Allowable (with safety factor)")
        if self.required_cycles is not None:
            draw_redline(axs[2], z, np.maximum(self.allowable_cycles, 0.5), self.required_cycles,
                         f"Required cycles ({self.required_cycles:g})", above=False,
                         mask=masks[f"Allowable cycles below required ({self.required_cycles:g})"], banner_slot=0)
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
        ws2 = wb.create_sheet("Yielding vs Fires")
        ws2.append(["Fires", "Peak plastic strain (-)", "z (in)", "Location", "Stations yielded", "Stations", "Ductility used (-)", "Fatigue life used (-)"])
        for r in self.yielding_history():
            ws2.append([r["fires"], r["peak_strain"], r["z"], LOCATION_LABELS[r["loc"]], r["stations_yielded"],
                        r["stations"], r["ductility_used"], r["fatigue_used"]])
        res = self.throat_ratcheting()
        if res is not None:
            ws3 = wb.create_sheet("Throat Ratcheting")
            ws3.append(["Regime", res["regime"], "Strain per fire (-)", res["d_eps_r"], "z (in)", res["z"]])
            ws3.append(["Fires", "Throat radius (m)", "Change (m)", "Throat area change (%)", "Chamber pressure change (%)"])
            ws3.append([0, res["r0"], 0.0, 0.0, 0.0])
            for r in res["rows"]:
                ws3.append([r["fires"], r["radius"], r["d_radius"], r["d_area_pct"], r["d_Pc_pct"]])
        wb.save(filename)
        print(f"Saved liner stress/life export: {filename}")
