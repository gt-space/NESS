"""
mdot_analysis.py - how well does THIS thrust chamber handle different propellant flow rates?

This is an off-design check, NOT a redesign: the contour, throat, channel geometry, wall thickness and
material are all held fixed. You give the oxidizer and fuel mass flows (or sweep them) and the code works out
what the engine does with them:

    1. Chamber pressure from the flows:   Pc = mdot_total * c*(Pc, MR) / At      (fixed-point on CEA c*)
    2. Mixture ratio:                     MR = mdot_ox / mdot_fuel
    3. Gas properties along the contour at that Pc and MR (CEA), via Engine.calcGasProperties
    4. Regen circuit re-solved with the fuel (or a flow you choose) as coolant, inlet T and P held at your inputs
    5. Liner stress, fatigue life, yielding and throat ratcheting re-evaluated
    6. Performance estimate (thrust, Isp) and flags for anything that would fail

Same convention as Engine.calcGasProperties: c* is the ideal CEA value, so the design flows reproduce the design Pc.
Thrust uses the ideal CEA thrust coefficient (no nozzle efficiency), so treat thrust/Isp as relative trends.

Main entry points
-----------------
OffDesignAnalyzer(engine, regen, solve_kwargs, ...)
    .operating_point(ox_mdot, fu_mdot)        one point (leaves engine + regen at that state for outputs())
    .throttle_sweep(ox_ref, fu_ref, fractions, MR_values=None)
    .restore()                                 put engine and regen back at the design point
plot_throttle_sweep(sweep, ...)               performance/thermal and structural dashboards
print_sweep_table(sweep), export_sweep_csv(sweep, path)
"""

import contextlib
import csv
import io

import numpy as np

from utils import _cea_obj_w_units_like_engine
from redlines import draw_redline, RED

PSI = 6894.76
LBM_PER_KG = 2.2046226218
SEPARATION_FRACTION = 0.4          # Summerfield: exit pressure below 0.4 * Pamb risks flow separation


def to_kg_per_s(value, units="kg/s"):
    """Convert a mass flow given in kg/s or lbm/s to kg/s (None passes through)."""
    if value is None:
        return None
    u = units.strip().lower().replace(" ", "")
    if u in ("kg/s", "kgs", "kg/sec"):
        return float(value)
    if u in ("lbm/s", "lb/s", "lbms", "lbm/sec"):
        return float(value) / LBM_PER_KG
    raise ValueError(f"Unknown mdot units '{units}'. Use 'kg/s' or 'lbm/s'.")


def chamber_pressure_from_mdots(engine, ox_mdot, fu_mdot, Pc_guess=None):
    """
    Chamber pressure [psia] that the throat will hold for these flows (kg/s), from Pc = mdot * c* / At with
    c* = CEA c*(Pc, MR) iterated to convergence. Returns (Pc_psia, MR, cstar [m/s]).
    """
    if ox_mdot <= 0 or fu_mdot <= 0:
        raise ValueError(f"Both flows must be positive (ox {ox_mdot}, fuel {fu_mdot} kg/s).")
    C = _cea_obj_w_units_like_engine(engine.oxName, engine.fuName)
    MR = ox_mdot / fu_mdot
    tot = ox_mdot + fu_mdot
    At_m2 = engine.At * 0.00064516
    Pc = float(Pc_guess if Pc_guess is not None else engine.Pc_psia)
    cstar = None
    for _ in range(60):
        cstar = float(C.get_Cstar(Pc, MR)) * engine.flow_cstar_factor()
        Pc_new = tot * cstar / At_m2 / PSI
        if abs(Pc_new - Pc) < 1e-7 * Pc_new:
            Pc = Pc_new
            break
        Pc = Pc_new
    return Pc, MR, cstar


def injector_stiffness(throttle, points):
    """
    Injector stiffness (injector pressure drop / chamber pressure) at a throttle fraction, from a linear schedule
    through `points` = ((throttle, stiffness), ...), e.g. ((0.40, 0.15), (1.00, 0.30)). Outside the given
    throttles the end segment is extended linearly (never below 0).
    """
    pts = sorted((float(a), float(b)) for a, b in points)
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])
    if len(pts) == 1:
        return float(max(ys[0], 0.0))
    if throttle < xs[0]:
        val = ys[0] + (throttle - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
    elif throttle > xs[-1]:
        val = ys[-1] + (throttle - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
    else:
        val = float(np.interp(throttle, xs, ys))
    return float(max(val, 0.0))


def _solve_regen_quiet(regen, coolant_mdot, solve_kwargs, inlet_pressure=None, quiet=True):
    """Solve the regen circuit at a coolant flow (and optionally a different inlet pressure) without limit exceptions."""
    regen.tot_coolant_mdot = float(coolant_mdot)
    regen.chan_mdot = regen.tot_coolant_mdot / regen.N
    kwargs = dict(solve_kwargs)
    if inlet_pressure is not None:
        kwargs["inlet_pressure"] = float(inlet_pressure)
    kwargs.update(show_progress=False, print_station_rows=False, enforce_melting_limit=False, enforce_boiling_limit=False)
    sink = io.StringIO()
    with (contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()):
        regen.solve_circuit(**kwargs)


def inlet_pressure_for_injector(regen, Pc_psia, stiffness, solve_kwargs, coolant_mdot, quiet=True, tol_psi=0.05):
    """
    Regen circuit inlet pressure that leaves the coolant at injector-inlet pressure Pc * (1 + stiffness) at the
    circuit outlet. The fuel flows through the channels first and then the injector, so

        P_regen_inlet = Pc * (1 + stiffness) + (pressure drop through the channels)

    and the channel pressure drop depends (weakly) on the inlet pressure, so it is iterated. The regen circuit is
    left solved at the converged inlet pressure. Returns a dict (pressures in Pa).
    """
    P_chamber = float(Pc_psia) * PSI
    P_target = P_chamber * (1.0 + float(stiffness))
    P_in = P_target
    for _ in range(10):
        _solve_regen_quiet(regen, coolant_mdot, solve_kwargs, inlet_pressure=P_in, quiet=quiet)
        dP = float(np.sum(regen.DP_arr))
        P_new = P_target + dP
        if abs(P_new - P_in) < tol_psi * PSI:
            P_in = P_new
            break
        P_in = P_new
    return dict(P_in_Pa=P_in, P_inj_in_Pa=P_target, dP_regen_Pa=float(np.sum(regen.DP_arr)),
                dP_inj_Pa=P_target - P_chamber, stiffness=float(stiffness))


def _nan():
    return float("nan")


class OffDesignAnalyzer:
    """
    Evaluate a fixed RegenCircuit/Engine at user-chosen propellant flows.

    Parameters
    ----------
    engine, regen : the design-point Engine and RegenCircuit (regen.solve_circuit already worked at design).
    solve_kwargs : dict passed to regen.solve_circuit (inlet_T_c, inlet_pressure, T_wg_bracket, ...).
    structural_kwargs : dict passed to regen.run_structural_analysis (Kt, life_safety_factor, required_cycles, t_jacket).
    T_hw_limit, dp_limit_psid : your own design limits (optional) used for failure flags.
    """

    def __init__(self, engine, regen, solve_kwargs, structural_kwargs=None, T_hw_limit=None, dp_limit_psid=None,
                 coolant_is_fuel=True, injector_stiffness_points=None):
        self.engine = engine
        self.regen = regen
        self.solve_kwargs = dict(solve_kwargs)
        self.structural_kwargs = dict(structural_kwargs or {})
        self.T_hw_limit = T_hw_limit
        self.dp_limit_psid = dp_limit_psid
        self.coolant_is_fuel = coolant_is_fuel
        # ((throttle, injector dP / Pc), ...): when given, the regen inlet pressure follows the throttle (see inlet_pressure_for_injector)
        self.injector_points = injector_stiffness_points
        self.required_cycles = self.structural_kwargs.get("required_cycles")
        # State the engine and regen circuit are in now (the "primary" run: the design point, or your actual flows if the
        # script ran the main solve at those). restore() puts everything back here.
        self.primary_Pc_psia = float(engine.Pc_psia)
        self.primary_MR = float(engine.MRcore)
        self.primary_fu_mdot = float(engine.fu_mdot)
        self.primary_ox_mdot = float(engine.ox_mdot)
        self.primary_coolant_mdot = float(regen.tot_coolant_mdot)
        # The point the hardware was actually designed for (marked on the sweep graphs)
        ds = getattr(engine, "design_state", None) or dict(
            Pc_psia=self.primary_Pc_psia, MR=self.primary_MR, ox_mdot=self.primary_ox_mdot, fu_mdot=self.primary_fu_mdot)
        self.design_ox_mdot = float(ds["ox_mdot"])
        self.design_fu_mdot = float(ds["fu_mdot"])

    # ------------------------------------------------------------------ helpers
    def _solve_regen(self, coolant_mdot, quiet=True, inlet_pressure=None):
        _solve_regen_quiet(self.regen, coolant_mdot, self.solve_kwargs, inlet_pressure=inlet_pressure, quiet=quiet)

    def _performance(self, Pc_psia, MR, ox_mdot, fu_mdot):
        eng = self.engine
        thrust_lbf = eng.cea_thrust_lbf(Pc_psia, MR)
        thrust_N = thrust_lbf * 4.4482216
        isp = thrust_N / ((ox_mdot + fu_mdot) * 9.80665)
        return thrust_lbf, isp, float(eng.P[-1])

    # ------------------------------------------------------------------ one operating point
    def operating_point(self, ox_mdot, fu_mdot, coolant_mdot=None, quiet=True, throttle=None):
        """
        Run one operating point (flows in kg/s). The engine and regen circuit are left at this state, so
        regen.outputs() / the structural analysis can be shown afterwards; call restore() to go back to design.
        Returns a dict of results (NaN where a stage could not run) with 'status' in {"ok", "FAIL", "error"}.
        """
        eng, regen = self.engine, self.regen
        res = dict(ox_mdot=float(ox_mdot), fu_mdot=float(fu_mdot), tot_mdot=float(ox_mdot + fu_mdot),
                   status="ok", error="", failures=[], warnings=[])
        for key in ("Pc_psia", "MR", "cstar", "thrust_lbf", "thrust_pct_design", "stiffness", "inlet_P_psia", "inj_inlet_P_psia", "inj_dP_psi", "isp_s", "Pe_psia", "T_chamber_K", "coolant_mdot",
                    "T_hw_max", "T_c_out", "DP_psid", "q_max_MW", "boil_margin_K", "vel_max", "Re_min",
                    "P_margin_psi", "melt_margin_K", "stress_over_yield", "cycles_to_failure", "allowable_cycles",
                    "plastic_1", "plastic_10", "plastic_20", "life_used_20_pct", "ratchet_regime", "mission_ok",
                    "ratchet_strain_pct", "throat_dr_um_20", "throat_dArea_pct_20", "throat_dPc_pct_20"):
            res[key] = _nan() if key != "ratchet_regime" else ""

        try:
            Pc, MR, cstar = chamber_pressure_from_mdots(eng, ox_mdot, fu_mdot)
            res.update(Pc_psia=Pc, MR=MR, cstar=cstar)
            eng.calcGasProperties(frozen=getattr(eng, "frozen", 1), Pc_psia_override=Pc, MR_override=MR)
            res["T_chamber_K"] = float(eng.Tcomb)
            res["thrust_lbf"], res["isp_s"], res["Pe_psia"] = self._performance(Pc, MR, ox_mdot, fu_mdot)
            res["thrust_pct_design"] = 100.0 * res["thrust_lbf"] / float(eng.thrust)

            cool = coolant_mdot if coolant_mdot is not None else (fu_mdot if self.coolant_is_fuel else ox_mdot)
            res["coolant_mdot"] = float(cool)
            if self.injector_points is not None and throttle is not None:
                s_inj = injector_stiffness(throttle, self.injector_points)
                feed = inlet_pressure_for_injector(self.regen, Pc, s_inj, self.solve_kwargs, cool, quiet=quiet)
                res.update(stiffness=s_inj, inlet_P_psia=feed["P_in_Pa"] / PSI, inj_inlet_P_psia=feed["P_inj_in_Pa"] / PSI,
                           inj_dP_psi=feed["dP_inj_Pa"] / PSI)
            else:
                self._solve_regen(cool, quiet=quiet)
                res["inlet_P_psia"] = float(self.solve_kwargs.get("inlet_pressure", float("nan"))) / PSI
        except Exception as ex:                                   # CEA out of range, non-physical pressure, ...
            res["status"] = "error"
            res["error"] = f"{type(ex).__name__}: {ex}"
            return res

        # ---- thermal / hydraulic metrics ----
        o, i = regen.circuit_outlet, regen.circuit_inlet
        P_gas = np.asarray(eng.P)[o:i + 1] * PSI
        res["T_hw_max"] = float(np.max(regen.T_hw_arr))
        res["T_c_out"] = float(np.max(regen.T_c_arr))
        res["DP_psid"] = float(np.sum(regen.DP_arr)) / PSI
        res["q_max_MW"] = float(np.max(regen.q_flux_arr)) / 1e6
        res["boil_margin_K"] = float(np.min(regen.coolant_boiling_arr - regen.T_c_arr))
        res["vel_max"] = float(np.max(regen.coolant_vel_arr))
        res["Re_min"] = float(np.min(regen.Re))
        res["P_margin_psi"] = float(np.min(regen.P_c_in_arr - P_gas)) / PSI
        melt = getattr(regen.material, "MeltingPoint", None)
        res["melt_margin_K"] = _nan() if melt is None else float(melt - res["T_hw_max"])

        # ---- structural metrics ----
        try:
            sa = regen.run_structural_analysis(**self.structural_kwargs)
            km = sa.key_metrics()
            res["stress_over_yield"] = km["max_stress_over_yield"]
            res["cycles_to_failure"] = km["cycles_to_failure"]
            res["allowable_cycles"] = km["allowable_cycles"]
            res["plastic_1"] = km["plastic_strain_pct_1_fire"]
            res["plastic_10"] = km["plastic_strain_pct_10_fires"]
            res["plastic_20"] = km["plastic_strain_pct_20_fires"]
            res["life_used_20_pct"] = km["fatigue_life_used_pct_20_fires"]
            res["ratchet_regime"] = km["throat_regime"]
            res["mission_ok"] = bool(km["survives_mission"])
            res["mission_fires"] = km["mission_fires"]
            if "throat_strain_per_fire_pct" in km:
                res["ratchet_strain_pct"] = km["throat_strain_per_fire_pct"]
                res["throat_dr_um_20"] = km["throat_radius_change_um_20_fires"]
                res["throat_dArea_pct_20"] = km["throat_area_change_pct_20_fires"]
                res["throat_dPc_pct_20"] = km["chamber_pressure_change_pct_20_fires"]
        except Exception as ex:
            res["warnings"].append(f"structural analysis skipped: {type(ex).__name__}: {ex}")

        self._flag(res)
        return res

    def _flag(self, res):
        """Failures (would hurt the engine) and warnings (worth knowing) for one result."""
        f, w = res["failures"], res["warnings"]
        melt = res["melt_margin_K"]
        if np.isfinite(melt) and melt < 0:
            f.append(f"hot wall {-melt:.0f} K above melting")
        if self.T_hw_limit is not None and res["T_hw_max"] > self.T_hw_limit:
            f.append(f"hot wall {res['T_hw_max'] - self.T_hw_limit:.0f} K above your {self.T_hw_limit:.0f} K limit")
        if res["boil_margin_K"] < 0:
            f.append("coolant above boiling limit")
        if res["P_margin_psi"] < 0:
            f.append("coolant pressure below gas pressure")
        if self.dp_limit_psid is not None and res["DP_psid"] > self.dp_limit_psid:
            f.append(f"channel dP {res['DP_psid']:.1f} psid above {self.dp_limit_psid:g}")
        if np.isfinite(res["stress_over_yield"]) and res["stress_over_yield"] > 1.0:
            w.append(f"stress {res['stress_over_yield']:.2f}x yield")
        if np.isfinite(res["allowable_cycles"]) and self.required_cycles is not None \
                and res["allowable_cycles"] < self.required_cycles:
            f.append(f"life {res['allowable_cycles']:.0f} cycles < required {self.required_cycles:g}")
        if res["ratchet_regime"] == "ratcheting":
            f.append("throat ratcheting (throat size drifts every fire)")
        elif res["ratchet_regime"] == "collapse":
            w.append("throat hoop stress above yield (plastic-collapse trend; the model is very conservative here)")
        if res["Pe_psia"] < SEPARATION_FRACTION * self.engine.Pamb:
            w.append("exit pressure below 0.4 x Pamb: flow separation risk")
        if res["Re_min"] < 2300:
            w.append("coolant laminar/transitional somewhere")
        if f:
            res["status"] = "FAIL"

    # ------------------------------------------------------------------ throttle sweep
    def throttle_sweep(self, ox_ref, fu_ref, fractions, MR_values=None, progress=True):
        """
        Sweep total mass flow as fractions of (ox_ref + fu_ref), flows in kg/s.
        Default holds the reference mixture ratio; pass MR_values to sweep several mixture ratios.
        Engine and regen are restored to the design point afterwards.
        """
        from tqdm import tqdm
        fractions = np.asarray(fractions, dtype=float)
        tot_ref = float(ox_ref + fu_ref)
        mr_ref = float(ox_ref / fu_ref)
        mr_list = [mr_ref] if MR_values is None else [float(m) for m in MR_values]
        results = []
        pts = [(m, f) for m in mr_list for f in fractions]
        design = None
        try:
            for m, f in tqdm(pts, desc="Throttle sweep", disable=not progress):
                tot = tot_ref * f
                fu = tot / (1.0 + m)
                ox = tot - fu
                r = self.operating_point(ox, fu, quiet=True, throttle=f)
                r["throttle"] = float(f)
                r["MR_case"] = float(m)
                results.append(r)
            design = self.operating_point(self.design_ox_mdot, self.design_fu_mdot, quiet=True,
                                          throttle=(self.design_ox_mdot + self.design_fu_mdot) / tot_ref)
            design["throttle"] = (self.design_ox_mdot + self.design_fu_mdot) / tot_ref
            design["MR_case"] = self.design_ox_mdot / self.design_fu_mdot
        finally:
            self.restore()
        return dict(fractions=fractions, MR_values=np.array(mr_list), results=results, design=design,
                    ox_ref=float(ox_ref), fu_ref=float(fu_ref), analyzer=self)

    def restore(self):
        """Put the engine and regen circuit back at the design point."""
        eng = self.engine
        eng.calcGasProperties(frozen=getattr(eng, "frozen", 1), Pc_psia_override=self.primary_Pc_psia,
                              MR_override=self.primary_MR)
        self._solve_regen(self.primary_coolant_mdot, quiet=True)
        try:
            self.regen.run_structural_analysis(**self.structural_kwargs)
        except Exception:
            pass


# ======================================================================================================
#  Reporting
# ======================================================================================================

def print_operating_point(res, analyzer=None):
    """Readable summary of one operating point."""
    print("\n==================== OFF-DESIGN OPERATING POINT ====================")
    print(f"  Flows            : ox {res['ox_mdot']:.4f} kg/s ({res['ox_mdot'] * LBM_PER_KG:.4f} lbm/s), "
          f"fuel {res['fu_mdot']:.4f} kg/s ({res['fu_mdot'] * LBM_PER_KG:.4f} lbm/s)")
    if res["status"] == "error":
        print(f"  COULD NOT RUN    : {res['error']}")
        return
    print(f"  Chamber pressure : {res['Pc_psia']:.1f} psia   Mixture ratio: {res['MR']:.3f}   c*: {res['cstar']:.0f} m/s")
    print(f"  Performance      : thrust {res['thrust_lbf']:.0f} lbf, Isp {res['isp_s']:.1f} s (ideal CF), "
          f"exit pressure {res['Pe_psia']:.2f} psia")
    print(f"  Coolant flow     : {res['coolant_mdot']:.4f} kg/s")
    print(f"  Hot wall (max)   : {res['T_hw_max']:.0f} K"
          + ("" if not np.isfinite(res["melt_margin_K"]) else
             (f"  ({res['melt_margin_K']:.0f} K below melting)" if res["melt_margin_K"] >= 0 else f"  ({-res['melt_margin_K']:.0f} K ABOVE melting)")))
    print(f"  Coolant          : outlet {res['T_c_out']:.0f} K, boiling margin {res['boil_margin_K']:.0f} K, "
          f"channel dP {res['DP_psid']:.1f} psid, max velocity {res['vel_max']:.1f} m/s, min Re {res['Re_min']:.0f}")
    print(f"  Heat flux (max)  : {res['q_max_MW']:.2f} MW/m^2")
    if np.isfinite(res["stress_over_yield"]):
        print(f"  Liner            : stress {res['stress_over_yield']:.2f}x yield, cycles to failure "
              f"{res['cycles_to_failure']:.3g} (allowable {res['allowable_cycles']:.3g}), "
              f"plastic strain {res['plastic_20']:.2f}% after 20 fires")
        print(f"  Throat           : {res['ratchet_regime']}, radius {res['throat_dr_um_20']:+.1f} um, "
              f"area {res['throat_dArea_pct_20']:+.2f}%, Pc {res['throat_dPc_pct_20']:+.2f}% after 20 fires")
    if res["failures"]:
        print("  !!!! WOULD FAIL  : " + "; ".join(res["failures"]))
    else:
        print("  No failure limits crossed.")
    for w in res["warnings"]:
        print(f"  note: {w}")
    print("====================================================================")


def _series(sweep, mr):
    rows = [r for r in sweep["results"] if np.isclose(r["MR_case"], mr)]
    rows.sort(key=lambda r: r["throttle"])
    return rows


def _col(rows, key):
    return np.array([r[key] for r in rows], dtype=float)


def _pct(x):
    return "-" if x is None or not np.isfinite(x) else f"{x * 100:.0f}%"


def _mission_text(r):
    m = r.get("mission_ok")
    return "-" if m is None or (isinstance(m, float) and not np.isfinite(m)) else ("YES" if m else "NO")


def print_sweep_table(sweep):
    """Compact table per mixture ratio: where it works and where it fails."""
    an = sweep["analyzer"]
    d = sweep.get("design")
    print("\n==== THROTTLE SWEEP REFERENCE ====")
    print(f"  100% = the flows you gave: ox {sweep['ox_ref']:.4f} kg/s + fuel {sweep['fu_ref']:.4f} kg/s "
          f"(total {sweep['ox_ref'] + sweep['fu_ref']:.4f} kg/s, MR {sweep['ox_ref'] / sweep['fu_ref']:.3f})")
    if d is not None and d["status"] != "error":
        print(f"  Design point : ox {d['ox_mdot']:.4f} + fuel {d['fu_mdot']:.4f} kg/s = {d['throttle'] * 100:.0f}% of that reference -> "
              f"Pc {d['Pc_psia']:.0f} psia, thrust {d['thrust_lbf']:.0f} lbf ({d['thrust_pct_design']:.0f}% of the {an.engine.thrust:.0f} lbf request)")
    for mr in sweep["MR_values"]:
        rows = _series(sweep, mr)
        print(f"\n---- THROTTLE SWEEP  (MR = {mr:.3f}) ----")
        print(f"{'Throttle':>8} | {'Pc [psia]':>9} | {'Thrust [lbf]':>12} | {'Isp [s]':>7} | {'T_hw [K]':>8} | "
              f"{'% design thrust':>15} | {'Inj stiff':>9} | {'Inlet P [psia]':>14} | {'Boil margin':>11} | {'dP [psid]':>9} | {'Life [cyc]':>10} | {'Survives':>8} | Status")
        for r in rows:
            status = r["status"] if r["status"] != "FAIL" else "FAIL: " + "; ".join(r["failures"])
            if r["status"] == "error":
                status = "ERROR: " + r["error"][:60]
            print(f"{r['throttle'] * 100:>7.0f}% | {r['Pc_psia']:>9.1f} | {r['thrust_lbf']:>12.0f} | {r['isp_s']:>7.1f} | "
                  f"{r['T_hw_max']:>8.0f} | {r['thrust_pct_design']:>14.0f}% | {_pct(r.get('stiffness')):>9} | {r['inlet_P_psia']:>14.0f} | {r['boil_margin_K']:>9.0f} K | {r['DP_psid']:>9.1f} | "
                  f"{r['allowable_cycles']:>10.3g} | {_mission_text(r):>8} | {status}")
        n_m = next((r["mission_fires"] for r in rows if np.isfinite(r.get("mission_fires", float("nan")))), None)
        surv = [r["throttle"] for r in rows if r.get("mission_ok") is True]
        if n_m is not None:
            if surv:
                print(f"  Survives {n_m:g} hotfires from {min(surv) * 100:.0f}% to {max(surv) * 100:.0f}% of the reference flow "
                      f"(lowest throttle that survives: {min(surv) * 100:.0f}%).")
            else:
                print(f"  No throttle point in this sweep survives {n_m:g} hotfires.")
        ok = [r["throttle"] for r in rows if r["status"] == "ok"]
        if ok:
            print(f"  Works without crossing any limit from {min(ok) * 100:.0f}% to {max(ok) * 100:.0f}% of the reference flow.")
        else:
            print("  No throttle point in this sweep stays inside every limit.")


def export_sweep_csv(sweep, path="throttle_sweep.csv"):
    keys = ["throttle", "MR_case", "ox_mdot", "fu_mdot", "coolant_mdot", "Pc_psia", "MR", "stiffness", "inlet_P_psia", "inj_dP_psi", "thrust_lbf", "isp_s", "Pe_psia",
            "T_hw_max", "T_c_out", "DP_psid", "q_max_MW", "boil_margin_K", "vel_max", "Re_min", "P_margin_psi",
            "melt_margin_K", "stress_over_yield", "cycles_to_failure", "allowable_cycles", "plastic_1", "plastic_10",
            "plastic_20", "life_used_20_pct", "ratchet_regime", "mission_ok", "mission_fires", "ratchet_strain_pct", "throat_dr_um_20",
            "throat_dArea_pct_20", "throat_dPc_pct_20", "status"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys + ["failures", "warnings"])
        for r in sweep["results"]:
            w.writerow([r.get(k, "") for k in keys] + ["; ".join(r["failures"]), "; ".join(r["warnings"] + ([r["error"]] if r["error"] else []))])
    print(f"Saved sweep CSV: {path}")
    return path


# ======================================================================================================
#  Graphs
# ======================================================================================================

def plot_throttle_sweep(sweep, show=True):
    """
    Two dashboards from a throttle sweep:
      1. Performance and thermal: Pc, thrust/Isp, exit pressure, hot wall T, boiling margin, heat flux,
         channel dP, coolant velocity/Re, coolant-over-gas pressure margin
      2. Structural: stress/yield, cycles to failure, plastic strain, fatigue life used, throat ratcheting
    Failures are drawn in red (limit line, X markers, shaded span, banner). Points that could not be run are gaps.
    """
    import matplotlib.pyplot as plt

    an = sweep["analyzer"]
    mrs = sweep["MR_values"]
    colors = [f"C{k}" for k in range(len(mrs))]
    multi = len(mrs) > 1
    melt = getattr(an.regen.material, "MeltingPoint", None)
    Pamb = an.engine.Pamb

    def lines(ax, key, ylabel, title, scale=1.0, fmt="o-", limit=None, limit_label=None, above=True, extra=None):
        """Plot one metric for every MR; draw a redline if a limit is given."""
        for slot, (mr, col) in enumerate(zip(mrs, colors)):
            rows = _series(sweep, mr)
            x = _col(rows, "throttle") * 100
            y = _col(rows, key) * scale
            ax.plot(x, y, fmt, color=col, markersize=4, label=f"MR {mr:.2f}" if multi else None)
            if limit is not None:
                draw_redline(ax, x, y, limit, limit_label, above=above, draw_limit=(slot == 0), banner_slot=slot,
                             x_name="throttle", x_unit="%")
            fail_mask = np.array([r["status"] in ("FAIL", "error") for r in rows])
            if extra is None and fail_mask.any():
                pass
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Throttle [% of reference flow]")
        ax.grid(True, alpha=0.3)
        if multi or limit is not None:
            ax.legend(fontsize=7, loc="best")

    tot_ref = sweep["ox_ref"] + sweep["fu_ref"]
    design = sweep.get("design")
    x_design = None if design is None or design["status"] == "error" else design["throttle"] * 100

    def decorate(fig_axes, label_top=True):
        """Mark the design flow and add a total-mdot axis along the top of every graph."""
        for k, ax in enumerate(np.atleast_1d(fig_axes).ravel()):
            if x_design is not None:
                lo, hi = ax.get_xlim()
                if lo <= x_design <= hi:
                    ax.axvline(x_design, color="k", linestyle="-.", linewidth=1.0, alpha=0.7)
                    ax.text(x_design, 0.02, " design flow", rotation=90, va="bottom", ha="right", fontsize=7,
                            transform=ax.get_xaxis_transform(), alpha=0.8)
            sec = ax.secondary_xaxis("top", functions=(lambda p: p / 100.0 * tot_ref, lambda m: m / tot_ref * 100.0))
            sec.tick_params(labelsize=7)
            if label_top and k < 3:
                sec.set_xlabel("Total mdot [kg/s]", fontsize=7)

    # ---------------------------------------------------------------- figure 1
    fig1, axs = plt.subplots(3, 3, figsize=(17, 12))
    lines(axs[0, 0], "Pc_psia", "Pc [psia]", "Chamber pressure (from the flows)")
    ax = axs[0, 1]
    lines(ax, "thrust_lbf", "Thrust [lbf]", "Thrust (left) and Isp (right), ideal CF")
    ax2 = ax.twinx()
    for mr, col in zip(mrs, colors):
        rows = _series(sweep, mr)
        ax2.plot(_col(rows, "throttle") * 100, _col(rows, "isp_s"), "s--", color=col, markersize=3, alpha=0.6)
    ax2.set_ylabel("Isp [s]")
    lines(axs[0, 2], "Pe_psia", "Exit pressure [psia]", f"Nozzle exit pressure (separation risk below 0.4 x Pamb)",
          limit=SEPARATION_FRACTION * Pamb, limit_label=f"Separation risk (0.4 x {Pamb:g} psia)", above=False)

    lines(axs[1, 0], "T_hw_max", "T_hw [K]", "Max hot wall temperature",
          limit=melt, limit_label=f"{an.regen.material.material_name} melting ({melt:.0f} K)" if melt else None)
    if an.T_hw_limit is not None:
        for slot, mr in enumerate(mrs):
            rows = _series(sweep, mr)
            draw_redline(axs[1, 0], _col(rows, "throttle") * 100, _col(rows, "T_hw_max"), an.T_hw_limit,
                         f"Your limit ({an.T_hw_limit:.0f} K)", draw_limit=(slot == 0), banner_slot=2 + slot,
                         x_name="throttle", x_unit="%")
        axs[1, 0].legend(fontsize=7)
    lines(axs[1, 1], "boil_margin_K", "Margin [K]", "Coolant boiling margin (min over the channel)",
          limit=0.0, limit_label="Coolant at its boiling limit", above=False)
    lines(axs[1, 2], "q_max_MW", "q [MW/m^2]", "Max gas-side heat flux")

    lines(axs[2, 0], "DP_psid", "dP [psid]", "Total channel pressure drop",
          limit=an.dp_limit_psid, limit_label=f"Max dP ({an.dp_limit_psid:g} psid)" if an.dp_limit_psid else None)
    ax = axs[2, 1]
    lines(ax, "vel_max", "Velocity [m/s]", "Max coolant velocity (left) and min Re (right)")
    ax2 = ax.twinx()
    for slot, (mr, col) in enumerate(zip(mrs, colors)):
        rows = _series(sweep, mr)
        x = _col(rows, "throttle") * 100
        y = _col(rows, "Re_min")
        ax2.plot(x, y, "s--", color=col, markersize=3, alpha=0.6)
        draw_redline(ax2, x, y, 2300.0, "Laminar below Re 2300", above=False, draw_limit=(slot == 0), banner_slot=slot,
                     x_name="throttle", x_unit="%")
    ax2.set_ylabel("Re_min")
    lines(axs[2, 2], "P_margin_psi", "P_coolant - P_gas [psi]", "Coolant pressure above chamber pressure (min)",
          limit=0.0, limit_label="Coolant below gas pressure", above=False)

    fig1.suptitle(f"Throttle sweep: performance and thermal  (geometry fixed; 100% = ox {sweep['ox_ref']:.3f} + fuel {sweep['fu_ref']:.3f} kg/s"
                  + ("" if x_design is None else f"; the {an.engine.thrust:.0f} lbf design flow is {x_design:.0f}%") + ")",
                  fontsize=13, fontweight="bold")
    decorate(axs)
    fig1.tight_layout()

    # ---------------------------------------------------------------- figure 2
    fig2, axs = plt.subplots(2, 3, figsize=(17, 8.5))
    lines(axs[0, 0], "stress_over_yield", "Stress / yield", "Worst elastic stress / yield (above 1 = yielding every fire)",
          limit=1.0, limit_label="Yield")
    lines(axs[0, 1], "allowable_cycles", "Cycles", "Allowable cycles to failure (with safety factor)",
          limit=an.required_cycles, limit_label=f"Required ({an.required_cycles:g})" if an.required_cycles else None, above=False)
    axs[0, 1].set_yscale("log")
    ax = axs[0, 2]
    for mr, col in zip(mrs, colors):
        rows = _series(sweep, mr)
        x = _col(rows, "throttle") * 100
        for key, ls, lab in (("plastic_1", ":", "after 1 fire"), ("plastic_10", "--", "after 10"), ("plastic_20", "-", "after 20")):
            ax.plot(x, _col(rows, key), ls, color=col, marker="o", markersize=3, label=(lab if mr == mrs[0] else None))
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_title("Peak plastic strain accumulated", fontsize=10)
    ax.set_ylabel("Plastic strain [%]")
    ax.set_xlabel("Throttle [% of reference flow]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    lines(axs[1, 0], "life_used_20_pct", "Fatigue life used [%]", "Fatigue life used after 20 fires",
          limit=100.0, limit_label="Life exhausted (100%)")
    lines(axs[1, 1], "throat_dr_um_20", "Radius change [um]", "Throat radius change after 20 fires (thermal ratcheting)")
    lines(axs[1, 2], "throat_dPc_pct_20", "Chamber pressure change [%]", "Chamber pressure change after 20 fires (throat closing)")

    fig2.suptitle("Throttle sweep: liner stress, life and ratcheting", fontsize=13, fontweight="bold")
    decorate(axs)
    fig2.tight_layout()

    # ---------------------------------------------------------------- figure 3: feed system
    fig3, axs3 = plt.subplots(1, 3, figsize=(17, 4.8))
    lines(axs3[0], "inlet_P_psia", "Pressure [psia]", "Regen circuit inlet pressure (what the feed system must supply)")
    if "inj_inlet_P_psia" in sweep["results"][0]:
        for mr, col in zip(mrs, colors):
            rows = _series(sweep, mr)
            axs3[0].plot(_col(rows, "throttle") * 100, _col(rows, "inj_inlet_P_psia"), "--", color=col, alpha=0.7)
            axs3[0].plot(_col(rows, "throttle") * 100, _col(rows, "Pc_psia"), ":", color=col, alpha=0.7)
        axs3[0].text(0.02, 0.97, "solid = regen inlet, dashed = injector inlet, dotted = chamber", transform=axs3[0].transAxes,
                     va="top", fontsize=7)
    lines(axs3[1], "stiffness", "Injector dP / Pc [%]", "Injector stiffness schedule", scale=100.0)
    lines(axs3[2], "inj_dP_psi", "Pressure drop [psid]", "Injector pressure drop (dashed: regen channels)")
    for mr, col in zip(mrs, colors):
        rows = _series(sweep, mr)
        axs3[2].plot(_col(rows, "throttle") * 100, _col(rows, "DP_psid"), "--", color=col, alpha=0.8)
    decorate(axs3)
    fig3.suptitle("Throttle sweep: feed system (regen inlet pressure follows the injector stiffness schedule)", fontsize=13, fontweight="bold")
    fig3.tight_layout()

    if show:
        plt.show()
    return fig1, fig2, fig3
