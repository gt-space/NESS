"""
redline_explain.py - plain-language explanations of what is physically happening where a redline is crossed.

explain_regen(regen, ...)   -> [(title, [lines])] for the thermal / pressure redlines of a solved RegenCircuit
explain_structural(sa)      -> [(title, [lines])] for the stress / life redlines of a LinerStructuralAnalysis
print_explanations(items)   -> pretty-prints either list
"""

import numpy as np

from redlines import span_text


def _section_name(regen, k):
    """Engine section for circuit station k (k = 0 at the injector end)."""
    eng = regen.engine
    i = regen.circuit_outlet + k
    if abs(i - eng.throat_ind) <= 2:
        return "the throat"
    if i <= eng.chamberEnd_ind:
        return "the combustion chamber"
    if i < eng.throat_ind:
        return "the converging section just upstream of the throat"
    return "the diverging nozzle"


def _fmt_span(regen, mask):
    return span_text(regen.circuit_z, mask)


def _wall_story(regen, k):
    """Where the temperature drop between hot gas and coolant is being spent at station k."""
    i = regen.circuit_outlet + k
    T_gas = float(regen.engine.T[i])
    T_hw, T_cw, T_c = float(regen.T_hw_arr[k]), float(regen.T_cw_arr[k]), float(regen.T_c_arr[k])
    d_gas, d_wall, d_cool = T_gas - T_hw, T_hw - T_cw, T_cw - T_c
    tot = max(d_gas + d_wall + d_cool, 1e-9)
    shares = {"gas film": d_gas / tot, "wall": d_wall / tot, "coolant film": d_cool / tot}
    dominant = max(shares, key=shares.get)
    return dict(
        T_gas=T_gas, T_hw=T_hw, T_cw=T_cw, T_c=T_c, d_gas=d_gas, d_wall=d_wall, d_cool=d_cool,
        shares=shares, dominant=dominant, q=float(regen.q_flux_arr[k]) / 1e6, h_hg=float(regen.h_hg_arr[k]),
        h_c=float(regen.h_c_arr[k]), Re=float(regen.Re[k]), v=float(regen.coolant_vel_arr[k]),
    )


def _wall_explanation(regen, k, what):
    w = _wall_story(regen, k)
    sec = _section_name(regen, k)
    lines = [
        f"At z = {regen.circuit_z[k]:.2f} in, in {sec}, the combustion gas ({w['T_gas']:.0f} K) pushes "
        f"{w['q']:.1f} MW/m^2 into the wall through a gas film with h = {w['h_hg']:.0f} W/m^2-K.",
        f"That heat has to cross three layers to reach the coolant: gas film {w['d_gas']:.0f} K, wall "
        f"{w['d_wall']:.0f} K, coolant film {w['d_cool']:.0f} K. The biggest temperature drop is in the "
        f"{w['dominant']} ({w['shares'][w['dominant']] * 100:.0f}%), which is what sets how hot {what} runs.",
    ]
    if w["dominant"] == "coolant film":
        regime = "laminar/transitional" if w["Re"] < 4000 else "turbulent"
        lines.append(
            f"The coolant is only carrying heat away at h = {w['h_c']:.0f} W/m^2-K (Re = {w['Re']:.0f}, {regime}, "
            f"{w['v']:.1f} m/s), so the wall surface has to get hot to push the heat into it. "
            "Fix: speed the coolant up here (narrower or shallower channel, or fewer channels), or add surface area."
        )
    elif w["dominant"] == "wall":
        lines.append(
            "The wall itself is the bottleneck: heat is struggling to conduct through it. "
            "Fix: thinner liner, a higher-conductivity material, or a lower heat load (e.g. film cooling)."
        )
    else:
        lines.append(
            "Most of the temperature drop is in the hot-gas boundary layer, i.e. the gas simply delivers a lot of "
            "heat here. Fix: lower the heat flux (lower Pc, different mixture ratio, a TBC) and make sure the "
            "coolant is not also the weak link."
        )
    return lines


def explain_regen(regen, T_hw_limit=None, T_cw_limit=None, dp_limit_psid=None):
    """
    Plain-language description of each crossed regen redline. Recomputes the crossings from the solved
    circuit, so it works whether or not plots were drawn.
    """
    out = []
    z = regen.circuit_z
    melt = getattr(regen.material, "MeltingPoint", None)

    def worst(mask, arr):
        idx = np.flatnonzero(mask)
        return int(idx[np.argmax(np.asarray(arr)[idx])])

    if melt is not None:
        m = regen.T_hw_arr > melt
        if m.any():
            k = worst(m, regen.T_hw_arr)
            out.append((
                f"Hot wall above {regen.material.material_name} melting point ({melt:.0f} K) at {_fmt_span(regen, m)}",
                [f"Peak hot wall {regen.T_hw_arr[k]:.0f} K. The metal would soften and melt through."]
                + _wall_explanation(regen, k, "the hot wall"),
            ))
    if T_hw_limit is not None:
        m = regen.T_hw_arr > T_hw_limit
        if m.any():
            k = worst(m, regen.T_hw_arr)
            out.append((
                f"Hot wall above your {T_hw_limit:.0f} K design limit at {_fmt_span(regen, m)}",
                [f"Peak hot wall {regen.T_hw_arr[k]:.0f} K ({regen.T_hw_arr[k] - T_hw_limit:.0f} K over)."]
                + _wall_explanation(regen, k, "the hot wall"),
            ))
    if T_cw_limit is not None:
        m = regen.T_cw_arr > T_cw_limit
        if m.any():
            k = worst(m, regen.T_cw_arr)
            out.append((
                f"Coolant-side wall above {T_cw_limit:.0f} K at {_fmt_span(regen, m)}",
                [f"Coolant-side wall reaches {regen.T_cw_arr[k]:.0f} K while the bulk coolant is only "
                 f"{regen.T_c_arr[k]:.0f} K. Fluid touching the wall can overheat, boil locally or coke even "
                 "though the bulk looks fine."]
                + _wall_explanation(regen, k, "the coolant-side wall"),
            ))

    m = regen.T_c_arr > regen.coolant_boiling_arr
    if m.any():
        k = worst(m, regen.T_c_arr - regen.coolant_boiling_arr)
        out.append((
            f"Bulk coolant above its boiling limit at {_fmt_span(regen, m)}",
            [f"By z = {z[k]:.2f} in the coolant has soaked up enough heat to reach {regen.T_c_arr[k]:.0f} K, above "
             f"the {regen.coolant_boiling_arr[k]:.0f} K boiling limit at its local pressure.",
             "Boiling forms vapor that insulates the wall and can starve the channel (flow instability), so the wall "
             "temperature can jump. Fix: more coolant flow, higher coolant pressure (raises the boiling point), or a "
             "lower total heat load."],
        ))

    P_c = regen.P_c_arr
    P_g = regen.engine.P[regen.circuit_outlet:regen.circuit_inlet + 1] * 6894.76
    m = P_c < P_g
    if m.any():
        k = worst(m, P_g - P_c)
        out.append((
            f"Coolant pressure below gas pressure at {_fmt_span(regen, m)}",
            [f"At z = {z[k]:.2f} in the coolant ({P_c[k] / 6894.76:.0f} psia) is lower than the chamber gas "
             f"({P_g[k] / 6894.76:.0f} psia), so the thin liner is pushed outward into the channels instead of being "
             "held inward, and hot gas could leak if it cracks. Fix: raise inlet pressure or cut the channel "
             "pressure drop."],
        ))

    if dp_limit_psid is not None:
        total = float(np.sum(regen.DP_arr)) / 6894.76
        if total > dp_limit_psid:
            k = int(np.argmax(regen.DP_arr))
            out.append((
                f"Total channel pressure drop {total:.1f} psid exceeds {dp_limit_psid:g} psid",
                [f"Friction in the channels eats {total:.1f} psid. The worst single section is at z = {z[k]:.2f} in "
                 f"({_section_name(regen, k)}), where the coolant moves at {regen.coolant_vel_arr[k]:.1f} m/s through "
                 f"a {regen.C_w_arr[k] * 1000:.2f} x {regen.C_h_arr[k] * 1000:.2f} mm channel.",
                 "Small channels and fast flow cost pressure. Fix: widen/deepen the channel where velocity peaks, add "
                 "channels (less flow each), or raise the supply pressure budget."],
            ))
    return out


def explain_structural(sa):
    """Plain-language description of each crossed stress/life redline of a LinerStructuralAnalysis."""
    from structural import LOCATIONS, LOCATION_LABELS, PA_TO_KSI
    out = []
    masks = sa._redline_masks()

    def worst_station(mask, score):
        idx = np.flatnonzero(mask)
        return int(idx[np.argmax(np.asarray(score)[idx])])

    ratio = np.max(np.vstack([sa.yield_ratio[loc] for loc in LOCATIONS]), axis=0)
    label = "Elastic stress above yield (plastic cycling)"
    if masks[label].any():
        k = worst_station(masks[label], ratio)
        loc = max(LOCATIONS, key=lambda name: sa.yield_ratio[name][k])
        th, hoop, bend = abs(sa.sigma_thermal_hot[k]), abs(sa.sigma_hoop_membrane[k]), abs(sa.sigma_bend_rib[k])
        dom = max({"thermal": th, "hoop": hoop, "bending": bend}.items(), key=lambda kv: kv[1])[0]
        expl = {
            "thermal": f"The gas side of the liner is {sa.T_hw[k]:.0f} K but the coolant side is {sa.T_cw[k]:.0f} K. "
                       "The hot face wants to expand but the cooler metal behind it holds it back, so it is squeezed "
                       "in compression.",
            "hoop": f"Coolant pressure is {sa.dP[k] / 6894.76:.0f} psi above the gas pressure, squeezing the thin liner "
                    "inward like a can being crushed around its circumference.",
            "bending": "Coolant pressure pushes the liner span between ribs like a clamped beam; a wide channel or "
                       "thin wall bows it and stresses the rib roots.",
        }[dom]
        out.append((
            f"Elastic stress above yield at {span_text(sa.z, masks[label])}",
            [f"Worst at z = {sa.z[k]:.2f} in ({LOCATION_LABELS[loc]}): {ratio[k]:.2f}x the yield strength "
             f"({sa.sigma_yield[loc][k] * PA_TO_KSI:.1f} ksi at that temperature). Main driver: {dom} "
             f"(thermal {th * PA_TO_KSI:.1f}, hoop {hoop * PA_TO_KSI:.1f}, bending {bend * PA_TO_KSI:.1f} ksi).",
             expl,
             "Yielding means the metal is deformed permanently every firing. It is not an instant failure, but it "
             "uses up fatigue life quickly. Fix: cut the wall temperature gradient (thinner wall, better cooling) or "
             "the pressure load (narrower channel / thicker wall)."],
        ))

    label = "Cycles to failure < 1 (fails on first cycle)"
    if masks[label].any():
        k = worst_station(masks[label], -sa.cycles_min)
        out.append((
            f"Liner fails on the first cycle at {span_text(sa.z, masks[label])}",
            [f"At z = {sa.z[k]:.2f} in the strain from one start-up is beyond what the material can take even once. "
             "Expect cracking or rupture on the first hot-fire."],
        ))

    req = sa.required_cycles
    if req is not None:
        label = f"Allowable cycles below required ({req:g})"
        if masks[label].any():
            k = worst_station(masks[label], -sa.cycles_min)
            loc = sa.governing_loc[k]
            out.append((
                f"Fatigue life shorter than the required {req:g} cycles at {span_text(sa.z, masks[label])}",
                [f"At z = {sa.z[k]:.2f} in the liner ({LOCATION_LABELS[loc]}) goes from relaxed to loaded and back "
                 f"each firing, with a strain swing of +/-{sa.eps_amp[loc][k] * 100:.3f}%. Repeating that wears the "
                 f"metal out after about {sa._fmt_cycles(sa.cycles_min[k])} cycles "
                 f"({sa._fmt_cycles(sa.allowable_cycles[k])} after the {sa.life_safety_factor:g}x safety factor).",
                 "Life is very sensitive to strain amplitude: lowering the wall temperature gradient or the pressure "
                 "bending a little buys a lot of cycles. Check the material fatigue data too (placeholders give "
                 "unreliable life)."],
            ))

    if sa.sigma_jacket_hoop is not None:
        label = "Jacket hoop stress above yield"
        if masks[label].any():
            k = worst_station(masks[label], sa.sigma_jacket_hoop / sa.jacket_yield)
            out.append((
                f"Outer jacket over-stressed at {span_text(sa.z, masks[label])}",
                [f"At z = {sa.z[k]:.2f} in the coolant pressure ({sa.P_c[k] / 6894.76:.0f} psia) inflates the outer "
                 f"jacket like a pressure vessel: {sa.sigma_jacket_hoop[k] * PA_TO_KSI:.1f} ksi against "
                 f"{sa.jacket_yield[k] * PA_TO_KSI:.1f} ksi yield. Fix: thicker jacket or lower coolant pressure."],
            ))
    return out


def print_explanations(items, header="WHAT IS HAPPENING AT EACH REDLINE"):
    """Pretty-print the output of explain_regen / explain_structural."""
    if not items:
        return
    print(f"\n---- {header} ----")
    for title, lines in items:
        print(f"\n* {title}")
        for ln in lines:
            print(f"    {ln}")
