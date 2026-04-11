"""
Throttle Sweep — Off-Design Regen Thermal Analysis
---------------------------------------------------
Fixed thrust-chamber + nozzle geometry and fixed regen channel geometry.
Sweeps throttle percentage (fraction of design thrust) and optionally MR.

At each operating point:
    1.  Pc  = Pc_design * throttle_fraction
    2.  Gas properties are recomputed via CEA at the new Pc (and MR).
    3.  Total mass flow = Pc * At / cstar  (constant cstar assumption — TODO: update)
    4.  Fuel / ox split from MR, coolant mdot from user-selected propellant.
    5.  Regen circuit is solved → T_hw_max, T_c_out, DP, q_flux_max, etc.

Failed runs appear as 0 on all plots.

Usage:
    python src/throttle_sweep.py
"""

import itertools
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from tqdm import tqdm

from engine import Engine
from materials import Material
from regen_circuit import RegenCircuit

warnings.filterwarnings("ignore")

# ------------------------------------------------------------------
#  ENGINE INPUTS  (design point — 100 % thrust)
# ------------------------------------------------------------------
name       = "Hopper SN1"
fuName     = "JetA"
oxName     = "LOX"
thrust     = 500        # lbf (design)
Pc_design  = 300        # psia (design, = 100 % throttle)
Pe         = 18.8       # psia
Pamb       = 14.8       # psia
MRcore     = 1.25        # design O/F
CR         = 6
Lstar      = 30         # in
cstarEff   = 0.85
pcentFFC   = 0.0
numPts     = 100
bell       = True
frozen     = True

# ------------------------------------------------------------------
#  REGEN CIRCUIT INPUTS  (fixed hardware)
# ------------------------------------------------------------------
t_w                     = 1 / 1000        # [mm] → [m]
N                       = 49
C_w                     = 1 / 1000        # [mm] → [m]
C_h                     = 2 / 1000        # [mm] → [m]
coolantName             = "n-Decane"
inlet_T_c               = 288             # [K]
inlet_P_c               = 410 * 6894.7    # [psia] → [Pa]
coolant_temp_limit      = 540             # [K] or None
tbc_thickness           = 0.0 / 1000    # [mm] → [m]
tbc_k                   = 1.0             # [W/m-K]
material                = Material("AlSi10Mg")

T_wg_bracket            = (300.0, 2500.0)
T_wg_residual_tol       = 1e-1
use_fixed_bartz_sigma   = False
fixed_bartz_sigma       = 1.0
sigma_clamp             = (0.5, 1.5)
coolant_htc_correlation = "Dittus-Boelter"

# Coolant source: "fuel" or "ox"
coolant_source = "fuel"

# ------------------------------------------------------------------
#  THROTTLE SWEEP CONFIGURATION
# ------------------------------------------------------------------
throttle_range  = (0.4, 1.0)       # fraction of design thrust
n_throttle      = 15               # number of throttle points

# MR mode: set to a single float for fixed MR, or a list/array for MR sweep.
MR_sweep = MRcore                  # fixed MR (use np.linspace(0.8, 1.4, 13) for sweep)

# ------------------------------------------------------------------


def _run_single(engine, coolant_mdot):
    """Run one regen solve. Returns dict of outputs or zeros on failure."""
    try:
        rc = RegenCircuit(
            t_w=t_w,
            material=material,
            tot_coolant_mdot=coolant_mdot,
            coolantName=coolantName,
            C_h=C_h,
            engine=engine,
            N=N,
            C_w=C_w,
            tbc_thickness=tbc_thickness,
            tbc_k=tbc_k,
        )
        rc.solve_circuit(
            inlet_T_c=inlet_T_c,
            inlet_pressure=inlet_P_c,
            T_wg_bracket=T_wg_bracket,
            T_wg_residual_tol=T_wg_residual_tol,
            fixed_sigma=(fixed_bartz_sigma if use_fixed_bartz_sigma else None),
            sigma_clamp=sigma_clamp,
            coolant_htc_correlation=coolant_htc_correlation,
            coolant_temp_limit=coolant_temp_limit,
            show_progress=False,
            print_station_rows=False,
        )
        return {
            "T_hw_max": float(np.max(rc.T_hw_arr)),
            "T_c_out": float(np.max(rc.T_c_arr)),
            "DP_psid": float(np.sum(rc.DP_arr)) / 6894.76,
            "q_flux_max_MW": float(np.max(rc.q_flux_arr)) / 1e6,
            "boiling_margin": float(np.min(rc.coolant_boiling_arr - rc.T_c_arr)),
            "vel_max": float(np.max(rc.coolant_vel_arr)),
        }
    except Exception:
        return {k: 0.0 for k in ("T_hw_max", "T_c_out", "DP_psid", "q_flux_max_MW", "boiling_margin", "vel_max")}


def main():
    mr_is_swept = not np.isscalar(MR_sweep)
    mr_axis = np.atleast_1d(MR_sweep).astype(float)
    throttle_axis = np.linspace(float(throttle_range[0]), float(throttle_range[1]), int(n_throttle))

    print("=" * 64)
    print("  NESS — Throttle Sweep (Off-Design Regen Analysis)")
    print("=" * 64)
    print(f"  Throttle range : {throttle_range[0]*100:.0f} % – {throttle_range[1]*100:.0f} %  ({n_throttle} pts)")
    if mr_is_swept:
        print(f"  MR range       : {mr_axis[0]:.2f} – {mr_axis[-1]:.2f}  ({len(mr_axis)} pts)")
    else:
        print(f"  MR (fixed)     : {mr_axis[0]:.2f}")
    print(f"  Coolant source : {coolant_source}")
    print(f"  Design Pc      : {Pc_design} psia")
    print(f"  Design thrust  : {thrust} lbf")
    print()

    # Build engine once at design point — locks geometry.
    print("Building engine at design point...")
    engine = Engine(
        thrust=thrust, Pc=Pc_design, Pe=Pe, MRcore=MRcore,
        oxName=oxName, fuName=fuName, CR=CR, Lstar=Lstar,
        cstarEff=cstarEff, pcentFFC=pcentFFC, frozen=frozen,
        name=name, numPts=numPts, verbose=False, bell=bell,
        Pamb=Pamb,
    )
    # Store design-point cstar (held constant across throttle — TODO: update per Pc).
    cstar_design = engine.cstar  # [m/s]

    output_keys = ("T_hw_max", "T_c_out", "DP_psid", "q_flux_max_MW", "boiling_margin", "vel_max")
    n_t = len(throttle_axis)
    n_mr = len(mr_axis)

    grids = {k: np.zeros((n_t, n_mr)) for k in output_keys}
    Pc_grid = np.zeros((n_t, n_mr))
    Pe_grid = np.zeros((n_t, n_mr))
    thrust_grid = np.zeros((n_t, n_mr))

    combos = list(itertools.product(range(n_t), range(n_mr)))
    print(f"Running {n_t} x {n_mr} = {len(combos)} cases...\n")
    t_start = time.time()

    for i_t, i_mr in tqdm(combos, desc="Throttle sweep"):
        throttle_frac = throttle_axis[i_t]
        mr = mr_axis[i_mr]

        Pc_throttled = Pc_design * throttle_frac  # psia
        Pc_grid[i_t, i_mr] = Pc_throttled

        # Recompute gas properties at this Pc and MR.
        engine.calcGasProperties(frozen=frozen, Pc_psia_override=Pc_throttled, MR_override=mr)

        # Derive mass flows using constant cstar assumption.
        At_m2 = engine.At * 0.00064516  # in^2 → m^2
        Pc_Pa = Pc_throttled * 6894.76
        tot_mdot = Pc_Pa * At_m2 / cstar_design
        fu_mdot = tot_mdot / (1.0 + mr)
        ox_mdot = tot_mdot - fu_mdot

        coolant_mdot = fu_mdot if coolant_source == "fuel" else ox_mdot

        thrust_actual = throttle_frac * thrust
        thrust_grid[i_t, i_mr] = thrust_actual

        # Nozzle exit pressure: Pe = P at last station (from isentropic expansion at this Pc).
        Pe_local = float(engine.P[-1])  # psia (engine.P is in psia before Pc conversion... 
        # Actually engine.P is built during calcGasProperties using P0=Pc_psia, so it's in psia.
        Pe_grid[i_t, i_mr] = Pe_local

        result = _run_single(engine, coolant_mdot)
        for k in output_keys:
            grids[k][i_t, i_mr] = result[k]

    elapsed = time.time() - t_start
    print(f"\nSweep complete in {elapsed:.1f} s")

    # Restore engine to design point.
    engine.calcGasProperties(frozen=frozen, Pc_psia_override=Pc_design, MR_override=MRcore)

    # ── Save CSV ──
    csv_path = "throttle_sweep_results.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("throttle_pct,MR,Pc_psia,thrust_lbf,Pe_psia," + ",".join(output_keys) + "\n")
        for i_t, i_mr in combos:
            vals = ",".join(f"{grids[k][i_t, i_mr]:.6f}" for k in output_keys)
            f.write(
                f"{throttle_axis[i_t]*100:.2f},{mr_axis[i_mr]:.4f},"
                f"{Pc_grid[i_t, i_mr]:.2f},{thrust_grid[i_t, i_mr]:.2f},"
                f"{Pe_grid[i_t, i_mr]:.4f},{vals}\n"
            )
    print(f"Results saved to {csv_path}")

    # ── Plotting ──
    throttle_pct = throttle_axis * 100.0
    sep_limit = 0.4 * Pamb  # psia

    if mr_is_swept:
        _plot_3d(throttle_pct, mr_axis, grids, Pe_grid, sep_limit)
    else:
        _plot_2d(throttle_pct, grids, Pe_grid, sep_limit)

    plt.show()


def _plot_2d(throttle_pct, grids, Pe_grid, sep_limit):
    """Line plots vs throttle % for fixed-MR case."""
    plot_specs = [
        ("Max Hot Wall Temp", grids["T_hw_max"][:, 0], "T_hw [K]"),
        ("Coolant Outlet Temp", grids["T_c_out"][:, 0], "T_c [K]"),
        ("Channel Pressure Drop", grids["DP_psid"][:, 0], "DP [psid]"),
        ("Max Heat Flux", grids["q_flux_max_MW"][:, 0], "q [MW/m²]"),
        ("Boiling Margin", grids["boiling_margin"][:, 0], "Margin [K]"),
        ("Max Coolant Velocity", grids["vel_max"][:, 0], "Vel [m/s]"),
    ]

    ncols = 2
    nrows = int(np.ceil(len(plot_specs) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(13, 4 * nrows), sharex=True)
    axs = axs.ravel()

    for ax, (title, ydata, ylabel) in zip(axs, plot_specs):
        ax.plot(throttle_pct, ydata, "o-", markersize=4)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    for ax in axs[len(plot_specs):]:
        ax.set_visible(False)
    for ax in axs:
        ax.set_xlabel("Throttle [%]")
    fig.suptitle(f"Throttle Sweep — MR = {MR_sweep}", fontsize=14, y=1.01)
    fig.tight_layout()

    # Pe plot with separation limit
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    ax2.plot(throttle_pct, Pe_grid[:, 0], "o-", color="C0", markersize=4, label="Nozzle Exit Pressure")
    ax2.axhline(sep_limit, color="r", linestyle="--", linewidth=1.5, label=f"Separation limit (0.4 × Pamb = {sep_limit:.2f} psia)")
    ax2.axhline(Pamb, color="gray", linestyle=":", linewidth=1, label=f"Pamb = {Pamb} psia")
    ax2.set_xlabel("Throttle [%]")
    ax2.set_ylabel("Exit Pressure [psia]")
    ax2.set_title("Nozzle Exit Pressure vs Throttle")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()


def _plot_3d(throttle_pct, mr_axis, grids, Pe_grid, sep_limit):
    """3D surface plots for throttle x MR sweep."""
    X, Y = np.meshgrid(mr_axis, throttle_pct)

    surface_specs = [
        ("Max Hot Wall Temp T_hw [K]", grids["T_hw_max"], cm.hot),
        ("Coolant Outlet Temp T_c [K]", grids["T_c_out"], cm.coolwarm),
        ("Channel DP [psid]", grids["DP_psid"], cm.viridis),
        ("Max Heat Flux [MW/m²]", grids["q_flux_max_MW"], cm.inferno),
        ("Boiling Margin [K]", grids["boiling_margin"], cm.RdYlGn),
        ("Max Coolant Velocity [m/s]", grids["vel_max"], cm.plasma),
    ]

    for title, Z, cmap in surface_specs:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")
        surf = ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor="k", linewidth=0.2, alpha=0.88)
        ax.set_xlabel("Mixture Ratio MR")
        ax.set_ylabel("Throttle [%]")
        ax.set_zlabel(title)
        ax.set_title(title)
        fig.colorbar(surf, shrink=0.55, pad=0.1)
        plt.tight_layout()

    # Pe surface + separation limit plane
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Pe_grid, cmap=cm.Blues, edgecolor="k", linewidth=0.2, alpha=0.85)
    sep_plane = np.full_like(Pe_grid, sep_limit)
    ax.plot_surface(X, Y, sep_plane, color="red", alpha=0.25, label=f"0.4×Pamb = {sep_limit:.2f} psia")
    ax.set_xlabel("Mixture Ratio MR")
    ax.set_ylabel("Throttle [%]")
    ax.set_zlabel("Exit Pressure [psia]")
    ax.set_title(f"Nozzle Exit Pressure (red = 0.4×Pamb = {sep_limit:.2f} psia)")
    plt.tight_layout()


if __name__ == "__main__":
    main()
