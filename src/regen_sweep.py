"""
Regen Channel Design Sweep
--------------------------
Sweep one or two independent variables from:
    - C_w [mm]
    - C_h [mm]
    - N   [count]
    - t_w [mm]
    - Pc  [psia]
    - MR  [O/F]

Set ``sweep_var_2 = None`` to sweep only ``sweep_var_1``; results are 1D and plots are
standard 2D line plots (metric vs swept variable). With two variables, 3D surface plots
are produced.

All non-swept variables are held at their nominal values.
Failed runs are recorded as zeros.
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
# Base engine settings
# ------------------------------------------------------------------
name = "Test Case"
fuName = "Ethanol"
oxName = "LOX"
thrust = 2500
Pe = 8
CR = 4
Lstar = 35
cstarEff = 1
pcentFFC = 0.0
numPts = 100
bell = True
frozen = True
heat_flux_knockdown = 1

# Conical fallback (used if bell=False)
chmbR = 1.35
chmbL = 4.066
contAngle = 45
throatR = 0.603
throatL = 1 / 8
exit_angle = 15
exitR = 1.151

# ------------------------------------------------------------------
# Regen fixed settings
# ------------------------------------------------------------------
coolantName = "Ethanol"
inlet_T_c = 288
inlet_P_c = 410 * 6894.7
coolant_temp_limit = 540
material = Material("AlSi10Mg")

T_wg_bracket = (200.0, 3500.0)
T_wg_residual_tol = 1e-1
use_fixed_bartz_sigma = True
fixed_bartz_sigma = 1.0
sigma_clamp = (0.3, 1.5)
coolant_htc_correlation = "Dittus-Boelter"

# ------------------------------------------------------------------
# Sweep configuration
# ------------------------------------------------------------------
# Nominals for all supported sweep variables.
nominal = {
    "C_w": 1.0,   # mm
    "C_h": 2.0,   # mm
    "N": 49,      # count
    "t_w": 1.0,   # mm
    "Pc": 350.0,  # psia    
    "MR": 1.4,    # O/F
}

# Sweep ranges (used for whichever variable(s) are swept).
ranges = {
    "C_w": np.linspace(1.0, 2.5, 20),             # mm
    "C_h": np.linspace(1.0, 4.0, 20),             # mm
    "N": np.arange(40, 190, 5).astype(float),      # count
    "t_w": np.linspace(0.7, 1.6, 15),             # mm
    "Pc": np.linspace(220.0, 420.0, 17),          # psia
    "MR": np.linspace(0.8, 1.4, 17),              # O/F
}

# Select one or two independent variables.
# 1D sweep: set sweep_var_2 = None (only sweep_var_1 is varied).
# 2D sweep: set both; they must differ.
sweep_var_1 = "C_w"
sweep_var_2 = "N"


def _axis_label(var_name):
    labels = {
        "C_w": "Channel Width C_w [mm]",
        "C_h": "Channel Height C_h [mm]",
        "N": "Number of Channels N",
        "t_w": "Wall Thickness t_w [mm]",
        "Pc": "Chamber Pressure Pc [psia]",
        "MR": "Mixture Ratio MR [O/F]",
    }
    return labels[var_name]


def _build_case(var1_name, val1, var2_name=None, val2=None):
    case = dict(nominal)
    case[var1_name] = float(val1)
    if var2_name is not None:
        case[var2_name] = float(val2)
    return case


def _build_engine(pc_psia, mr):
    return Engine(
        thrust=thrust,
        Pc=float(pc_psia),
        Pe=Pe,
        MRcore=float(mr),
        oxName=oxName,
        fuName=fuName,
        CR=CR,
        Lstar=Lstar,
        cstarEff=cstarEff,
        pcentFFC=pcentFFC,
        frozen=frozen,
        name=name,
        numPts=numPts,
        verbose=False,
        bell=bell,
        chmbR=chmbR,
        chmbL=chmbL,
        contAngle=contAngle,
        throatR=throatR,
        throatL=throatL,
        expAngle=exit_angle,
        exitR=exitR,
    )


def run_single(case, engine_cache, fail_log=None):
    """
    Return tuple:
        (T_hw_max, T_c_out, DP_psid)
    Failed run -> (0, 0, 0) and the error is appended to *fail_log*.
    """
    try:
        engine_key = (float(case["Pc"]), float(case["MR"]))
        if engine_key not in engine_cache:
            engine_cache[engine_key] = _build_engine(engine_key[0], engine_key[1])
        engine = engine_cache[engine_key]

        rc = RegenCircuit(
            t_w=float(case["t_w"]) / 1000.0,
            material=material,
            tot_coolant_mdot=engine.fu_mdot,
            coolantName=coolantName,
            C_h=float(case["C_h"]) / 1000.0,
            engine=engine,
            heat_flux_knockdown=heat_flux_knockdown,
            N=int(round(case["N"])),
            C_w=float(case["C_w"]) / 1000.0,
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

        t_hw_max = float(np.max(rc.T_hw_arr))
        t_c_out = float(np.max(rc.T_c_arr))
        dp_psid = float(np.sum(rc.DP_arr)) / 6894.76
        return t_hw_max, t_c_out, dp_psid
    except Exception as exc:
        case_str = ", ".join(f"{k}={v}" for k, v in case.items())
        entry = f"FAIL [{case_str}]: {type(exc).__name__}: {exc}"
        if fail_log is not None:
            fail_log.append(entry)
        return 0.0, 0.0, 0.0


def main():
    supported = set(nominal.keys())
    one_d = sweep_var_2 is None

    if sweep_var_1 not in supported:
        raise ValueError(f"Supported sweep vars: {sorted(supported)}")
    if not one_d:
        if sweep_var_2 not in supported:
            raise ValueError(f"Supported sweep vars: {sorted(supported)}")
        if sweep_var_1 == sweep_var_2:
            raise ValueError("sweep_var_1 and sweep_var_2 must be different.")

    axis1 = ranges[sweep_var_1]
    label1 = _axis_label(sweep_var_1)
    if not one_d:
        axis2 = ranges[sweep_var_2]
        label2 = _axis_label(sweep_var_2)

    swept = (sweep_var_1,) if one_d else (sweep_var_1, sweep_var_2)
    fixed_vars = [v for v in ("C_w", "C_h", "N", "t_w", "Pc", "MR") if v not in swept]

    print("=" * 64)
    print("  NESS — Regen Design Sweep")
    print("=" * 64)
    if one_d:
        print(f"  Sweep axis : {sweep_var_1} (1D)")
    else:
        print(f"  Sweep axes : {sweep_var_1} x {sweep_var_2}")
    print("  Fixed vars : " + ", ".join([f"{v}={nominal[v]}" for v in fixed_vars]))
    print()

    n1 = len(axis1)
    engine_cache = {}

    if one_d:
        t_hw_grid = np.zeros(n1)
        t_c_grid = np.zeros(n1)
        dp_grid = np.zeros(n1)
        combos = [(i,) for i in range(n1)]
        n_cases = n1
    else:
        n2 = len(axis2)
        t_hw_grid = np.zeros((n1, n2))
        t_c_grid = np.zeros((n1, n2))
        dp_grid = np.zeros((n1, n2))
        combos = list(itertools.product(range(n1), range(n2)))
        n_cases = n1 * n2

    fail_log = []
    print(f"Running {n_cases} cases...")
    t_start = time.time()
    for idx in tqdm(combos, desc="Sweep"):
        if one_d:
            i1 = idx[0]
            case = _build_case(sweep_var_1, axis1[i1])
            t_hw, t_c, dp = run_single(case, engine_cache, fail_log)
            t_hw_grid[i1] = t_hw
            t_c_grid[i1] = t_c
            dp_grid[i1] = dp
        else:
            i1, i2 = idx
            case = _build_case(sweep_var_1, axis1[i1], sweep_var_2, axis2[i2])
            t_hw, t_c, dp = run_single(case, engine_cache, fail_log)
            t_hw_grid[i1, i2] = t_hw
            t_c_grid[i1, i2] = t_c
            dp_grid[i1, i2] = dp

    elapsed = time.time() - t_start
    n_fail = len(fail_log)
    print(f"Sweep complete in {elapsed:.1f} s  ({n_cases - n_fail}/{n_cases} converged, {n_fail} failed)")
    if fail_log:
        print(f"\n{'='*64}")
        print(f"  FAILED CASES ({n_fail})")
        print(f"{'='*64}")
        for entry in fail_log:
            print(f"  {entry}")
        print()

    csv_path = "regen_sweep_results.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("C_w_mm,C_h_mm,N,t_w_mm,Pc_psia,MR,T_hw_max_K,T_c_out_K,DP_psid\n")
        for idx in combos:
            if one_d:
                i1 = idx[0]
                case = _build_case(sweep_var_1, axis1[i1])
                th, tc, dp = t_hw_grid[i1], t_c_grid[i1], dp_grid[i1]
            else:
                i1, i2 = idx
                case = _build_case(sweep_var_1, axis1[i1], sweep_var_2, axis2[i2])
                th, tc, dp = t_hw_grid[i1, i2], t_c_grid[i1, i2], dp_grid[i1, i2]
            f.write(
                f"{case['C_w']:.6f},{case['C_h']:.6f},{int(round(case['N']))},"
                f"{case['t_w']:.6f},{case['Pc']:.6f},{case['MR']:.6f},"
                f"{th:.6f},{tc:.6f},{dp:.6f}\n"
            )
    print(f"Results saved to {csv_path}")

    fixed_text = ", ".join([f"{v}={nominal[v]}" for v in fixed_vars])

    if one_d:
        fig, axs = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
        plot_rows = [
            ("Max Hot Wall Temp T_hw [K]", t_hw_grid, "T_hw [K]"),
            ("Coolant Outlet Temp T_c [K]", t_c_grid, "T_c [K]"),
            ("Channel Pressure Drop DP [psid]", dp_grid, "DP [psid]"),
        ]
        for ax, (title, ydata, ylab) in zip(axs, plot_rows):
            ax.plot(axis1, ydata, "o-", markersize=4)
            ax.set_ylabel(ylab)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
        axs[-1].set_xlabel(label1)
        fig.suptitle(f"Regen design sweep — {label1}\n(Fixed: {fixed_text})", fontsize=11, y=1.02)
        plt.tight_layout()
    else:
        x, y = np.meshgrid(axis2, axis1)
        plot_data = [
            ("Max Hot Wall Temp T_hw [K]", t_hw_grid, cm.hot),
            ("Coolant Outlet Temp T_c [K]", t_c_grid, cm.coolwarm),
            ("Channel Pressure Drop DP [psid]", dp_grid, cm.viridis),
        ]
        for title, z, cmap in plot_data:
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111, projection="3d")
            surf = ax.plot_surface(x, y, z, cmap=cmap, edgecolor="k", linewidth=0.25, alpha=0.88)
            ax.set_xlabel(label2)
            ax.set_ylabel(label1)
            ax.set_zlabel(title)
            ax.set_title(f"{title}\n(Fixed: {fixed_text})")
            fig.colorbar(surf, shrink=0.6, pad=0.08)
            plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()