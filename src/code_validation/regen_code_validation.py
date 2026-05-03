from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# --- USER OPTIONS ---
plot_all_outputs = True
selected_test_case = "Test Case 3"  # "Test Case 1", "Test Case 2", "Test Case 3"
save_plot = False
save_plot_filename = "regen_validation_comparison.png"


BASE_DIR = Path(__file__).resolve().parent
NESS_RESULTS_PATH = BASE_DIR / "NESS Test Case Results.xlsx"
RPA_RESULTS_PATH = BASE_DIR / "RPA Test Case Results.xlsx"
TEST_CASES_PATH = BASE_DIR / "Regen Test Cases.xlsx"


def load_case_sheet(excel_path, sheet_name):
    xls = pd.ExcelFile(excel_path)
    if sheet_name not in xls.sheet_names:
        raise ValueError(
            f"Sheet '{sheet_name}' was not found in '{excel_path.name}'. "
            f"Available sheets: {xls.sheet_names}"
        )
    return pd.read_excel(excel_path, sheet_name=sheet_name)


def get_unit_normalized_series(ness_df, rpa_df):
    """
    Build normalized comparison series with common units:
    - axial: inches
    - temperature: K
    - velocity: m/s
    - pressure: psia
    - heat flux: W/m^2
    - HTC: W/m^2-K
    """
    mpa_to_psia = 145.0377377
    rpa_location_mm = rpa_df["Location (mm)"].to_numpy(dtype=float)
    rpa_radius_mm = rpa_df[" Radius, mm"].to_numpy(dtype=float)
    rpa_throat_idx = int(np.argmin(rpa_radius_mm))
    rpa_throat_location_mm = float(rpa_location_mm[rpa_throat_idx])
    rpa_axial_centered_in = (rpa_location_mm - rpa_throat_location_mm) / 25.4

    series_map = {
        "Axial Location": {
            "ness": ness_df["Axial Position (in)"].to_numpy(dtype=float),
            "rpa": rpa_axial_centered_in,
            "unit": "in",
        },
        "Hot Gas Heat Transfer Coefficient": {
            "ness": ness_df["Hot Gas HTC (W/m^2-K)"].to_numpy(dtype=float),
            "rpa": (rpa_df["Hot Gas HTC (kW/m^2-K)"].to_numpy(dtype=float) * 1000.0),
            "unit": "W/m^2-K",
        },
        "Heat Flux": {
            "ness": ness_df["Heat Flux (W/m^2)"].to_numpy(dtype=float),
            "rpa": (rpa_df["Heat Flux (kW/m^2)"].to_numpy(dtype=float) * 1000.0),
            "unit": "W/m^2",
        },
        "Twg (Hot Wall Temp)": {
            "ness": ness_df["Hot Wall Temp (K)"].to_numpy(dtype=float),
            "rpa": rpa_df["     Twg, K"].to_numpy(dtype=float),
            "unit": "K",
        },
        "Twc (Cold Wall Temp)": {
            "ness": ness_df["Cold Wall Temp (K)"].to_numpy(dtype=float),
            "rpa": rpa_df["     Twc, K"].to_numpy(dtype=float),
            "unit": "K",
        },
        "Tc (Coolant Temp)": {
            "ness": ness_df["Coolant Bulk Temp (K)"].to_numpy(dtype=float),
            "rpa": rpa_df["      Tc, K"].to_numpy(dtype=float),
            "unit": "K",
        },
        "Pc (Coolant Pressure)": {
            "ness": ness_df["Coolant Pressure (psia)"].to_numpy(dtype=float),
            "rpa": (rpa_df["    pc, MPa"].to_numpy(dtype=float) * mpa_to_psia),
            "unit": "psia",
        },
        "wc (Coolant Velocity)": {
            "ness": ness_df["Coolant Velocity (m/s)"].to_numpy(dtype=float),
            "rpa": rpa_df["    wc, m/s"].to_numpy(dtype=float),
            "unit": "m/s",
        },
        "rho (Coolant Density)": {
            "ness": ness_df["Coolant Density (kg/m^3)"].to_numpy(dtype=float),
            "rpa": rpa_df[" rho, kg/m^3"].to_numpy(dtype=float),
            "unit": "kg/m^3",
        },
    }
    return series_map


def plot_all_case_outputs(test_case_name):
    ness_df = load_case_sheet(NESS_RESULTS_PATH, test_case_name)
    rpa_df = load_case_sheet(RPA_RESULTS_PATH, test_case_name)

    series_map = get_unit_normalized_series(ness_df, rpa_df)
    x_ness = series_map["Axial Location"]["ness"]
    x_rpa = series_map["Axial Location"]["rpa"]

    y_keys = [k for k in series_map.keys() if k != "Axial Location"]
    fig, axes = plt.subplots(3, 3, figsize=(18, 12), sharex=False)
    axes = axes.ravel()

    for i, key in enumerate(y_keys):
        ax = axes[i]
        y_ness = series_map[key]["ness"]
        y_rpa = series_map[key]["rpa"]
        unit = series_map[key]["unit"]

        ax.plot(x_ness, y_ness, linewidth=2.0, label="NESS")
        ax.plot(x_rpa, y_rpa, "--", linewidth=2.0, label="RPA")
        ax.set_title(key)
        ax.set_xlabel("Axial Location [in]")
        ax.set_ylabel(f"{key} [{unit}]")
        ax.grid(True, alpha=0.3)
        ax.legend()

    for j in range(len(y_keys), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Regen Code Validation: NESS vs RPA ({test_case_name})", fontsize=15)
    fig.tight_layout()

    if save_plot:
        out_path = BASE_DIR / save_plot_filename
        fig.savefig(out_path, dpi=200, bbox_inches="tight")

    plt.show()


def print_selected_test_case_table(test_case_name):
    if not TEST_CASES_PATH.exists():
        return
    test_cases_df = pd.read_excel(TEST_CASES_PATH)
    if test_case_name not in test_cases_df.columns:
        raise ValueError(
            f"'{test_case_name}' not found in '{TEST_CASES_PATH.name}'. "
            f"Available: {list(test_cases_df.columns[1:])}"
        )
    table_df = test_cases_df[["Parameter", test_case_name]].copy()
    table_df.columns = ["Input Parameter", "Value"]
    table_df = table_df.dropna(how="any")
    print()
    print(f"{test_case_name} Inputs")
    print(table_df.to_string(index=False))


if __name__ == "__main__":
    print_selected_test_case_table(selected_test_case)
    if plot_all_outputs:
        plot_all_case_outputs(selected_test_case)
