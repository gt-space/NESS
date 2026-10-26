import os
from pathlib import Path


from materials import Material
from regen_circuit import RegenCircuit
import numpy as np
from engine import Engine
from utils import cea_chamber_tcomb_mr_sweep
import time
import matplotlib.pyplot as plt

### --- RUNTIME TRACKER --- ###
start_time = time.time()

### --- DESIGN & PLOT OPTIONS --- ###
design_engine = True                 # Design engine contour and performance calculations
design_regen = True                  # Solve regen circuit and performance analysis
design_regen = True                  # Solve regen circuit and performance analysis

# REGEN OUTPUT PLOT OPTIONS
display_regen_contour_plot = True   # Output Plot of Regen Channels
display_regen_outputs = True        # Show plots of regen circuit
show_regen_temps = True             # Regen Circuit Temps plot
show_cold_temps = False              # Cold Temps plot
show_dp = True                       # DP plot
show_pressures = True               # Pressures plot
show_qdot = False                    # Qdot plot
show_re = False                      # Re plot
show_pr = False                     # Pr plot
show_coolant_density = False         # Coolant density plot
show_coolant_velocity = False        # Coolant velocity plot
show_wall_thermal_gradient = False   # Wall thermal gradient plot
show_tangential_stresses = False     # Tangential stresses plot
show_longitudinal_stresses = False   # Longitudinal stresses plot
show_htc = True                     # HTC plot
show_coolant_htc_both = True       # If True, overlay Gnielinski and Dittus-Boelter coolant HTCs
show_bartz_sigma = False             # Bartz sigma subplot
show_heat_flux = False               # Heat flux (MW/m^2) plot

# BARTZ T_wg ROOT SOLVER OPTIONS
T_wg_bracket = (50.0, 3500.0)      # [K] global bracket for root solve
T_wg_residual_tol = 1e-1            # [K] residual tolerance for T_wg - T_hw
use_fixed_bartz_sigma = False       # If True, override Bartz sigma with fixed value
fixed_bartz_sigma = 1.0             # Used only when use_fixed_bartz_sigma=True
sigma_clamp = (0.5, 1.5)                  # If set to (min_sigma, max_sigma), clamps sigma in Bartz station solve
coolant_htc_correlation = "Gnielinski"    # "Gnielinski" or "Dittus-Boelter"
use_recovery_temp = False          # If True, drive heat flux with the adiabatic wall (recovery) temp instead of static gas temp

# STRESS & LIFECYCLE (structural.py)
run_structural_analysis = True      # Thermal + hoop/pressure stress, Neuber-corrected LCF life
stress_concentration_Kt = 1.0       # Local Kt at the liner/rib fillet (1.0 = none)
life_safety_factor = 4.0            # Allowable cycles = predicted cycles / this
required_cycles = 10                 # Hotfires the hardware must survive (pass/fail + the 'will it survive' verdict). Nothing is reported past 20.
t_jacket = None                     # [m] outer jacket thickness for the jacket hoop check (None to skip)
show_lifecycle_plots = True        # Yielding / fatigue build-up vs number of fires
show_ratcheting_plots = True       # Thermal ratcheting map + throat size vs fires
export_structural_csv = True       # Plain CSVs: key metrics, yielding vs fires, throat vs fires, per-station
show_structural_plots = True       # Stress components / von Mises / life vs position
export_structural_excel = True     # Per-station stress & life table
structural_excel_filename = "Liner Stress Life.xlsx"


# OFF-DESIGN FLOW ANALYSIS (mdot_analysis.py) - NOT a redesign: the geometry, channels, wall and material stay fixed
# and you see how this thrust chamber copes with different propellant flows. Inlet coolant T and P stay as set above.
use_injector_stiffness = True        # Regen inlet pressure follows throttle through the injector stiffness (replaces the fixed inlet_P_c)
injector_stiffness_points = ((0.40, 0.15), (1.00, 0.30))   # (throttle fraction, injector dP / Pc); linear between and beyond.
                                    # Inlet pressure = Pc*(1 + stiffness) + regen channel dP. 100% = the flows below (your actual flows).
primary_run_at_actual_flows = True  # True: the MAIN run (scorecard, "will it survive", redlines, all plots and exports) uses the
                                    # actual flows below instead of the design flows. Needs at least one actual flow to be set.
mdot_units = "kg/s"                 # units of the flows below: "kg/s" or "lbm/s"
run_mdot_analysis = True           # Evaluate ONE operating point at the flows below
analysis_ox_mdot = 0.4664             # oxidizer flow;  None = the design value
analysis_fu_mdot = 0.3731             # fuel flow;      None = the design value
analysis_coolant_mdot = None        # coolant flow;   None = all of the fuel flow
show_mdot_analysis_plots = True     # Regen + stress/life graphs for that operating point
run_throttle_sweep = True          # Sweep total flow (throttle) at fixed mixture ratio and graph the results
throttle_sweep_range = (0.4, 1.2)   # fractions of the reference flow (reference = the flows above, or the design flows)
n_throttle_points = 13
throttle_sweep_MR_values = None     # None = hold the reference MR; or a list like [1.6, 2.0, 2.4] to overlay several MRs
export_throttle_csv = False         # Save the sweep to throttle_sweep.csv
print_regen_station_rows = True   # If True, print per-station lines like "T_hw : … K" (tqdm still runs)

# ENGINE / PERFORMANCE OUTPUTS
display_nozzle_mesh = True         # Display mesh of the contour
show_nozzle_plot = True            # Show engine contour plot
show_bartz_plot = True             # Show Plot of Bartz HTC
show_gas_temp_plot = True          # Show Plot of Gas Temperatures along the contour
show_engine_perf_outputs = False    # Show Engine Performance Outputs
show_gas_props = False              # Generate plots of ALL CEA gas properties along the contour            
# Chamber Tcomb vs MR: RocketCEA-only sweep via utils (fast). MR_sweep_range None => prompt at run.
show_chamber_T_vs_MR = True
MR_sweep_range = (0.5, 4)         # (MR_min, MR_max); set None to prompt for min,max when plot enabled
num_MR_sweep_points = 25

# FILE EXPORT
export_nozzle = False                # Export inner (ID) nozzle contour as .txt to CAD
export_nozzle_outer = False          # Export outer (OD) contour using nozzle_wall_od_thickness_in
export_regen_chans = False           # Export Regen Circuit Channelas as a .txt to CAD
export_bartz_data = False            # Export Bartz HTC to a Excel file
export_gas_temps = False             # Export Gas Temps to an Excel file
export_pressures = False             # Export Pressures to an Excel file
export_regen_performance_excel = False   # Export axial regen performance outputs to Excel (.xlsx)
regen_performance_excel_filename = "Test Case 3.xlsx"

### --- ENGINE PERFORMANCE INPUTS --- ###
name = "Hopper Regen"
fuName = "RP-1"  # Jet-A
oxName = "LOX"
thrust = 550 # [lbf] 
Pc = 300 # psia
Pe = 18.8 # psia
MRcore = 1.25
CR = 5.994
Lstar = 30 # in
cstarEff = 0.9
pcentFFC = 0.0 #Film Cooling
numPts = 100 
bell = False
frozen = True
cstar_efficiency_in_flows = False  # True: flows = Pc*At/(cstarEff*c*), i.e. real flow at your combustion efficiency (changes the design coolant flow)

### --- Conical Engine Inputs --- ###
chmbR = 1.42 # in
chmbL = 4 # in, cylindrical portion of chamber
contAngle = 28 # degrees
chmbFilletR = 1.5 # in, fillet radius where the chamber meets the converging cone (0 = sharp corner)
throatFilletR = 0.9 # in, fillet radius where the converging cone meets the throat (0 = sharp corner); also the Bartz throat curvature radius
throatR = 0.58 # in; set to None to size the throat automatically for the thrust above (see the DESIGN CHECK printout)
throatL = 1/8 # in
exit_angle = 20 # degrees
exitR = .9925 # in

### --- MAIN CODE --- ###

if design_engine:
    # Create Engine
    engine = Engine(
        thrust = thrust,
        Pc = Pc,
        Pe = Pe,
        MRcore = MRcore,
        oxName=oxName,
        fuName=fuName,
        CR = CR,
        Lstar = Lstar,
        cstarEff=cstarEff,
        pcentFFC=pcentFFC,
        frozen=frozen,
        name=name,
        numPts=numPts,
        verbose=display_nozzle_mesh,
        bell=bell, 
        chmbR=chmbR, 
        chmbL=chmbL, 
        contAngle=contAngle, 
        throatR=throatR, 
        throatL=throatL, 
        expAngle=exit_angle, 
        exitR=exitR,
        chmbFilletR=chmbFilletR,
        throatFilletR=throatFilletR,
        cstarEffInFlow=cstar_efficiency_in_flows)
    engine.print_design_check()

    # Move the engine from the design point to your actual operating flows (geometry is unchanged).
    if primary_run_at_actual_flows and (analysis_ox_mdot is not None or analysis_fu_mdot is not None):
        from mdot_analysis import chamber_pressure_from_mdots, to_kg_per_s
        actual_ox = to_kg_per_s(analysis_ox_mdot, mdot_units) or engine.ox_mdot
        actual_fu = to_kg_per_s(analysis_fu_mdot, mdot_units) or engine.fu_mdot
        Pc_actual, MR_actual, _ = chamber_pressure_from_mdots(engine, actual_ox, actual_fu)
        engine.calcGasProperties(frozen=frozen, Pc_psia_override=Pc_actual, MR_override=MR_actual)
        print(f"\n==== MAIN RUN AT YOUR ACTUAL FLOWS (not the design flows) ====")
        print(f"  ox {actual_ox:.4f} kg/s + fuel {actual_fu:.4f} kg/s -> Pc {Pc_actual:.1f} psia (design {engine.design_state['Pc_psia']:g}), "
              f"MR {MR_actual:.3f} (design {engine.design_state['MR']:g}), ideal thrust {engine.cea_thrust_lbf(Pc_actual, MR_actual):.0f} lbf")

# Chamber combustion temperature vs mixture ratio (RocketCEA Tcomb only; no per-MR Engine)
if design_engine and show_chamber_T_vs_MR:
    _mr_range = MR_sweep_range
    if _mr_range is None:
        _raw = input("Enter MR sweep min and max (comma-separated, e.g. 1.0,2.5): ").strip()
        _parts = [float(x.strip()) for x in _raw.split(",")]
        if len(_parts) != 2:
            raise ValueError("Expected two values: MR_min, MR_max")
        _mr_range = (min(_parts), max(_parts))
    mr_lo, mr_hi = float(_mr_range[0]), float(_mr_range[1])
    mr_axis, tcomb_axis, tcomb_at_mrcore = cea_chamber_tcomb_mr_sweep(
        oxName,
        fuName,
        Pc,
        mr_lo,
        mr_hi,
        num_MR_sweep_points,
        mr_highlight=MRcore,
    )
    plt.figure()
    plt.plot(mr_axis, tcomb_axis, color="C0", label="Tcomb (CEA)")
    plt.scatter(
        [MRcore],
        [tcomb_at_mrcore],
        color="C3",
        s=80,
        zorder=5,
        edgecolors="black",
        linewidths=0.8,
        label=f"MRcore = {MRcore}",
    )
    plt.xlabel("Mixture ratio O/F")
    plt.ylabel("Chamber temperature Tcomb (K)")
    plt.title(f"Tcomb vs MR ({oxName} / {fuName}, Pc = {Pc} psia)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

# Gas Props
if show_gas_props:
    engine.plot_gas_props()

### --- REGEN CIRCUIT INPUTS --- ###
t_w = 1 / 1000 # [mm] --> [m]
nozzle_wall_od_thickness_in = 0.25  # [in] ID→OD structural thickness for export_nozzle_outer (CAD outer contour)
N = 50                               # number of channels

# Channel geometry is defined at four stations along the circuit (coolant flows inlet -> outlet).
# C_w = channel width, C_h = channel height, both in meters. "z" is the axial position in inches
# (same axis as engine.Contour_z); leave it out to use the contour start / end of chamber / throat / contour end.
# Dimensions are linearly interpolated in z between stations.
regen_stations = {
    "inlet":   {"C_w": 1.0 / 1000, "C_h": 2.0 / 1000},   # nozzle end, where coolant enters  [mm] --> [m]
    "throat":  {"C_w": 1.0 / 1000, "C_h": 1.0 / 1000},   # at the throat
    "chamber": {"C_w": 1.0 / 1000, "C_h": 2.0 / 1000},   # along the chamber
    "outlet":  {"C_w": 1.0 / 1000, "C_h": 2.0 / 1000},   # injector end, where coolant leaves
}
min_fin_width = 0.5 / 1000      # [m] smallest rib allowed anywhere on the circuit
heat_flux_knockdown = 1
coolantName = "n-Dodecane"
tot_coolant_mdot = engine.fu_mdot # kg/s
if primary_run_at_actual_flows and analysis_coolant_mdot is not None:
    from mdot_analysis import to_kg_per_s as _to_kg
    tot_coolant_mdot = _to_kg(analysis_coolant_mdot, mdot_units)   # your coolant flow instead of all of the fuel
inlet_T_c = 300 # [K]
inlet_P_c = 410 * 6894.7 # [psia] --> [Pa]
coolant_temp_limit = None # [K]; set a float to override fluid saturation/pseudocritical limit, set as JetA
tbc_thickness = 0.0 / 1000  # [mm] → [m]; set 0 for no coating
tbc_k = 1.0                 # [W/m-K] thermal conductivity of TBC (only used when tbc_thickness > 0)
# circuit_inlet = -100 # Not used for now
material = Material("AlSi10Mg") # Material for the regen circuit (AlSi10Mg, Inconel 718, etc.)

#print(f"Fuel Mdot : {round(tot_coolant_mdot, 3)} kg/s")
#print(f"OX Mdot : {round(engine.ox_mdot, 3)} kg/s")
#print(f"Throat Diameter : {engine.Dt * 39.3701} in")

### --- REGEN CIRCUIT CONSTRAINTS --- ###
T_hw_max = 1050 # K
channel_dp_max = 100 # [psid] --> [Pa]

# Create Regen Circuit
if design_regen:

    regen_circuit = RegenCircuit(
        t_w=t_w,
        material=material,
        tot_coolant_mdot=tot_coolant_mdot,
        coolantName=coolantName,
        engine=engine,
        heat_flux_knockdown=heat_flux_knockdown,
        N=N,
        stations=regen_stations,
        min_fin_width=min_fin_width,
        tbc_thickness=tbc_thickness,
        tbc_k=tbc_k,
    )

    if display_regen_contour_plot:
        regen_circuit.plot_regen_geometry(engine)

    # Regen inlet pressure: fixed input, or set by the injector stiffness at 100% throttle (= the main run's flows)
    inlet_P_c_used = inlet_P_c
    if use_injector_stiffness:
        from mdot_analysis import injector_stiffness, inlet_pressure_for_injector
        stiff_main = injector_stiffness(1.0, injector_stiffness_points)
        feed = inlet_pressure_for_injector(
            regen_circuit, engine.Pc_psia, stiff_main,
            dict(inlet_T_c=inlet_T_c, inlet_pressure=inlet_P_c, T_wg_bracket=T_wg_bracket, T_wg_residual_tol=T_wg_residual_tol,
                 fixed_sigma=(fixed_bartz_sigma if use_fixed_bartz_sigma else None), sigma_clamp=sigma_clamp,
                 coolant_htc_correlation=coolant_htc_correlation, coolant_temp_limit=coolant_temp_limit,
                 use_recovery_temp=use_recovery_temp),
            coolant_mdot=regen_circuit.tot_coolant_mdot,
        )
        inlet_P_c_used = feed["P_in_Pa"]
        print(f"\n==== REGEN INLET PRESSURE FROM INJECTOR STIFFNESS ====")
        print(f"  Pc {engine.Pc_psia:.1f} psia, stiffness {stiff_main * 100:.0f}% -> injector dP {feed['dP_inj_Pa'] / 6894.7:.1f} psid, "
              f"injector inlet {feed['P_inj_in_Pa'] / 6894.7:.1f} psia, channel dP {feed['dP_regen_Pa'] / 6894.7:.1f} psid "
              f"-> regen inlet {inlet_P_c_used / 6894.7:.1f} psia (the fixed input was {inlet_P_c / 6894.7:.0f} psia)")

    # Solve Regen Circuit
    regen_circuit.solve_circuit(
        inlet_T_c=inlet_T_c,
        inlet_pressure=inlet_P_c_used,
        T_wg_bracket=T_wg_bracket,
        T_wg_residual_tol=T_wg_residual_tol,
        fixed_sigma=(fixed_bartz_sigma if use_fixed_bartz_sigma else None),
        sigma_clamp=sigma_clamp,
        coolant_htc_correlation=coolant_htc_correlation,
        coolant_temp_limit=coolant_temp_limit,
        print_station_rows=print_regen_station_rows,
        use_recovery_temp=use_recovery_temp,
    )

    if run_structural_analysis:
        structural = regen_circuit.run_structural_analysis(
            Kt=stress_concentration_Kt,
            life_safety_factor=life_safety_factor,
            required_cycles=required_cycles,
            t_jacket=t_jacket,
        )
        structural.summary()
        if show_structural_plots:
            structural.plot()
        if show_lifecycle_plots:
            structural.plot_lifecycle()
        if show_ratcheting_plots:
            structural.plot_ratcheting()
        if export_structural_csv:
            structural.export_csv("Liner")
        if export_structural_excel:
            structural.export_to_excel(structural_excel_filename)

    if display_regen_outputs:
        regen_circuit.outputs(
            show_regen_temps=show_regen_temps,
            show_cold_temps=show_cold_temps,
            show_dp=show_dp,
            show_pressures=show_pressures,
            show_qdot=show_qdot,
            show_re=show_re,
            show_pr=show_pr,
            show_coolant_density=show_coolant_density,
            show_coolant_velocity=show_coolant_velocity,
            show_wall_thermal_gradient=show_wall_thermal_gradient,
            show_tangential_stresses=show_tangential_stresses,
            show_longitudinal_stresses=show_longitudinal_stresses,
            show_htc=show_htc,
            show_coolant_htc_both=show_coolant_htc_both,
            show_bartz_sigma=show_bartz_sigma,
            show_heat_flux=show_heat_flux,
            T_hw_limit=T_hw_max,          # redline: your hot wall design limit [K]
            dp_limit_psid=channel_dp_max, # redline: max total channel pressure drop [psid]
        ) # Plot regen circuit values (crossed failure limits are drawn in red)

    if export_regen_chans:
        regen_circuit.generate_single_channel_curves()

    if export_regen_performance_excel:
        regen_circuit.export_performance_to_excel(
            filename=regen_performance_excel_filename,
        )

if show_bartz_plot:
    print(regen_circuit.h_hg_arr)
    plt.figure()
    plt.title("Bartz Testing")
    plt.xlabel("Axial Position (in)")
    plt.ylabel("HTC (W/m^2-K)")
    plt.plot(engine.Contour_z, regen_circuit.h_hg_arr)
    plt.show()

if export_bartz_data:
    data = np.column_stack((engine.Contour_z, regen_circuit.h_hg_arr))
    np.savetxt("BartzHeatsinkMax.csv", data, delimiter=",", header="X Position (in), HTC (W/m^2-K)", comments="")

if show_gas_temp_plot:
    plt.figure()
    plt.title("Temps vs Position")
    plt.xlabel("Axial Position (in)")
    plt.ylabel("Temperature (K)")
    plt.plot(engine.Contour_z, engine.T)
    plt.show()

if export_gas_temps:
    import numpy as np
    data = np.column_stack((engine.Contour_z, engine.T))
    np.savetxt("TempsHeatsinkMax.csv", data, delimiter=",", header="X Position (in), Gas Temperature (K)", comments="")

### --- RUNTIME ----
end_time = time.time()
print(f"Runtime: {round(end_time - start_time, 2)} sec")

### --- PLOT OUTPUTS ---

# Plot Engine Contour
if show_nozzle_plot and design_regen:
    # Contour with the regen channel cutaway alongside it
    regen_circuit.plot_contour_with_channels(
        title=f'Hopper Engine Profile - {thrust} lbf' + (' (Bell)' if bell else ' (Conical)'),
        t_jacket=t_jacket,
    )
elif show_nozzle_plot:
    if bell:
        engine.R.geomObj.plot_geometry(title=f'Hopper Engine Profile - {thrust} lbf', show_grid=True)
    else:
        plt.figure(figsize=(10, 5))
        plt.plot(engine.Contour_z, engine.Contour_r, color='blue')
        plt.plot(engine.Contour_z, -np.array(engine.Contour_r), color='blue')
        plt.plot([engine.Contour_z[0], engine.Contour_z[-1]], [0, 0], 'k--', linewidth=1)
        plt.title(f'Hopper Engine Profile - {thrust} lbf (Conical)')
        plt.xlabel('Axial Distance (in)')
        plt.ylabel('Radial Distance (in)')
        plt.grid()
        plt.axis('equal')
        plt.show()

if export_nozzle:
    engine.exportGeometry(filename="Hopper Engine Contour 550 lbf 4_11_26")
if export_nozzle_outer:
    engine.exportGeometry(
        filename="Hopper Engine Contour 550 lbf 4_11_26 OD",
        wall_thickness_in=nozzle_wall_od_thickness_in,
    )


if export_pressures:
    import numpy as np
    data = np.column_stack((engine.Contour_z, engine.P))
    np.savetxt("PressuresHeatsinkMax.csv", data, delimiter=",", header="X Position (in), Pressure (Psia)", comments="")

# Show Engine Performance Outputs
if show_engine_perf_outputs:
    print(engine.R.get_summ_str())


### --- OFF-DESIGN FLOW ANALYSIS / THROTTLE SWEEP ---
if (run_mdot_analysis or run_throttle_sweep) and design_regen:
    from mdot_analysis import (OffDesignAnalyzer, to_kg_per_s, print_operating_point, plot_throttle_sweep,
                               print_sweep_table, export_sweep_csv)

    analyzer = OffDesignAnalyzer(
        engine, regen_circuit,
        solve_kwargs=dict(
            inlet_T_c=inlet_T_c,
            inlet_pressure=inlet_P_c_used,
            T_wg_bracket=T_wg_bracket,
            T_wg_residual_tol=T_wg_residual_tol,
            fixed_sigma=(fixed_bartz_sigma if use_fixed_bartz_sigma else None),
            sigma_clamp=sigma_clamp,
            coolant_htc_correlation=coolant_htc_correlation,
            coolant_temp_limit=coolant_temp_limit,
            use_recovery_temp=use_recovery_temp,
        ),
        structural_kwargs=dict(
            Kt=stress_concentration_Kt,
            life_safety_factor=life_safety_factor,
            required_cycles=required_cycles,
            t_jacket=t_jacket,
        ),
        T_hw_limit=T_hw_max,
        dp_limit_psid=channel_dp_max,
        injector_stiffness_points=(injector_stiffness_points if use_injector_stiffness else None),
    )
    ox_ref = to_kg_per_s(analysis_ox_mdot, mdot_units) or analyzer.primary_ox_mdot
    if analysis_ox_mdot is None and analysis_fu_mdot is None:
        print("NOTE: no actual operating flows were given (analysis_ox_mdot / analysis_fu_mdot), so the design flows are the 100% reference.")
    fu_ref = to_kg_per_s(analysis_fu_mdot, mdot_units) or analyzer.primary_fu_mdot

    if run_mdot_analysis:
        point = analyzer.operating_point(ox_ref, fu_ref, coolant_mdot=to_kg_per_s(analysis_coolant_mdot, mdot_units), throttle=1.0)
        print_operating_point(point)
        if show_mdot_analysis_plots and point["status"] != "error":
            regen_circuit.outputs(
                show_regen_temps=True, show_cold_temps=True, show_dp=True, show_pressures=True,
                show_qdot=False, show_re=True, show_pr=False, show_coolant_density=False, show_coolant_velocity=True,
                show_wall_thermal_gradient=False, show_tangential_stresses=False, show_longitudinal_stresses=False,
                show_htc=True, show_bartz_sigma=False, show_heat_flux=True, show_fin_efficiency=False,
                T_hw_limit=T_hw_max, dp_limit_psid=channel_dp_max,
            )
            structural_pt = regen_circuit.structural
            structural_pt.summary()
            structural_pt.plot()
            structural_pt.plot_lifecycle()
            structural_pt.plot_ratcheting()
        analyzer.restore()

    if run_throttle_sweep:
        sweep = analyzer.throttle_sweep(
            ox_ref, fu_ref,
            fractions=np.linspace(throttle_sweep_range[0], throttle_sweep_range[1], int(n_throttle_points)),
            MR_values=throttle_sweep_MR_values,
        )
        print_sweep_table(sweep)
        if export_throttle_csv:
            export_sweep_csv(sweep, "throttle_sweep.csv")
        plot_throttle_sweep(sweep)
