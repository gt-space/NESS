import os
from pathlib import Path

print("CWD =", os.getcwd())
print("CSV will be at =", Path("Hopper_Bartz_HTC_550lbf_2_2_26.csv").resolve())

#import regen_circuit as ta
from materials import Material
from regen_circuit import RegenCircuit
import numpy as np
from engine import Engine
import time
import matplotlib.pyplot as plt

### --- RUNTIME TRACKER --- ###
start_time = time.time()

### --- DESIGN & PLOT OPTIONS --- ###
design_engine = True                 # Design engine contour and performance calculations
design_regen = True                  # Solve regen circuit and performance analysis
display_regen_contour_plot = True   # Output Plot of Regen Channels
display_regen_outputs = True        # Show plots of regen circuit
display_nozzle_mesh = False           # Display mesh of the contour
show_nozzle_plot = False              # Show engine contour plot
export_nozzle = False                # Export Nozzle as .txt to CAD    
export_regen_chans = False           # Export Regen Circuit Channelas as a .txt to CAD
show_bartz_plot = False               # Show Plot of Bartz HTC
export_bartz_data = False            # Export Bartz HTC to a Excel file
show_gas_temp_plot = False            # Show Plot of Gas Temperatures along the contour
export_gas_temps = False             # Export Gas Temps to an Excel file
show_engine_perf_outputs = False     # Show Engine Performance Outputs

### --- ENGINE PERFORMANCE INPUTS --- ###
name = "Hopper SN1"
fuName = "Kerosene"  # Jet-A
oxName = "LOX"
thrust = 500 # [lbf] --> [N]
Pc = 300 # psia
Pe = 14.8 # psia
MRcore = 1.8
CR = 5
Lstar = 40 # in
cstarEff = 0.85
numPts = 100

### --- Conical Engine Inputs --- ###
chmbR = 1.35 # in
chmbL = 4.066 # in, cylindrical portion of chamber
contAngle = 45 # degrees
throatR = 0.603 # in
throatL = 1/8 # in
exit_angle = 15 # degrees
exitR = 1.151 # in

# 40 in L* is what Vespula Heatsink uses
# L* = 43 in
# CR = 6.5

# L* = 30 in
# CR = 5 

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
        exitR=exitR)
              

# Gas Props Testing
#engine.plot_gas_props()

### --- REGEN CIRCUIT INPUTS --- ###
t_w = 1 / 1000 # [mm] --> [m]
N = 30
C_w = 1 / 1000 # [mm] --> [m]
C_h = 2.5 / 1000 # [mm] --> [m]
coolantName = "n-Dodecane"
tot_coolant_mdot = engine.fu_mdot # kg/s
inlet_T_c = 293 # [K]
inlet_P_c = 420 * 6894.7 # [psia] --> [Pa]
# circuit_inlet = -100 # Not used for now
material = Material("Pure Copper")

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
        C_h=C_h,
        engine=engine,
        N =N,
        C_w=C_w,
    )

    if display_regen_contour_plot:
        regen_circuit.plot_regen_geometry(engine)

    # Solve Regen Circuit
    regen_circuit.solve_circuit(
        inlet_T_c=inlet_T_c,
        inlet_pressure=inlet_P_c,
    )

    if display_regen_outputs:
        regen_circuit.outputs() # Plot regen circuit values

    if export_regen_chans:
        regen_circuit.generate_single_channel_curves()

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
    np.savetxt("HTC_550lbf_FIXED_ConicalV3.csv", data, delimiter=",", header="Y Position (in), HTC (W/m^2-K)", comments="")

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
    np.savetxt("Hopper_Gas_Temps_550lbf_ConicalV3.csv", data, delimiter=",", header="Y Position (in), Gas Temperature (K)", comments="")

### --- PLOT OUTPUTS ---

# Plot Engine Contour
# Plot Engine Contour
if show_nozzle_plot:
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
    engine.exportGeometry(filename="Hopper Engine Contour 550 lbf 2_6_26")


#if export_HTC_hg:

# Show Engine Performance Outputs
if show_engine_perf_outputs:
    print(engine.R.get_summ_str())

### --- OUTPUTS --- ###

#regen_circuit.outputs()

end_time = time.time()
print(f"Runtime: {round(end_time - start_time, 2)} sec")