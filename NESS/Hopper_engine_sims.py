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
from conical import conicalContour

### --- RUNTIME TRACKER --- ###
start_time = time.time()

### --- DESIGN & PLOT OPTIONS --- ###
design_engine = True
bell = False
design_regen = True
display_regen_contour_plot = False
display_regen_outputs = False
display_nozzle_mesh = False
show_nozzle_plot = True
export_nozzle = True
export_regen_chans = False
show_bartz_plot = True 
export_bartz_data = True 
show_gas_temp_plot = True
export_gas_temps = True

### --- ENGINE PERFORMANCE INPUTS --- ###
name = "Hopper SN1"
fuName = "Kerosene"  # Jet-A
oxName = "LOX"
thrust = 550 # [lbf] --> [N]
Pc = 300 # psia
Pe = 14.8 # psia
MRcore = 2
CR = 5
Lstar = 30 # in
cstarEff = 0.85
numPts = 100

### --- Conical Engine Inputs --- ###
chmbR = 1.35 # in
chmbL = 4.066 # in, cylindrical portion of chamber
contAngle = 30 # degrees
throatR = 0.603 # in
throatL = 0.5 # in
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
    conical_engine = conicalContour(chmbR=chmbR, chmbL=chmbL, contAngle=contAngle, 
    throatR=throatR, throatL=throatL, expAngle=exit_angle, exitR=exitR, numpts=numPts)
    conical_engine.makeContour()
              

# Gas Props Testing
#engine.plot_gas_props()

### --- REGEN CIRCUIT INPUTS --- ###
t_w = 1 / 1000 # [mm] --> [m]
N = 60
C_w = 1.5 / 1000 # [mm] --> [m]
C_h = 2.5 / 1000 # [mm] --> [m]
coolantName = "Ethanol"
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
    np.savetxt("Hopper_Bartz_HTC_550lbf_ConicalV1.csv", data, delimiter=",", header="Y Position (in), HTC (W/m^2-K)", comments="")

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
    np.savetxt("Hopper_Gas_Temps_550lbf_ConicalV1.csv", data, delimiter=",", header="Y Position (in), Gas Temperature (K)", comments="")

### --- PLOT OUTPUTS ---

# Plot Engine Contour
if show_nozzle_plot:
    if bell:
        engine.R.geomObj.plot_geometry( title=f'Hopper Engine Bell Nozzle Profile - {thrust} lbf', show_grid=True )
        engine.R.geomObj.plot_geometry( title=f'Hopper Engine Profile - {thrust} lbf', show_grid=True )
    else:
        conical_engine.plotContour()

if export_nozzle:
    if bell:
        engine.exportGeometry(filename="Hopper Engine Contour 550 lbf 2_6_26")
    else:
        conical_engine.saveContour(filename="Hopper Engine Contour 550 lbf ConicalV1.csv")

#if export_HTC_hg:


print(engine.R.get_summ_str())

#plt.show()

### --- OUTPUTS --- ###

#regen_circuit.outputs()

end_time = time.time()
print(f"Runtime: {round(end_time - start_time, 2)} sec")