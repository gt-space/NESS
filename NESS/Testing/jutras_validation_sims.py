
#import regen_circuit as ta
from materials import Material
from regen_circuit import RegenCircuit

from engine import Engine

import time
import matplotlib.pyplot as plt

start_time = time.time()

### --- ENGINE PERFORMANCE INPUTS --- ###
name = "Jutras Engine PDR"
fuName = "Kerosene"  # IPA
oxName = "LOX"
thrust = 786 # [lbf] --> [N]
Pc = 261 # psia
Pe = 8 # psia
MRcore = 1.8
CR = 5
Lstar = 41.33 # in
cstarEff = 0.85
numPts = 300

### --- REGEN CIRCUIT INPUTS --- ###
t_w = 1 / 1000 # [mm] --> [m]
N = 72
C_w = 1 / 1000 # [mm] --> [m]
C_h = 4.88 / 1000 # [mm] --> [m]
coolantName = "Dodecane"
tot_coolant_mdot = 0.5 # kg/s
inlet_T_c = 293 # [K]
inlet_P_c = 460 * 6894.7 # [psia] --> [Pa]
material = Material("Pure Copper")

### --- REGEN CIRCUIT CONSTRAINTS --- ###
T_hw_max = 1050 # K
channel_dp_max = 100 # [psid] --> [Pa]

### --- MAIN CODE --- ###

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
    verbose=False
)

# Gas Props Testing
#engine.plot_gas_props()

# Create Regen Circuit
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

# Solve Regen Circuit
regen_circuit.solve_circuit(
    inlet_T_c=inlet_T_c,
    inlet_pressure=inlet_P_c
)

print(regen_circuit.h_hg_arr)
plt.figure()
plt.title("Bartz Testing")
plt.xlabel("Axial Position (in)")
plt.ylabel("HTC (W/m^2-K)")
plt.plot(engine.Contour_z, regen_circuit.h_hg_arr)
plt.show()

#regen_circuit.outputs() # Plotting function for regen circuit

#print(engine.AreaRatio)
#print(engine.M)

### --- PLOT OUTPUTS ---
'''
#engine.exportGeometry(filename="Hopper Engine Contour 10_23_25")
engine.R.geomObj.plot_geometry( title=f'Hopper Engine Profile 1023 - {thrust} lbf', show_grid=True )
print(engine.R.get_summ_str())

plt.figure()
plt.plot(engine.Contour_z, engine.AreaRatio)
plt.xlabel("Axial Position (m)")
plt.ylabel("Area Ratio")
plt.title("Area Ratio vs. Axial Position")

plt.figure()
plt.plot(engine.Contour_z, engine.M)
plt.xlabel("Axial Position (m)")
plt.ylabel("Mach Number")
plt.title("Mach Number vs. Axial Position")

plt.show()
'''
# CHAMBER INPUTS FROM PERFORMANCE ANALYSIS
# Chamber Radius
# Chamber Length
# Nozzle Geometry (Converging Ratio, Throat Diameter, Expansion Ratio, Bell Parameters)

#regen_geo_out = ta.regen_geometry(
#    r_ch=0.1016,
#    t_w= 1/1000,
#    N=50,
#    C_w=1/1000
#)
#print(regen_geo_out)


### --- SOLUTION --- ###



### --- OUTPUTS --- ###
#regen_circuit.outputs()

end_time = time.time()
print(f"Runtime: {round(end_time - start_time, 2)} sec")