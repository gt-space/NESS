import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import math

from CoolProp.CoolProp import PropsSI
from tqdm import tqdm
from rocketcea.cea_obj import CEA_Obj

# High-level layout

### --- PROPELLANTS --- ###
oxName = "N2O"
fuName = "Ethanol"

### --- PRESSURES & TEMPS --- ###
ox_tank_press = 500 # psia
fu_tank_press = 500 # psia
FU_temp = 293 # K
Pc = np.array([200, 350])
MR = np.array([2, 4])

### --- CdA NETWORK --- ###
OTV_CdA_range = np.array([9, 30.24]) # mm^2
FTV_CdA_range = np.array([2, 10.4]) # mm^2

OINJ_CdA = 32.3 # mm^2
FINJ_CdA = 11.2 # mm^2

### --- ENGINE SPECS --- ###
nominal_thrust = 500 # lbf
nominal_Pc = 350 # psia
nominal_MR = 3
thrust_range = [40, 100] # Thrust Range
cstarEff = 0.85
ideal_cstar = 1513 # m/s
eps = 3.79
throat_area = 0.877 # in^2
Pa = 14.8

# Plots
# Thrust over PcMR Range
# OX Valve % Open over PcMR Range
# FU Valve % Open over PcMR Range


def create_engine_map(
    ox_tank_press,
    fu_tank_press,
    OTV_CdA_range,
    FTV_CdA_range,
    OINJ_CdA,
    FINJ_CdA,
    throat_area,
    fu_temp,
    ideal_cstar,
    eps,
    Pa
):  
    # Unit Conversions
    ox_tank_press = ox_tank_press * 6894.7 # psia --> Pa
    fu_tank_press = fu_tank_press * 6894.7 # psia --> Pa
    OTV_CdA_range = OTV_CdA_range * 1E-6 # mm^2 --> m^2
    FTV_CdA_range = FTV_CdA_range * 1E-6 # mm^2 --> m^2
    OINJ_CdA = OINJ_CdA * 1E-6 # mm^2 --> m^2
    FINJ_CdA = FINJ_CdA * 1E-6 # mm^2 --> m^2
    At = throat_area * 0.00064516 # in^2 --> m^2
    
    # Create Valve CdA Grid
    OTV_CdA_map = np.linspace(OTV_CdA_range[0], OTV_CdA_range[1], 10)
    FTV_CdA_map = np.linspace(FTV_CdA_range[0], FTV_CdA_range[1], 10)
    OTV_CdA_vals, FTV_CdA_vals = np.meshgrid(OTV_CdA_map, FTV_CdA_map)

    # OTV/FTV Grid Print
    #print("OTV values:\n", OTV_CdA_vals)
    #print("\nFTV values:\n", FTV_CdA_vals)
    '''
    OTV_CdA = 32 / 1000000 # m^2
    FTV_CdA = 11 / 1000000 # m^2
    ox_temp = 293
    fu_temp = 293
    ideal_cstar = 1513 # m/s

    # calc_mdot testing
    calc_mdot(
        ox_tank_press,
        fu_tank_press,
        ox_temp,
        fu_temp,
        OTV_CdA,
        FTV_CdA,
        OINJ_CdA,
        FINJ_CdA,
        oxName,
        fuName,
        cstarEff,
        ideal_cstar,
        At
    )
    '''
    MR_grid = np.zeros_like(OTV_CdA_vals)
    Pc_grid = np.zeros_like(OTV_CdA_vals)
    thrust_grid = np.zeros_like(OTV_CdA_vals)

    for i in tqdm(range(len(FTV_CdA_map))):
        FTV_CdA = FTV_CdA_map[i]
        for j in range(len(OTV_CdA_map)):
            OTV_CdA = OTV_CdA_map[j]
            
            OX_mdot, FU_mdot, tot_mdot, Pc = calc_mdot(
                ox_tank_press,
                fu_tank_press,
                fu_temp,
                OTV_CdA,
                FTV_CdA,
                OINJ_CdA,
                FINJ_CdA,
                oxName,
                fuName,
                cstarEff,
                ideal_cstar,
                At
            )

            MR = OX_mdot / FU_mdot
            print(f"MR : {MR}")

            thrust = calc_thrust(
                oxName,
                fuName,
                Pc,
                MR,
                eps,
                At,
                Pa,
                tot_mdot
            )
            
            MR_grid[i, j] = MR
            Pc_grid[i, j] = Pc
            thrust_grid[i, j] = thrust
            print(f'Thrust : {thrust} N')

    # Normalize Valve CdAs & Pc Unit Conversion
    OTV_CdA_vals = OTV_CdA_vals / OTV_CdA_range[1]
    FTV_CdA_vals = FTV_CdA_vals / FTV_CdA_range[1]
    Pc_grid = Pc_grid / 6894.7 # Pa --> psia
    thrust_grid = thrust_grid * 0.224809 # N --> lbf
    ox_tank_press = ox_tank_press / 6894.7 # Pa --> psia
    fu_tank_press = fu_tank_press / 6894.7 # Pa --> psia

    # MR Heat Map
    plot_2d_heatmap(x_grid=OTV_CdA_vals, y_grid=FTV_CdA_vals, z_values=MR_grid, 
                    title=f"MR over Throttle Range | OX Tank = {ox_tank_press} psia | FU Tank = {fu_tank_press} psia",
                    xlabel="OX Valve Open %", ylabel="FU Valve Open %", cbar_label="MR")
    # Pc Heat Map
    plot_2d_heatmap(x_grid=OTV_CdA_vals, y_grid=FTV_CdA_vals, z_values=Pc_grid, 
                    title=f"Pc over Throttle Range | OX Tank = {ox_tank_press} psia | FU Tank = {fu_tank_press} psia",
                    xlabel="OX Valve Open %", ylabel="FU Valve Open %", cbar_label="Pc (psia)")
    
    # Thrust Heat Map
    plot_2d_heatmap(x_grid=OTV_CdA_vals, y_grid=FTV_CdA_vals, z_values=thrust_grid, 
                    title=f"Thrust over Throttle Range | OX Tank = {ox_tank_press} psia | FU Tank = {fu_tank_press} psia",
                    xlabel="OX Valve Open %", ylabel="FU Valve Open %", cbar_label="Thrust (lbf)")

    plt.show()



def plot_2d_heatmap(x_grid, y_grid, z_values, title="2D Heatmap", 
                    xlabel="X", ylabel="Y", cbar_label="Value",
                    cmap='viridis', figsize=(10, 8), show_values=False):
    """
    Plot a 2D heatmap from grid coordinates and values.
    
    Parameters:
    -----------
    x_grid : 2D numpy array
        X-coordinates of grid points (same shape as z_values)
    y_grid : 2D numpy array
        Y-coordinates of grid points (same shape as z_values)
    z_values : 2D numpy array
        Values to plot as heatmap
    title : str
        Plot title
    xlabel, ylabel : str
        Axis labels
    cbar_label : str
        Colorbar label
    cmap : str
        Colormap name (e.g., 'viridis', 'plasma', 'hot', 'coolwarm', 'jet')
    figsize : tuple
        Figure size (width, height)
    show_values : bool
        If True, display numeric values on each cell
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the heatmap using pcolormesh
    im = ax.pcolormesh(x_grid, y_grid, z_values, cmap=cmap, shading='auto')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label=cbar_label)
    
    # Set labels and title
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Optionally show values on cells
    if show_values:
        for i in tqdm(range(z_values.shape[0])):
            for j in range(z_values.shape[1]):
                text = ax.text(x_grid[i, j], y_grid[i, j], 
                             f'{z_values[i, j]:.2f}',
                             ha="center", va="center", 
                             color="white" if z_values[i, j] < np.mean(z_values) else "black",
                             fontsize=8)
    
    ax.set_aspect('equal')
    plt.tight_layout()
    
    return fig, ax

    
    

def calc_mdot(
    ox_tank_press,
    fu_tank_press,
    fu_temp,
    OTV_CdA,
    FTV_CdA,
    OINJ_CdA,
    FINJ_CdA,
    oxName,
    fuName,
    cstarEff,
    ideal_cstar,
    At
):
    OX_sys_CdA = np.sqrt((1 / ((1 / (OTV_CdA)**2) + (1 / (OINJ_CdA)**2))))
    FU_sys_CdA = np.sqrt((1 / ((1 / (FTV_CdA)**2) + (1 / (FINJ_CdA)**2))))
    
    ox_rho = PropsSI("D", "Q", 0, "P", ox_tank_press, oxName)
    fu_rho = PropsSI("D", "T", fu_temp, "P", fu_tank_press, fuName)

    actual_cstar = cstarEff * ideal_cstar # m/s

    # System of Equations Solver for Pc and mdot
    OX_mdot, FU_mdot, tot_mdot, Pc = sp.symbols('ox_mdot fu_mdot tot_mdot Pc')

    eq1 = sp.Eq(OX_sys_CdA * sp.sqrt(2 * ox_rho * (ox_tank_press - Pc)), OX_mdot)
    eq2 = sp.Eq(FU_sys_CdA * sp.sqrt(2 * fu_rho * (fu_tank_press - Pc)), FU_mdot)
    eq3 = sp.Eq(OX_mdot + FU_mdot, tot_mdot)
    eq4 = sp.Eq((Pc * At) / tot_mdot, actual_cstar)

    solution = sp.solve([eq1, eq2, eq3, eq4], (OX_mdot, FU_mdot, tot_mdot, Pc))

    # Extract the tuple
    sol_tuple = solution[0]

    # Map to a dictionary
    sol_dict = dict(zip((OX_mdot, FU_mdot, tot_mdot, Pc), sol_tuple))

    # Extract real parts and convert to float (removes tiny imaginary noise)
    numeric_solution = {var: complex(expr).real 
                    for var, expr in sol_dict.items()}

    # Output
    print(f"OX Tank Press : {ox_tank_press / 6894.7} psia")
    print(f"FU Tank Press : {fu_tank_press / 6894.7} psia")
    print(f"Pc : {numeric_solution[Pc] / 6894.7} psia")
    print(f"OX Mdot : {numeric_solution[OX_mdot]} kg/s")
    print(f"FU Mdot : {numeric_solution[FU_mdot]} kg/s")
    print(f"OX Rho : {ox_rho} kg/m^3")
    print(f"FU Rho : {fu_rho} kg/m^3")
    print(f"At : {At} m^2")

    OX_mdot = numeric_solution[OX_mdot]
    FU_mdot = numeric_solution[FU_mdot]
    tot_mdot = numeric_solution[tot_mdot]
    Pc = numeric_solution[Pc]

    return OX_mdot, FU_mdot, tot_mdot, Pc

    #return OX_mdot, FU_mdot

def calc_thrust(
    oxName,
    fuName,
    Pc,
    MR,
    eps,
    At,
    Pa,
    tot_mdot
):  
    # Unit Conversion
    Pa = Pa * 6894.7 # psia --> Pa

    C = CEA_Obj(oxName=oxName, fuelName=fuName)
    Pe = Pc / C.get_PcOvPe(Pc=Pc, MR=MR, eps=eps)
    M = C.get_MachNumber(Pc=Pc, MR=MR, eps=eps)
    _, _, a_exit = C.get_SonicVelocities(Pc=Pc, MR=MR, eps=eps)
    a_exit = a_exit * 0.3048 # ft/s --> m/s
    u_e = M * a_exit

    momentum_thrust = tot_mdot * u_e
    pressure_thrust = At * eps * (Pe - Pa)

    print(f"Ue : {u_e} m/s")
    print(f"M : {M}")
    print(f"a : {a_exit} m/s")
    print(f"Pe : {Pe / 6894.7} psia")
    print(f"Momentum Thrust : {momentum_thrust} N")
    print(f"Pressure Thrust : {pressure_thrust} N")

    thrust = momentum_thrust + pressure_thrust

    return thrust


# Main Code
create_engine_map(
    ox_tank_press=ox_tank_press,
    fu_tank_press=fu_tank_press,
    OTV_CdA_range=OTV_CdA_range,
    FTV_CdA_range=FTV_CdA_range,
    OINJ_CdA=OINJ_CdA,
    FINJ_CdA=FINJ_CdA,
    throat_area=throat_area,
    fu_temp=FU_temp,
    ideal_cstar=ideal_cstar,
    eps=eps,
    Pa=Pa
)
