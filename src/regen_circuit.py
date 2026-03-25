import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

from scipy.optimize import fsolve 
from utils import solve_system, check_defined_vars
from CoolProp.CoolProp import PropsSI
from fluid import Fluid
from tqdm import tqdm

'''
All functions, equations, and examples are documented here:

https://www.notion.so/Regen-Code-Overview-27a27990d78e806b8886ec2d1a555117?source=copy_link

'''

# Start with ONLY fixed width and height channels.

class RegenCircuit:

    def __init__(
        self, 
        t_w, 
        material,
        tot_coolant_mdot,
        coolantName,
        C_h,
        engine,
        N = None,
        f_w = None,
        C_w = None,
    ):
        self.t_w = t_w
        self.material = material
        self.tot_coolant_mdot = tot_coolant_mdot
        self.N = N
        self.chan_mdot = self.tot_coolant_mdot / self.N
        self.coolantName = coolantName
        self.f_w = f_w
        self.C_w = C_w
        self.C_h = C_h
        self.ch_AR = self.C_h / self.C_w
        self.engine = engine # Engine object from engine.py

        # Check that the inputted regen geometry is compatible with the Contour.
        self.regen_geometry(
            r_throat=self.engine.r_throat,
            t_w=self.t_w,
            N=self.N,
            f_w=self.f_w,
            C_w=self.C_w,
        )

        # User Output - signifies completion of circuit creation
        print(" === Regen Circuit has been successfully created! ===")

    '''
    def define_constraints(
        T_hw_MAX,


    )
    '''

    def bartz(
        self,
        engine,
        Cp,
        Pr,
        mu,
        Pc,
        AreaRatio,
        gamma,
        M,
        T_hg

    ) -> float:
        '''
        Calculate the Heat Transfer Coefficient for hot gas in the chamber and nozzle from combustion gas properties using the well known Bartz Correlation.

        See Eqn. [] for more documentation.

        Inputs:

        Outputs:

        Raises:

        '''

        omega = 0.6 # From Bartz Paper

        recovery_factor = math.sqrt(Pr)
        T_wg = T_hg * (1 + (recovery_factor / 2)*(gamma - 1)*M**2)
        # Wall gas should be slightly colder than the freestream gas - got this eqn from PK Week 10 Aero

        h_gA = 0.026 / (engine.Dt**0.2)
        h_gB = (Cp * mu**0.2) / (Pr**0.6)
        h_gC = (Pc / engine.cstar)**0.8
        h_gD = (engine.Dt / engine.r_up_throat)**0.1
        h_gE = (1/AreaRatio)**0.9
        
        sigma_A = (0.5 * (T_wg / T_hg) * (1 + 0.5*(gamma - 1)*M**2) + 0.5)**(0.8 - 0.2*omega)
        sigma_B = (1 + 0.5*(gamma - 1)*M**2)**(0.2*omega)
        sigma = 1 / (sigma_A * sigma_B)

        h_hg = h_gA * h_gB * h_gC * h_gD * h_gE * sigma

        # Bartz Debug
        '''
        print(f"T_wg : {T_wg}")
        print(f"h_gA : {h_gA}")
        print(f"h_gB : {h_gB}")
        print(f"h_gC : {h_gC}")
        print(f"h_gD : {h_gD}")
        print(f"h_gE : {h_gE}")
        print(f"sigma_A : {sigma_A}")
        print(f"sigma_B : {sigma_B}")
        print(f"sigma : {sigma}")
        '''

        return h_hg

    def gnielinski(
        self,
        fluid,
        eps = None,
        debug = False

    ) -> float:
        '''
        Calculate the Coolant Heat Transfer Coefficient using Gnielinski's correlation.
        
        '''
        
        P = 2 * (self.C_w + self.C_h) # Channel Wetted Perimeter
        A_c = self.C_w * self.C_h # Channel Flow Area
        D_h = (4 * A_c) / P # Channel Hydraulic Diameter
        Re = (self.chan_mdot * D_h) / (A_c * fluid.mu) # Channel Reynold's Number
        Pr = fluid.Pr # Prandlt Number

        f_d = self.calc_friction_factor(
            D_h=D_h,
            Re=Re,
            eps=eps
        )

        # Gnielinksi
        Nu = ((f_d / 8) * (Re - 1000) * Pr) / (1 + 12.7 * ((f_d/8)**(1/2)) * (Pr**(2/3) - 1)) # Nusselt Number
        h_l = (Nu * fluid.k) / (D_h)  # Coolant HTC

        # Dittus-Boelter
        Nu_db = 0.023 * (Re**(0.8)) * (Pr**0.4)
        h_l_db = (Nu_db * fluid.k) / D_h

        if debug:
            print(f"Channel Wetted Perimeter : {P} m")
            print(f"Flow Area : {A_c} m^2")
            print(f"Hydraulic Diameter : {D_h} m")
            print(f"Re : {Re}")
            print(f"Pr : {Pr}")
            print(f"f_d : {f_d}")
            print(f"Nu : {Nu}")
            print(f"Gnielinksi HTC : {h_l} W/m^2-K")
            print(f"DB HTC : {h_l_db} W/m^2-K")
            print(f"Nusselt DB : {Nu_db}")

        return h_l
    
    def calculate_stress(
        self,
        T_hw,
        T_cw,
        material,
        qdot,
        P_c,
        P_hg,
    ):
        '''
        Calculate the stress in the liner due to thermal gradient between the hot gas and coolant.
        '''
        # Update temperature-dependent material properties (if available).
        # For hardcoded-only materials, this call leaves values unchanged.
        mat_props = material.update_material_properties(T_hw)

        E = mat_props.get("E", getattr(material, "E", None))
        alpha = mat_props.get("alpha", getattr(material, "alpha", None))
        v = mat_props.get("v", getattr(material, "v", None))
        k = mat_props.get("k", getattr(material, "k", None))

        if any(val is None for val in [E, alpha, v, k]):
            raise ValueError(
                "Material is missing required properties for stress calculation. "
                "Need E, alpha, v, and k."
            )

        # Tangential Stress
        sigma_pt = ((P_c - P_hg) / 2) * (self.C_w / self.t_w)**2
        sigma_tt = (E * alpha * qdot * self.t_w) / (2 * (1 - v) * k)
        sigma_t = sigma_pt + sigma_tt

        # Longitudinal Stress
        # pressure ignored maybe update later
        sigma_tl = E * alpha * (T_hw - T_cw)
        sigma_l = sigma_tl
        
        return sigma_pt, sigma_tt, sigma_t, sigma_tl, sigma_l

    def thermal_network(
        self,
        r_eng,
        t_w, 
        c_w,
        dz,
        material,
        h_hg,
        h_c,
        qdot = None,
        T_hg = None,
        T_hw = None,
        T_cw = None,
        T_c = None,

    ) -> float:
        '''
        Specify 2 of the following only:
            - qdot, T_hg, T_hw, T_cw, T_c
        '''

        qdot_new = qdot
        T_hg_new = T_hg
        T_hw_new = T_hw
        T_cw_new = T_cw
        T_c_new = T_c

        qdot, T_hg, T_hw, T_cw, T_c = sp.symbols('qdot T_hg T_hw T_cw T_c')

        inputs = {
            qdot : qdot_new,
            T_hg : T_hg_new,
            T_hw : T_hw_new,
            T_cw : T_cw_new,
            T_c : T_c_new
        }

        # this not working??
        if check_defined_vars(inputs, min_required=2):
            pass
        else:
            TypeError("1 Input must be defined for the function thermal_network.")

        ## Fin Effective Area Calculations ##
        # Single Fin #
        m = math.sqrt((h_c * (2 * self.f_w + 2 * dz)) / (material.k * self.f_w * dz))
        eta_f = np.tanh(m * self.C_h)

        # Multi-Fin #
        A_b = self.N * c_w * dz
        A_hwall_ch = (2 * np.pi * r_eng * dz) / self.N
        A_f = 2 * self.N * self.C_h * dz
        A_t = A_b + A_f
        eta_o = (A_b + eta_f * A_f) / A_t
        A_ch_eff = (eta_o * A_t) / self.N

        # Thermal Resistances
        R_hghw = 1 / (h_hg * A_hwall_ch)
        R_wall = t_w / (material.k * A_hwall_ch)
        R_cwc = 1 / (h_c * A_ch_eff)
        R_tot = R_hghw + R_wall + R_cwc

        '''
        # Debug Prints
        print(f"T_hg : {T_hg_new} K")
        print(f"T_c : {T_c_new} K")
        print(f"A_ch_eff : {A_ch_eff} m^2")
        print(f"R_hghw : {R_hghw}")
        print(f"R_wall : {R_wall}")
        print(f"R_cwc : {R_cwc}")
        print(f"R_tot : {R_tot}")
        print(f"Material k : {material.k} W/mK")
        print(f"H_hg : {h_hg} W/m^2-K")
        print(f"H_c : {h_c} W/m^2-K")
        print(f"Wall Thickness : {t_w} m")
        '''
        
        #eq1 = sp.Eq(qdot, (T_hg - T_c) / R_tot)
        eq2 = sp.Eq(qdot, (T_hg - T_hw) / R_hghw)
        eq3 = sp.Eq(qdot, (T_hw - T_cw) / R_wall)
        eq4 = sp.Eq(qdot, (T_cw - T_c) / R_cwc)

        eqns = [eq2, eq3, eq4]
        vars_all = [qdot, T_hg, T_hw, T_cw, T_c]

        output = solve_system(eqns, vars_all, inputs)

        # Assign Output
        T_hg = output[T_hg]
        T_c = output[T_c]
        T_cw = output[T_cw]
        T_hw = output[T_hw]
        qdot = output[qdot]

        return T_hg, T_c, T_cw, T_hw, qdot, A_ch_eff
    
    def solve_circuit(
        self,
        inlet_T_c,
        inlet_pressure,
        circuit_inlet = None, # [m]
        
    ):
        # Set fluid inlet position by default to bottom of nozzle
        if circuit_inlet is None:
            circuit_inlet = len(self.engine.Contour_z) - 1
            circuit_inlet_pos = self.engine.Contour_z[circuit_inlet] * 0.0254 # [in] -- > [m]
            circuit_outlet = 0

            print(circuit_inlet)
            print(circuit_outlet)

        
        # can input a min num of Pts, then find velocity at max velocity in channel, then set timestep to dz / max velocity
        #dz = round(self.Contour.z[-1] / self.numPts, 0)
        time_tracker = 0

        # Calculate dz - use the arc length based on the specified dz NOT the dz since that will underestimate volume
        # will lead to underconservative heating
        dz = (circuit_inlet_pos - (-self.engine.Contour_z[0] * 0.0254)) / self.engine.numPts
        dr = np.diff(self.engine.Contour_r * 0.0254)
        ds = np.sqrt(dz**2 + dr**2)

        # repeat end value so arrays have same shapes
        ds_end = ds[-1]
        ds = np.append(ds, ds_end)

        #print(dr)
        #print(dz)
        #print(ds)

        # Initialize Arrays (preallocate once for faster writes).
        n_stations = circuit_inlet - circuit_outlet + 1
        self.T_hw_arr = np.empty(n_stations)
        self.T_cw_arr = np.empty(n_stations)
        self.T_c_arr = np.empty(n_stations)
        self.Qdot_arr = np.empty(n_stations)
        self.coolant_vel_arr = np.empty(n_stations)
        self.P_c_arr = np.empty(n_stations)
        self.h_hg_arr = np.empty(n_stations)
        self.h_c_arr = np.empty(n_stations)
        self.DP_arr = np.empty(n_stations)
        self.Re = np.empty(n_stations)
        self.sigma_pt = np.empty(n_stations)
        self.sigma_tt = np.empty(n_stations)
        self.sigma_t = np.empty(n_stations)
        self.sigma_pl = np.empty(n_stations)
        self.sigma_l = np.empty(n_stations)

        # T and P Tracker
        P_c = inlet_pressure
        T_c = inlet_T_c
        
        for i in range(len(self.engine.T)):
            h_hg = self.bartz(
                engine=self.engine,
                Cp=self.engine.Cp[i],
                Pr=self.engine.Pr[i],
                mu=self.engine.mu[i],
                Pc=self.engine.Pc,
                AreaRatio=self.engine.AreaRatio[i],
                gamma=self.engine.gamma[i],
                M=self.engine.M[i],
                T_hg=self.engine.T[i]
            )
        
        # Get Boiling Point
        self.coolant_boiling = self.calculate_boiling_point(P_c, T_c, self.coolantName, plot_cp_curve=False)
        print(f" Coolant BP : {self.coolant_boiling} K")
        
        # Regen Circuit Loop
        for i in tqdm(range(circuit_inlet, circuit_outlet - 1, -1)):
            # calculate coolant fluid properties from (inlet pressure and T_c)
            # calculate coolant volume from channel geometry (C_w * C_h * dz)
            # solve thermal network with inputted T_c and T_hg --> outputs T_hw, T_cw, and Qdot
            # calculate Q from Qdot and timestep and use Q = mCpdT to calculate new T_c
            # Use coolant velocity to calculate new z using timestep
            print(f"T_c : {T_c} K")
            print(f"P_c : {P_c / 6894.7} psia")
            print(f"Ind : {i}")
            coolant = Fluid("T", T_c, "P", P_c, self.coolantName)

            coolantVolume = self.C_w * self.C_h * ds[i]
            m = coolantVolume * coolant.rho

            h_hg = self.bartz(
                engine=self.engine,
                Cp=self.engine.Cp[i],
                Pr=self.engine.Pr[i],
                mu=self.engine.mu[i],
                Pc=self.engine.Pc,
                AreaRatio=self.engine.AreaRatio[i],
                gamma=self.engine.gamma[i],
                M=self.engine.M[i],
                T_hg=self.engine.T[i]
            )

            print(f"Bartz HTC : {h_hg} W/m^2-K")

            h_c = self.gnielinski(fluid=coolant, debug=True)
            print(f"Coolant HTC : {h_c} W/m^2-K")


            T_hg, T_c, T_cw, T_hw, qdot, A_ch_eff = self.thermal_network(
                r_eng=self.engine.Contour_r[i] * 0.0254,
                t_w=self.t_w,
                c_w=self.C_w,
                dz=ds[i],
                material=self.material,
                h_hg=h_hg,
                h_c=h_c,
                T_hg=self.engine.T[i],
                T_c=T_c
                )
            q_flux = qdot / ds[i] * ((2 * np.pi * self.engine.r_ch) / self.N) # W/m^2
            
            # Stress Calculations
            sigma_pt, sigma_tt, sigma_t, sigma_pl, sigma_l = self.calculate_stress(
                T_hw=T_hw,
                T_cw=T_cw,
                material=self.material,
                qdot=q_flux,
                P_c=P_c,
                P_hg=self.engine.P[i] * 6894.76,
            )
            
            '''
            print("Thermal Network Outputs")
            print(f"qdot: {qdot} W")
            print(f"T_hg : {T_hg} K")
            print(f"T_hw : {T_hw} K")
            print(f"T_cw : {T_cw} K")
            print(f"T_c : {T_c} K")
            '''

            coolant_vel = self.chan_mdot / (coolant.rho * self.C_w * self.C_h)
            D_h = (4*(self.C_w * self.C_h)) / (2*(self.C_w + self.C_h))
            Re = (coolant.rho * coolant_vel * D_h) / coolant.mu
            DP = self.calc_DP(coolant, coolant_vel, ds[i], D_h, Re)
            P_c_new = P_c - DP
            P_c = P_c_new

            dt = ds[i] / coolant_vel
            T_c_new = ((qdot * dt) / (m * coolant.Cp)) + T_c
            T_c = T_c_new

            # Debug Print
            #print(f"T_c : {T_c} K")
            #print(f"P_c : {P_c} psia")
            print(f"Coolant Velocity : {coolant_vel} m/s")
            print(f"Re : {Re}")
            print(f"DP : {DP / 6894.7} psid")

            # Boiling Check
            if T_c > self.coolant_boiling:
                raise ValueError(f" === Coolant boiled at {round(self.coolant_boiling, 2)} K!")
                break
            
            # Indexed update (faster than repeated np.append in loop).
            out_idx = i - circuit_outlet
            self.T_hw_arr[out_idx] = T_hw
            self.T_cw_arr[out_idx] = T_cw
            self.T_c_arr[out_idx] = T_c_new
            self.Qdot_arr[out_idx] = qdot
            self.coolant_vel_arr[out_idx] = coolant_vel
            self.P_c_arr[out_idx] = P_c_new
            self.h_hg_arr[out_idx] = h_hg
            self.h_c_arr[out_idx] = h_c
            self.DP_arr[out_idx] = DP
            self.Re[out_idx] = Re
            self.sigma_pt[out_idx] = sigma_pt
            self.sigma_tt[out_idx] = sigma_tt
            self.sigma_t[out_idx] = sigma_t
            self.sigma_pl[out_idx] = sigma_pl
            self.sigma_l[out_idx] = sigma_l

            #break
    
    def calculate_boiling_point(
        self,
        P_c,
        T_c,
        fluid,
        plot_cp_curve=False
    ):
        """
        Determine a reference "boiling limit" temperature.

        - If saturation at the given pressure exists, returns the saturated boiling point (Q=0)
          and treats the current state as subcooled vs superheated based on `T_c`.
        - If the state is supercritical (or saturation temperature cannot be computed), returns
          the temperature at which Cp(T) peaks at this pressure (pseudocritical Cp peak).
        - Reference Paper: https://www.osti.gov/servlets/purl/4159664
        """

        # Try to compute subcritical saturation boiling temperature at this pressure.
        try:
            subcritical_boiling_point = PropsSI("T", "P", P_c, "Q", 0, fluid)  # [K]
            # Subcritical: reference limit is the saturation boiling point.
            if T_c <= subcritical_boiling_point:
                return subcritical_boiling_point
        except Exception:
            # If saturation data is unavailable, assume we're in (or near) supercritical regime.
            subcritical_boiling_point = None

        # Supercritical (or unable to compute saturation): use pseudocritical temperature
        # where Cp(T) reaches its maximum at this pressure.
        T_min = 200.0  # [K] broad sweep guard
        T_max = 2000.0 # [K] broad sweep guard
        n_samples = 200
        T_candidates = np.linspace(T_min, T_max, n_samples)

        cp_values = []
        valid_T = []

        for T in T_candidates:
            try:
                cp = PropsSI("C", "T", float(T), "P", float(P_c), self.coolantName)  # [J/kg/K]
                if np.isfinite(cp):
                    cp_values.append(float(cp))
                    valid_T.append(float(T))
            except Exception:
                # CoolProp can fail outside the fluid's valid domain for a given P.
                continue

        if not cp_values:
            raise ValueError(
                f"Could not compute Cp(T) sweep for fluid '{self.coolantName}' at "
                f"P={P_c:.3e} Pa. No valid CoolProp states were found in the "
                f"temperature sweep range."
            )

        peak_idx = int(np.argmax(cp_values))
        supercritical_boiling_point = valid_T[peak_idx]

        if plot_cp_curve:
            plt.figure()
            plt.plot(valid_T, cp_values, label=f"{self.coolantName} Cp(T)")
            plt.axvline(
                supercritical_boiling_point,
                linestyle="--",
                label=f"Cp peak T = {supercritical_boiling_point:.2f} K"
            )
            plt.scatter(
                [supercritical_boiling_point],
                [cp_values[peak_idx]],
                zorder=3
            )
            plt.title(f"Cp vs Temperature at P={P_c/6894.7:.2f} psia")
            plt.xlabel("Temperature [K]")
            plt.ylabel("Cp [J/kg-K]")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.show()

        return supercritical_boiling_point


    def calc_DP(
        self,
        fluid,
        velocity,
        L,
        D_h,
        Re,
        eps = None
    ):
        
        f_D = self.calc_friction_factor(D_h, Re, eps)
        print(f"Friction Factor : {f_D}")
        print(f"D_h : {D_h} m")
        print(f"Rho : {fluid.rho} kg/m^3")
        dP = 0.5 * fluid.rho * velocity**2 * (L / D_h) * f_D # Pa

        #print(f"f_d : {f_D}")
        print(f"Length : {L} m")

        return dP

    def calc_friction_factor(
        self,
        D_h, 
        Re,
        eps = None
    ):
        '''
        Use the Colebrook-White Equation to calculate the Darcy Friction Factor.

        See Eqn. [] for more documentation.

        Inputs:

        Outputs:

        Raises:

        '''

        def colebrook_white(f, eps, D_h, Re):

            return (1 / math.sqrt(f)) + 2*math.log((eps / (3.8 * D_h)) + (2.51 / (Re * math.sqrt(f))), 10) 

        f_guess = 0.04
        
        if eps is None:
            eps = 0.015 / 1000 # [mm] UPDATE - ALSO UPDATE FOR LPBF SURFACES

        f = fsolve(colebrook_white, f_guess, args=(eps, D_h, Re))[0]

        return f


    def regen_geometry(
        self,
        r_throat,
        t_w,
        N = None,
        f_w = None,
        C_w = None,

    ) -> float:
        N_new = N
        f_w_new = f_w
        C_w_new = C_w

        N, f_w, C_w = sp.symbols('N f_w C_w')

        inputs = {
            N : N_new,
            f_w : f_w_new,
            C_w : C_w_new
        }

        if check_defined_vars(inputs, min_required=2):
            pass
        else:
            TypeError("2 Inputs must be defined for the function regen_geometry.")

        r_throat = r_throat * 0.0254 # [in] --> [m]
        eq1 = sp.Eq(f_w, ((2*np.pi*(r_throat + t_w)) / N) - C_w)

        eqns = [eq1]
        vars_all = [N, f_w, C_w]

        output = solve_system(eqns, vars_all, inputs)
        if output[f_w] <= 0:
            raise ValueError(f"ERROR : Fin Width ({round(output[f_w] * 1000, 3)} mm) is negative or 0. Reduce N ({output[N]}) OR C_w ({output[C_w] * 1000} mm)")

        # Reset all variables
        self.f_w = output[f_w]
        self.N = output[N]
        self.C_w = output[C_w]

        print(f"Fin Width: {self.f_w}")
        print(f"Throat Radius: {r_throat}")

        return output

    def plot_regen_geometry(
        self,
        engine,
        show_dimensions = True    
    ):
        '''
        Plot regen channels.
        '''

        # Unit Conversions
        engine.r_ch = engine.r_ch * 25.4

        outer_wall_thickness = 5 # [mm]
        inner_wall_thickness = self.t_w * 1000  # [mm]

        inner_radius = engine.r_ch
        hot_wall_inner_radius = inner_radius - inner_wall_thickness
        channel_outer_radius = inner_radius + self.C_h * 1000
        outer_radius = channel_outer_radius + outer_wall_thickness
        
        # Calculate angular spacing - distribute evenly around full circle
        total_angle_per_channel = 2 * np.pi / self.N  # Total angle per channel+rib pair
        
        # Calculate angular widths
        channel_angle = (self.C_w * 1000) / inner_radius  # Angular width of channel
        wall_angle = (self.t_w * 1000) / inner_radius    # Angular width of rib
        
        # Verify geometry fits
        if (channel_angle + wall_angle) > total_angle_per_channel:
            print(f"WARNING: Channel + rib width ({channel_angle + wall_angle:.4f} rad) exceeds available space ({total_angle_per_channel:.4f} rad)")
            print(f"Adjusting channel width to fit...")
            channel_angle = total_angle_per_channel - wall_angle
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Draw hot wall annulus (red filled region) - inner wall thickness
        theta = np.linspace(0, 2*np.pi, 200)
        
        # Outer boundary of hot wall (inner radius)
        x_hot_outer = inner_radius * np.cos(theta)
        y_hot_outer = inner_radius * np.sin(theta)
        
        # Inner boundary of hot wall
        x_hot_inner = hot_wall_inner_radius * np.cos(theta[::-1])
        y_hot_inner = hot_wall_inner_radius * np.sin(theta[::-1])
        
        # Create hot wall annulus polygon
        hot_wall_annulus = np.vstack([
            np.column_stack([x_hot_outer, y_hot_outer]),
            np.column_stack([x_hot_inner, y_hot_inner])
        ])
        
        hot_wall_poly = Polygon(hot_wall_annulus, facecolor='red', 
                               edgecolor='darkred', linewidth=2, alpha=0.7, zorder=10)
        ax.add_patch(hot_wall_poly)
        
        # Draw combustion chamber (inside hot wall)
        ax.fill(x_hot_inner, y_hot_inner, color='lightcoral', alpha=0.3, zorder=9)
        ax.plot(x_hot_inner, y_hot_inner, 'r--', linewidth=1.5, alpha=0.5, zorder=9)
        
        # Draw outer wall (solid shell) - background layer
        theta_outer = np.linspace(0, 2*np.pi, 200)
        x_outer = outer_radius * np.cos(theta_outer)
        y_outer = outer_radius * np.sin(theta_outer)
        
        # Fill entire region from channel_outer_radius to outer_radius as solid outer wall
        outer_annulus_outer = np.column_stack([x_outer, y_outer])
        x_ch_bound = channel_outer_radius * np.cos(theta_outer[::-1])
        y_ch_bound = channel_outer_radius * np.sin(theta_outer[::-1])
        outer_annulus_inner = np.column_stack([x_ch_bound, y_ch_bound])
        outer_annulus = np.vstack([outer_annulus_outer, outer_annulus_inner])
        
        outer_wall_poly = Polygon(outer_annulus, facecolor='darkgray', 
                                 edgecolor='black', linewidth=2, alpha=0.4, zorder=1)
        ax.add_patch(outer_wall_poly)
        
        # Draw channel outer boundary
        ax.plot(x_ch_bound, y_ch_bound, 'b--', linewidth=1.5, alpha=0.5, 
                label='Channel depth', zorder=2)
        
        # Draw the annular region between inner_radius and channel_outer_radius
        # This will be filled with channels and ribs
        x_annulus_outer = channel_outer_radius * np.cos(theta_outer)
        y_annulus_outer = channel_outer_radius * np.sin(theta_outer)
        x_annulus_inner = inner_radius * np.cos(theta_outer[::-1])
        y_annulus_inner = inner_radius * np.sin(theta_outer[::-1])
        annulus_full = np.vstack([
            np.column_stack([x_annulus_outer, y_annulus_outer]),
            np.column_stack([x_annulus_inner, y_annulus_inner])
        ])
        
        # Fill the annulus with rib material (gray background)
        annulus_poly = Polygon(annulus_full, facecolor='gray', 
                              edgecolor='none', alpha=0.6, zorder=2)
        ax.add_patch(annulus_poly)
        
        # Draw channels (cut out from rib material)
        channels = []
        
        current_angle = 0
        for i in range(self.N):
            # Channel start and end angles
            theta_start = current_angle
            theta_end = current_angle + channel_angle
            
            # Create channel polygon (cuts through the rib material)
            n_points = 20
            theta_channel = np.linspace(theta_start, theta_end, n_points)
            
            # Inner arc
            x_ch_inner = inner_radius * np.cos(theta_channel)
            y_ch_inner = inner_radius * np.sin(theta_channel)
            
            # Outer arc
            x_ch_outer = channel_outer_radius * np.cos(theta_channel[::-1])
            y_ch_outer = channel_outer_radius * np.sin(theta_channel[::-1])
            
            # Combine into polygon
            x_channel = np.concatenate([x_ch_inner, x_ch_outer])
            y_channel = np.concatenate([y_ch_inner, y_ch_outer])
            
            channel_poly = Polygon(np.column_stack([x_channel, y_channel]), 
                                facecolor='lightblue', edgecolor='blue', 
                                linewidth=1, alpha=0.8, zorder=5)
            ax.add_patch(channel_poly)
            channels.append(channel_poly)
            
            # Move to next channel (use total angle per channel to distribute evenly)
            current_angle = (i + 1) * total_angle_per_channel
        
        # Add dimensions if requested
        if show_dimensions:
            # Inner radius only
            ax.plot([0, inner_radius], [0, 0], 'k-', linewidth=0.5, zorder=12)
            ax.annotate('', xy=(inner_radius, 0), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1.5), zorder=12)
            ax.text(inner_radius/2, -0.05*outer_radius, f'R_inner = {inner_radius:.2f} mm',
                ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='white'), zorder=12)
        
        # Formatting
        ax.set_aspect('equal')
        ax.set_xlim(-outer_radius*1.3, outer_radius*1.3)
        ax.set_ylim(-outer_radius*1.3, outer_radius*1.3)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('X [mm]', fontsize=12)
        ax.set_ylabel('Y [mm]', fontsize=12)
        ax.set_title(f'Regenerative Cooling Channel Cross-Section\n{self.N} Channels', 
                    fontsize=14, fontweight='bold')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', edgecolor='darkred', label='Hot wall', alpha=0.7),
            Patch(facecolor='lightcoral', edgecolor='r', label='Combustion chamber', alpha=0.3),
            Patch(facecolor='lightblue', edgecolor='b', label='Coolant channels'),
            Patch(facecolor='gray', edgecolor='k', label='Structural ribs'),
            Patch(facecolor='darkgray', edgecolor='k', label='Outer wall', alpha=0.4),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
        
        # Add text summary (top left)
        textstr = f'Chamber ID: {engine.r_ch:.2f} mm\n'
        textstr += f'Hot wall thickness: {inner_wall_thickness:.2f} mm\n'
        textstr += f'Channel height: {self.C_h * 1000:.2f} mm\n'
        textstr += f'Channel width: {self.C_w * 1000:.2f} mm\n'
        textstr += f'Rib thickness: {self.t_w * 1000:.2f} mm\n'
        textstr += f'Outer wall: {outer_wall_thickness:.2f} mm\n'
        textstr += f'Number of channels: {self.N}\n'
        textstr += f'Total OD: {2*outer_radius:.2f} mm'
        
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.show()

    def generate_single_channel_curves(
    self,
    output_dir=".",
    plot=True
):
        """
        Generate 4 corner curves for a single rectangular cooling channel for SolidWorks.
        
        The 4 curves represent the corners of one rectangular channel:
        - Corner 1: Inner radius, channel start (hot wall side, one edge)
        - Corner 2: Inner radius, channel end (hot wall side, other edge)
        - Corner 3: Outer radius, channel end (cold side, other edge)
        - Corner 4: Outer radius, channel start (cold side, one edge)
        
        Parameters from self:
        ---------------------
        self.C_h : float
            Channel height (radial) in meters
        self.C_w : float
            Channel width (circumferential) in meters
        self.t_w : float
            Hot wall thickness in meters
        self.N : int
            Number of channels (for reference)
        self.f_w : float
            Fin/rib width in meters
        self.engine.Contour_z : array
            Axial positions in inches
        self.engine.Contour_r : array
            Radial positions in inches
            
        Output Files:
        -------------
        - corner1_inner_start.txt : Hot wall, channel start edge
        - corner2_inner_end.txt : Hot wall, channel end edge
        - corner3_outer_end.txt : Channel outer, end edge
        - corner4_outer_start.txt : Channel outer, start edge
        
        Each file format: X, Y, Z (Cartesian coordinates)
        
        SolidWorks Import:
        ------------------
        1. Insert -> Curve -> Curve Through XYZ Points
        2. Import all 4 corner curves
        3. Use Loft or Boundary Surface to create channel
        4. Circular pattern around axis (N times)
        """
        
        # Unit conversions - keep everything in inches
        # Copy to local variables so we don't modify self attributes
        z_contour = np.asarray(self.engine.Contour_z).copy()  # inches
        r_contour = np.asarray(self.engine.Contour_r).copy()  # inches
        
        # Convert to inches for local use only
        C_h = self.C_h * 39.3701  # meters to inches (local variable)
        C_w = self.C_w * 39.3701  # meters to inches (local variable)
        t_w = self.t_w * 39.3701  # meters to inches (local variable)
        f_w = self.f_w * 39.3701  # meters to inches (local variable)
        
        # Calculate radii
        r_inner = r_contour + t_w  # Inner surface (where channel starts radially)
        r_outer = r_contour + C_h  # Outer surface (channel outer boundary)
        
        # Calculate angular positions for one channel
        # Channel starts at angle = 0
        # Channel width in angular terms
        theta_channel_start = 0.0
        theta_channel_end = C_w / r_contour[0]  # Approximate angle at first station
        
        # For each axial station, calculate the 4 corners in 3D space
        n_points = len(z_contour)
        
        # Initialize corner arrays (X, Y, Z coordinates)
        corner1 = np.zeros((n_points, 3))  # Inner radius, start angle
        corner2 = np.zeros((n_points, 3))  # Inner radius, end angle
        corner3 = np.zeros((n_points, 3))  # Outer radius, end angle
        corner4 = np.zeros((n_points, 3))  # Outer radius, start angle
        
        for i in range(n_points):
            z = z_contour[i]
            r_in = r_inner[i]
            r_out = r_outer[i]
            
            # Recalculate angular width at each station (accounts for changing radius)
            theta_end = C_w / r_in
            theta_start = 0.0
            
            # Corner 1: Inner radius, start angle (θ=0)
            corner1[i, 0] = r_in * np.cos(theta_start)  # X
            corner1[i, 1] = r_in * np.sin(theta_start)  # Y
            corner1[i, 2] = z  # Z
            
            # Corner 2: Inner radius, end angle
            corner2[i, 0] = r_in * np.cos(theta_end)
            corner2[i, 1] = r_in * np.sin(theta_end)
            corner2[i, 2] = z
            
            # Corner 3: Outer radius, end angle
            corner3[i, 0] = r_out * np.cos(theta_end)
            corner3[i, 1] = r_out * np.sin(theta_end)
            corner3[i, 2] = z
            
            # Corner 4: Outer radius, start angle
            corner4[i, 0] = r_out * np.cos(theta_start)
            corner4[i, 1] = r_out * np.sin(theta_start)
            corner4[i, 2] = z
        
        # Store in dictionary
        corners = {
            'corner1_inner_start': corner1,
            'corner2_inner_end': corner2,
            'corner3_outer_end': corner3,
            'corner4_outer_start': corner4
        }
        
        # Write output files
        filenames = {
            'corner1_inner_start': f'{output_dir}/corner1_inner_start.txt',
            'corner2_inner_end': f'{output_dir}/corner2_inner_end.txt',
            'corner3_outer_end': f'{output_dir}/corner3_outer_end.txt',
            'corner4_outer_start': f'{output_dir}/corner4_outer_start.txt'
        }
        
        for key, filename in filenames.items():
            coords = corners[key]
            
            # Write XYZ format for SolidWorks with "in" suffix
            with open(filename, 'w') as f:
                for i in range(len(coords)):
                    # Format: X\tY\tZ with "in" suffix, tab-separated
                    f.write(f"{coords[i,0]:.6f}in\t{-coords[i,2]:.6f}in\t{coords[i,1]:.6f}in\n")
            
            print(f"✓ Saved: {filename}")
        
        # Visualization
        if plot:
            fig = plt.figure(figsize=(16, 6))
            
            # Plot 1: 3D view of channel corners
            ax1 = fig.add_subplot(121, projection='3d')
            
            ax1.plot(corner1[:,0], corner1[:,1], corner1[:,2], 'r-', linewidth=2, label='Corner 1 (Inner-Start)')
            ax1.plot(corner2[:,0], corner2[:,1], corner2[:,2], 'b-', linewidth=2, label='Corner 2 (Inner-End)')
            ax1.plot(corner3[:,0], corner3[:,1], corner3[:,2], 'g-', linewidth=2, label='Corner 3 (Outer-End)')
            ax1.plot(corner4[:,0], corner4[:,1], corner4[:,2], 'orange', linewidth=2, label='Corner 4 (Outer-Start)')
            
            # Draw connecting lines at a few stations to show rectangle
            indices = [0, len(z_contour)//4, len(z_contour)//2, 3*len(z_contour)//4, -1]
            for idx in indices:
                # Draw rectangle at this axial station
                rect_x = [corner1[idx,0], corner2[idx,0], corner3[idx,0], corner4[idx,0], corner1[idx,0]]
                rect_y = [corner1[idx,1], corner2[idx,1], corner3[idx,1], corner4[idx,1], corner1[idx,1]]
                rect_z = [corner1[idx,2], corner2[idx,2], corner3[idx,2], corner4[idx,2], corner1[idx,2]]
                ax1.plot(rect_x, rect_y, rect_z, 'k-', alpha=0.3, linewidth=1)
            
            ax1.set_xlabel('X (inches)', fontsize=10)
            ax1.set_ylabel('Y (inches)', fontsize=10)
            ax1.set_zlabel('Z (inches)', fontsize=10)
            ax1.set_title('Single Rectangular Channel - 4 Corners', fontsize=12, fontweight='bold')
            ax1.legend(loc='best', fontsize=8)
            ax1.grid(True, alpha=0.3)
            
            # Plot 2: 2D projection showing channel cross-section
            ax2 = fig.add_subplot(122)
            
            # Plot profile view (R vs Z)
            ax2.plot(z_contour, r_inner, 'r-', linewidth=2, label='Inner radius (hot wall)')
            ax2.plot(z_contour, r_outer, 'b-', linewidth=2, label='Outer radius (channel outer)')
            ax2.fill_between(z_contour, r_inner, r_outer, alpha=0.3, color='lightblue', label='Channel height')
            
            ax2.set_xlabel('Axial Position Z (inches)', fontsize=10)
            ax2.set_ylabel('Radius R (inches)', fontsize=10)
            ax2.set_title('Channel Profile View', fontsize=12, fontweight='bold')
            ax2.legend(loc='best')
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f'{output_dir}/channel_corners_visualization.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # Print summary
        print(f"\n{'='*70}")
        print("Rectangular Channel Corner Curves Generated")
        print(f"{'='*70}")
        print(f"Number of axial points: {len(z_contour)}")
        print(f"Channel height (radial): {C_h:.4f} inches")
        print(f"Channel width (circumferential): {C_w:.4f} inches")
        print(f"Hot wall thickness: {t_w:.4f} inches")
        print(f"Fin/rib width: {f_w:.4f} inches")
        print(f"Number of channels: {self.N}")
        print(f"\nAxial extent: Z = {z_contour[0]:.2f} to {z_contour[-1]:.2f} inches")
        print(f"Radial extent: R = {r_inner[0]:.2f} to {r_outer[-1]:.2f} inches")
        print(f"\nOutput files saved to: {output_dir}/")
        print(f"{'='*70}\n")
        
        print("SolidWorks Import Instructions:")
        print("=" * 70)
        print("1. Insert → Curve → Curve Through XYZ Points")
        print("2. Import all 4 corner files:")
        print("   - corner1_inner_start.txt")
        print("   - corner2_inner_end.txt")
        print("   - corner3_outer_end.txt")
        print("   - corner4_outer_start.txt")
        print("\n3. Create channel surface:")
        print("   - Use Loft or Boundary Surface between the 4 corner curves")
        print("   - Or create 4 ruled surfaces and join")
        print("\n4. Pattern the channel:")
        print(f"   - Circular Pattern: {self.N} instances around Z-axis")
        print(f"   - This creates all {self.N} channels")
        print("\n5. Create solid bodies from surfaces as needed")
        print("=" * 70)
        
        return corners

    def outputs(
        self,
        show_regen_temps=True,
        show_cold_temps=True,
        show_dp=True,
        show_pressures=True,
        show_qdot=True,
        show_re=True,
        show_coolant_density=True,
        show_coolant_velocity=True,
        show_wall_thermal_gradient=True,
        show_tangential_stresses=True,
        show_longitudinal_stresses=True,
        show_htc=True,
    ):
        #### ---- PLOTS ---- ####
        pa_to_ksi = 1.0 / 6894757.293168
        boiling_limit = self.coolant_boiling
        coolant_rho_arr = np.array([
            PropsSI("D", "T", float(T_val), "P", float(P_val), self.coolantName)
            for T_val, P_val in zip(self.T_c_arr, self.P_c_arr)
        ])

        if show_regen_temps or show_cold_temps:
            fig, axs = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
            if show_regen_temps:
                axs[0].set_title("Regen Circuit Temperatures")
                axs[0].set_ylabel("Temperature [K]")
                axs[0].plot(self.engine.Contour_z, self.engine.T, color='orange', label='Combustion Gas Temp')
                axs[0].plot(self.engine.Contour_z, self.T_hw_arr, color='r', label='Hot Wall Temp')
                axs[0].plot(self.engine.Contour_z, self.T_cw_arr, color='b', label='Cold Wall Temp')
                axs[0].plot(self.engine.Contour_z, self.T_c_arr, color='y', label='Bulk Coolant Temp')
                axs[0].axhline(boiling_limit, color='k', linestyle='--', label=f'{self.coolantName} Boiling Temp')
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_cold_temps:
                axs[1].set_title("Regen Circuit Cold Side Temperatures")
                axs[1].set_xlabel("Axial Position [m]")
                axs[1].set_ylabel("Temperature [K]")
                axs[1].plot(self.engine.Contour_z, self.T_cw_arr, color='b', label='Cold Wall Temp')
                axs[1].plot(self.engine.Contour_z, self.T_c_arr, color='y', label='Bulk Coolant Temp')
                axs[1].axhline(boiling_limit, color='k', linestyle='--', label=f'{self.coolantName} Boiling Temp')
                axs[1].legend()
            else:
                axs[1].set_visible(False)
            fig.tight_layout()

        if show_dp or show_pressures:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            if show_dp:
                axs[0].plot(self.engine.Contour_z, self.DP_arr / 6894.7, label='Pressure Drop')
                axs[0].set_title("Pressure Drop along Engine")
                axs[0].set_ylabel("DP [psia]")
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_pressures:
                axs[1].plot(self.engine.Contour_z, self.engine.P, label='Chamber Static Pressure')
                axs[1].plot(self.engine.Contour_z, self.P_c_arr / 6894.7, label='Coolant Pressure')
                axs[1].set_title("Pressures along Engine")
                axs[1].set_xlabel("Axial Position [m]")
                axs[1].set_ylabel("Pressure [psia]")
                axs[1].legend()
            else:
                axs[1].set_visible(False)
            fig.tight_layout()

        if show_qdot or show_re or show_coolant_density or show_coolant_velocity:
            fig, axs = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
            if show_qdot:
                axs[0, 0].plot(self.engine.Contour_z, self.Qdot_arr, label='Qdot')
                axs[0, 0].set_title("Radial Heat Flow")
                axs[0, 0].set_ylabel("Qdot [J]")
                axs[0, 0].legend()
            else:
                axs[0, 0].set_visible(False)

            if show_re:
                axs[0, 1].plot(self.engine.Contour_z, self.Re, label='Reynolds Number')
                axs[0, 1].set_title("Coolant Re")
                axs[0, 1].set_ylabel("Re [---]")
                axs[0, 1].legend()
            else:
                axs[0, 1].set_visible(False)

            if show_coolant_density:
                axs[1, 0].plot(self.engine.Contour_z, coolant_rho_arr, label='Coolant Density')
                axs[1, 0].set_title("Coolant Density")
                axs[1, 0].set_xlabel("Axial Position [m]")
                axs[1, 0].set_ylabel("Density [kg/m^3]")
                axs[1, 0].legend()
            else:
                axs[1, 0].set_visible(False)

            if show_coolant_velocity:
                axs[1, 1].plot(self.engine.Contour_z, self.coolant_vel_arr, label='Coolant Velocity')
                axs[1, 1].set_title("Coolant Velocity")
                axs[1, 1].set_xlabel("Axial Position [m]")
                axs[1, 1].set_ylabel("Velocity [m/s]")
                axs[1, 1].legend()
            else:
                axs[1, 1].set_visible(False)
            fig.tight_layout()

        if show_wall_thermal_gradient or show_tangential_stresses or show_longitudinal_stresses:
            fig, axs = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
            if show_wall_thermal_gradient:
                axs[0].plot(self.engine.Contour_z, self.T_hw_arr - self.T_cw_arr, label='Wall Thermal Gradient')
                axs[0].set_title("Wall Thermal Gradient")
                axs[0].set_ylabel("T_hw - T_cw [K]")
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_tangential_stresses:
                axs[1].plot(self.engine.Contour_z, self.sigma_t * pa_to_ksi, label='Total Tangential Stress')
                axs[1].plot(self.engine.Contour_z, self.sigma_pt * pa_to_ksi, label='Pressure Tangential Stress')
                axs[1].plot(self.engine.Contour_z, self.sigma_tt * pa_to_ksi, label='Thermal Tangential Stress')
                axs[1].set_title("Tangential Stresses")
                axs[1].set_ylabel("Stress [ksi]")
                axs[1].legend()
            else:
                axs[1].set_visible(False)

            if show_longitudinal_stresses:
                axs[2].plot(self.engine.Contour_z, self.sigma_l * pa_to_ksi, label='Total Longitudinal Stress')
                axs[2].plot(self.engine.Contour_z, self.sigma_pl * pa_to_ksi, label='Thermal Longitudinal Stress')
                axs[2].set_title("Longitudinal Stresses")
                axs[2].set_xlabel("Axial Position [m]")
                axs[2].set_ylabel("Stress [ksi]")
                axs[2].legend()
            else:
                axs[2].set_visible(False)
            fig.tight_layout()

        if show_htc:
            print(f"Hg Length: {len(self.h_hg_arr)}")
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(self.engine.Contour_z, self.h_hg_arr, label='Hot Gas HTC')
            ax.plot(self.engine.Contour_z, self.h_c_arr, label='Coolant HTC')
            ax.set_title("HTCs along Engine")
            ax.set_xlabel("Axial Position [m]")
            ax.set_ylabel("HTC [W/m^2-K]")
            ax.legend()
            fig.tight_layout()

        # Unit conversions and calculations prior to printing
        self.channel_DP = (self.P_c_arr[-1] - self.P_c_arr[0]) / 6894.76 # Pa
        boiling_temp_margin = self.coolant_boiling - max(self.T_c_arr)
        coolant = Fluid("T", self.T_c_arr[-1], "P", self.T_c_arr[-1], self.coolantName)
        rho1 = coolant.rho
        coolant = Fluid("T", self.T_c_arr[0], "P", self.T_c_arr[0], self.coolantName)
        rho2 = coolant.rho
        avg_rho = (rho1 + rho2) / 2
        self.chan_CdA = (self.chan_mdot / (np.sqrt(2 * avg_rho * (self.P_c_arr[-1] - self.P_c_arr[0])))) * 1000000

        #### ---- SCALAR OUTPUTS ---- ####
        print("\n#### ---- REGEN CIRCUIT OUTPUTS ---- ####")
        print(f"Channel CdA [mm^2]: {round(self.chan_CdA, 3)}")
        print(f"Channel Coolant Mdot [kg/s]: {round(self.chan_mdot, 4)}")
        print(f"Total Coolant Mdot [kg/s]: {round(self.tot_coolant_mdot, 3)}")
        print(f"Channel DP [psid]: {round(self.channel_DP, 2)}")
        print(f"# of Channels: {self.N}")
        print(f"Channel Width [mm]: {round(self.C_w * 1000, 2)}")
        print(f"Channel Height [mm]: {round(self.C_h * 1000, 2)}")
        print(f"Channel AR: {round(self.ch_AR, 2)}")
        print(f"Max Hot Wall Temp [K]: {round(max(self.T_hw_arr), 2)}")
        print(f"Coolant Outlet Temp [K]: {round(max(self.T_c_arr), 2)}")
        print(f"Coolant Boiling Temp Margin [K]: {round(boiling_temp_margin, 2)}")

        plt.show()
        
    

