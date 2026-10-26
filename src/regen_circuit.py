import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
from pathlib import Path

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

from scipy.optimize import fsolve, root_scalar 
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

    # Channel geometry is defined at these four stations, listed in the order they occur along the
    # engine axis (injector end -> nozzle end). Coolant flows from "inlet" to "outlet".
    STATION_NAMES = ("outlet", "chamber", "throat", "inlet")

    def __init__(
        self,
        t_w,
        material,
        tot_coolant_mdot,
        coolantName,
        engine,
        heat_flux_knockdown,
        N,
        stations = None,
        C_w = None,
        C_h = None,
        tbc_thickness = 0.0,
        tbc_k = None,
        min_fin_width = 0.0,
    ):
        '''
        Regen circuit with channel geometry fully defined by the user.

        Inputs:
            t_w : hot wall (liner) thickness [m]
            N : number of channels [-]
            stations : dict with the four keys "outlet", "chamber", "throat", "inlet".
                Each value is a dict {"C_w": [m], "C_h": [m], "z": [in] (optional)}.
                    - "inlet"   : where coolant enters the circuit (nozzle end)
                    - "chamber" : channel dimensions along the chamber
                    - "throat"  : channel dimensions at the throat
                    - "outlet"  : where coolant leaves the circuit (injector end)
                "z" is an axial position in the same coordinates as engine.Contour_z [in]. If omitted
                it defaults to the contour start (outlet), end of the cylindrical chamber (chamber),
                the throat (throat) and the contour end (inlet). Channel width/height are linearly
                interpolated in z between consecutive stations.
            C_w, C_h : legacy shortcut - constant width/height [m] used when `stations` is None.
            min_fin_width : smallest rib width allowed anywhere on the circuit [m].
        '''
        self.t_w = t_w
        self.material = material
        self.tot_coolant_mdot = tot_coolant_mdot
        self.N = int(round(float(N)))
        if self.N < 1:
            raise ValueError(f"N must be a positive integer, got {N}.")
        self.chan_mdot = self.tot_coolant_mdot / self.N
        self.coolantName = coolantName
        self.engine = engine # Engine object from engine.py
        self.heat_flux_knockdown = heat_flux_knockdown
        self.tbc_thickness = float(tbc_thickness)
        self.tbc_k = float(tbc_k) if tbc_k is not None else None
        self.min_fin_width = float(min_fin_width)

        # Build the per-station channel geometry and check it against the Contour.
        self.build_channel_geometry(stations=stations, C_w=C_w, C_h=C_h)

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
        T_0,
        T_wg,

    ) -> float:
        '''
        Calculate the Heat Transfer Coefficient for hot gas in the chamber and nozzle from combustion gas properties using the well known Bartz Correlation.

        h = [0.026 / Dt^0.2 * (mu^0.2 Cp / Pr^0.6) * (Pc / c*)^0.8 * (Dt / r_c)^0.1] * (At / A)^0.9 * sigma

        sigma = 1 / { [0.5 (T_wg / T_0)(1 + (gamma-1)/2 M^2) + 0.5]^(0.8 - omega/5) * (1 + (gamma-1)/2 M^2)^(omega/5) }

        T_0 is the gas STAGNATION (chamber) temperature, not the local static temperature. Bartz's
        sigma already contains the (1 + (gamma-1)/2 M^2) term that converts T_0 to static, so passing the
        static temperature here would double count the Mach effect.

        Inputs:
            Cp, Pr, mu : local gas properties [J/kg-K, -, Pa-s]
            Pc : chamber pressure [Pa]
            AreaRatio : local A / At [-]
            gamma, M : local gas gamma and Mach number
            T_0 : gas stagnation temperature [K]
            T_wg : gas-side wall temperature [K]

        Outputs:
            h_hg [W/m^2-K], sigma [-]
        '''

        omega = 0.6 # From Bartz Paper


        h_gA = 0.026 / (engine.Dt**0.2)
        h_gB = (Cp * mu**0.2) / (Pr**0.6)   
        h_gC = (Pc / engine.cstar)**0.8
        h_gD = (engine.Dt / engine.r_up_throat)**0.1
        h_gE = (1/AreaRatio)**0.9
        
        sigma_A = (0.5 * (T_wg / T_0) * (1 + 0.5*(gamma - 1)*M**2) + 0.5)**(0.8 - 0.2*omega)
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

        return h_hg, sigma

    def solve_bartz_twg_station(
        self,
        i,
        T_c,
        h_c,
        dz,
        T_wg_bracket,
        residual_tol=1e-2,
        fixed_sigma=None,
        sigma_clamp=None,
        use_recovery_temp=False,
    ):
        """
        Solve for T_wg such that residual = T_wg_input - T_tbc_surface = 0 at one station.
        When tbc_thickness == 0, T_tbc_surface == T_hw so this reduces to the original formulation.
        """
        if T_wg_bracket is None or len(T_wg_bracket) != 2:
            raise ValueError("T_wg_bracket must be a 2-element tuple/list: (T_wg_min, T_wg_max).")

        twg_min = float(T_wg_bracket[0])
        twg_max = float(T_wg_bracket[1])

        if sigma_clamp is not None:
            sigma_min = float(sigma_clamp[0])
            sigma_max = float(sigma_clamp[1])
            if sigma_min > sigma_max:
                raise ValueError("sigma_clamp must be (min_sigma, max_sigma).")

        def _resolve_sigma_and_h(
            T_wg_trial,
        ):
            h_hg_model, sigma_model = self.bartz(
                engine=self.engine,
                Cp=self.engine.Cp[i],
                Pr=self.engine.Pr[i],
                mu=self.engine.mu[i],
                Pc=self.engine.Pc,
                AreaRatio=self.engine.AreaRatio[i],
                gamma=self.engine.gamma[i],
                M=self.engine.M[i],
                T_0=self.engine.Tcomb,
                T_wg=float(T_wg_trial),
            )

            if fixed_sigma is None:
                sigma_use = sigma_model
            else:
                sigma_use = float(fixed_sigma)

            # Optional clamp to constrain solution space.
            if sigma_clamp is not None:
                sigma_use = min(max(sigma_use, sigma_min), sigma_max)

            # Bartz predicts h_hg = (other_factors) * sigma.
            # If sigma_use differs from sigma_model, scale h_hg accordingly.
            if abs(sigma_model) < 1e-12 and abs(sigma_use - sigma_model) > 0:
                raise ValueError(
                    f"Bartz sigma ~0 at station {i} for T_wg_trial={T_wg_trial}; "
                    "cannot safely scale h_hg."
                )

            if sigma_use == sigma_model:
                h_hg_use = h_hg_model
            else:
                h_hg_use = h_hg_model * (sigma_use / sigma_model)

            return h_hg_use, sigma_use, sigma_model

        # Gas temperature driving the heat flux: local static temperature by default, or the adiabatic
        # wall (recovery) temperature T_aw = T_0 (1 + r (g-1)/2 M^2) / (1 + (g-1)/2 M^2), r = Pr^(1/3).
        if use_recovery_temp:
            _g = self.engine.gamma[i]
            _M = self.engine.M[i]
            _r = self.engine.Pr[i] ** (1.0 / 3.0)
            T_drive = self.engine.Tcomb * (1 + _r * 0.5 * (_g - 1) * _M**2) / (1 + 0.5 * (_g - 1) * _M**2)
        else:
            T_drive = self.engine.T[i]

        def _station_eval(T_wg_trial):
            h_hg_use, sigma_use, _ = _resolve_sigma_and_h(T_wg_trial)

            T_hg_out, T_c_out, T_cw_out, T_hw_out, qdot_out, A_ch_eff_out, T_tbc_out, eta_f_out = self.thermal_network(
                r_eng=self.engine.Contour_r[i] * 0.0254,
                t_w=self.t_w,
                c_w=self.C_w_contour[i],
                C_h=self.C_h_contour[i],
                f_w=self.f_w_contour[i],
                dz=dz,
                material=self.material,
                h_hg=h_hg_use,
                h_c=h_c,
                T_hg=T_drive,
                T_c=T_c,
            )

            # Residual: T_wg trial must equal the TBC outer surface temperature.
            # When tbc_thickness == 0, T_tbc_out == T_hw_out.
            residual = float(T_wg_trial) - float(T_tbc_out)
            return residual, h_hg_use, sigma_use, T_hg_out, T_c_out, T_cw_out, T_hw_out, qdot_out, A_ch_eff_out, T_tbc_out, eta_f_out

        # --- Bracket endpoints only, then brentq on the sub-interval with a sign change ---
        twg_grid = np.linspace(twg_min, twg_max, 2)
        f_vals = []
        for twg_trial in twg_grid:
            f_vals.append(_station_eval(float(twg_trial))[0])

        # Identify sign-change intervals
        intervals = []
        anchor = 0.5 * (twg_min + twg_max)
        for j in range(len(twg_grid) - 1):
            f1 = f_vals[j]
            f2 = f_vals[j + 1]
            if abs(f1) <= residual_tol:
                # Already within tolerance at sample point.
                return float(twg_grid[j]), *_station_eval(float(twg_grid[j]))[1:]
            if f1 * f2 < 0:
                a = float(twg_grid[j])
                b = float(twg_grid[j + 1])
                mid = 0.5 * (a + b)
                intervals.append((a, b, mid))

        if not intervals:
            raise ValueError(
                f"No sign change found during T_wg scan at station {i} "
                f"in bracket [{twg_min:.3f}, {twg_max:.3f}]. "
                "Try expanding/adjusting the bracket or clamp sigma."
            )

        # Choose the sub-interval closest to bracket midpoint.
        a_sel, b_sel, _ = min(intervals, key=lambda x: abs(x[2] - anchor))

        sol = root_scalar(
            lambda t: _station_eval(t)[0],
            bracket=(a_sel, b_sel),
            method="brentq",
            xtol=residual_tol,
        )

        if not sol.converged:
            raise ValueError(f"T_wg root solve did not converge at station {i}.")

        station = _station_eval(sol.root)
        residual = station[0]
        if abs(residual) > residual_tol:
            raise ValueError(
                f"T_wg root residual too high at station {i}: |residual|={abs(residual):.6f} > {residual_tol}."
            )

        _, h_hg_use, sigma_use, T_hg_out, T_c_out, T_cw_out, T_hw_out, qdot_out, A_ch_eff_out, T_tbc_out, eta_f_out = station
        return sol.root, h_hg_use, sigma_use, T_hg_out, T_c_out, T_cw_out, T_hw_out, qdot_out, A_ch_eff_out, T_tbc_out, eta_f_out

    @staticmethod
    def laminar_rect_duct_nu(C_w, C_h):
        """
        Fully developed laminar Nusselt number for a rectangular duct, constant heat flux (H1).
        Shah & London polynomial fit in aspect ratio alpha = min(C_w, C_h) / max(C_w, C_h).
        """
        a = min(C_w, C_h) / max(C_w, C_h)
        return 8.235 * (1 - 2.0421*a + 3.0853*a**2 - 2.4765*a**3 + 1.0578*a**4 - 0.1861*a**5)

    def coolant_htc(
        self,
        fluid,
        C_w,
        C_h,
        correlation="Gnielinski",
        eps=None,
        debug=False,
        enforce_validity=True,
    ):
        """
        Calculate coolant-side HTC at one station, for a channel of width C_w and height C_h [m].

        The default (and intended) correlation is Gnielinski, using the Darcy friction factor with wall
        roughness from calc_friction_factor. Gnielinski is valid for Re > 2300; with
        enforce_validity=True, below that the fully developed laminar rectangular-duct Nusselt number
        is used instead (and Gnielinski is floored at that value so Nu never drops below laminar).
        enforce_validity=False applies the Gnielinski expression at every Re (extrapolating it).

        Parameters
        ----------
        correlation : str
            "Gnielinski" or "Dittus-Boelter" (Dittus-Boelter is also always returned for comparison)
        """
        P = 2 * (C_w + C_h)             # Channel wetted perimeter
        A_c = C_w * C_h                 # Channel flow area
        D_h = (4 * A_c) / P             # Hydraulic diameter
        Re = (self.chan_mdot * D_h) / (A_c * fluid.mu)
        Pr = fluid.Pr

        f_d = self.calc_friction_factor(
            D_h=D_h,
            Re=Re,
            eps=eps,
        )

        # Gnielinski
        Nu_lam = self.laminar_rect_duct_nu(C_w, C_h)
        if enforce_validity and Re < 2300.0:
            Nu_gn = Nu_lam
        else:
            Nu_gn = ((f_d / 8) * (Re - 1000) * Pr) / (1 + 12.7 * ((f_d / 8) ** 0.5) * (Pr ** (2 / 3) - 1))
            if enforce_validity:
                Nu_gn = max(Nu_gn, Nu_lam)
        h_gn = (Nu_gn * fluid.k) / D_h

        # Dittus-Boelter
        Nu_db = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
        h_db = (Nu_db * fluid.k) / D_h

        corr_key = str(correlation).strip().lower()
        if corr_key in ("gnielinski", "gn"):
            h_sel = h_gn
        elif corr_key in ("dittus-boelter", "dittus_boelter", "db"):
            h_sel = h_db
        else:
            raise ValueError(
                f"Unknown coolant HTC correlation '{correlation}'. "
                "Use 'Gnielinski' or 'Dittus-Boelter'."
            )

        if debug:
            print(f"Coolant HTC correlation : {correlation}")
            print(f"Hydraulic Diameter : {D_h} m")
            print(f"Re : {Re}")
            print(f"Pr : {Pr}")
            print(f"f_d : {f_d}")
            print(f"Nu_Gnielinski : {Nu_gn}")
            print(f"Nu_DittusBoelter : {Nu_db}")
            print(f"h_c_Gnielinski : {h_gn} W/m^2-K")
            print(f"h_c_DittusBoelter : {h_db} W/m^2-K")
            print(f"h_c_Selected : {h_sel} W/m^2-K")

        return h_sel, h_gn, h_db
    
    def calculate_stress(
        self,
        T_hw,
        T_cw,
        material,
        qdot,
        P_c,
        P_hg,
        C_w,
    ):
        '''
        Quick per-station stress estimate used inside the circuit march (C_w = local channel width [m]).
        The full stress/lifecycle analysis lives in structural.py.
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
        sigma_pt = ((P_c - P_hg) / 2) * (C_w / self.t_w)**2
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
        C_h,
        f_w,
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
        m = math.sqrt((h_c * (2 * dz)) / (material.k * f_w * dz))
        eta_f = np.tanh(m * C_h) / (m * C_h)

        # Multi-Fin #
        A_b = self.N * c_w * dz
        A_hwall_ch = (2 * np.pi * r_eng * dz) / self.N
        A_f = 2 * self.N * C_h * dz
        A_t = A_b + A_f
        eta_o = (A_b + eta_f * A_f) / A_t
        A_ch_eff = (eta_o * A_t) / self.N

        # Thermal Resistances
        R_hghw = 1 / (self.heat_flux_knockdown * h_hg * A_hwall_ch)
        R_tbc = 0.0
        if self.tbc_thickness > 0.0 and self.tbc_k is not None:
            R_tbc = self.tbc_thickness / (self.tbc_k * A_hwall_ch)
        R_wall = t_w / (material.k * A_hwall_ch)
        R_cwc = 1 / (h_c * A_ch_eff)
        R_tot = R_hghw + R_tbc + R_wall + R_cwc

        # Fast path: only T_hg and T_c specified (regen station solve). Series network is linear;
        # SymPy is O(100–1000×) slower per call than closed form.
        if (
            qdot_new is None
            and T_hw_new is None
            and T_cw_new is None
            and T_hg_new is not None
            and T_c_new is not None
        ):
            qdot_val = (T_hg_new - T_c_new) / R_tot
            T_tbc_val = T_hg_new - qdot_val * R_hghw           # TBC outer surface (gas-side)
            T_hw_val = T_tbc_val - qdot_val * R_tbc             # metal hot wall under coating
            T_cw_val = T_hw_val - qdot_val * R_wall
            return T_hg_new, T_c_new, T_cw_val, T_hw_val, qdot_val, A_ch_eff, T_tbc_val, eta_f

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
        eq2 = sp.Eq(qdot, (T_hg - T_hw) / (R_hghw + R_tbc))
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

        T_tbc_val = float(T_hg) - float(qdot) * R_hghw

        return T_hg, T_c, T_cw, T_hw, qdot, A_ch_eff, T_tbc_val, eta_f
    
    def solve_circuit(
        self,
        inlet_T_c,
        inlet_pressure,
        T_wg_bracket=(300.0, 4500.0), # [K]
        T_wg_residual_tol=1e-2,
        fixed_sigma=None,
        sigma_clamp=None,
        coolant_htc_correlation="Gnielinski",
        coolant_temp_limit=None,
        show_progress=True,
        print_station_rows=False,
        use_recovery_temp=False,
        enforce_gnielinski_validity=True,
    ):
        # Circuit extent comes from the "inlet" and "outlet" station positions (see build_channel_geometry).
        circuit_inlet = self.circuit_inlet
        circuit_outlet = self.circuit_outlet

        # can input a min num of Pts, then find velocity at max velocity in channel, then set timestep to dz / max velocity
        #dz = round(self.Contour.z[-1] / self.numPts, 0)
        time_tracker = 0

        # Calculate dz - use the arc length based on the specified dz NOT the dz since that will underestimate volume
        # will lead to underconservative heating
        dz = np.abs(np.diff(self.engine.Contour_z * 0.0254))
        dr = np.diff(self.engine.Contour_r * 0.0254)
        ds = np.sqrt(dz**2 + dr**2)

        # repeat end value so arrays have same shapes
        ds_end = ds[-1]
        ds = np.append(ds, ds_end)
        self.regen_curve_length_m = float(np.sum(ds[circuit_outlet:circuit_inlet + 1]))

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
        self.P_c_in_arr = np.empty(n_stations)
        self.h_hg_arr = np.empty(n_stations)
        self.sigma = np.empty(n_stations)
        self.h_c_arr = np.empty(n_stations)
        self.h_c_gn_arr = np.empty(n_stations)
        self.h_c_db_arr = np.empty(n_stations)
        self.DP_arr = np.empty(n_stations)
        self.f_D_arr = np.empty(n_stations)
        self.Re = np.empty(n_stations)
        self.Pr_arr = np.empty(n_stations)
        self.sigma_pt = np.empty(n_stations)
        self.sigma_tt = np.empty(n_stations)
        self.sigma_t = np.empty(n_stations)
        self.sigma_pl = np.empty(n_stations)
        self.sigma_l = np.empty(n_stations)
        self.coolant_boiling_arr = np.empty(n_stations)
        self.T_tbc_arr = np.empty(n_stations)
        self.q_flux_arr = np.empty(n_stations)
        self.eta_f_arr = np.empty(n_stations)

        # T and P Tracker
        P_c = inlet_pressure
        T_c = inlet_T_c
        
        # Inlet reference limit. If user provides coolant_temp_limit, use it everywhere.
        if coolant_temp_limit is None:
            self.coolant_boiling_inlet = self.calculate_boiling_point(P_c, T_c, self.coolantName, plot_cp_curve=False)
            print(f" Coolant BP (inlet P): {round(self.coolant_boiling_inlet, 2)} K")
        else:
            self.coolant_boiling_inlet = float(coolant_temp_limit)
            print(f" Coolant temperature limit (user): {round(self.coolant_boiling_inlet, 2)} K")

        station_range = range(circuit_inlet, circuit_outlet - 1, -1)
        for i in tqdm(station_range, disable=(not show_progress)):
            T_c_in = T_c
            P_c_in = P_c
            coolant = Fluid("T", T_c, "P", P_c, self.coolantName)

            C_w_i = self.C_w_contour[i]
            C_h_i = self.C_h_contour[i]

            coolantVolume = C_w_i * C_h_i * ds[i]
            m = coolantVolume * coolant.rho

            h_c, h_c_gn, h_c_db = self.coolant_htc(
                fluid=coolant,
                C_w=C_w_i,
                C_h=C_h_i,
                correlation=coolant_htc_correlation,
                debug=False,
                enforce_validity=enforce_gnielinski_validity,
            )

            T_wg_sol, h_hg, sigma, T_hg, T_c, T_cw, T_hw, qdot, A_ch_eff, T_tbc, eta_f = self.solve_bartz_twg_station(
                i=i,
                T_c=T_c,
                h_c=h_c,
                dz=ds[i],
                T_wg_bracket=T_wg_bracket,
                residual_tol=T_wg_residual_tol,
                fixed_sigma=fixed_sigma,
                sigma_clamp=sigma_clamp,
                use_recovery_temp=use_recovery_temp,
            )
            q_flux = qdot / (ds[i] * ((2 * np.pi * self.engine.Contour_r[i] * 0.0254) / self.N)) # W/m^2
            
            # Stress Calculations
            sigma_pt, sigma_tt, sigma_t, sigma_pl, sigma_l = self.calculate_stress(
                T_hw=T_hw,
                T_cw=T_cw,
                material=self.material,
                qdot=q_flux,
                P_c=P_c,
                P_hg=self.engine.P[i] * 6894.76,
                C_w=C_w_i,
            )
            
            '''
            print("Thermal Network Outputs")
            print(f"qdot: {qdot} W")
            print(f"T_hg : {T_hg} K")
            print(f"T_hw : {T_hw} K")
            print(f"T_cw : {T_cw} K")
            print(f"T_c : {T_c} K")
            '''

            coolant_vel = self.chan_mdot / (coolant.rho * C_w_i * C_h_i)
            D_h = (4*(C_w_i * C_h_i)) / (2*(C_w_i + C_h_i))
            Re = (coolant.rho * coolant_vel * D_h) / coolant.mu
            f_D = self.calc_friction_factor(D_h, Re)
            DP = self.calc_DP(coolant, coolant_vel, ds[i], D_h, Re, f_D=f_D)
            P_c_new = P_c - DP
            if P_c_new <= 0.0:
                raise ValueError(
                    "Coolant pressure became non-physical (<= 0 Pa) during pressure-drop integration. "
                    f"station={i}, P_old={P_c/6894.76:.2f} psia, dP={DP/6894.76:.2f} psid, "
                    f"P_new={P_c_new/6894.76:.2f} psia, Re={Re:.0f}, f_D={f_D:.4f}, "
                    f"vel={coolant_vel:.2f} m/s, Dh={D_h*1000:.3f} mm, L={ds[i]*1000:.3f} mm. "
                    "Reduce per-channel mass flow (increase N) or increase channel area (C_w/C_h)."
                )
            P_c = P_c_new

            dt = ds[i] / coolant_vel
            T_c_new = ((qdot * dt) / (m * coolant.Cp)) + T_c
            T_c = T_c_new

            # Boiling limit at local (P, T): T_sat(P) drops as pressure drops along the channel.
            if coolant_temp_limit is None:
                coolant_boiling_local = self.calculate_boiling_point(P_c, T_c, self.coolantName, plot_cp_curve=False)
            else:
                coolant_boiling_local = float(coolant_temp_limit)
            if T_c > coolant_boiling_local:
                raise ValueError(
                    f" === Coolant above local boiling limit at z index {i}: "
                    f"T_c={round(T_c, 2)} K vs limit {round(coolant_boiling_local, 2)} K at P={round(P_c / 6894.76, 2)} psia!"
                )

            # Indexed update (faster than repeated np.append in loop).
            out_idx = i - circuit_outlet
            self.T_hw_arr[out_idx] = T_hw
            self.T_cw_arr[out_idx] = T_cw
            self.T_c_arr[out_idx] = T_c_new
            self.Qdot_arr[out_idx] = qdot
            self.coolant_vel_arr[out_idx] = coolant_vel
            self.P_c_arr[out_idx] = P_c_new
            self.P_c_in_arr[out_idx] = P_c_in
            self.h_hg_arr[out_idx] = h_hg
            self.sigma[out_idx] = sigma
            self.h_c_arr[out_idx] = h_c
            self.h_c_gn_arr[out_idx] = h_c_gn
            self.h_c_db_arr[out_idx] = h_c_db
            self.DP_arr[out_idx] = DP
            self.f_D_arr[out_idx] = f_D
            self.Re[out_idx] = Re
            self.Pr_arr[out_idx] = coolant.Pr
            self.sigma_pt[out_idx] = sigma_pt
            self.sigma_tt[out_idx] = sigma_tt
            self.sigma_t[out_idx] = sigma_t
            self.sigma_pl[out_idx] = sigma_pl
            self.sigma_l[out_idx] = sigma_l
            self.coolant_boiling_arr[out_idx] = coolant_boiling_local
            self.T_tbc_arr[out_idx] = T_tbc
            self.q_flux_arr[out_idx] = q_flux
            self.eta_f_arr[out_idx] = eta_f

            if print_station_rows:
                z_in = float(self.engine.Contour_z[i])
                print(f"index : {i}  z : {z_in:.4f} in")
                print(f"T_c in : {T_c_in:.2f} K")
                print(f"T_c out : {T_c_new:.2f} K")
                print(f"P_c : {P_c_new / 6894.76:.2f} psia")
                print(f"T_hg : {T_hg:.2f} K")
                print(f"T_tbc : {T_tbc:.2f} K")
                print(f"T_hw : {T_hw:.2f} K")
                print(f"T_cw : {T_cw:.2f} K")
                print(f"T_wg : {T_wg_sol:.2f} K")
                print(f"h_hg : {h_hg:.0f} W/m^2-K")
                print(f"h_c : {h_c:.0f} W/m^2-K")
                print(f"h_c_Gnielinski : {h_c_gn:.0f} W/m^2-K")
                print(f"h_c_DittusBoelter : {h_c_db:.0f} W/m^2-K")
                print(f"qdot : {qdot:.4g} W")
                print(f"coolant_vel : {coolant_vel:.3f} m/s")
                print(f"Re : {Re:.1f}")
                print(f"f_D : {f_D:.4f}")
                print(f"Pr : {coolant.Pr:.3f}")
                print(f"dp : {DP / 6894.76:.3f} psid")
                print(f"sigma : {sigma:.4f}")
                print()

        n_low_re = int(np.sum(self.Re < 2300.0))
        if n_low_re > 0:
            print(
                f" WARNING: coolant Re < 2300 (laminar/transitional) at {n_low_re} of {n_stations} stations "
                f"(min Re = {float(np.min(self.Re)):.0f}); Gnielinski is not valid there"
                + (" - laminar Nu used." if enforce_gnielinski_validity else " - extrapolated as requested.")
                + " Consider more flow per channel (fewer/smaller channels)."
            )

        # Minimum local boiling limit along the channel (tightest vs bulk T_c for scalar summaries).
        self.coolant_boiling = float(np.min(self.coolant_boiling_arr))

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
        eps = None,
        f_D = None,
    ):
        if f_D is None:
            f_D = self.calc_friction_factor(D_h, Re, eps)
        #print(f"Friction Factor : {f_D}")
        #print(f"D_h : {D_h} m")
        #print(f"Rho : {fluid.rho} kg/m^3")
        dP = 0.5 * fluid.rho * velocity**2 * (L / D_h) * f_D # Pa

        #print(f"f_d : {f_D}")
        #print(f"Length : {L} m")

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

        if eps is None:
            eps = 0.020 / 1000 # [mm] UPDATE - ALSO UPDATE FOR LPBF SURFACES

        if Re <= 0:
            raise ValueError(f"Re must be > 0 to compute friction factor, got Re={Re}.")

        # Laminar regime (Darcy factor).
        if Re < 2300.0:
            return 64.0 / Re

        # Churchill correlation (Darcy factor): valid across laminar/transitional/turbulent
        # and avoids iterative Colebrook root-finding instability.
        rel_rough = eps / D_h
        A = (2.457 * np.log(1.0 / (((7.0 / Re) ** 0.9) + (0.27 * rel_rough)))) ** 16
        B = (37530.0 / Re) ** 16
        f = 8.0 * (((8.0 / Re) ** 12) + (1.0 / ((A + B) ** 1.5))) ** (1.0 / 12.0)
        return float(f)


    def build_channel_geometry(
        self,
        stations = None,
        C_w = None,
        C_h = None,
    ):
        '''
        Build the channel geometry from user inputs at four stations (outlet, chamber, throat, inlet).

        Channel width and height are linearly interpolated in axial position between consecutive
        stations. Rib (fin) width follows from the local channel pitch:

            f_w(z) = 2 pi (r(z) + t_w) / N - C_w(z)

        Sets (all per-contour-point arrays are indexed like engine.Contour_z, [m]):
            self.stations, self.station_z                 - resolved station inputs
            self.circuit_outlet, self.circuit_inlet       - contour indices of the circuit ends
            self.C_w_contour, self.C_h_contour, self.f_w_contour
            self.circuit_z                                - Contour_z over the circuit [in]
            self.C_w_arr, self.C_h_arr, self.f_w_arr, self.ch_AR_arr - same, over the circuit only
        '''
        z_all = np.asarray(self.engine.Contour_z, dtype=float)    # [in]
        r_all = np.asarray(self.engine.Contour_r, dtype=float)    # [in]

        # Legacy shortcut: constant channel over the whole contour.
        if stations is None:
            if C_w is None or C_h is None:
                raise ValueError(
                    "Define the channel geometry with `stations` (outlet, chamber, throat, inlet), "
                    "or pass both C_w and C_h for a constant-section channel."
                )
            stations = {name: {"C_w": C_w, "C_h": C_h} for name in self.STATION_NAMES}

        missing = [n for n in self.STATION_NAMES if n not in stations]
        extra = [n for n in stations if n not in self.STATION_NAMES]
        if missing or extra:
            raise ValueError(
                f"stations must have exactly the keys {self.STATION_NAMES}. Missing: {missing}, unexpected: {extra}."
            )

        # Default axial positions for any station that does not specify "z".
        default_z = {
            "outlet": z_all[0],
            "chamber": z_all[self.engine.chamberEnd_ind],
            "throat": z_all[self.engine.throat_ind],
            "inlet": z_all[-1],
        }
        resolved = {}
        for name in self.STATION_NAMES:
            st = stations[name]
            for key in ("C_w", "C_h"):
                if key not in st or st[key] is None or float(st[key]) <= 0.0:
                    raise ValueError(f"stations['{name}']['{key}'] must be a positive value in meters.")
            z_st = st.get("z", None)
            resolved[name] = {
                "z": float(default_z[name] if z_st is None else z_st),
                "C_w": float(st["C_w"]),
                "C_h": float(st["C_h"]),
            }

        station_z = np.array([resolved[n]["z"] for n in self.STATION_NAMES])
        if np.any(np.diff(station_z) <= 0.0):
            raise ValueError(
                "Station z positions must strictly increase outlet < chamber < throat < inlet, got "
                + ", ".join(f"{n}={z:.4f} in" for n, z in zip(self.STATION_NAMES, station_z))
            )
        if station_z[0] < z_all[0] - 1e-9 or station_z[-1] > z_all[-1] + 1e-9:
            raise ValueError(
                f"Outlet/inlet z ({station_z[0]:.4f}, {station_z[-1]:.4f} in) must lie within the engine contour "
                f"({z_all[0]:.4f} to {z_all[-1]:.4f} in)."
            )

        self.stations = resolved
        self.station_z = station_z
        self.circuit_outlet = int(np.argmin(np.abs(z_all - station_z[0])))
        self.circuit_inlet = int(np.argmin(np.abs(z_all - station_z[-1])))

        # Linear interpolation of channel dimensions between stations (flat beyond the end stations).
        cw_st = np.array([resolved[n]["C_w"] for n in self.STATION_NAMES])
        ch_st = np.array([resolved[n]["C_h"] for n in self.STATION_NAMES])
        self.C_w_contour = np.interp(z_all, station_z, cw_st)
        self.C_h_contour = np.interp(z_all, station_z, ch_st)
        self.f_w_contour = (2.0 * np.pi * (r_all * 0.0254 + self.t_w)) / self.N - self.C_w_contour

        o, i = self.circuit_outlet, self.circuit_inlet
        self.circuit_z = z_all[o:i + 1]
        self.C_w_arr = self.C_w_contour[o:i + 1]
        self.C_h_arr = self.C_h_contour[o:i + 1]
        self.f_w_arr = self.f_w_contour[o:i + 1]
        self.ch_AR_arr = self.C_h_arr / self.C_w_arr

        # Rib width must be positive (and above the manufacturing minimum) everywhere on the circuit.
        worst = int(np.argmin(self.f_w_arr))
        if self.f_w_arr[worst] <= max(self.min_fin_width, 0.0):
            raise ValueError(
                f"Fin width {self.f_w_arr[worst] * 1000:.3f} mm at z={self.circuit_z[worst]:.3f} in "
                f"(r={r_all[o + worst]:.3f} in) is below the minimum {max(self.min_fin_width, 0.0) * 1000:.3f} mm. "
                f"Reduce N ({self.N}) or the channel width C_w ({self.C_w_arr[worst] * 1000:.3f} mm) there."
            )

        print(f"Channels: N={self.N}, circuit z = {station_z[0]:.3f} -> {station_z[-1]:.3f} in (coolant flows inlet -> outlet)")
        for n in self.STATION_NAMES:
            st = resolved[n]
            print(f"  {n:8s} z={st['z']:8.3f} in  C_w={st['C_w'] * 1000:6.3f} mm  C_h={st['C_h'] * 1000:6.3f} mm")
        print(f"Min fin width: {self.f_w_arr[worst] * 1000:.3f} mm at z={self.circuit_z[worst]:.3f} in")

    def plot_regen_geometry(
        self,
        engine,
        show_dimensions = True    
    ):
        '''
        Plot regen channel cross-sections at the outlet, throat, and inlet stations
        in a single window using subplots (channel size is taken from the local station geometry).
        '''
        # Unit conversions
        contour_r_mm = np.asarray(engine.Contour_r) * 25.4
        contour_z_in = np.asarray(engine.Contour_z)
        throat_idx = int(np.argmin(contour_r_mm))

        section_indices = [self.circuit_outlet, throat_idx, self.circuit_inlet]
        section_titles = ["Regen Outlet Plane", "Throat Plane", "Regen Inlet Plane"]

        outer_wall_thickness = 5.0 # [mm]
        inner_wall_thickness = self.t_w * 1000.0

        # Common angular spacing for channels/ribs around circumference.
        # Cast SymPy numeric types to native Python floats/ints for NumPy compatibility.
        n_channels = int(round(float(self.N)))
        total_angle_per_channel = float(2.0 * np.pi / n_channels)

        fig, axs = plt.subplots(1, 3, figsize=(18, 6))

        for ax, idx, title in zip(axs, section_indices, section_titles):
            # Engine contour is gas-side wall ID. Build layers radially outward:
            # gas wall ID -> hot wall OD -> channel OD -> outer shell OD.
            channel_height = self.C_h_contour[idx] * 1000.0
            channel_width = self.C_w_contour[idx] * 1000.0
            wall_width = self.f_w_contour[idx] * 1000.0

            gas_wall_inner_radius = contour_r_mm[idx]
            hot_wall_outer_radius = gas_wall_inner_radius + inner_wall_thickness
            channel_outer_radius = hot_wall_outer_radius + channel_height
            outer_radius = channel_outer_radius + outer_wall_thickness

            # Local channel/rib angular widths at this station
            channel_angle = float(channel_width / hot_wall_outer_radius)
            wall_angle = float(wall_width / hot_wall_outer_radius)
            if (channel_angle + wall_angle) > total_angle_per_channel:
                channel_angle = float(total_angle_per_channel - wall_angle)

            theta = np.linspace(0, 2*np.pi, 200)

            # Combustion chamber core
            ax.fill(
                gas_wall_inner_radius * np.cos(theta),
                gas_wall_inner_radius * np.sin(theta),
                color='lightcoral', alpha=0.3
            )

            # Hot wall annulus
            hot_wall_annulus = np.vstack([
                np.column_stack([hot_wall_outer_radius * np.cos(theta), hot_wall_outer_radius * np.sin(theta)]),
                np.column_stack([gas_wall_inner_radius * np.cos(theta[::-1]), gas_wall_inner_radius * np.sin(theta[::-1])])
            ])
            ax.add_patch(Polygon(hot_wall_annulus, facecolor='red', edgecolor='darkred', linewidth=1.5, alpha=0.7))

            # Outer wall annulus
            outer_annulus = np.vstack([
                np.column_stack([outer_radius * np.cos(theta), outer_radius * np.sin(theta)]),
                np.column_stack([channel_outer_radius * np.cos(theta[::-1]), channel_outer_radius * np.sin(theta[::-1])])
            ])
            ax.add_patch(Polygon(outer_annulus, facecolor='darkgray', edgecolor='black', linewidth=1.2, alpha=0.4))

            # Rib region annulus background
            rib_annulus = np.vstack([
                np.column_stack([channel_outer_radius * np.cos(theta), channel_outer_radius * np.sin(theta)]),
                np.column_stack([hot_wall_outer_radius * np.cos(theta[::-1]), hot_wall_outer_radius * np.sin(theta[::-1])])
            ])
            ax.add_patch(Polygon(rib_annulus, facecolor='gray', edgecolor='none', alpha=0.6))

            # Channels
            current_angle = 0.0
            for _ in range(n_channels):
                theta_start = current_angle
                theta_end = current_angle + channel_angle
                theta_channel = np.linspace(theta_start, theta_end, 20)

                x_ch_inner = hot_wall_outer_radius * np.cos(theta_channel)
                y_ch_inner = hot_wall_outer_radius * np.sin(theta_channel)
                x_ch_outer = channel_outer_radius * np.cos(theta_channel[::-1])
                y_ch_outer = channel_outer_radius * np.sin(theta_channel[::-1])

                x_channel = np.concatenate([x_ch_inner, x_ch_outer])
                y_channel = np.concatenate([y_ch_inner, y_ch_outer])
                ax.add_patch(
                    Polygon(
                        np.column_stack([x_channel, y_channel]),
                        facecolor='lightblue', edgecolor='blue', linewidth=0.8, alpha=0.85
                    )
                )

                current_angle = float(current_angle + total_angle_per_channel)

            if show_dimensions:
                ax.text(
                    0.02, 0.98,
                    f"R_gas_ID = {gas_wall_inner_radius:.2f} mm\nOD = {2*outer_radius:.2f} mm\nZ = {contour_z_in[idx]:.3f} in",
                    transform=ax.transAxes,
                    verticalalignment='top',
                    fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )

            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.set_xlabel('X [mm]')
            ax.set_ylabel('Y [mm]')
            ax.set_aspect('equal')
            ax.set_xlim(-outer_radius*1.2, outer_radius*1.2)
            ax.set_ylim(-outer_radius*1.2, outer_radius*1.2)
            ax.grid(True, alpha=0.3)

        fig.suptitle(f"Regenerative Cooling Channel Cross-Sections ({self.N} Channels)", fontsize=14, fontweight='bold')
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
        self.C_h_arr : array
            Channel height (radial) in meters, per station over the circuit
        self.C_w_arr : array
            Channel width (circumferential) in meters, per station over the circuit
        self.t_w : float
            Hot wall thickness in meters
        self.N : int
            Number of channels (for reference)
        self.f_w_arr : array
            Fin/rib width in meters, per station over the circuit
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
        o, inl = self.circuit_outlet, self.circuit_inlet
        z_contour = np.asarray(self.engine.Contour_z).copy()[o:inl + 1]  # inches
        r_contour = np.asarray(self.engine.Contour_r).copy()[o:inl + 1]  # inches
        
        # Convert to inches for local use only (arrays: one value per axial point)
        C_h = self.C_h_arr * 39.3701  # meters to inches (local variable)
        C_w = self.C_w_arr * 39.3701  # meters to inches (local variable)
        t_w = self.t_w * 39.3701  # meters to inches (local variable)
        f_w = self.f_w_arr * 39.3701  # meters to inches (local variable)
        
        # Radial extent: gas-side contour r_contour, then liner t_w, then channel depth C_h.
        # Channel inner (hot side of coolant path) and outer (cold side) — both in inches.
        r_inner = r_contour + t_w
        r_outer = r_inner + C_h  # was wrongly r_contour + C_h (ignored t_w; radial span was C_h - t_w)

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
            theta_end = C_w[i] / r_in
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
        print(f"Channel height (radial): {C_h.min():.4f} to {C_h.max():.4f} inches")
        print(f"Channel width (circumferential): {C_w.min():.4f} to {C_w.max():.4f} inches")
        print(f"Hot wall thickness: {t_w:.4f} inches")
        print(f"Fin/rib width: {f_w.min():.4f} to {f_w.max():.4f} inches")
        print(f"Number of channels: {self.N}")
        print(f"\nAxial extent: Z = {z_contour[0]:.2f} to {z_contour[-1]:.2f} inches")
        print(f"Radial extent: R = {r_inner[0]:.2f} to {r_outer.max():.2f} inches")
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

    def export_performance_to_excel(self, filename="regen_performance_outputs.xlsx"):
        """
        Export per-station regen and gas-side performance outputs to an Excel workbook.
        """
        try:
            from openpyxl import Workbook
        except ImportError as exc:
            raise ImportError(
                "openpyxl is required for Excel export. Install it with: pip install openpyxl"
            ) from exc

        z_in = np.asarray(self.circuit_z, dtype=float)
        z_m = z_in * 0.0254
        coolant_rho_arr = np.array([
            PropsSI("D", "T", float(T_val), "P", float(P_val), self.coolantName)
            for T_val, P_val in zip(self.T_c_arr, self.P_c_arr)
        ])

        wb = Workbook()
        ws = wb.active
        ws.title = "Axial Performance"

        headers = [
            "Axial Position (in)",
            "Axial Position (m)",
            "Gas Temperature (K)",
            "Gas Pressure (psia)",
            "Gas Density (kg/m^3)",
            "Hot Wall Temp (K)",
            "Cold Wall Temp (K)",
            "Coolant Bulk Temp (K)",
            "Coolant Pressure (psia)",
            "Coolant Density (kg/m^3)",
            "Coolant Velocity (m/s)",
            "Qdot (W)",
            "Heat Flux (W/m^2)",
            "Heat Flux (MW/m^2)",
            "Hot Gas HTC (W/m^2-K)",
            "Coolant HTC Selected (W/m^2-K)",
            "Coolant HTC Gnielinski (W/m^2-K)",
            "Coolant HTC Dittus-Boelter (W/m^2-K)",
            "Bartz Sigma (-)",
            "Reynolds Number (-)",
            "Prandtl Number (-)",
            "Darcy Friction Factor (-)",
            "Pressure Drop dP (Pa)",
            "Pressure Drop dP (psid)",
            "Coolant Boiling Limit (K)",
            "TBC Surface Temp (K)",
            "Tangential Stress Total (Pa)",
            "Tangential Stress Pressure (Pa)",
            "Tangential Stress Thermal (Pa)",
            "Longitudinal Stress Total (Pa)",
            "Longitudinal Stress Thermal (Pa)",
        ]
        ws.append(headers)

        o = self.circuit_outlet
        for i in range(len(z_in)):
            ws.append([
                float(z_in[i]),
                float(z_m[i]),
                float(self.engine.T[o + i]),
                float(self.engine.P[o + i]),
                float(self.engine.rho[o + i]),
                float(self.T_hw_arr[i]),
                float(self.T_cw_arr[i]),
                float(self.T_c_arr[i]),
                float(self.P_c_arr[i] / 6894.76),
                float(coolant_rho_arr[i]),
                float(self.coolant_vel_arr[i]),
                float(self.Qdot_arr[i]),
                float(self.q_flux_arr[i]),
                float(self.q_flux_arr[i] / 1e6),
                float(self.h_hg_arr[i]),
                float(self.h_c_arr[i]),
                float(self.h_c_gn_arr[i]),
                float(self.h_c_db_arr[i]),
                float(self.sigma[i]),
                float(self.Re[i]),
                float(self.Pr_arr[i]),
                float(self.f_D_arr[i]),
                float(self.DP_arr[i]),
                float(self.DP_arr[i] / 6894.76),
                float(self.coolant_boiling_arr[i]),
                float(self.T_tbc_arr[i]),
                float(self.sigma_t[i]),
                float(self.sigma_pt[i]),
                float(self.sigma_tt[i]),
                float(self.sigma_l[i]),
                float(self.sigma_pl[i]),
            ])

        output_path = Path(filename)
        if output_path.suffix.lower() != ".xlsx":
            output_path = output_path.with_suffix(".xlsx")
        wb.save(output_path)
        print(f"Saved regen axial performance export: {output_path.resolve()}")

    def run_structural_analysis(self, **kwargs):
        """
        Thermal + mechanical (hoop/pressure bending) stress and Neuber-corrected low-cycle-fatigue life
        along the circuit. Keyword arguments go to structural.LinerStructuralAnalysis (Kt,
        life_safety_factor, required_cycles, include_liner_hoop, t_jacket). Call after solve_circuit.
        """
        from structural import LinerStructuralAnalysis
        self.structural = LinerStructuralAnalysis(self, **kwargs).solve()
        return self.structural

    def outputs(
        self,
        show_regen_temps=True,
        show_cold_temps=True,
        show_dp=True,
        show_pressures=True,
        show_qdot=True,
        show_re=True,
        show_pr=False,
        show_coolant_density=True,
        show_coolant_velocity=True,
        show_wall_thermal_gradient=True,
        show_tangential_stresses=True,
        show_longitudinal_stresses=True,
        show_htc=True,
        show_coolant_htc_both=False,
        show_bartz_sigma=True,
        show_heat_flux=True,
        show_fin_efficiency=True,
    ):
        #### ---- PLOTS ---- ####
        pa_to_ksi = 1.0 / 6894757.293168
        coolant_rho_arr = np.array([
            PropsSI("D", "T", float(T_val), "P", float(P_val), self.coolantName)
            for T_val, P_val in zip(self.T_c_arr, self.P_c_arr)
        ])

        if show_regen_temps or show_cold_temps:
            fig, axs = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
            if show_regen_temps:
                axs[0].set_title("Regen Circuit Temperatures")
                axs[0].set_ylabel("Temperature [K]")
                axs[0].plot(self.circuit_z, self.engine.T[self.circuit_outlet:self.circuit_inlet + 1], color='orange', label='Combustion Gas Temp')
                if self.tbc_thickness > 0.0:
                    axs[0].plot(self.circuit_z, self.T_tbc_arr, color='m', linestyle='--', label='TBC Surface Temp')
                axs[0].plot(self.circuit_z, self.T_hw_arr, color='r', label='Hot Wall Temp')
                axs[0].plot(self.circuit_z, self.T_cw_arr, color='b', label='Cold Wall Temp')
                axs[0].plot(self.circuit_z, self.T_c_arr, color='y', label='Bulk Coolant Temp')
                axs[0].plot(
                    self.circuit_z,
                    self.coolant_boiling_arr,
                    color='k',
                    linestyle='--',
                    label=f'{self.coolantName} Boiling limit (local P)',
                )
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_cold_temps:
                axs[1].set_title("Regen Circuit Cold Side Temperatures")
                axs[1].set_xlabel("Axial Position [m]")
                axs[1].set_ylabel("Temperature [K]")
                axs[1].plot(self.circuit_z, self.T_cw_arr, color='b', label='Cold Wall Temp')
                axs[1].plot(self.circuit_z, self.T_c_arr, color='y', label='Bulk Coolant Temp')
                axs[1].plot(
                    self.circuit_z,
                    self.coolant_boiling_arr,
                    color='k',
                    linestyle='--',
                    label=f'{self.coolantName} Boiling limit (local P)',
                )
                axs[1].legend()
            else:
                axs[1].set_visible(False)
            fig.tight_layout()

        if show_dp or show_pressures:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            if show_dp:
                axs[0].plot(self.circuit_z, self.DP_arr / 6894.7, label='Pressure Drop')
                axs[0].set_title("Pressure Drop along Engine")
                axs[0].set_ylabel("DP [psia]")
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_pressures:
                axs[1].plot(self.circuit_z, self.engine.P[self.circuit_outlet:self.circuit_inlet + 1], label='Chamber Static Pressure')
                axs[1].plot(self.circuit_z, self.P_c_arr / 6894.7, label='Coolant Pressure')
                axs[1].set_title("Pressures along Engine")
                axs[1].set_xlabel("Axial Position [m]")
                axs[1].set_ylabel("Pressure [psia]")
                axs[1].legend()
            else:
                axs[1].set_visible(False)
            fig.tight_layout()

        if show_qdot or show_re or show_pr or show_coolant_density or show_coolant_velocity or show_fin_efficiency:
            plot_specs = []
            if show_qdot:
                plot_specs.append(("Radial Heat Flow", self.Qdot_arr, "Qdot", "Qdot [J]"))
            if show_re:
                plot_specs.append(("Coolant Re", self.Re, "Reynolds Number", "Re [---]"))
            if show_pr:
                plot_specs.append(("Coolant Pr", self.Pr_arr, "Prandtl Number", "Pr [---]"))
            if show_coolant_density:
                plot_specs.append(("Coolant Density", coolant_rho_arr, "Coolant Density", "Density [kg/m^3]"))
            if show_coolant_velocity:
                plot_specs.append(("Coolant Velocity", self.coolant_vel_arr, "Coolant Velocity", "Velocity [m/s]"))
            if show_fin_efficiency:
                plot_specs.append(("Fin Efficiency", self.eta_f_arr, "Fin Efficiency", "η_f [---]"))

            nplots = len(plot_specs)
            ncols = 2
            nrows = int(np.ceil(nplots / ncols))
            fig, axs = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), sharex=True)
            axs = np.atleast_1d(axs).ravel()

            for ax, (title, yvals, label, ylab) in zip(axs, plot_specs):
                ax.plot(self.circuit_z, yvals, label=label)
                ax.set_title(title)
                ax.set_ylabel(ylab)
                ax.legend()

            for ax in axs[nplots:]:
                ax.set_visible(False)

            for ax in axs[:nplots]:
                ax.set_xlabel("Axial Position [m]")

            fig.tight_layout()

        if show_wall_thermal_gradient or show_tangential_stresses or show_longitudinal_stresses:
            fig, axs = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
            if show_wall_thermal_gradient:
                axs[0].plot(self.circuit_z, self.T_hw_arr - self.T_cw_arr, label='Wall Thermal Gradient')
                axs[0].set_title("Wall Thermal Gradient")
                axs[0].set_ylabel("T_hw - T_cw [K]")
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_tangential_stresses:
                axs[1].plot(self.circuit_z, self.sigma_t * pa_to_ksi, label='Total Tangential Stress')
                axs[1].plot(self.circuit_z, self.sigma_pt * pa_to_ksi, label='Pressure Tangential Stress')
                axs[1].plot(self.circuit_z, self.sigma_tt * pa_to_ksi, label='Thermal Tangential Stress')
                axs[1].set_title("Tangential Stresses")
                axs[1].set_ylabel("Stress [ksi]")
                axs[1].legend()
            else:
                axs[1].set_visible(False)

            if show_longitudinal_stresses:
                axs[2].plot(self.circuit_z, self.sigma_l * pa_to_ksi, label='Total Longitudinal Stress')
                axs[2].plot(self.circuit_z, self.sigma_pl * pa_to_ksi, label='Thermal Longitudinal Stress')
                axs[2].set_title("Longitudinal Stresses")
                axs[2].set_xlabel("Axial Position [m]")
                axs[2].set_ylabel("Stress [ksi]")
                axs[2].legend()
            else:
                axs[2].set_visible(False)
            fig.tight_layout()

        if show_htc or show_bartz_sigma:
            print(f"Hg Length: {len(self.h_hg_arr)}")
            fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            if show_htc:
                axs[0].plot(self.circuit_z, self.h_hg_arr, label='Hot Gas HTC')
                axs[0].plot(self.circuit_z, self.h_c_arr, label='Coolant HTC (Selected)')
                if show_coolant_htc_both:
                    axs[0].plot(self.circuit_z, self.h_c_gn_arr, '--', label='Coolant HTC (Gnielinski)')
                    axs[0].plot(self.circuit_z, self.h_c_db_arr, '--', label='Coolant HTC (Dittus-Boelter)')
                axs[0].set_title("HTCs along Engine")
                axs[0].set_ylabel("HTC [W/m^2-K]")
                axs[0].legend()
            else:
                axs[0].set_visible(False)

            if show_bartz_sigma:
                axs[1].plot(self.circuit_z, self.sigma, label='Bartz Sigma')
                axs[1].set_title("Bartz Sigma along Engine")
                axs[1].set_xlabel("Axial Position [m]")
                axs[1].set_ylabel("Sigma [---]")
                axs[1].legend()
            else:
                axs[1].set_visible(False)
            fig.tight_layout()

        if show_heat_flux:
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            ax.plot(self.circuit_z, self.q_flux_arr / 1e6, color='C3', label='Heat Flux')
            ax.set_title("Gas-Side Heat Flux into Chamber Wall")
            ax.set_xlabel("Axial Position [in]")
            ax.set_ylabel("Heat Flux [MW/m²]")
            ax.legend()
            ax.grid(True, alpha=0.3)
            fig.tight_layout()

        # Unit conversions and calculations prior to printing
        self.channel_DP = float(np.sum(self.DP_arr)) / 6894.76 # [psid]
        boiling_temp_margin = float(np.min(self.coolant_boiling_arr - self.T_c_arr))
        coolant = Fluid("T", self.T_c_arr[-1], "P", self.P_c_arr[-1], self.coolantName)
        rho1 = coolant.rho
        coolant = Fluid("T", self.T_c_arr[0], "P", self.P_c_arr[0], self.coolantName)
        rho2 = coolant.rho
        avg_rho = (rho1 + rho2) / 2
        self.chan_CdA = (self.chan_mdot / (np.sqrt(2 * avg_rho * (self.P_c_arr[-1] - self.P_c_arr[0])))) * 1000000

        #### ---- SCALAR OUTPUTS ---- ####
        print("\n#### ---- REGEN CIRCUIT OUTPUTS ---- ####")
        print(f"Channel CdA [mm^2]: {round(self.chan_CdA, 3)}")
        print(f"Channel Coolant Mdot [kg/s]: {round(self.chan_mdot, 4)}")
        print(f"Total Coolant Mdot [kg/s]: {round(self.tot_coolant_mdot, 3)}")
        print(f"Channel DP [psid]: {round(self.channel_DP, 2)}")
        print(f"Min Coolant Pressure [psia]: {round(float(np.min(self.P_c_arr)) / 6894.76, 2)}")
        print(f"Max Coolant Velocity [m/s]: {round(float(np.max(self.coolant_vel_arr)), 2)}")
        print(f"Max Friction Factor f_D [---]: {round(float(np.max(self.f_D_arr)), 4)}")
        print(f"Regen Curve Length [m]: {round(self.regen_curve_length_m, 4)}")
        print(f"Regen Curve Length [in]: {round(self.regen_curve_length_m / 0.0254, 2)}")
        print(f"# of Channels: {self.N}")
        for name in self.STATION_NAMES:
            st = self.stations[name]
            print(f"  {name:8s} z={st['z']:.3f} in  C_w={st['C_w'] * 1000:.2f} mm  C_h={st['C_h'] * 1000:.2f} mm  AR={st['C_h'] / st['C_w']:.2f}")
        print(f"Min Fin Width [mm]: {round(float(np.min(self.f_w_arr)) * 1000, 3)}")
        if self.tbc_thickness > 0.0:
            print(f"TBC Thickness [mm]: {round(self.tbc_thickness * 1000, 3)}")
            print(f"TBC k [W/m-K]: {self.tbc_k}")
            print(f"Max TBC Surface Temp [K]: {round(float(np.max(self.T_tbc_arr)), 2)}")
        print(f"Max Hot Wall Temp [K]: {round(max(self.T_hw_arr), 2)}")
        print(f"Coolant Outlet Temp [K]: {round(max(self.T_c_arr), 2)}")
        print(f"Max Heat Flux [MW/m^2]: {round(float(np.max(self.q_flux_arr)) / 1e6, 3)}")
        print(f"Coolant Boiling Temp Margin [K]: {round(boiling_temp_margin, 2)}")

        plt.show()
        
    

