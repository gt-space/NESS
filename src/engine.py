## --- ENGINE.PY --- ###
# Define an engine with the class Engine to size the contour and calculate performance parameters.

# Brief Description - Engine Class holds all geometry, performance parameters, and gas properties along the 
# engine. RocketIsp used for all performance calculations.

import matplotlib.pyplot as plt
import numpy as np

from rocketisp.geometry import Geometry
from rocketisp.efficiencies import Efficiencies
from rocketisp.stream_tubes import CoreStream
from rocketisp.rocket_isp import RocketThruster
from rocketcea.cea_obj_w_units import CEA_Obj
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d

### --- ENGINE AXIAL PROPERTIES --- ###
# 1. Prandtl Number
# 2. Dynamic Viscosity (mu)
# 3. Static Pressure (Pa)
# 4. Area Ratio (A / A*)
# 5. Mach Number
# 6. Static Temperature (T)
# 7.  Static Pressure (P)
# 8. Static Density (Rho)
# 9. Specific Heat Constant pressure (Cp)
# 10. Gamma
# 11. Molecular Weight
# 12. Thermal Conductivity (k)

class Engine:

    def __init__(self, thrust, Pc, Pe, MRcore, oxName, fuName, name, CR, Lstar, cstarEff, frozen, pcentFFC = 0.0, Pamb = 14.8, numPts = 200, verbose = False,
                 bell = True, chmbR=None, chmbL=None, contAngle=None, throatR=None, throatL=None, expAngle=None, exitR=None,
                 chmbFilletR=0.0, throatFilletR=0.0, cstarEffInFlow=False):
        self.thrust = thrust
        self.Pc_psia = Pc
        self.Pc = Pc
        self.Pe = Pe
        self.Pamb = Pamb
        self.MRcore = MRcore
        self.oxName = oxName
        self.fuName = fuName
        self.pcentFFC = pcentFFC
        self.CR = CR
        self.cstarEff = cstarEff
        self.Lstar = Lstar
        self.name = name
        self.Lprime = self.Lstar / self.CR
        self.numPts = numPts
        self.bell = bell
        self.chmbR = chmbR
        self.chmbL = chmbL
        self.contAngle = contAngle
        self.throatR = throatR
        self.throatL = throatL
        self.expAngle = expAngle
        self.exitR = exitR
        self.chmbFilletR = 0.0 if chmbFilletR is None else float(chmbFilletR)   # [in] chamber-to-cone fillet (conical only)
        self.throatFilletR = 0.0 if throatFilletR is None else float(throatFilletR)   # [in] cone-to-throat fillet (conical only)
        self.frozen = frozen
        # False (default, original behavior): flows follow Pc*At/c*_ideal. True: flows follow Pc*At/(cstarEff*c*_ideal),
        # i.e. real propellant flow for a given Pc when the combustion efficiency is below 1.
        self.cstarEffInFlow = bool(cstarEffInFlow)
        self.createEngine(verbose=verbose)
        self.calcGasProperties(frozen)
        # The point the hardware was designed for (Pc, MR and flows). Off-design analysis can move the engine away from it.
        self.design_state = dict(Pc_psia=float(self.Pc_psia), MR=float(self.MRcore),
                                 ox_mdot=float(self.ox_mdot), fu_mdot=float(self.fu_mdot))

    def createEngine(self, verbose=False):
        G = Geometry(CR=self.CR, LchamberInp=self.Lprime, RchmConv=2.5, pcentBell=80)
        E = Efficiencies( ERE = self.cstarEff )
        C = CoreStream(geomObj=G, effObj=E,
                       oxName=self.oxName, fuelName=self.fuName, MRcore=self.MRcore,
                       Pc=self.Pc, Pamb=self.Pamb, adjCstarODE=self.cstarEff,
                       pcentFFC=self.pcentFFC, ko=0.035)
        
        self.R = RocketThruster(name=self.name, 
                           coreObj=C, injObj=None
                           )
        
        self.R.set_eps_to_equal_pexit( Pexit_psia=self.Pe )
        self.R.scale_Rt_to_Thrust( self.thrust, Pamb=self.Pamb )
        # What RocketIsp sized for the requested thrust (kept for the design consistency check)
        self.rocketisp_throatR = float(np.sqrt(self.R.geomObj.At / np.pi))                     # [in]
        self.rocketisp_tot_mdot = float((self.R.coreObj.wdotOx + self.R.coreObj.wdotFl) * 0.453592)   # [kg/s]
        if (not self.bell) and self.throatR is None:
            self.throatR = self.rocketisp_throatR      # conical + throatR=None: size the throat for the requested thrust

        self.R.geomObj.getNozObj()
        self.Contour_z = [-self.R.geomObj.Lcham] + self.R.geomObj.nozObj.abs_zContour
        self.Contour_r = [self.R.geomObj.Rinj] + self.R.geomObj.nozObj.abs_rContour

        if self.bell == True:
            self.Contour_z = [-self.R.geomObj.Lcham] + self.R.geomObj.nozObj.abs_zContour
            self.Contour_r = [self.R.geomObj.Rinj] + self.R.geomObj.nozObj.abs_rContour
        else:
            self.Contour_z, self.Contour_r = self.makeConicalContour()
        
        self.Contour_z, self.Contour_r, self.throat_ind, self.chamberEnd_ind = self.resample_nozzle_contour(
        z_contour = self.Contour_z, r_contour = self.Contour_r, numPts = self.numPts, verbose=verbose)
        if not self.bell:
            # The cubic-spline resample can undershoot the throat radius by ~1e-5 in; the conical contour never goes below it.
            self.Contour_r = np.maximum(np.asarray(self.Contour_r, dtype=float), self.throatR)        
        
        print(f"Chamber End Ind: {self.chamberEnd_ind}")
        print(f"Throat Ind: {self.throat_ind}")
        
        if self.bell == True:
            self.At = self.R.geomObj.At
        else:
            self.At = np.pi * (self.throatR ** 2)

        self.Dt = np.sqrt((4 / np.pi) * self.At) * 0.0254 # [in] --> [m]
        self.eps = self.R.geomObj.eps
        self.r_ch = self.Contour_r[0]
        self.r_throat = self.Contour_r[self.throat_ind]
        self.r_up_throat = self.R.geomObj.RupThroat * 0.0254 # [in] --> [m]
        if not self.bell and self.throatFilletR > 0.0:
            # Conical contour: the cone-to-throat fillet IS the upstream throat radius of curvature used by Bartz.
            self.r_up_throat = self.throatFilletR * 0.0254 # [in] --> [m]
        self.fu_mdot = self.R.coreObj.wdotFl * 0.453592
        self.ox_mdot = self.R.coreObj.wdotOx * 0.453592

    def conicalGeometry(self):
        """
        Key points and fillets of the conical contour (all lengths in inches, axial origin at the throat centre).

        Layout, left to right:  cylindrical chamber -> [chamber fillet] -> converging cone -> [throat fillet]
        -> flat throat (length throatL, between tangent points) -> diverging cone.

        chmbFilletR rounds the chamber/cone corner and throatFilletR rounds the cone/throat corner. Each fillet
        is a circular arc tangent to both walls it joins. Because the throat arc needs room beyond the flat
        throat, the converging side is shifted upstream by its tangent length; the flat throat, the diverging
        cone and the exit stay where they are.

        Returns a dict with x_c1 / x_c2 (the two virtual sharp corners), x_start (chamber start), x_fs / x_fe (flat
        throat start/end), x_exit, and fillet dicts f1 / f2 (None when the radius is 0). A fillet dict has the
        arc centre (xc, rc), radius R, tangent length T and axial extent x_start..x_end.
        """
        th = np.radians(self.contAngle)
        al = np.radians(self.expAngle)
        R1 = max(float(self.chmbFilletR), 0.0)
        R2 = max(float(self.throatFilletR), 0.0)
        T1 = R1 * np.tan(th / 2.0)
        T2 = R2 * np.tan(th / 2.0)

        x_fs = -self.throatL / 2.0
        x_fe = self.throatL / 2.0
        x_c2 = x_fs - T2                                           # cone / throat virtual corner
        x_c1 = x_c2 - (self.chmbR - self.throatR) / np.tan(th)     # chamber / cone virtual corner
        x_start = x_c1 - self.chmbL
        x_exit = x_fe + (self.exitR - self.throatR) / np.tan(al)

        cone_length = (self.chmbR - self.throatR) / np.sin(th)     # along the cone between the virtual corners
        if T1 + T2 >= 0.999 * cone_length:
            raise ValueError(
                f"The fillets need {T1 + T2:.3f} in of converging cone (chmbFilletR {R1} in + throatFilletR {R2} in) "
                f"but only {cone_length:.3f} in is available. Reduce a fillet radius or lengthen the cone."
            )
        if T1 > self.chmbL:
            raise ValueError(
                f"chmbFilletR = {R1} in needs {T1:.3f} in of straight chamber but chmbL is only {self.chmbL} in. "
                "Reduce chmbFilletR or lengthen the chamber."
            )

        f1 = f2 = None
        if R1 > 0.0:
            f1 = dict(xc=x_c1 - T1, rc=self.chmbR - R1, R=R1, T=T1,
                      x_start=x_c1 - T1, x_end=x_c1 - T1 + R1 * np.sin(th))
        if R2 > 0.0:
            f2 = dict(xc=x_fs, rc=self.throatR + R2, R=R2, T=T2,
                      x_start=x_fs - R2 * np.sin(th), x_end=x_fs)
        return dict(x_c1=x_c1, x_c2=x_c2, x_start=x_start, x_fs=x_fs, x_fe=x_fe, x_exit=x_exit, f1=f1, f2=f2)

    def conicalContourRadius(self, x, g=None):
        """Gas-side wall radius [in] of the conical contour at axial position x [in]."""
        g = g or self.conicalGeometry()
        th = np.radians(self.contAngle)
        al = np.radians(self.expAngle)
        f1, f2 = g["f1"], g["f2"]
        x_cyl_end = f1["x_start"] if f1 else g["x_c1"]
        x_cone_end = f2["x_start"] if f2 else g["x_fs"]
        if x <= x_cyl_end:
            return self.chmbR
        if f1 is not None and x < f1["x_end"]:
            return f1["rc"] + np.sqrt(max(f1["R"] ** 2 - (x - f1["xc"]) ** 2, 0.0))
        if x <= x_cone_end:
            return self.throatR + np.tan(th) * (g["x_c2"] - x)
        if f2 is not None and x < f2["x_end"]:
            return f2["rc"] - np.sqrt(max(f2["R"] ** 2 - (x - f2["xc"]) ** 2, 0.0))
        if x <= g["x_fe"]:
            return self.throatR
        return self.throatR + np.tan(al) * (x - g["x_fe"])

    def makeConicalContour(self):
        if self.chmbR is None or self.throatR is None:
            raise ValueError("You set bell=False but did not provide the conical dimensions")

        g = self.conicalGeometry()
        # A fillet needs a finer raw contour so the arc is smooth before it is resampled to numPts.
        has_fillet = g["f1"] is not None or g["f2"] is not None
        n_raw = self.numPts if not has_fillet else max(self.numPts, 800)
        z_coords = np.linspace(g["x_start"], g["x_exit"], n_raw)
        temp_r_coords = [self.conicalContourRadius(x, g) for x in z_coords]

        return z_coords, np.array(temp_r_coords)

    def calcGasProperties(self, frozen=1, Pc_psia_override=None, MR_override=None):
        # Optionally override Pc and MR for off-design analysis (throttle sweeps, etc.)
        if Pc_psia_override is not None:
            self.Pc_psia = float(Pc_psia_override)
        if MR_override is not None:
            self.MRcore = float(MR_override)

        # Reset Pc to psia for CEA calls (safe for re-entry).
        self.Pc = self.Pc_psia

        # Initialize some Arrays
        self.AreaRatio = np.array([])
        self.M = np.array([])
        self.P = np.array([])
        self.T = np.array([])
        self.rho = np.array([])

        # Setup CEA Object, general calculations
        C = CEA_Obj(
            oxName=self.oxName, fuelName=self.fuName,
            isp_units='sec',
            cstar_units='m/s',
            pressure_units='psia',
            temperature_units='K',
            sonic_velocity_units='m/s',
            enthalpy_units='J/kg',
            density_units='kg/m^3',
            specific_heat_units='J/kg-K',
            viscosity_units='poise',
            thermal_cond_units="W/cm-degC",
            )
        T0 = C.get_Tcomb(self.Pc_psia, self.MRcore) * self.cstarEff**2 # Combustion temp knockdown
        self.Tcomb = float(T0)  # Chamber (combustion) temperature from CEA [K]
        P0 = self.Pc_psia
        rho0 = C.get_Chamber_Density(self.Pc_psia, self.MRcore)
        self.cstar = C.get_Cstar(self.Pc_psia, self.MRcore)

        # Calculate Area Ratio
        for i in range(len(self.Contour_r)):
            # Calculate Propertiese
            AreaRat = (np.pi * self.Contour_r[i]**2) / self.At
            self.AreaRatio = np.append(self.AreaRatio, AreaRat)

        # --- Chamber Gas Properties Extrapolation --- 
        Cp_chamb, mu_chamb, k_chamb, Pr_chamb = C.get_Chamber_Transport(self.Pc_psia, self.MRcore, frozen=frozen)
        MW_chamb, gamma_chamb = C.get_Chamber_MolWt_gamma(self.Pc_psia, self.MRcore)

        # Assign Properties
        self.Cp = np.full(self.chamberEnd_ind + 1, Cp_chamb)
        self.mu = np.full(self.chamberEnd_ind + 1, mu_chamb)
        self.k = np.full(self.chamberEnd_ind + 1, k_chamb)
        self.Pr = np.full(self.chamberEnd_ind + 1, Pr_chamb)
        self.MW = np.full(self.chamberEnd_ind + 1, MW_chamb)
        self.gamma = np.full(self.chamberEnd_ind + 1, gamma_chamb)

        # --- Converging Section Gas Properties (Linear Correlation) --- 
        # Throat Properties
        Cp_t, mu_t, k_t, Pr_t = C.get_Throat_Transport(self.Pc_psia, self.MRcore, frozen=frozen)
        MW_t, gamma_t = C.get_Throat_MolWt_gamma(self.Pc_psia, self.MRcore, frozen=frozen)
        
        # Property Z-Positions
        z_prop1 = self.Contour_z[self.chamberEnd_ind]
        z_prop2 = self.Contour_z[self.throat_ind]
        z_arr = self.Contour_z[self.chamberEnd_ind + 1 : self.throat_ind]

        # Property Extrapolations
        Cp = self.extrapolate_properties(Cp_chamb, z_prop1, Cp_t, z_prop2, z_arr)
        mu = self.extrapolate_properties(mu_chamb, z_prop1, mu_t, z_prop2, z_arr)
        k = self.extrapolate_properties(k_chamb, z_prop1,k_t, z_prop2, z_arr)
        Pr = self.extrapolate_properties(Pr_chamb, z_prop1, Pr_t, z_prop2, z_arr)
        MW = self.extrapolate_properties(MW_chamb, z_prop1, MW_t, z_prop2, z_arr)
        gamma = self.extrapolate_properties(gamma_chamb, z_prop1, gamma_t, z_prop2, z_arr)

        # Assign Properties
        self.Cp = np.concatenate((self.Cp, Cp))
        self.mu = np.concatenate((self.mu, mu))
        self.k = np.concatenate((self.k, k))
        self.Pr = np.concatenate((self.Pr, Pr))
        self.MW = np.concatenate((self.MW, MW))
        self.gamma = np.concatenate((self.gamma, gamma))

        # --- Throat ---
        Cp, mu, k, Pr = C.get_Throat_Transport(self.Pc_psia, self.MRcore, frozen=frozen)
        MW, gamma = C.get_Throat_MolWt_gamma(self.Pc_psia, self.MRcore, frozen=frozen)
        
        # Assign Properties
        self.Cp = np.append(self.Cp, Cp)
        self.mu = np.append(self.mu, mu)
        self.k = np.append(self.k, k)
        self.Pr = np.append(self.Pr, Pr)
        self.MW = np.append(self.MW, MW)
        self.gamma = np.append(self.gamma, gamma)

        # --- Diverging Section Gas Properties --- 
        for eps in self.AreaRatio[self.throat_ind + 1 :]:
            Cp, mu, k, Pr = C.get_Exit_Transport(self.Pc_psia, self.MRcore, eps=eps, frozen=frozen)
            MW, gamma = C.get_exit_MolWt_gamma(self.Pc_psia, self.MRcore, eps=eps, frozen=frozen)

            # Assign Properties
            self.Cp = np.append(self.Cp, Cp)
            self.mu = np.append(self.mu, mu)
            self.k = np.append(self.k, k)
            self.Pr = np.append(self.Pr, Pr)
            self.MW = np.append(self.MW, MW)
            self.gamma = np.append(self.gamma, gamma)
            
        self.solveAreaMach() # Solves for M from Area Ratio

        for i in range(len(self.Contour_z)):
            
            # Calculate Properties
            P = P0 / ((1 + 0.5*(self.gamma[i] - 1)*self.M[i]**2)**(self.gamma[i] / (self.gamma[i] - 1)))
            T = T0 / (1 + 0.5*(self.gamma[i] - 1)*self.M[i]**2)
            rho = rho0 / ((1 + 0.5*(self.gamma[i] - 1)*self.M[i]**2)**(1 / (self.gamma[i] - 1)))

            # Assign Properties
            self.P = np.append(self.P, P)
            self.T = np.append(self.T, T) 
            self.rho = np.append(self.rho, rho)
        
        # Unit Conversions (for Bartz)
        self.mu = self.mu / 10 # [Poise] --> [Pa * s]
        self.Pc = self.Pc_psia * 6894.76 # [psia] --> [Pa]

        # Compute mass flows from Pc, cstar, At for off-design re-calls.
        # At is in in^2; cstar in m/s; Pc in Pa.
        At_m2 = self.At * 0.00064516  # in^2 → m^2
        tot_mdot = (self.Pc * At_m2) / (self.cstar * self.flow_cstar_factor())  # kg/s
        self.fu_mdot = tot_mdot / (1.0 + self.MRcore)
        self.ox_mdot = tot_mdot - self.fu_mdot

    def flow_cstar_factor(self):
        """Multiplier on the CEA c* used in the Pc <-> mass-flow relation (cstarEff if cstarEffInFlow else 1)."""
        return float(self.cstarEff) if getattr(self, "cstarEffInFlow", False) else 1.0

    def cea_thrust_lbf(self, Pc_psia, MR):
        """Ideal CEA thrust [lbf] of THIS geometry (throat At, exit area ratio) at chamber pressure and mixture ratio."""
        from utils import _cea_obj_w_units_like_engine
        C = _cea_obj_w_units_like_engine(self.oxName, self.fuName)
        eps = float(self.AreaRatio[-1])
        _, CFamb, _ = C.get_PambCf(Pamb=self.Pamb, Pc=float(Pc_psia), MR=float(MR), eps=eps)
        return float(CFamb) * float(Pc_psia) * self.At          # psia * in^2 = lbf

    def design_consistency(self):
        """
        Does the geometry you entered actually make the thrust you asked for at the design Pc and MR?
        Returns a dict (all lengths in inches, thrust in lbf, flows in kg/s).
        """
        thrust_geom = self.cea_thrust_lbf(self.Pc_psia, self.MRcore)
        return dict(
            requested_thrust=float(self.thrust),
            geometry_thrust=thrust_geom,
            ratio=thrust_geom / float(self.thrust),
            throatR=float(self.throatR) if self.throatR is not None else float(np.sqrt(self.At / np.pi)),
            rocketisp_throatR=self.rocketisp_throatR,
            design_flow=float(self.fu_mdot + self.ox_mdot),
            rocketisp_flow=self.rocketisp_tot_mdot,
        )

    def print_design_check(self):
        d = self.design_consistency()
        print("\n---- DESIGN CHECK: thrust requested vs thrust of this geometry ----")
        print(f"  Requested thrust            : {d['requested_thrust']:.0f} lbf at Pc {self.Pc_psia:g} psia, MR {self.MRcore:g}")
        print(f"  This geometry (ideal CEA CF): {d['geometry_thrust']:.0f} lbf  ({d['ratio'] * 100:.0f}% of requested)")
        print(f"  Throat radius               : {d['throatR']:.3f} in used; RocketIsp sizes {d['rocketisp_throatR']:.3f} in for {d['requested_thrust']:.0f} lbf")
        print(f"  Design flow used by this code: {d['design_flow']:.3f} kg/s; RocketIsp's {d['requested_thrust']:.0f} lbf flow: {d['rocketisp_flow']:.3f} kg/s"
              + ("" if self.cstarEffInFlow else f"  (c* efficiency {self.cstarEff:g} is not applied to flows; set cstarEffInFlow=True to apply it)"))
        if abs(d["ratio"] - 1.0) > 0.05:
            print(f"  NOTE: the throat is {'smaller' if d['ratio'] < 1 else 'larger'} than a {d['requested_thrust']:.0f} lbf design needs. "
                  f"Set throatR = {d['rocketisp_throatR']:.3f} in (or None to size it automatically) to match.")

    def exportGeometry(self, filename, wall_thickness_in=0.0):
        """
        Export the nozzle contour as [X, Y] coordinates in a .txt file for import into CAD.

        Parameters
        ----------
        filename : str
            Base path (``.txt`` is appended).
        wall_thickness_in : float, optional
            Inches added to the inner (gas-side) radius at each axial station. Use ``0`` for the
            inner contour (ID). Set to the structural wall thickness ID→OD to export the outer
            contour (OD).
        """
        dr = float(wall_thickness_in)
        # Export the nozzle contour as [X, Y] coordinates in a .txt file for import into CAD.
        if self.bell:
            x = [float(r) + dr for r in self.Contour_r]
            y = [-n for n in self.Contour_z]

            with open(filename + '.txt', "w") as f:
                for xi, yi in zip(x, y):
                    f.write(f"{xi:.6f}in\t{yi:.6f}in\t0in\n")
        else:
            # Conical nozzle: export the critical points, with sampled arcs where a fillet is used
            g = self.conicalGeometry()
            th = np.radians(self.contAngle)
            f1, f2 = g["f1"], g["f2"]

            # x_crit = radius (r), y_crit = axial position (z)
            x_crit = [self.chmbR]
            y_crit = [g["x_start"]]
            if f1 is None:
                x_crit.append(self.chmbR)
                y_crit.append(g["x_c1"])
            else:
                ang = np.linspace(0.0, th, 13)
                y_crit += list(f1["xc"] + f1["R"] * np.sin(ang))
                x_crit += list(f1["rc"] + f1["R"] * np.cos(ang))
            if f2 is None:
                x_crit.append(self.throatR)
                y_crit.append(g["x_fs"])
            else:
                ang = np.linspace(-th, 0.0, 13)
                y_crit += list(f2["xc"] + f2["R"] * np.sin(ang))
                x_crit += list(f2["rc"] - f2["R"] * np.cos(ang))
            x_crit += [self.throatR, self.exitR]
            y_crit += [g["x_fe"], g["x_exit"]]

            # Export exact critical points
            with open(filename + '.txt', "w") as f:
                for xi, yi in zip(x_crit, y_crit):
                    f.write(f"{xi:.6f}in\t{yi:.6f}in\t0in\n")
        
    def solveAreaMach(self, tol=1e-6):
        """
        Solve for Mach number with Area Ratio.
        
        Returns:
            np.array of Mach numbers
        """
        for i in range(len(self.Contour_r)):
            # Handle M = 1 case
            if abs(self.AreaRatio[i] - 1.0) < tol:
                self.M = np.append(self.M, 1)
                continue
            
            area_mach = lambda M: ((1/M) * ((2/(self.gamma[i]+1)*(1 + (self.gamma[i]-1)/2 * M**2)))**((self.gamma[i]+1)/(2*(self.gamma[i]-1)))) - self.AreaRatio[i]

            # Subsonic root (if exists)
            try:
                sol_sub = root_scalar(area_mach, bracket=(1e-8, 1.0), method='brentq')
                M_sub = sol_sub.root
            except ValueError:
                M_sub = None

            # Supersonic root (if exists)
            try:
                sol_sup = root_scalar(area_mach, bracket=(1.0, 20.0), method='brentq')
                M_sup = sol_sup.root
            except ValueError:
                M_sup = None

            # Determine if we're in converging or diverging section
            # Check if area is increasing or decreasing
            if not self.bell:
                is_diverging = i > self.throat_ind
            else:
                if i > 0:
                    dA = self.AreaRatio[i] - self.AreaRatio[i-1]
                    is_diverging = dA > 0
                else:
                    is_diverging = False
            
            # Pick branch based on nozzle section
            if is_diverging and M_sup is not None:
                self.M = np.append(self.M, M_sup)
            elif M_sub is not None:
                self.M = np.append(self.M, M_sub)
            else:
                raise ValueError(f"No valid Mach number found for A/A*={self.AreaRatio[i]:.6f}")

    def extrapolate_properties(self, prop1, z_pos1, prop2, z_pos2, z_arr):
        m = (prop2 - prop1) / (z_pos2 - z_pos1)
        b = prop1 - m * z_pos1

        extrap_props = []

        for val in z_arr:
            new_prop = m * val + b
            extrap_props.append(new_prop)

        return np.array(extrap_props)

    def plot_gas_props(self):
        # Gamma
        plt.figure()
        plt.plot(self.Contour_z, self.gamma)
        plt.title("Gamma")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Gamma")

        # Cp
        plt.figure()
        plt.plot(self.Contour_z, self.Cp)
        plt.title("Specific Heat Constant Pressure (Cp)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Specific Heat Constant Pressure (Cp) (J/kg-K)")

        # Mu
        plt.figure()
        plt.plot(self.Contour_z, self.mu)
        plt.title("Dynamic Viscosity (Mu)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Dynamic Viscosity (Mu) (Pa * s)")

        # k
        plt.figure()
        plt.plot(self.Contour_z, self.k)
        plt.title("Thermal Conductivity (k)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Thermal Conductivity (k) (W/cm-C)")

        # Pr
        plt.figure()
        plt.plot(self.Contour_z, self.Pr)
        plt.title("Prandtl Number (Pr)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Prandtl Number (Pr)")

        # P
        print(len(self.P))
        plt.figure()
        plt.plot(self.Contour_z, self.P)
        plt.title("Static Pressure (P)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Static Pressure (P) (psia)")

        # T
        plt.figure()
        plt.plot(self.Contour_z, self.T)
        plt.title("Static Temperature (T)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Static Temperature (T) (K)")

        # Rho
        plt.figure()
        plt.plot(self.Contour_z, self.rho)
        plt.title("Static Density (Rho)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Static Density (Rho) (kg/m^3)")

        # MW
        plt.figure()
        plt.plot(self.Contour_z, self.MW)
        plt.title("Molecular Weight (MW)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Molecular Weight (MW)")

        # M
        plt.figure()
        plt.plot(self.Contour_z, self.M)
        plt.title("Mach Number (M)")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Mach Number (M)")

        # Area Ratio
        plt.figure()
        plt.plot(self.Contour_z, self.AreaRatio)
        plt.title("Area Ratio")
        plt.xlabel("Axial Position (in)")
        plt.ylabel("Area Ratio")

        plt.show()
    
    def resample_nozzle_contour(self, z_contour, r_contour, numPts, verbose=False):
        """
        Resample nozzle contour with equal spacing along the arc length.
        First two points define a horizontal chamber section.
        Throat location is preserved exactly from original contour.
        
        Parameters:
        -----------
        z_contour : array-like
            Original axial positions (first 2 points = chamber start/end)
        r_contour : array-like
            Original radial positions (first 2 points should have same r)
        numPts : int
            Desired number of points along the contour
        verbose : bool, optional
            If True, print diagnostics and show plot. Default is False.

        Returns:
        --------
        z_new : np.array
            Resampled axial positions
        r_new : np.array
            Resampled radial positions
        throat_idx : int
            Index of the throat (minimum radius) in the resampled contour
        chamber_end_idx : int
            Index of the last point in the chamber section
        """
        from scipy.interpolate import CubicSpline
        
        z_contour = np.asarray(z_contour, dtype=float).flatten()
        r_contour = np.asarray(r_contour, dtype=float).flatten()
        
        # Find original throat location (minimum radius in nozzle section)
        throat_idx_orig = np.argmin(r_contour)
        z_throat_exact = z_contour[throat_idx_orig]
        r_throat_exact = r_contour[throat_idx_orig]
        
        # Find the end of the flat chamber dynamically
        # It's the last point that has the exact same radius as the injector face
        r_inj = r_contour[0]
        chamber_mask = np.isclose(r_contour, r_inj, rtol=1e-5, atol=1e-5)
        
        # Find the last contiguous True value from the start
        chamber_end_idx_orig = 0
        for i in range(len(chamber_mask)):
            if chamber_mask[i]:
                chamber_end_idx_orig = i
            else:
                break
                
        # Separate chamber from nozzle using the true dynamic index
        z_chamber = [z_contour[0], z_contour[chamber_end_idx_orig]]
        r_chamber = [r_contour[0], r_contour[chamber_end_idx_orig]]
        r_chamber_avg = np.mean(r_chamber)

        z_nozzle = z_contour[chamber_end_idx_orig:]
        r_nozzle = r_contour[chamber_end_idx_orig:]

        
        # Calculate arc lengths
        chamber_length = abs(z_chamber[1] - z_chamber[0])
        
        # Nozzle: use parametric curve
        dz_noz = np.diff(z_nozzle)
        dr_noz = np.diff(r_nozzle)
        chord_noz = np.sqrt(dz_noz**2 + dr_noz**2)
        t_noz = np.concatenate([[0], np.cumsum(chord_noz)])
        t_noz = t_noz / t_noz[-1] if t_noz[-1] > 0 else t_noz
        
        # Create splines for nozzle section
        if len(z_nozzle) < 4:
            z_noz_spline = interp1d(t_noz, z_nozzle, kind='linear')
            r_noz_spline = interp1d(t_noz, r_nozzle, kind='linear')
        else:
            z_noz_spline = CubicSpline(t_noz, z_nozzle, bc_type='natural')
            r_noz_spline = CubicSpline(t_noz, r_nozzle, bc_type='natural')
        
        # High-resolution nozzle curve for arc length
        t_noz_fine = np.linspace(0, 1, 1000)
        z_noz_fine = z_noz_spline(t_noz_fine)
        r_noz_fine = r_noz_spline(t_noz_fine)
        
        dz_noz_fine = np.diff(z_noz_fine)
        dr_noz_fine = np.diff(r_noz_fine)
        ds_noz_fine = np.sqrt(dz_noz_fine**2 + dr_noz_fine**2)
        s_noz_fine = np.concatenate([[0], np.cumsum(ds_noz_fine)])
        nozzle_length = s_noz_fine[-1]
        
        # Find arc length to throat in the fine resolution curve
        throat_idx_fine = np.argmin(np.abs(z_noz_fine - z_throat_exact))
        s_throat = s_noz_fine[throat_idx_fine]
        
        # Total arc length
        total_length = chamber_length + nozzle_length
        
        # Distribute points proportionally
        chamber_fraction = chamber_length / total_length
        num_chamber_pts = max(2, int(numPts * chamber_fraction))
        num_nozzle_pts = numPts - num_chamber_pts
        
        # Generate chamber points (horizontal line)
        z_chamber_new = np.linspace(z_chamber[0], z_chamber[1], num_chamber_pts)
        r_chamber_new = np.full(num_chamber_pts, r_chamber_avg)
        
        # Generate nozzle points with throat constraint
        if num_nozzle_pts > 0:
            # Calculate how many points before and after throat
            throat_fraction = s_throat / nozzle_length
            num_before_throat = int(num_nozzle_pts * throat_fraction)
            num_after_throat = num_nozzle_pts - num_before_throat - 1  # -1 for throat point itself
            
            # Ensure at least one point on each side
            num_before_throat = max(1, num_before_throat)
            num_after_throat = max(1, num_after_throat)
            
            # Generate arc length values
            s_before = np.linspace(0, s_throat, num_before_throat, endpoint=False)
            s_throat_point = np.array([s_throat])
            s_after = np.linspace(s_throat, nozzle_length, num_after_throat + 1)[1:]  # Skip first (throat)
            
            s_noz_new = np.concatenate([s_before, s_throat_point, s_after])
            
            # Map arc length to parameter t
            t_noz_from_s = interp1d(s_noz_fine, t_noz_fine, kind='linear', 
                                    bounds_error=False, fill_value='extrapolate')
            t_noz_new = t_noz_from_s(s_noz_new)
            
            z_nozzle_new = np.asarray(z_noz_spline(t_noz_new)).flatten()
            r_nozzle_new = np.asarray(r_noz_spline(t_noz_new)).flatten()
            
            # Force throat point to exact original values
            throat_idx_in_nozzle = num_before_throat
            z_nozzle_new[throat_idx_in_nozzle] = z_throat_exact
            r_nozzle_new[throat_idx_in_nozzle] = r_throat_exact
            
            # Combine chamber and nozzle (remove duplicate point at junction)
            z_new = np.concatenate([z_chamber_new[:-1], z_nozzle_new])
            r_new = np.concatenate([r_chamber_new[:-1], r_nozzle_new])
        else:
            z_new = z_chamber_new
            r_new = r_chamber_new
        
        # Find throat index in new array
        throat_idx = int(np.argmin(r_new))
        
        # Verify throat is exact
        if verbose:
            print(f"Original throat: z={z_throat_exact:.6f}, r={r_throat_exact:.6f}")
            print(f"Resampled throat: z={z_new[throat_idx]:.6f}, r={r_new[throat_idx]:.6f}")
            print(f"Difference: dz={abs(z_new[throat_idx]-z_throat_exact):.9f}, dr={abs(r_new[throat_idx]-r_throat_exact):.9f}")
        
        # Chamber end index
        chamber_end_idx = int(num_chamber_pts - 1)
        
        # Verbose output
        if verbose:
            print(f"\nChamber length: {chamber_length:.4f}")
            print(f"Nozzle length: {nozzle_length:.4f}")
            print(f"Total length: {total_length:.4f}")
            print(f"Chamber points: {num_chamber_pts}, Nozzle points: {num_nozzle_pts}")
            print(f"Throat index: {throat_idx}")
            print(f"Chamber end index: {chamber_end_idx}")
            
            # Visualization
            plt.figure(figsize=(14, 7))
            
            # Plot original points
            plt.scatter(z_contour, r_contour, c='red', s=150, zorder=5, 
                        label='Original points', edgecolors='black', linewidths=2)
            
            # Highlight original throat
            plt.scatter([z_throat_exact], [r_throat_exact], c='darkgreen', s=400, 
                        marker='X', zorder=7, edgecolors='black', linewidths=2,
                        label='Original throat (exact)')
            
            # Plot resampled points
            plt.plot(z_new, r_new, 'b-', label='Resampled contour', 
                    linewidth=2, marker='o', markersize=4)
            
            # Highlight chamber section
            plt.plot(z_chamber_new, r_chamber_new, 'cyan', linewidth=3, 
                    label='Chamber (resampled)', alpha=0.7)
            
            # Mark special locations
            plt.axvline(x=z_new[throat_idx], color='green', linestyle='--', 
                        linewidth=2, label=f'Throat (idx={throat_idx})')
            plt.axvline(x=z_new[chamber_end_idx], color='orange', linestyle='--', 
                        linewidth=2, label=f'Chamber end (idx={chamber_end_idx})')
            
            # Mark resampled throat point
            plt.scatter([z_new[throat_idx]], [r_new[throat_idx]], c='lime', 
                        s=300, marker='*', zorder=6, edgecolors='black', linewidths=1.5,
                        label='Resampled throat')
            
            plt.xlabel('Axial position (z)', fontsize=12)
            plt.ylabel('Radius (r)', fontsize=12)
            plt.title(f'Nozzle Contour Resampling - Throat Preserved Exactly', fontsize=14)
            plt.legend(fontsize=10, loc='best')
            plt.grid(True, alpha=0.3)
            plt.axis('equal')
            plt.tight_layout()
            plt.show()
        
        return z_new, r_new, throat_idx, chamber_end_idx