
## --- ENGINE.PY --- ###
# Define an engine with the class Engine to size the contour and calculate performance parameters.

# Brief Description - Engine Class holds all geometry, performance parameters, and gas properties along the 
# engine. RocketIsp used for all performance calculations.

import matplotlib.pyplot as plt
import numpy as np
from conical import conicalContour

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

    def __init__(self, thrust, Pc, Pe, MRcore, oxName, fuName, name, CR, Lstar, cstarEff, pcentFFC = 0.0, Pamb = 14.8, numPts = 200, verbose = False,
                 bell = True, chmbR=None, chmbL=None, contAngle=None, throatR=None, throatL=None, expAngle=None, exitR=None):
        self.thrust = thrust
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
        self.createEngine(verbose=verbose)
        self.calcGasProperties()

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

        self.R.geomObj.getNozObj()
        self.Contour_z = [-self.R.geomObj.Lcham] + self.R.geomObj.nozObj.abs_zContour
        self.Contour_r = [self.R.geomObj.Rinj] + self.R.geomObj.nozObj.abs_rContour

        if self.bell == False:
            my_conical = conicalContour(self.chmbR, self.chmbL, self.contAngle, self.throatR, self.throatL, self.expAngle, self.exitR, self.numPts)
            self.Contour_z, self.Contour_r = my_conical.makeContour()
        
        self.Contour_z, self.Contour_r, self.throat_ind, self.chamberEnd_ind = self.resample_nozzle_contour(
        z_contour = self.Contour_z, r_contour = self.Contour_r, numPts = self.numPts, verbose=verbose)        
        
        print(f"Chamber End Ind: {self.chamberEnd_ind}")
        print(f"Throat Ind: {self.throat_ind}")
        
        self.At = self.R.geomObj.At     
        self.Dt = np.sqrt((4 / np.pi) * self.At) * 0.0254 # [in] --> [m]
        self.eps = self.R.geomObj.eps
        self.r_ch = self.Contour_r[0]
        self.r_up_throat = self.R.geomObj.RupThroat * 0.0254 # [in] --> [m]
        self.fu_mdot = self.R.coreObj.wdotFl * 0.453592
        self.ox_mdot = self.R.coreObj.wdotOx * 0.453592

    def calcGasProperties(self, frozen=1):
        # Area Ratio (A / A*) - DONE
        # Mach Number - DONE
        # P0 - DONE
        # P - DONE
        # T0 - DONE
        # T - DONE
        # Rho0 - DONE
        # Rho - DONE
        # Specific Enthalpy (h)
        # Specific Heat constant Pressure (Cp) - DONE
        # Specific Heat constant Volume (Cv)
        # Gamma - DONE
        # MW - DONE
        # Speed of Sound (a)
        # Viscosity (mu) - DONE
        # Thermal Conductivity (k) - DONE
        # Prandtl Number (Pr) - DONE

        # Initialize some Arrays
        self.AreaRatio = np.array([]) # - DONE
        self.M = np.array([]) # - DONE
        self.P = np.array([]) # - DONE
        self.T = np.array([]) # - DONE
        self.rho = np.array([]) # - DONE
        #self.h = np.array([])
        #self.Cv = np.array([])
        #self.a = np.array([])

        # Setup CEA Object, general calculations
        C = CEA_Obj(
            oxName=self.oxName, fuelName=self.fuName,
            isp_units='sec',      # Isp in seconds (default)
            cstar_units='m/s',    # C* in m/s (default is ft/sec)
            pressure_units='psia',  # Chamber pressure units
            temperature_units='K', # Temperature units
            sonic_velocity_units='m/s',
            enthalpy_units='J/kg',
            density_units='kg/m^3',
            specific_heat_units='J/kg-K',
            viscosity_units='poise',
            thermal_cond_units="W/cm-degC",
            )
        T0 = C.get_Tcomb(self.Pc, self.MRcore)
        P0 = self.Pc
        rho0 = C.get_Chamber_Density(self.Pc, self.MRcore)
        self.cstar = C.get_Cstar(self.Pc, self.MRcore)

        # Calculate Area Ratio
        for i in range(len(self.Contour_r)):
            # Calculate Propertiese
            AreaRat = (np.pi * self.Contour_r[i]**2) / self.At
            self.AreaRatio = np.append(self.AreaRatio, AreaRat)

        # --- Chamber Gas Properties Extrapolation --- 
        Cp_chamb, mu_chamb, k_chamb, Pr_chamb = C.get_Chamber_Transport(self.Pc, self.MRcore)
        MW_chamb, gamma_chamb = C.get_Chamber_MolWt_gamma(self.Pc, self.MRcore)

        # Assign Properties
        self.Cp = np.full(self.chamberEnd_ind + 1, Cp_chamb)
        self.mu = np.full(self.chamberEnd_ind + 1, mu_chamb)
        self.k = np.full(self.chamberEnd_ind + 1, k_chamb)
        self.Pr = np.full(self.chamberEnd_ind + 1, Pr_chamb)
        self.MW = np.full(self.chamberEnd_ind + 1, MW_chamb)
        self.gamma = np.full(self.chamberEnd_ind + 1, gamma_chamb)

        # --- Converging Section Gas Properties (Linear Correlation) --- 
        # Throat Properties
        Cp_t, mu_t, k_t, Pr_t = C.get_Throat_Transport(self.Pc, self.MRcore)
        MW_t, gamma_t = C.get_Throat_MolWt_gamma(self.Pc, self.MRcore)
        
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
        Cp, mu, k, Pr = C.get_Throat_Transport(self.Pc, self.MRcore)
        MW, gamma = C.get_Throat_MolWt_gamma(self.Pc, self.MRcore)
        
        # Assign Properties
        self.Cp = np.append(self.Cp, Cp)
        self.mu = np.append(self.mu, mu)
        self.k = np.append(self.k, k)
        self.Pr = np.append(self.Pr, Pr)
        self.MW = np.append(self.MW, MW)
        self.gamma = np.append(self.gamma, gamma)

        # --- Diverging Section Gas Properties --- 
        for eps in self.AreaRatio[self.throat_ind + 1 :]:
            Cp, mu, k, Pr = C.get_Exit_Transport(self.Pc, self.MRcore, eps=eps)
            MW, gamma = C.get_exit_MolWt_gamma(self.Pc, self.MRcore, eps=eps)

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
        self.Pc = self.Pc * 6894.76 # [psia] --> [Pa]

    def exportGeometry(self, filename):
        # Export the nozzle contour as [X, Y] coordinates in a .txt file for import into CAD.
        
        x = self.Contour_r
        y = [-n for n in self.Contour_z]

        with open(filename + '.txt', "w") as f:
            for xi, yi in zip(x, y):
                f.write(f"{xi:.6f}in\t{yi:.6f}in\t0in\n")

    def solveAreaMach(self, tol=1e-8):
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
            if i > 0:
                dA = self.AreaRatio[i] - self.AreaRatio[i-1]
                is_diverging = dA > 0
            else:
                is_diverging = False  # First point assumed converging
            
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
        
        # Separate chamber (first 2 points) from nozzle (remaining points)
        z_chamber = z_contour[:2]
        r_chamber = r_contour[:2]
        r_chamber_avg = np.mean(r_chamber)
        
        z_nozzle = z_contour[1:]
        r_nozzle = r_contour[1:]
        
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