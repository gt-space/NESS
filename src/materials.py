import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
class Material:

    def __init__(self, material_name):
        self.material_name = material_name

        if self.material_name == "Inconel 718":
            self.rho = 8190 # kg/m^3
            self.E = 160E9 # Pa
            self.k = 11.4 # W/mK
            self.Cp = 435 # J/kgK
            self.MeltingPoint = 1533 # K
            self.v = 0.34
            self.alpha = 24.8E-6 # m/m/C
            #self.YTS = 
            #self.UTS = 

            # --- Structural / fatigue data (used by structural.py) ---
            # PLACEHOLDER: approximate wrought, aged Inconel 718 values. Replace with your own test or
            # vendor data before trusting any life number.
            self.set_strength_tables(
                T=[293.0, 811.0, 922.0, 1033.0, 1144.0],           # K
                sigma_y_MPa=[1100.0, 1000.0, 960.0, 900.0, 700.0],
                sigma_uts_MPa=[1300.0, 1170.0, 1100.0, 1000.0, 800.0],
            )
            self.eps_f = 0.35        # true fracture ductility ln(1/(1-RA)) [-]
            self.n_cyc = 0.10        # cyclic strain hardening exponent n' [-]
            self.cyc_yield_ratio = 1.0   # cyclic yield / monotonic yield [-]
            self.structural_data_verified = False

        elif self.material_name == "Pure Copper":
            self.rho = 8960 # kg/m^3
            self.E = 110E9 # Pa
            self.k = 401 # W/mK
            self.Cp = 385 # J/kgK
            self.MeltingPoint = 1350 # K
            self.v = 0.34
            self.alpha = 24.8E-6 # m/m/C

        elif self.material_name == "AlSi10Mg":
            self.rho = 2670 # kg/m^3
            self.k = 140 # W/mK
            self.E = 70E9 # Pa
            self.alpha = 20.E-6 # m/m/C
            self.v = 0.35
            
            # Yield Stress vs. Temperature Lookup Table
            T = np.array([298.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15]) # K
            sigma_y = np.array([29.587752, 28.717524, 26.251878, 26.396916, 22.916004, 19.145016, 10.15266, 4.35114, 1.740456]) # ksi
            self.stress_interp = PchipInterpolator(T, sigma_y)

        elif self.material_name == "GRCop-42":
            self.rho = 8960 # kg/m^3 - UPDATE
            self.E = 110E9 # Pa - UPDATE
            self.k = 344 # W/mK - UPDATE
            self.Cp = 385 # J/kgK - UPDATE
            self.MeltingPoint = 1350 # K - UPDATE
            self.v = 0.34 # - UPDATE
            self.alpha = 24.8E-6 # m/m/C - UPDATE

            # --- Structural / fatigue data (used by structural.py) ---
            # PLACEHOLDER: rough GRCop-42 strength-vs-temperature trend. Replace with your own test or
            # NASA/vendor data (as-built vs HIP'd LPBF differ) before trusting any life number.
            self.set_strength_tables(
                T=[293.0, 573.0, 773.0, 973.0, 1073.0],            # K
                sigma_y_MPa=[230.0, 190.0, 150.0, 90.0, 60.0],
                sigma_uts_MPa=[390.0, 330.0, 250.0, 150.0, 100.0],
            )
            self.eps_f = 0.5         # true fracture ductility ln(1/(1-RA)) [-]
            self.n_cyc = 0.15        # cyclic strain hardening exponent n' [-]
            self.cyc_yield_ratio = 1.0   # cyclic yield / monotonic yield [-]
            self.structural_data_verified = False

        '''
        elif self.material_name == "316L SS":
            self.rho = 
            self.k = 
            self.YS = 
            self.E = 

        elif self.material_name == "CuCrZr":
            self.rho = 
            self.k = 
            self.YS = 
            self.E = 

        elif self.material_name == "GrCop-42":
            self.rho = 
            self.k = 
            self.YS = 
            self.E = 
        Also GrCop-84, 
    '''

    def set_strength_tables(self, T, sigma_y_MPa, sigma_uts_MPa):
        '''
        Temperature-dependent monotonic yield and ultimate strength tables (T in K, strengths in MPa).
        '''
        self._strength_T = np.asarray(T, dtype=float)
        self._sigma_y_tab = np.asarray(sigma_y_MPa, dtype=float) * 1e6    # [Pa]
        self._sigma_uts_tab = np.asarray(sigma_uts_MPa, dtype=float) * 1e6  # [Pa]

    def structural_props(self, T):
        '''
        Properties needed by the stress / lifecycle analysis at temperature T [K].

        Strengths are linearly interpolated from the tables set by set_strength_tables (held flat beyond
        the ends of the table); E, alpha, v are the constants for this material. Anything the material
        does not define comes back as None so the caller can report exactly what is missing.

        Optional Coffin-Manson-Basquin overrides (sigma_f, b, eps_f_prime, c) may be set on the material
        as a dict `fatigue`; otherwise structural.py derives them from UTS, E and eps_f with Manson's
        Universal Slopes.
        '''
        props = {
            "E": getattr(self, "E", None),
            "alpha": getattr(self, "alpha", None),
            "v": getattr(self, "v", None),
            "k": getattr(self, "k", None),
            "eps_f": getattr(self, "eps_f", None),
            "n_cyc": getattr(self, "n_cyc", None),
            "cyc_yield_ratio": getattr(self, "cyc_yield_ratio", 1.0),
            "fatigue": getattr(self, "fatigue", None),
            "sigma_y": None,
            "sigma_uts": None,
            "in_table": True,
        }
        if hasattr(self, "_strength_T"):
            props["sigma_y"] = float(np.interp(T, self._strength_T, self._sigma_y_tab))
            props["sigma_uts"] = float(np.interp(T, self._strength_T, self._sigma_uts_tab))
            props["in_table"] = bool(self._strength_T[0] <= T <= self._strength_T[-1])
        elif hasattr(self, "stress_interp"):
            props["sigma_y"] = float(self.stress_interp(T)) * 6894757.293168
        return props

    def update_material_properties(self, T):
        '''
        Update all material properties based on material name and temperature.
        '''
        # Temperature-dependent properties are optional. If a material only has
        # hardcoded constants, leave values unchanged and simply return current properties.
        if hasattr(self, "stress_interp"):
            self.sigma_y_ksi = float(self.stress_interp(T))
            self.sigma_y = self.sigma_y_ksi * 6894757.293168  # [ksi] -> [Pa]

        props = {}
        for attr in [
            "rho", "k", "Cp", "E", "v", "alpha",
            "sigma_y", "sigma_y_ksi", "MeltingPoint"
        ]:
            if hasattr(self, attr):
                props[attr] = getattr(self, attr)

        return props

        
    def get_yield_temp(self, temp):
        
       return stress_interp(temp)

    def plot_material_props(self):
        
        # Yield Stress vs. Temp
        T_fine = np.linspace(T.min(), T.max(), 500)
        stress_fine = self.stress_interp(T_fine)

        # Yield Stress vs. Temp Plot
        plt.figure()
        plt.plot(T_fine, stress_fine)
        plt.xlabel("Temperature (K)")
        plt.ylabel("Yield Stress (ksi)")
        plt.title("AlSi10Mg Yield Stress vs. Temp")
        plt.show()


'''
T_fine = np.linspace(T.min(), T.max(), 500)
stress_fine = stress_interp(T_fine)

stress_at_500K = stress_interp(500)
print(f"Stress @ 500 K : {stress_at_500K}")

plt.figure()
plt.plot(T_fine, stress_fine)
plt.xlabel("Temperature (K)")
plt.ylabel("Yield Stress (ksi)")
plt.title("AlSi10Mg Yield Stress vs. Temp")
plt.show()
'''