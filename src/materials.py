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
            self.k = 130 # W/mK
            self.E = 70E9 # Pa
            
            # Yield Stress vs. Temperature Lookup Table
            T = np.array([298.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15]) # K
            sigma_y = np.array([29.587752, 28.717524, 26.251878, 26.396916, 22.916004, 19.145016, 10.15266, 4.35114, 1.740456]) # ksi
            self.stress_interp = PchipInterpolator(T, sigma_y)

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