class Material:

    def __init__(self, material_name):
        self.material_name = material_name

        if self.material_name == "Inconel 718":
            self.rho = 8190 # kg/m^3
            self.k = 11.4 # W/mK
            self.Cp = 0.435 # J/kgK
            self.MeltingPoint = 1533 # K
            #self.YTS =
            #self.UTS = 
            #self.E = 
        elif self.material_name == "Pure Copper":
            self.rho = 8960 # kg/m^3
            self.k = 401 # W/mK
            self.Cp = 385 # J/kgK
            self.MeltingPoint = 1350 # K

'''
        elif self.material_name == "316L SS":
            self.rho = 
            self.k = 
            self.YS = 
            self.E = 
        
        elif self.material_name == "AlSi10Mg":
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