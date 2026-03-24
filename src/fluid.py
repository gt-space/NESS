from CoolProp.CoolProp import PropsSI

class Fluid:
  
    def __init__(
        self,
        prop1,
        val1,
        prop2,
        val2,
        fluid
    ):
      self.prop1 = prop1
      self.val1 = val1
      self.prop2 = prop2
      self.val2 = val2
      self.fluid = fluid

      self.calcProps()

    def calcProps(self):
      self.rho = PropsSI("D", self.prop1, self.val1, self.prop2, self.val2, self.fluid) # [kg/m^3]
      self.Cp = PropsSI("C", self.prop1, self.val1, self.prop2, self.val2, self.fluid) # [J/kg/K]
      self.mu = PropsSI("V", self.prop1, self.val1, self.prop2, self.val2, self.fluid) # [Pa*s]
      self.k = PropsSI("L", self.prop1, self.val1, self.prop2, self.val2, self.fluid) # [W/m/K]
      self.Pr = PropsSI("PRANDTL", self.prop1, self.val1, self.prop2, self.val2, self.fluid)
    

