import numpy as np

# Capabilities
# 1) 3D Position Simulation of TVC
# 2) Plots
    # Gimbal Angle
    # TVC1/2 Angle + Force
    # Angular Acceleration
    # Angular Velocity

# Position Inputs
gimbal = np.array([0, 0]) # in
actuator = np.array([4, 0]) # in

# Load Inputs
thrust = 500 # lbf
mu = 0.6

# Static Loads
F_feedline = 30 # lbf
F_fric = 