import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Constants & Specs ---
THRUST = 500.0        # N (~3300 lbs)
I_LOCAL = np.diag([15.0, 15.0, 5.0]) # Engine Inertia Tensor
DT = 0.01

# Parasitic Loads
FEEDLINE_FORCE_MAG = 133.4  # N (~30 lbs)
FEEDLINE_LOCAL_PT = np.array([0.1, 0.0, -0.2]) # Feedline attach point
MU_GIMBAL = 0.15            # Friction coefficient
R_PIN = 0.025               # 25mm gimbal pin radius

# Actuator Geometry (with offsets)
# Points where actuators attach to the rocket body
ACT_BODY_PTS = np.array([
    [0.4, 0.1, 0.2],   # Act 1 body mount
    [-0.1, 0.4, 0.2]   # Act 2 body mount
])
# Points where actuators attach to the engine frame (local)
ACT_ENG_LOCAL_PTS = np.array([
    [0.3, 0.1, -0.6],  # Act 1 engine lug
    [-0.1, 0.3, -0.6]  # Act 2 engine lug
])

def get_R(a, b):
    """Full 3D Rotation Matrix (Pitch alpha, Yaw beta)."""
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    return Ry @ Rx

def dynamics(t, state):
    # state = [alpha, beta, wx, wy, wz]
    alpha, beta, wx, wy, wz = state
    omega = np.array([wx, wy, wz])
    R = get_R(alpha, beta)
    I_body = R @ I_LOCAL @ R.T
    
    # 1. Thrust Torque
    r_cp = R @ np.array([0, 0, -1.5])
    f_thrust = R @ np.array([0, 0, -THRUST])
    tau_thrust = np.cross(r_cp, f_thrust)
    
    # 2. Feedline Resistance (Horizontal always in Body Frame)
    # Acting at the transformed connection point
    r_feed = R @ FEEDLINE_LOCAL_PT
    f_feed = np.array([-FEEDLINE_FORCE_MAG, 0, 0]) # Resisting X-axis tilt
    tau_feed = np.cross(r_feed, f_feed)
    
    # 3. Gimbal Friction Torque
    # Opposes angular velocity. magnitude = T * mu * r_pin
    if np.linalg.norm(omega) > 1e-5:
        tau_fric = - (omega / np.linalg.norm(omega)) * (THRUST * MU_GIMBAL * R_PIN)
    else:
        tau_fric = np.zeros(3)
        
    # 4. Actuator Torque (PD Control example)
    # Note: In a real sim, you'd calculate actual linear force vectors here
    kp, kd = 2500.0, 600.0
    tau_ctrl = np.array([
        kp * (np.radians(10) - alpha) - kd * wx, # Target 10 deg Pitch
        kp * (0 - beta) - kd * wy,
        0
    ])

    # Total Torque
    tau_total = tau_thrust + tau_feed + tau_fric + tau_ctrl
    
    # Newton-Euler: dw/dt = I^-1 * (tau - w x Iw)
    dw_dt = np.linalg.inv(I_body) @ (tau_total - np.cross(omega, I_body @ omega))
    
    return [wx, wy, dw_dt[0], dw_dt[1], dw_dt[2]]

# Run Simulation
sol = solve_ivp(dynamics, (0, 3), [0,0,0,0,0], t_eval=np.linspace(0, 3, 300))

# Plot Results
plt.figure(figsize=(10, 5))
plt.plot(sol.t, np.degrees(sol.y[0]), label='Pitch (deg)')
plt.plot(sol.t, np.degrees(sol.y[1]), label='Yaw (deg)')
plt.axhline(10, color='r', linestyle='--', alpha=0.5, label='Target')
plt.title("TVC Response with Feedline Force & Gimbal Friction")
plt.xlabel("Time (s)"); plt.ylabel("Angle (deg)"); plt.legend(); plt.grid(True)
plt.show()