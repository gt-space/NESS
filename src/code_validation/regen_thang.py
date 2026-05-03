import os
import re
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import CEA_Obj as CEA_Obj_Imperial

# ============================================================
# COMBINED ENGINE PERFORMANCE + REGENERATIVE COOLING MODEL
#
# FLOW:
#   1. solve_engine_performance() computes chamber pressure,
#      thrust, burn duration, scaled combustion temperature, etc.
#   2. run_regen_model() uses that solved chamber pressure and
#      shared inputs for the thermal model
#
# NOTE:
#   m_dot_fuel and m_dot_c are taken to be the same
# ============================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# -----------------------------
# SHARED DEFAULT INPUTS
# -----------------------------
DEFAULTS = {
    # propellants
    'ox_name': 'LOX',
    'fuel_name': 'RP1',

    # main engine operating inputs
    'm_dot_fuel': 1.6,                 # kg/s, also used as coolant mass flow
    'MR': 2.0,
    'eps': 4.8,
    'eta_cstar': 0.90,

    # throat / nozzle geometry
    'D_t': 3.00 * 0.0254,              # m

    # tank / burn sizing
    'num_fuel_tanks': 2,
    'num_LOX_tanks': 2,
    'fuel_tank_volume_L_each': 60.0,   # L
    'LOX_tank_volume_L_each': 60.0,    # L
    'rho_fuel': 815.0,                 # kg/m^3
    'rho_LOX': 1141.0,                 # kg/m^3

    # ambient
    'p_amb_psia': 14.696,

    # coolant
    'fluid': 'n-Dodecane',
    'T_c_in': 300.0,                   # K
    'P_c_in': 400 * 6895.0,            # Pa

    # chamber geometry
    'L_chamber': 0.204,                # m

    # contour files
    'inner_contour_file': os.path.join(BASE_DIR, 'inner_contour.sldcrv'),
    'outer_contour_file': os.path.join(BASE_DIR, 'outer_contour.sldcrv'),

    # wall fallback
    't_wall': 0.001,                   # m
    'k_wall': 344.0,                   # W/m-K

    # channel geometry
    'N_channels': 126,
    'channel_width': 1.0e-3,           # m
    'channel_height': 3.0e-3,          # m

    # roughness
    'roughness': 0.0,                  # m

    # coolant correlation
    'coolant_method': 'GN',            # 'GN' or 'DB'

    # soot / deposit model
    'use_soot_deposit': True,

    # mesh
    'N': 500,

    # boundary layer sizing
    'target_y_plus': 1.0,
    'N_blayers': 18,
    'growth_rate': 1.15,
}

# -----------------------------
# CONSTANTS
# -----------------------------
psi_to_Pa = 6894.75729
ft_to_m = 0.3048
lbf_to_N = 4.4482216152605
R_to_K = 5.0 / 9.0


# ============================================================
# PERFORMANCE MODEL
# ============================================================
def solve_engine_performance(p, max_iter=100, tol_psia=1e-4):
    m_dot_fuel = p['m_dot_fuel']
    MR = p['MR']
    D_t = p['D_t']
    eps = p['eps']
    eta_cstar = p['eta_cstar']

    if m_dot_fuel <= 0.0 or MR <= 0.0:
        raise ValueError('m_dot_fuel and MR must be > 0.')

    A_t = np.pi * (D_t / 2.0) ** 2

    m_dot_ox = MR * m_dot_fuel
    m_dot_total = m_dot_fuel + m_dot_ox

    m_fuel_tank = (
        p['fuel_tank_volume_L_each'] * p['num_fuel_tanks'] / 1000.0 * p['rho_fuel']
    )
    m_LOX_tank = (
        p['LOX_tank_volume_L_each'] * p['num_LOX_tanks'] / 1000.0 * p['rho_LOX']
    )

    burn_duration = min(m_fuel_tank / m_dot_fuel, m_LOX_tank / m_dot_ox)

    p_amb_psia = p['p_amb_psia']
    p_amb_Pa = p_amb_psia * psi_to_Pa

    ispObj = CEA_Obj_Imperial(oxName=p['ox_name'], fuelName=p['fuel_name'])

    pc_psia = 300.0

    for _ in range(max_iter):
        c_star_ideal = ispObj.get_Cstar(Pc=pc_psia, MR=MR) * ft_to_m
        c_star_real = eta_cstar * c_star_ideal

        pc_new_Pa = c_star_real * m_dot_total / A_t
        pc_new_psia = pc_new_Pa / psi_to_Pa

        if abs(pc_new_psia - pc_psia) < tol_psia:
            pc_psia = pc_new_psia
            break

        pc_psia = pc_new_psia

    pc_Pa = pc_psia * psi_to_Pa

    T_comb_ideal_R = ispObj.get_Tcomb(Pc=pc_psia, MR=MR)
    T_comb_ideal_K = T_comb_ideal_R * R_to_K
    T_comb_scaled_K = T_comb_ideal_K * eta_cstar**2

    M_exit = ispObj.get_MachNumber(
        Pc=pc_psia, MR=MR, eps=eps, frozen=1, frozenAtThroat=1
    )

    _, gam_throat = ispObj.get_Throat_MolWt_gamma(Pc=pc_psia, MR=MR, eps=eps)

    p_exit_psia = pc_psia * (
        1.0 + (gam_throat - 1.0) / 2.0 * M_exit**2
    ) ** (-gam_throat / (gam_throat - 1.0))
    p_exit_Pa = p_exit_psia * psi_to_Pa

    a_exit_ideal = ispObj.get_SonicVelocities(
        Pc=pc_psia, MR=MR, eps=eps, frozen=1, frozenAtThroat=1
    )[-1] * ft_to_m

    a_exit_scaled = a_exit_ideal * np.sqrt(T_comb_scaled_K / T_comb_ideal_K)
    v_exit = M_exit * a_exit_scaled

    A_exit = eps * A_t
    F_momentum = m_dot_total * v_exit
    F_pressure = (p_exit_Pa - p_amb_Pa) * A_exit
    F_total = F_momentum + F_pressure

    total_impulse = F_total * burn_duration

    return {
        'm_dot_fuel': m_dot_fuel,
        'm_dot_c': m_dot_fuel,
        'm_dot_ox': m_dot_ox,
        'm_dot_total': m_dot_total,
        'MR': MR,
        'eps': eps,
        'A_t': A_t,
        'Pc_psia': pc_psia,
        'Pc_Pa': pc_Pa,
        'Pc_psig': pc_psia - p_amb_psia,
        'p_exit_psia': p_exit_psia,
        'p_exit_Pa': p_exit_Pa,
        'T_comb_ideal_K': T_comb_ideal_K,
        'T_comb_scaled_K': T_comb_scaled_K,
        'v_exit_m_s': v_exit,
        'thrust_N': F_total,
        'thrust_lbf': F_total / lbf_to_N,
        'burn_duration_s': burn_duration,
        'total_impulse_Ns': total_impulse,
        'c_star_ideal_m_s': c_star_ideal,
        'c_star_real_m_s': c_star_real,
    }


# ============================================================
# REGEN MODEL HELPER FUNCTIONS
# ============================================================
def bartz_sigma(M, gamma, Tw, T0_ref):
    term = 1.0 + 0.5 * (gamma - 1.0) * M**2
    return (0.5 * (Tw / T0_ref) * term + 0.5) ** (-0.68) * term ** (-0.12)


def friction_factor(Re, roughness, D_h):
    if Re <= 0.0:
        return np.nan
    if Re < 2300.0:
        return 64.0 / Re
    term = (roughness / (3.7 * D_h)) ** 1.11 + 6.9 / Re
    return (-1.8 * np.log10(term)) ** (-2)


def fin_efficiency(h, k, t_fin, H_fin):
    m = np.sqrt(2.0 * h / (k * t_fin))
    return np.tanh(m * H_fin) / (m * H_fin)


def first_cell_height_from_yplus(y_plus, mu, rho, u_tau):
    return y_plus * mu / (rho * u_tau)


def friction_velocity_from_darcy_f(f_darcy, V_bulk):
    return V_bulk * np.sqrt(f_darcy / 8.0)


def inflation_thickness(y1, n_layers, r):
    if abs(r - 1.0) < 1.0e-12:
        return y1 * n_layers
    return y1 * (r**n_layers - 1.0) / (r - 1.0)


def soot_resistance_Rd(area_ratio):
    """
    Carbon deposit thermal resistance from the user-provided RP-1 soot graph.
    Returns Rd in m^2-K/W.
    area_ratio = A/A*
    """

    # convert (in^2-sec-degF)/BTU to (m^2-K/W)
    return 700 * 0.00064516*5/9/1055.06


def corrected_bartz_hg(hg_bartz, area_ratio, use_soot=True):
    """
    Apply the carbon deposit correction:
        h_gc = 1 / (1/h_g + R_d)

    Returns:
        h_gc, R_d
    """
    if not use_soot:
        return hg_bartz, 0.0

    Rd = soot_resistance_Rd(area_ratio)
    h_gc = 1.0 / (1.0 / hg_bartz + Rd)
    return h_gc, Rd


def _read_sldcrv(path):
    r_vals = []
    x_vals = []

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            nums = re.findall(r'([-+]?\d*\.?\d+)in', line)
            if len(nums) >= 2:
                r_vals.append(float(nums[0]) * 0.0254)
                x_vals.append(float(nums[1]) * 0.0254)

    if len(r_vals) < 2:
        raise ValueError(f'Could not parse contour file: {path}')

    r = np.asarray(r_vals, dtype=float)
    x = np.asarray(x_vals, dtype=float)

    order = np.argsort(x)
    x = x[order]
    r = r[order]

    return x, r


def _find_chamber_end_index(r_inner):
    r0 = r_inner[0]
    const_idx = np.where(np.isclose(r_inner, r0, rtol=0.0, atol=1.0e-9))[0]
    if len(const_idx) == 0:
        raise ValueError('Could not identify the straight chamber section in the inner contour.')
    return int(const_idx[-1])


def build_full_geometry(p):
    x_inner_raw, r_inner_raw = _read_sldcrv(p['inner_contour_file'])

    outer_file = p.get('outer_contour_file')
    if outer_file is not None and os.path.exists(outer_file):
        x_outer_raw, r_outer_raw = _read_sldcrv(outer_file)
        r_outer_raw_on_inner = np.interp(x_inner_raw, x_outer_raw, r_outer_raw)
        t_wall_raw = r_outer_raw_on_inner - r_inner_raw
    else:
        t_wall_raw = np.full_like(r_inner_raw, p['t_wall'])

    if np.any(t_wall_raw <= 0.0):
        raise ValueError('Outer contour must stay outside the inner contour everywhere.')

    i_chamber_end = _find_chamber_end_index(r_inner_raw)
    L_chamber_file = x_inner_raw[i_chamber_end]
    L_chamber_user = p['L_chamber']
    delta_L = L_chamber_user - L_chamber_file

    x_new = x_inner_raw.copy()
    x_new[:i_chamber_end + 1] = np.linspace(0.0, L_chamber_user, i_chamber_end + 1)
    x_new[i_chamber_end + 1:] = x_inner_raw[i_chamber_end + 1:] + delta_L

    if np.any(np.diff(x_new) <= 0.0):
        raise ValueError('Requested L_chamber is too short for this contour resolution.')

    x_uniform = np.linspace(x_new[0], x_new[-1], p['N'])
    r_inner = np.interp(x_uniform, x_new, r_inner_raw)
    t_wall = np.interp(x_uniform, x_new, t_wall_raw)
    r_outer = r_inner + t_wall

    D_inner = 2.0 * r_inner
    circumference_inner = 2.0 * np.pi * r_inner

    return {
        'x_raw': x_new,
        'r_inner_raw': r_inner_raw,
        't_wall_raw': t_wall_raw,
        'x': x_uniform,
        'r_inner': r_inner,
        'r_outer': r_outer,
        't_wall': t_wall,
        'D_inner': D_inner,
        'circumference_inner': circumference_inner,
        'L_total': x_uniform[-1],
        'L_chamber_file': L_chamber_file,
        'i_chamber_end_raw': i_chamber_end,
    }


# ============================================================
# REGEN MODEL
# ============================================================
def run_regen_model(overrides=None):
    p = DEFAULTS.copy()
    if overrides is not None:
        p.update(overrides)

    perf = solve_engine_performance(p)

    Pc_gas = perf['Pc_Pa']
    MR = perf['MR']
    eps = perf['eps']
    T_aw_override = perf['T_comb_scaled_K']

    fluid = p['fluid']
    m_dot_c = p['m_dot_fuel']
    T_c_in = p['T_c_in']
    P_c_in = p['P_c_in']

    D_t = p['D_t']
    k_wall = p['k_wall']
    N_channels = p['N_channels']
    channel_width = p['channel_width']
    channel_height = p['channel_height']
    roughness = p['roughness']
    coolant_method = p['coolant_method']
    use_soot = p['use_soot_deposit']
    target_y_plus = p['target_y_plus']
    N_blayers = p['N_blayers']
    growth_rate = p['growth_rate']

    geom = build_full_geometry(p)
    x = geom['x']
    dx = x[1] - x[0]
    D_local = geom['D_inner']
    t_wall_local = geom['t_wall']
    circumference_local = geom['circumference_inner']
    D_chamber_ref = D_local[0]

    pitch_local = circumference_local / N_channels
    fin_thickness_local = pitch_local - channel_width

    if np.any(fin_thickness_local <= 0.0):
        i_bad = int(np.argmin(fin_thickness_local))
        raise ValueError(
            f'Computed fin thickness is <= 0 somewhere. First failing station x={x[i_bad]:.6f} m'
        )

    cea = CEA_Obj(
        oxName=p['ox_name'],
        fuelName=p['fuel_name'],
        pressure_units='Pa',
        temperature_units='K',
        cstar_units='m/s',
        specific_heat_units='J/kg-K',
        viscosity_units='poise',
        thermal_cond_units='W/cm-degC'
    )

    c_star_ideal = cea.get_Cstar(Pc=Pc_gas, MR=MR)
    T0_ideal = cea.get_Tcomb(Pc=Pc_gas, MR=MR)
    mw, gamma_g = cea.get_Chamber_MolWt_gamma(Pc=Pc_gas, MR=MR, eps=eps)
    cp_g, mu_g_poise, k_g_W_per_cmK, Pr_g = cea.get_Chamber_Transport(Pc=Pc_gas, MR=MR, eps=eps)

    mu_g = mu_g_poise * 0.1
    k_g = k_g_W_per_cmK * 100.0
    T_aw = T_aw_override

    A_flow_per_channel = channel_width * channel_height
    A_flow_total = N_channels * A_flow_per_channel
    P_channel = channel_width + 2.0 * channel_height
    D_h = 4.0 * A_flow_per_channel / P_channel
    A_star = np.pi * (D_t / 2.0) ** 2

    def bartz_hg_from_diameter(Tw, D_here):
        M = 0.0
        sigma = bartz_sigma(M, gamma_g, Tw, T_aw)
        hg = (
            0.026
            / (D_t**0.2)
            * ((mu_g**0.2 * cp_g) / (Pr_g**0.6))
            * ((Pc_gas / c_star_ideal)**0.8)
            * ((D_t / D_here)**1.8)
            * sigma
        )
        return hg

    N_pts = len(x)

    T_c = np.zeros(N_pts)
    P_c = np.zeros(N_pts)

    coolant_inlet_index = N_pts - 1
    coolant_outlet_index = 0

    h_g = np.zeros(N_pts)
    h_g_bartz = np.zeros(N_pts)
    h_c = np.zeros(N_pts)
    qpp = np.zeros(N_pts)
    T_w_hot = np.zeros(N_pts)
    T_w_cold = np.zeros(N_pts)
    Rd_arr = np.zeros(N_pts)
    area_ratio_arr = np.zeros(N_pts)

    Re_arr = np.zeros(N_pts)
    Pr_arr = np.zeros(N_pts)
    Nu_arr = np.zeros(N_pts)
    f_arr = np.zeros(N_pts)
    V_arr = np.zeros(N_pts)
    eta_f_arr = np.zeros(N_pts)
    P_eff_arr = np.zeros(N_pts)

    rho_arr = np.zeros(N_pts)
    mu_arr = np.zeros(N_pts)
    u_tau_arr = np.zeros(N_pts)
    y1_arr = np.zeros(N_pts)
    inflation_thickness_arr = np.zeros(N_pts)
    delta_y_wall_unit_arr = np.zeros(N_pts)

    T_c[coolant_inlet_index] = T_c_in
    P_c[coolant_inlet_index] = P_c_in

    Tw_guess = 700.0

    for march_idx in range(N_pts - 1, 0, -1):
        i = march_idx
        i_next = march_idx - 1

        Tc = T_c[i]
        Pc = P_c[i]

        rho_c = CP.PropsSI('D', 'T', Tc, 'P', Pc, fluid)
        mu_c = CP.PropsSI('V', 'T', Tc, 'P', Pc, fluid)
        cp_c = CP.PropsSI('C', 'T', Tc, 'P', Pc, fluid)
        k_c = CP.PropsSI('L', 'T', Tc, 'P', Pc, fluid)

        V_c = m_dot_c / (rho_c * A_flow_total)
        Re = rho_c * V_c * D_h / mu_c
        Pr = cp_c * mu_c / k_c

        Re_arr[i] = Re
        Pr_arr[i] = Pr
        V_arr[i] = V_c
        rho_arr[i] = rho_c
        mu_arr[i] = mu_c

        f = friction_factor(Re, roughness, D_h)
        f_arr[i] = f

        if coolant_method == 'DB':
            Nu = 0.023 * Re**0.8 * Pr**0.4
        else:
            Nu = ((f / 8.0) * (Re - 1000.0) * Pr) / (
                1.0 + 12.7 * np.sqrt(f / 8.0) * (Pr**(2.0 / 3.0) - 1.0)
            )

        Nu_arr[i] = Nu
        h_c_local = Nu * k_c / D_h
        h_c[i] = h_c_local

        eta_f = fin_efficiency(h_c_local, k_wall, fin_thickness_local[i], channel_height)
        eta_f_arr[i] = eta_f

        P_base_per_channel = channel_width
        P_fin_per_channel = 2.0 * channel_height
        P_eff = N_channels * (P_base_per_channel + eta_f * P_fin_per_channel)
        P_eff_arr[i] = P_eff

        Tw_hot_local = Tw_guess
        for _ in range(5):
            A_local = np.pi * (D_local[i] / 2.0) ** 2
            area_ratio = A_local / A_star

            h_g_bartz_local = bartz_hg_from_diameter(Tw_hot_local, D_local[i])
            h_g_corr_local, Rd_local = corrected_bartz_hg(
                h_g_bartz_local,
                area_ratio,
                use_soot=use_soot,
            )

            R_tot = 1.0 / h_g_corr_local + t_wall_local[i] / k_wall + 1.0 / h_c_local
            q_local = (T_aw - Tc) / R_tot
            Tw_hot_new = T_aw - q_local / h_g_corr_local

            if abs(Tw_hot_new - Tw_hot_local) < 1.0e-6:
                Tw_hot_local = Tw_hot_new
                break
            Tw_hot_local = Tw_hot_new

        h_g_bartz[i] = h_g_bartz_local
        h_g[i] = h_g_corr_local
        Rd_arr[i] = Rd_local
        area_ratio_arr[i] = area_ratio
        qpp[i] = q_local
        T_w_hot[i] = Tw_hot_local
        T_w_cold[i] = Tc + q_local / h_c_local

        u_tau = friction_velocity_from_darcy_f(f, V_c)
        y1 = first_cell_height_from_yplus(target_y_plus, mu_c, rho_c, u_tau)
        t_infl = inflation_thickness(y1, N_blayers, growth_rate)

        u_tau_arr[i] = u_tau
        y1_arr[i] = y1
        inflation_thickness_arr[i] = t_infl
        delta_y_wall_unit_arr[i] = first_cell_height_from_yplus(1.0, mu_c, rho_c, u_tau)

        dQ = q_local * P_eff * dx
        dT = dQ / (m_dot_c * cp_c)
        T_c[i_next] = Tc + dT

        dPdx = -f * (rho_c * V_c**2) / (2.0 * D_h)
        dP = dPdx * dx
        P_c[i_next] = max(Pc + dP, 1.0e5)

        Tw_guess = T_w_hot[i]

    for arr in [
        h_g, h_g_bartz, h_c, qpp, T_w_hot, T_w_cold, Rd_arr, area_ratio_arr,
        Re_arr, Pr_arr, Nu_arr, f_arr, V_arr, eta_f_arr, P_eff_arr,
        rho_arr, mu_arr, u_tau_arr, y1_arr, inflation_thickness_arr,
        delta_y_wall_unit_arr
    ]:
        arr[0] = arr[1]

    y1_min = np.min(y1_arr)
    y1_max = np.max(y1_arr)
    inflation_min = np.min(inflation_thickness_arr)
    inflation_max = np.max(inflation_thickness_arr)
    i_worst = np.argmin(y1_arr)

    return {
        'params': p,
        'performance': perf,
        'geometry': geom,
        'x': x,
        'T_c': T_c,
        'P_c': P_c,
        'h_g': h_g,
        'h_g_bartz': h_g_bartz,
        'h_c': h_c,
        'qpp': qpp,
        'T_w_hot': T_w_hot,
        'T_w_cold': T_w_cold,
        'T_aw': T_aw,
        'Rd': Rd_arr,
        'area_ratio': area_ratio_arr,
        'c_star_ideal': c_star_ideal,
        'T0_ideal': T0_ideal,
        'gamma_g': gamma_g,
        'mw': mw,
        'cp_g': cp_g,
        'mu_g': mu_g,
        'k_g': k_g,
        'Pr_g': Pr_g,
        'circumference': circumference_local[0],
        'pitch': pitch_local[0],
        'pitch_local': pitch_local,
        'fin_thickness': fin_thickness_local[0],
        'fin_thickness_local': fin_thickness_local,
        'D_h': D_h,
        'A_flow_total': A_flow_total,
        'D_local': D_local,
        't_wall_local': t_wall_local,
        'y1_min': y1_min,
        'y1_max': y1_max,
        'inflation_min': inflation_min,
        'inflation_max': inflation_max,
        'i_worst': i_worst,
        'Re_arr': Re_arr,
        'V_arr': V_arr,
        'f_arr': f_arr,
        'eta_f_arr': eta_f_arr,
        'max_hot_wall_temp': np.max(T_w_hot),
        'D_chamber_reference': D_chamber_ref,
        'coolant_inlet_index': coolant_inlet_index,
        'coolant_outlet_index': coolant_outlet_index,
    }


# ============================================================
# PRINTING
# ============================================================
def print_case_results(result):
    p = result['params']
    perf = result['performance']
    i_worst = result['i_worst']
    x = result['x']
    geom = result['geometry']

    print('========== PERFORMANCE ==========')
    print('m_dot_fuel [kg/s]               =', perf['m_dot_fuel'])
    print('m_dot_c [kg/s]                  =', perf['m_dot_c'])
    print('m_dot_ox [kg/s]                 =', perf['m_dot_ox'])
    print('m_dot_total [kg/s]              =', perf['m_dot_total'])
    print('MR                              =', perf['MR'])
    print('eps                             =', perf['eps'])
    print('Chamber Pressure [psia]         =', perf['Pc_psia'])
    print('Chamber Pressure [Pa]           =', perf['Pc_Pa'])
    print('Scaled Tcomb [K]                =', perf['T_comb_scaled_K'])
    print('Exit pressure [psia]            =', perf['p_exit_psia'])
    print('Exit velocity [m/s]             =', perf['v_exit_m_s'])
    print('Thrust [N]                      =', perf['thrust_N'])
    print('Thrust [lbf]                    =', perf['thrust_lbf'])
    print('Burn duration [s]               =', perf['burn_duration_s'])
    print('Total impulse [N-s]             =', perf['total_impulse_Ns'])

    print('\n========== CEA GAS-SIDE INPUTS ==========')
    print('Ideal c_star [m/s]              =', result['c_star_ideal'])
    print('eta_cstar                       =', p['eta_cstar'])
    print('Ideal T0_gas [K]                =', result['T0_ideal'])
    print('Effective T0_gas [K]            =', result['T_aw'])
    print('gamma_g                         =', result['gamma_g'])
    print('MW [kg/kmol]                    =', result['mw'])
    print('cp_g [J/kg-K]                   =', result['cp_g'])
    print('mu_g [Pa-s]                     =', result['mu_g'])
    print('k_g [W/m-K]                     =', result['k_g'])
    print('Pr_g                            =', result['Pr_g'])

    print('\n========== GEOMETRY ==========')
    print('Contour chamber length [m]      =', geom['L_chamber_file'])
    print('Requested chamber length [m]    =', p['L_chamber'])
    print('Total cooled length [m]         =', geom['L_total'])
    print('Number of channels              =', p['N_channels'])
    print('Chamber circumference [m]       =', result['circumference'])
    print('Chamber pitch [m]               =', result['pitch'])
    print('Channel width [m]               =', p['channel_width'])
    print('Chamber fin thickness [m]       =', result['fin_thickness'])
    print('Minimum fin thickness [m]       =', np.min(result['fin_thickness_local']))
    print('Maximum wall thickness [m]      =', np.max(result['t_wall_local']))
    print('Minimum wall thickness [m]      =', np.min(result['t_wall_local']))
    print('Hydraulic diameter [m]          =', result['D_h'])
    print('Total flow area [m^2]           =', result['A_flow_total'])

    print('\n========== THERMAL RESULTS ==========')
    i_in = result['coolant_inlet_index']
    i_out = result['coolant_outlet_index']
    print('Coolant inlet temperature [K]   =', result['T_c'][i_in])
    print('Coolant outlet temperature [K]  =', result['T_c'][i_out])
    print('Coolant inlet pressure [Pa]     =', result['P_c'][i_in])
    print('Coolant outlet pressure [Pa]    =', result['P_c'][i_out])
    print('Total pressure drop [Pa]        =', result['P_c'][i_in] - result['P_c'][i_out])
    print('Average coolant velocity [m/s]  =', np.mean(result['V_arr']))
    print('Average Reynolds number         =', np.mean(result['Re_arr']))
    print('Average fin efficiency          =', np.mean(result['eta_f_arr']))
    print('Average Bartz h_g [W/m^2-K]     =', np.mean(result['h_g_bartz']))
    print('Average corrected h_g [W/m^2-K] =', np.mean(result['h_g']))
    print('Average soot Rd [m^2-K/W]       =', np.mean(result['Rd']))
    print('Max heat flux [MW/m^2]          =', np.max(result['qpp']) / 1e6)
    print('Max hot wall temp [K]           =', np.max(result['T_w_hot']))
    print('Max cold wall temp [K]          =', np.max(result['T_w_cold']))

    print('\n========== BOUNDARY LAYER SIZING ==========')
    print('Target y+                       =', p['target_y_plus'])
    print('Number of BL layers             =', p['N_blayers'])
    print('Growth rate                     =', p['growth_rate'])
    print('Worst-case axial location [m]   =', x[i_worst])
    print('Worst-case Reynolds number      =', result['Re_arr'][i_worst])
    print('Worst-case bulk velocity [m/s]  =', result['V_arr'][i_worst])
    print('Worst-case friction factor      =', result['f_arr'][i_worst])

    print('\n--- Recommended first-layer height ---')
    print('Minimum required y1 [m]         =', result['y1_min'])
    print('Minimum required y1 [micron]    =', result['y1_min'] * 1e6)
    print('Maximum local y1 [m]            =', result['y1_max'])
    print('Maximum local y1 [micron]       =', result['y1_max'] * 1e6)

    print('\n--- Inflation stack thickness ---')
    print('Minimum local stack thickness [m]   =', result['inflation_min'])
    print('Minimum local stack thickness [mm]  =', result['inflation_min'] * 1e3)
    print('Maximum local stack thickness [m]   =', result['inflation_max'])
    print('Maximum local stack thickness [mm]  =', result['inflation_max'] * 1e3)

    print('\n--- SUMMARY ---')
    print(f"Coolant Outlet Temp: {result['T_c'][i_out]:.1f} K")
    print(f"Thrust: {perf['thrust_lbf']:.1f} lbf")
    print(f"Chamber Pressure: {perf['Pc_psia']:.1f} psia")

# ============================================================
# BASELINE PLOTTING
# ============================================================
def plot_temperature_vs_axial_length(result):
    x = result['x']
    plt.figure()
    plt.plot(x, result['T_c'], label='Coolant temperature')
    plt.plot(x, result['T_w_hot'], label='Hot wall temperature')
    plt.plot(x, result['T_w_cold'], label='Cold wall temperature')
    plt.axhline(result['T_aw'], linestyle='--', label='Adiabatic wall temperature')
    plt.xlabel('Axial distance [m]')
    plt.ylabel('Temperature [K]')
    plt.title('Temperature vs Axial Length')
    plt.legend()
    plt.grid(True)


def plot_heat_transfer_coefficients_vs_axial_length(result):
    x = result['x']
    plt.figure()
    plt.plot(x, result['h_g_bartz'], label='Gas-side h (Bartz)')
    plt.plot(x, result['h_g'], label='Gas-side h (corrected)')
    plt.plot(x, result['h_c'], label='Coolant-side h')
    plt.xlabel('Axial distance [m]')
    plt.ylabel('Heat transfer coefficient [W/m^2-K]')
    plt.title('Heat Transfer Coefficients vs Axial Length')
    plt.legend()
    plt.grid(True)


def plot_full_solution(result):
    x = result['x']

    plot_temperature_vs_axial_length(result)

    plt.figure()
    plt.plot(x, result['qpp'] / 1e6)
    plt.xlabel('Axial distance [m]')
    plt.ylabel('Heat flux [MW/m^2]')
    plt.title('Heat Flux vs Axial Length')
    plt.grid(True)

    plot_heat_transfer_coefficients_vs_axial_length(result)

    plt.figure()
    plt.plot(x, result['Rd'])
    plt.xlabel('Axial distance [m]')
    plt.ylabel('Soot resistance Rd [m^2-K/W]')
    plt.title('Carbon Deposit Thermal Resistance vs Axial Length')
    plt.grid(True)

    plt.figure()
    plt.plot(x, result['geometry']['r_inner'], label='Inner contour radius')
    plt.plot(x, result['geometry']['r_outer'], label='Outer contour radius')
    plt.xlabel('Axial distance [m]')
    plt.ylabel('Radius [m]')
    plt.title('Contour Radius vs Axial Length')
    plt.legend()
    plt.grid(True)

    # Channel geometry vs axial position
    channel_width_arr = result['pitch_local'] - result['fin_thickness_local']
    channel_height_arr = np.full_like(x, result['params']['channel_height'])
    fig, ax1 = plt.subplots()
    ax1.plot(x, channel_width_arr * 1e3, color='tab:blue', label='Channel width')
    ax1.plot(x, channel_height_arr * 1e3, color='tab:orange', label='Channel height')
    ax1.set_xlabel('Axial distance [m]')
    ax1.set_ylabel('Dimension [mm]')
    ax1.set_title('Channel Width and Height vs Axial Position')
    ax1.legend()
    ax1.grid(True)


# ============================================================
# SWEEP / PLOTTING FUNCTIONS
# ============================================================
def sweep_mdot_for_multiple_MR(m_dot_values, MR_values, base_overrides=None):
    results = {}

    for MR in MR_values:
        total_impulse_arr = []
        coolant_outlet_temp_arr = []
        max_hot_wall_temp_arr = []

        for m_dot in m_dot_values:
            overrides = {} if base_overrides is None else dict(base_overrides)
            overrides['m_dot_fuel'] = m_dot
            overrides['MR'] = MR

            result = run_regen_model(overrides)

            x = result['x']
            T_w_hot = result['T_w_hot']
            L = result['params']['L_chamber']

            mask = (x >= 0.0) & (x <= L)
            max_hot_wall_temp_0_to_L = np.max(T_w_hot[mask])

            total_impulse_arr.append(result['performance']['total_impulse_Ns'])
            coolant_outlet_temp_arr.append(result['T_c'][result['coolant_outlet_index']])
            max_hot_wall_temp_arr.append(max_hot_wall_temp_0_to_L)

        results[MR] = {
            'm_dot': np.array(m_dot_values, dtype=float),
            'total_impulse_Ns': np.array(total_impulse_arr, dtype=float),
            'coolant_outlet_temp_K': np.array(coolant_outlet_temp_arr, dtype=float),
            'max_hot_wall_temp_K': np.array(max_hot_wall_temp_arr, dtype=float),
        }

    return results


def plot_total_impulse_vs_mdot(results_by_MR):
    plt.figure()
    for MR, data in results_by_MR.items():
        plt.plot(data['m_dot'], data['total_impulse_Ns'], marker='o', label=f'MR = {MR:.2f}')
    plt.xlabel('Fuel / coolant mass flow rate [kg/s]')
    plt.ylabel('Total impulse [N-s]')
    plt.title('Total Impulse vs Mass Flow Rate for Various MR')
    plt.grid(True)
    plt.legend()


def plot_coolant_outlet_temp_vs_mdot(results_by_MR):
    plt.figure()
    for MR, data in results_by_MR.items():
        plt.plot(data['m_dot'], data['coolant_outlet_temp_K'], marker='o', label=f'MR = {MR:.2f}')
    plt.xlabel('Fuel / coolant mass flow rate [kg/s]')
    plt.ylabel('Coolant outlet temperature [K]')
    plt.title('Coolant Outlet Temperature vs Mass Flow Rate for Various MR')
    plt.grid(True)
    plt.legend()


def plot_max_wall_temp_vs_mdot(results_by_MR):
    plt.figure()
    for MR, data in results_by_MR.items():
        plt.plot(data['m_dot'], data['max_hot_wall_temp_K'], marker='o', label=f'MR = {MR:.2f}')
    plt.xlabel('Fuel / coolant mass flow rate [kg/s]')
    plt.ylabel('Maximum hot wall temperature [K]')
    plt.title('Maximum Hot Wall Temperature vs Mass Flow Rate for Various MR')
    plt.grid(True)
    plt.legend()


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    # -------------------------
    # BASELINE CASE
    # -------------------------
    baseline_result = run_regen_model()
    print_case_results(baseline_result)
    plot_full_solution(baseline_result)

    # -------------------------
    # SWEEP SETTINGS
    # -------------------------
    m_dot_values = np.linspace(1.4, 2.6, 20)
    MR_values = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

    # optional extra overrides applied to every sweep case
    base_overrides = None
    # example:
    # base_overrides = {
    #     'eta_cstar': 0.85,
    #     'T_c_in': 300.0,
    #     'L_chamber': 0.204,
    #     'use_soot_deposit': True,
    # }

    # -------------------------
    # RUN SWEEP
    # -------------------------
    sweep_results = sweep_mdot_for_multiple_MR(
        m_dot_values=m_dot_values,
        MR_values=MR_values,
        base_overrides=base_overrides
    )

    # -------------------------
    # SWEEP PLOTS
    # -------------------------
    plot_total_impulse_vs_mdot(sweep_results)
    plot_coolant_outlet_temp_vs_mdot(sweep_results)
    plot_max_wall_temp_vs_mdot(sweep_results)

    plt.show()
