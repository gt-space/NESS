import sympy as sp
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj


def _cea_obj_w_units_like_engine(oxName, fuelName):
    """
    RocketCEA CEA_Obj with unit settings aligned to Engine.calcGasProperties.
    Pc is interpreted as psia; Tcomb returned in K.
    """
    return CEA_Obj(
        oxName=oxName,
        fuelName=fuelName,
        isp_units="sec",
        cstar_units="m/s",
        pressure_units="psia",
        temperature_units="K",
        sonic_velocity_units="m/s",
        enthalpy_units="J/kg",
        density_units="kg/m^3",
        specific_heat_units="J/kg-K",
        viscosity_units="poise",
        thermal_cond_units="W/cm-degC",
    )


def cea_chamber_tcomb(oxName, fuelName, Pc_psia, MR):
    """Single-point chamber (combustion) temperature from CEA [K]."""
    C = _cea_obj_w_units_like_engine(oxName, fuelName)
    return float(C.get_Tcomb(float(Pc_psia), float(MR)))


def cea_chamber_tcomb_mr_sweep(oxName, fuelName, Pc_psia, mr_min, mr_max, n_points, mr_highlight=None):
    """
    Sweep Tcomb vs mixture ratio using RocketCEA only (no Engine / RocketIsp).

    Parameters
    ----------
    oxName, fuelName : str
        Propellant names for CEA (same as Engine).
    Pc_psia : float
        Chamber pressure [psia].
    mr_min, mr_max : float
        O/F sweep endpoints (inclusive linspace endpoints).
    n_points : int
        Number of MR samples.
    mr_highlight : float, optional
        If given, also returns Tcomb at this MR (e.g. design MRcore) using the same CEA object.

    Returns
    -------
    mr_axis : ndarray
    tcomb_axis : ndarray
    tcomb_highlight : float or None
        Tcomb at mr_highlight, or None if mr_highlight is None.
    """
    if mr_min >= mr_max:
        raise ValueError("mr_min must be < mr_max")
    n_points = int(n_points)
    if n_points < 2:
        raise ValueError("n_points must be >= 2")

    C = _cea_obj_w_units_like_engine(oxName, fuelName)
    Pc_psia = float(Pc_psia)
    mr_axis = np.linspace(float(mr_min), float(mr_max), n_points)
    tcomb_axis = np.array([float(C.get_Tcomb(Pc_psia, float(mr))) for mr in mr_axis], dtype=float)

    if mr_highlight is None:
        return mr_axis, tcomb_axis, None
    tcomb_highlight = float(C.get_Tcomb(Pc_psia, float(mr_highlight)))
    return mr_axis, tcomb_axis, tcomb_highlight


def solve_system(equations, variables, inputs):
        """
        Solve a system of equations with flexible inputs.

        equations : list of sympy Eq
        variables : list of sympy symbols
        inputs : dict {symbol: value or None}
                    - If value is None, it's treated as unknown
                    - If value is numeric, it's treated as given

        Returns: list of dict solutions (each includes all variables)
        """
        # Drop None values (unknowns)
        inputs = {k: v for k, v in inputs.items() if v is not None}

        # Figure out unknowns
        unknowns = [v for v in variables if v not in inputs]

        # Solve system for unknowns
        sol = sp.solve(equations, unknowns, dict=True)

        if not sol:
            raise ValueError("No solution found.")

        results = []
        for s in sol:
            # Substitute known values into the symbolic solution
            evaluated = {k: v.subs(inputs) for k, v in s.items()}
            # Merge known + solved variables
            results.append({**inputs, **evaluated})

        results_r = results[0]
        return results_r

def check_defined_vars(
    inputs : dict,
    min_required : int,
        
) -> bool:
    # Checks if there are the required number of defined variables.
    
    defined_count = sum(value is not None for value in inputs.values())

    return defined_count >= min_required