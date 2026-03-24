import sympy as sp

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