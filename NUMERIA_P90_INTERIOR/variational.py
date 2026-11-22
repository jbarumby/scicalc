"""P90 exact local interior variational derivatives for one movable point.

Only the one-node interior function is studied, not every term of the
original Preissmann Lemma5 objective. The derivative changes sign: blindly
'straightening' or equalizing a gap does NOT monotonically change this term.
All identities are exact Fraction statements for rational positions.
"""
from fractions import Fraction as F

def one_point_profile(t):
    t=F(t)
    if not 0<t<1:
        raise ValueError("interior point must be between endpoints 0 and 1")
    if t<=F(1,2):
        value=F(1)/(t*(1-t)**2)
        # d/dt value = (3t-1)/[t²(1-t)³]
        derivative=(3*t-1)/(t*t*(1-t)**3)
        branch="left_gap_is_nearest"
    else:
        value=F(1)/(t*t*(1-t))
        # d/dt value = (3t-2)/[t³(1-t)²]
        derivative=(3*t-2)/(t**3*(1-t)**2)
        branch="right_gap_is_nearest"
    return {"t":str(t),"kernel_exact":str(value),
            "piecewise_derivative_exact":str(derivative),
            "branch":branch,
            "derivative_sign":1 if derivative>0 else -1 if derivative<0 else 0}

def compare_profiles():
    pts=[F(1,6),F(1,4),F(1,3),F(2,5),
         F(1,2),F(3,5),F(2,3),F(3,4),F(5,6)]
    data=[one_point_profile(t) for t in pts]
    signs=[x["derivative_sign"] for x in data]
    if signs!=[-1,-1,0,1,1,-1,0,1,1]:
        raise ArithmeticError("the one-node variational derivative pattern changed")
    # In this single-term objective, 1/3 and 2/3 are local minima;
    # 1/2 is a sharp local maximum, not a universally monotone direction.
    return {
        "piecewise_exact_profiles":data,
        "one_interior_derivative_changes_sign":True,
        "stationary_positions":["1/3","2/3"],
        "cusp_maximum_position":"1/2",
        "full_Preissmann_Lemma5_variational_objective_reconstructed":False,
        "all_multinode_interior_variations_reproved":False
    }
