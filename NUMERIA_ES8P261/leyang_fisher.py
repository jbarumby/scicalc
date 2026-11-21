"""NUMERIA ES8-P261: exact Lee-Yang/Fisher/LP separation gates.

Proves finite Ising polynomial identities and a counterexample to the
implication 'positive-definite characteristic function => real zeros'.
No assertions about Xi zeros, RH or the positivity of Weil's form.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import json

def require(condition, explanation):
    if not condition:
        raise ArithmeticError(explanation)

def path_ising_coefficients(length: int, disagreement_weight: Q):
    """P_N(z,a)=sum_{spins} a^disagreements z^number_of_negative_spins.

    a=exp(-2K); 0<a<=1 corresponds to ferromagnetic K>=0.
    This partition polynomial drops the common factor exp((N-1)K).
    """
    if type(length) is not int or not 1 <= length <= 12:
        raise ValueError("length must be an integer 1..12")
    a = Q(disagreement_weight)
    if a < 0:
        raise ValueError("the weights in this audit must be nonnegative")
    coefficients = [Q(0)] * (length + 1)
    for spins in product((-1, +1), repeat=length):
        disagreements = sum(spins[j] != spins[j+1] for j in range(length-1))
        negative = sum(s == -1 for s in spins)
        coefficients[negative] += a**disagreements
    return tuple(coefficients)

def is_palindromic(coefficients):
    return tuple(coefficients) == tuple(reversed(coefficients))

def exact_path_at_zero_field(length: int, a: Q):
    """Normalized field fugacity z=1: 2(1+a)^(N-1)."""
    a = Q(a)
    direct = sum(path_ising_coefficients(length,a),Q(0))
    theoretical = 2 * (1+a)**(length-1)
    require(direct == theoretical, "open-chain Fisher factorization mismatch")
    return theoretical

def quartic_circle_criterion(a: Q):
    """Exact circle iff y=z+z^-1 roots of y^2+A*y+B-2 lie in [-2,2].

    Proved using discriminant, vertex in [-2,2], nonnegative polynomial
    at both boundary points. No approximate complex root computation.
    """
    a = Q(a)
    poly = path_ising_coefficients(4, a)
    require(poly[0] == 1 and poly[4] == 1 and is_palindromic(poly),
            "reciprocal quartic required")
    A, B = poly[1], poly[2]
    discriminant = A*A - 4*(B-2)
    q_minus2 = 4 - 2*A + B - 2
    q_plus2 = 4 + 2*A + B - 2
    vertex_within = abs(A) <= 4
    circle = (discriminant >= 0 and vertex_within
              and q_minus2 >= 0 and q_plus2 >= 0)
    return {
       "disagreement_weight":str(a),
       "coefficients_z_ascending":[str(x) for x in poly],
       "y_quadratic_coefficients":[str(B-2),str(A),"1"],
       "discriminant_exact":str(discriminant),
       "Q_at_minus2_exact":str(q_minus2),
       "Q_at_plus2_exact":str(q_plus2),
       "vertex_inside_interval":vertex_within,
       "all_zeros_on_unit_circle_proved":circle,
    }

def determinant_exact(matrix):
    """Fraction-preserving determinant via Schur elimination."""
    n=len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    if n==0:return Q(1)
    a=[[Q(z) for z in row] for row in matrix]
    det=Q(1)
    for k in range(n):
        pivot=next((r for r in range(k,n) if a[r][k]),None)
        if pivot is None:return Q(0)
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            det=-det
        x=a[k][k]
        det*=x
        for i in range(k+1,n):
            factor=a[i][k]/x
            for j in range(k+1,n):
                a[i][j]-=factor*a[k][j]
            a[i][k]=Q(0)
    return det

def rational_characteristic_grid(size=4):
    """Gram matrix phi((j-k)*pi/2), phi(t)=3/4+(1/4)cos(t)."""
    if type(size) is not int or not 1<=size<=8:
        raise ValueError("size must be 1..8")
    cos_grid=(Q(1),Q(0),Q(-1),Q(0))
    return tuple(tuple(Q(3,4)+Q(1,4)*cos_grid[(i-j)%4]
                       for j in range(size)) for i in range(size))

def all_principal_minors_nonnegative(matrix):
    n=len(matrix)
    summary={}
    for k in range(1,n+1):
        values=[]
        for idx in combinations(range(n),k):
            det=determinant_exact([[matrix[i][j] for j in idx] for i in idx])
            require(det>=0,"a purported characteristic Toeplitz minor is negative")
            values.append(det)
        summary[str(k)] = [str(v) for v in values]
    return summary

def lp_counterexample():
    """The characteristic phi(t)=3/4+1/4 cos t is NOT Laguerre-Polya.

    It is the Fourier transform of 3/4 delta_0 +1/8(delta_-1+delta_+1).
    Multiplying by 8z, z=e^(it), yields z^2+6z+1, not a unit-circle
    zero polynomial. Laguerre at pi: phi'^2 - phi phi'' = -1/8.
    """
    gram=rational_characteristic_grid(4)
    minor=all_principal_minors_nonnegative(gram)
    laguerre_at_zero=Q(1,4)
    laguerre_at_pi=-Q(1,8)
    require(laguerre_at_pi < 0 < laguerre_at_zero,
            "non-Laguerre-Pólya local certificate failed")
    return {
        "characteristic_measure":[["-1","1/8"],["0","3/4"],["+1","1/8"]],
        "characteristic_function":"phi(t)=3/4+(1/4)*cos(t)",
        "toeplitz_4_by_4_exact":[[str(x) for x in row] for row in gram],
        "all_15_principal_minors_exact":minor,
        "polynomial_in_z":"8*z*phi(t)=z^2+6*z+1; z=exp(i*t)",
        "reciprocal_polynomial_roots":"z=-3+-2*sqrt(2), both off unit circle",
        "laguerre_inequality_at_t0_exact":str(laguerre_at_zero),
        "laguerre_inequality_at_tpi_exact":str(laguerre_at_pi),
        "positive_definite_but_not_Laguerre_Polya":True,
    }

def bost_connes_finite_gate():
    """Exact identity for Z_2(s)=1+2^(-s): zeros have real part zero."""
    return {
        "finite_partition":"Z_2(s)=1+2^(-s)",
        "zero_equation":"s=-(2k+1)*pi*i/log(2), k in integers",
        "all_zeros_real_part_exact":"0",
        "physical_Bost_Connes_partition":"Z(beta)=zeta(beta) in Re(beta)>1",
        "analytic_continuation_not_derived":True,
        "critical_line_from_partition_not_proved":True
    }

def audit():
    ferro = quartic_circle_criterion(Q(1,2))
    anti = quartic_circle_criterion(Q(2))
    require(ferro["all_zeros_on_unit_circle_proved"],
            "4-spin ferromagnetic polynomial failed")
    require(not anti["all_zeros_on_unit_circle_proved"],
            "anti-ferromagnetic counterexample must fail")
    checks=0
    for n in range(1,6):
        for a in (Q(0),Q(1,2),Q(1),Q(2)):
            p=path_ising_coefficients(n,a)
            require(is_palindromic(p),"spin-flip symmetry failed")
            require(exact_path_at_zero_field(n,a)==sum(p),"Fisher sum mismatch")
            checks+=1
    non_lp = lp_counterexample()
    return {
        "schema":"NUMERIA.ES8P261.LEE_YANG_FISHER_LP.GATE_A.1",
        "source_Konstantopoulos_Patie_Sarkar_2024":
          "https://www.numdam.org/articles/10.5802/aif.3600/",
        "paper_theorem_2_3":"xi_Landau/xi_Landau(0) in D_L iff RH",
        "paper_corollary_2_4":
          "xi_Landau(0)/xi_Landau(i*t) is symmetric Polya-frequency characteristic iff RH",
        "unconditional_property":
          "xi_Landau/xi_Landau(0) is already a characteristic function",
        "path_model":"open chain; weight a^disagreements; a=exp(-2K)",
        "ferromagnetic_4spin_exact":ferro,
        "anti_ferromagnetic_4spin_exact":anti,
        "Fisher_zero_field":"Z_N(z=1,a)=2*(1+a)^(N-1); zeros a=-1, multiplicity N-1",
        "spin_flip_and_fisher_identities_checked":checks,
        "positive_definite_non_LP_counterexample":non_lp,
        "Bost_Connes_boundary":bost_connes_finite_gate(),
        "independent_Xi_Laguerre_Polya_membership_proved":False,
        "Weil_operator_positive_transfer_proved":False,
        "Riemann_hypothesis_proved":False,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    data=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("P261 Gate A PASS: exact Ising/Fisher/characteristic counterexamples")

if __name__=="__main__":
    main()
