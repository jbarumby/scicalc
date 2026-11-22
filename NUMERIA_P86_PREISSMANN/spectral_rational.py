"""P86: finite spectral audit of Preissmann's weighted Hilbert inequality.

For a finite rational ordered configuration x_i, let
    delta_i = min_{j!=i}|x_i-x_j|,
    K_ij = 1/(x_i-x_j) for i!=j; K_ii=0.
    B=D^(1/2) K D^(1/2), so B is real skew-symmetric.
Then norm(B)^2 < rho is EQUIVALENT to positive definiteness of
    H = rho*D^(-1) + K*D*K,
which has RATIONAL entries. We certify H by exact Fraction-based LDL.
This is rigorous for EACH explicit finite point set, NOT a uniform
proof for arbitrary N (the missing content of Preissmann's lemmas).

Additionally, norm(B)<=N-1 is proved uniformly for every real node
configuration of size N<=5; since 4<pi*m0, this supplies a genuine
but low-dimensional theorem independent of Preissmann.

Original source: E. Preissmann, L'Enseignement Mathématique 30
(1984), 95--113 (E-Periodica), lemmas 1--6 and pp. 111--113.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse
import json

SCHEMA="NUMERIA.P86.PREISSMANN.SPECTRAL.RATIONAL.1"
SOURCES={
 "Preissmann_1984":"https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29",
 "Ramaré_Lemma6_5":"https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf",
 "Carneiro_Littmann_2023":"https://ems.press/content/serial-article-files/44548",
}
PI_LOWER=F(157,50)
M0_LOWER=F(263,200)
RHO=(PI_LOWER*M0_LOWER)**2
PAIRS={
 "two_nodes":["0","1"],
 "three_uniform":["0","1","2"],
 "five_uniform":["0","1","2","3","4"],
 "nine_triangular":[str(i*(i+1)//2) for i in range(9)],
 "twelve_uniform":[str(i) for i in range(12)],
 "twelve_squares":[str(i*i) for i in range(12)],
 "clustered_six":["0","1/1000","1","100","1000001/10000","1000"],
 "irregular_seven":["0","1/7","2/3","10/3","100/9","99/2","300"],
}

def source_m0_pi():
    """Exact lower bounds for pi and m0, no floating point."""
    # Machin identity + alternating arctan: pi > 16(1/5-1/375)-4/239.
    pi_lb=16*(F(1,5)-F(1,375))-F(4,239)
    if not pi_lb>PI_LOWER:
        raise ArithmeticError("pi lower bound incorrect")
    lo=(F(3,2)*(M0_LOWER**2-1))
    if not (lo>0 and lo**2<F(6,5) and 4<PI_LOWER*M0_LOWER):
        raise ArithmeticError("m0 lower bound or small-N inequality incorrect")
    return {"pi_lower_strict":str(PI_LOWER),
            "pi_Machin_rational_lower":str(pi_lb),
            "m0_lower_strict":str(M0_LOWER),
            "pi_m0_lower_rational":str(PI_LOWER*M0_LOWER),
            "rho_strict_below_(pi_m0)^2":str(RHO)}

def rational_nodes(values):
    nodes=tuple(F(x) for x in values)
    if len(nodes)<2 or len(set(nodes))!=len(nodes):
        raise ValueError("at least two distinct nodes required")
    return nodes

def gaps(nodes):
    return tuple(min(abs(x-y) for y in nodes if y!=x) for x in nodes)

def mat_mul(A,B):
    if not A or not B or len(A[0])!=len(B):
        raise ValueError("invalid matrix dimensions")
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0))
             for j in range(len(B[0]))] for i in range(len(A))]

def build(nodes,threshold=RHO):
    x=rational_nodes(nodes)
    n=len(x)
    d=gaps(x)
    K=[[F(0) if i==j else F(1)/(x[i]-x[j])
        for j in range(n)] for i in range(n)]
    if any(K[i][j]!=-K[j][i] for i in range(n) for j in range(n)):
        raise ArithmeticError("K must be skew-symmetric")
    # A=K*D is similar to B = sqrt(D)*K*sqrt(D).
    A=[[K[i][j]*d[j] for j in range(n)] for i in range(n)]
    A2=mat_mul(A,A)
    # H = (rho I+A^2)*D^(-1) = rho D^(-1)+K D K,
    # because A^2=K D K D. This is symmetric.
    H=[[ (F(threshold,d[i]) if i==j else F(0)) +
         sum((K[i][k]*d[k]*K[k][j] for k in range(n)),F(0))
        for j in range(n)] for i in range(n)]
    if any(H[i][j]!=H[j][i] for i in range(n) for j in range(n)):
        raise ArithmeticError("H must be symmetric")
    for i in range(n):
        for j in range(n):
            if H[i][j]*d[j] != (F(threshold) if i==j else F(0))+A2[i][j]:
                raise ArithmeticError("rational B^2 similarity identity failed")
    return x,d,K,A,H

def exact_ldl_positive(H):
    """No machine floats, no numerical tolerance; negative/zero pivots fail."""
    n=len(H)
    L=[[F(0)]*n for _ in range(n)]
    d=[F(0)]*n
    for i in range(n):
        L[i][i]=F(1)
        for j in range(i):
            L[i][j]=(H[i][j]-sum((L[i][k]*d[k]*L[j][k]
                                   for k in range(j)),F(0)))/d[j]
        d[i]=H[i][i]-sum((L[i][k]**2*d[k] for k in range(i)),F(0))
        if d[i]<=0:
            raise ArithmeticError(f"non-positive exact pivot at {i}")
    for i in range(n):
        for j in range(n):
            v=sum((L[i][k]*d[k]*L[j][k]
                   for k in range(min(i,j)+1)),F(0))
            if H[i][j]!=v:
                raise ArithmeticError("LDL exact reconstruction failed")
    return d

def row_energies(nodes):
    """Exact rational values Q_i=sum_(j!=i)delta_i delta_j/(xi-xj)^2.

    Uniform proof, for ALL REAL nodes in README:
    on the right, nearest-neighbor contribution <=1, and every
    further term telescopes <=1; likewise on the left. Hence Q_i<=4.
    """
    x=rational_nodes(nodes)
    d=gaps(x)
    total=[]
    for i in range(len(x)):
        q=sum((d[i]*d[j]/(x[i]-x[j])**2
               for j in range(len(x)) if j!=i),F(0))
        if not 0<q<=4:
            raise ArithmeticError("general row energy <=4 failed for rational sample")
        total.append(q)
    return total

def example(name,values):
    nodes,d,K,A,H=build(values)
    pivots=exact_ldl_positive(H)
    energies=row_energies(nodes)
    # These are exactly verified witnesses of a finite-configuration
    # inequality stronger than norm(B)<pi*m0.
    raw=";".join(str(x) for x in nodes).encode()
    return {
        "name":name,
        "n":len(nodes),
        "nodes_exact_rational":[str(x) for x in nodes],
        "delta_exact_rational":[str(x) for x in d],
        "nodes_sha256":sha256(raw).hexdigest(),
        "rho_lower_than_pi2_m02_exact":str(RHO),
        "ldl_pivots_exact":[str(x) for x in pivots],
        "row_energies_exact":[str(x) for x in energies],
        "row_energy_le_4_every_index":all(q<=4 for q in energies),
        "min_pivot_sign_verified":all(v>0 for v in pivots),
        "spectral_norm_strictly_below_pi_m0_for_these_nodes":True,
    }

def run():
    constants=source_m0_pi()
    samples=[example(name,pts) for name,pts in PAIRS.items()]
    if len(samples)!=len(PAIRS):
        raise ArithmeticError("incomplete finite sample audit")
    return {
      "schema":SCHEMA,
      "sources":SOURCES,
      "constants":constants,
      "finite_samples":samples,
      "samples_count":len(samples),
      "verified_sizes":[v["n"] for v in samples],
      "spectral_equivalence":
          "||sqrt(D) K sqrt(D)||_2^2 < rho iff rho*D^(-1)+K D K is positive definite",
      "proof_of_equivalence":"similarity A=KD and congruence by sqrt(D)",
      "general_all_N_row_energy_le_4":
          "For ANY distinct real x and ANY N>=2, sum_(j!=i) delta_i*delta_j/(xi-xj)^2 <= 4 by two telescoping one-sided series",
      "general_all_N_energy_le_4_proved_independently":True,
      "Preissmann_sharper_pi2_over3_row_energy_proved_here":False,
      "general_uniform_small_n_le_5":
          "for any distinct REAL x_1..x_N with N<=5: |B_ij|<=1, ||B||_2 <= N-1 <= 4 < pi*m0",
      "general_all_N_Preissmann_proved_here":False,
      "Preissmann_Lemma1_general_spacing_maximizer_proved_here":False,
      "Preissmann_Lemmas2_to_5_global_spacing_perturbations_proved_here":False,
      "Preissmann_Lemma6_three_point_denominator_bound_proved_here":False,
      "Preissmann_final_uniform_spectral_estimate_proved_here":False,
      "finite_sample_certificates_are_not_a_proof_for_all_N":True,
      "RZ_double_sum_bound_reproved":False,
      "X_lower_3e10_3e12_covered":False,
      "RH_proved":False,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    opts=parser.parse_args()
    data=run()
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":"FINITE_RATIONAL_LDL_CLOSED",
                      "sample_sizes":data["verified_sizes"],
                      "general_Preissmann_proved":False,
                      "small_N_le5_all_configurations_proved":True},indent=2))

if __name__=="__main__":
    main()
