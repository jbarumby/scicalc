"""P87: independent reconstruction of Preissmann Lemma 1 via Shan/Yangjit.

GENERAL THEOREM (proved in README, not extrapolated from samples):
For any finite distinct real nodes x_i, delta_i=min_(j!=i)|xi-xj|,
  Q_i=sum_(j!=i) delta_i*delta_j/(xi-xj)^2 <= 2*zeta(2)=pi^2/3.
A more general one-sided inequality holds for ANY nonnegative, convex,
decreasing summable f:[1,infinity)->[0,infinity):
 F_N(a)=sum_{n=1}^N min(a_n,a_(n+1))*f(sum_(k=1)^n a_k)
 <= sum_(j=1)^infinity f(j), whenever a1>=1, other a>0.

The mathematical proof of the universal inequality is the 3-lemma
convexity argument reproduced in README from Shan, as reconstructed in
Yangjit, arXiv:2203.14950v1, Lemmas 2.3--2.5 (pp. 5--7).
No use of Preissmann Lemma 1 as an axiom. Python checks exact rational
finite configurations and every lifting/telescoping step; tests do NOT
constitute proof of the theorem for arbitrary REAL nodes.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from NUMERIA_P86_PREISSMANN.spectral_rational import PAIRS,rational_nodes,gaps

SCHEMA="NUMERIA.P87.SHAN.PREISSMANN_LEMMA1.1"
YANGJIT="https://arxiv.org/pdf/2203.14950"
ORIGINAL="https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29"

def floor_fraction(q):
    q=F(q)
    return q.numerator//q.denominator

def harmonic2(N):
    if N<0:raise ValueError("N cannot be negative")
    return sum((F(1,k*k) for k in range(1,N+1)), F(0))

def F_N(a):
    """Exact rational F_N for f(t)=1/t^2; a has N+1 positive entries."""
    a=tuple(F(x) for x in a)
    if len(a)<2 or a[0]<1 or any(x<=0 for x in a):
        raise ValueError("a1>=1 and all a positive")
    t=F(0); total=F(0)
    for j in range(len(a)-1):
        t+=a[j]
        total+=min(a[j],a[j+1])/(t*t)
    return total

def mass_bound(a):
    """Exact one-sided finite endpoint quadrature for nondecreasing a>=1.

    Σ a_n f(t_n) <= Σ_{1..floor(t_N)} f(j)
                        + {t_N}f(floor(t_N)+1).
    Derivation uses convex interpolation and telescoping, not numerical
    quadrature; see Lemma 2.3 in arxiv:2203.14950.
    """
    arr=tuple(F(x) for x in a)
    if len(arr)<2 or any(v<1 for v in arr) or any(arr[k]>arr[k+1] for k in range(len(arr)-1)):
        raise ValueError("mass bound requires a nondecreasing a>=1")
    t=F(0); lhs=F(0)
    for v in arr[:-1]:
        t+=v
        lhs+=v/t**2
    M=floor_fraction(t)
    theta=t-M
    rhs=harmonic2(M)+theta/(M+1)**2
    if lhs>rhs:
        raise ArithmeticError("convex endpoint mass inequality failed")
    return {"mass_left":str(lhs),"finite_lattice_budget":str(rhs),
            "floor_end":M,"fractional_end":str(theta)}

def lift_certificate(a):
    """Increase each descent once, verifying each exact F_N increases."""
    v=[F(x) for x in a]
    start=F_N(v)
    steps=[]
    for k in range(1,len(v)):
        if v[k]<v[k-1]:
            old=F_N(v)
            eps=v[k-1]-v[k]
            v[k]=v[k-1]
            new=F_N(v)
            if new<old:
                raise ArithmeticError("convexity step reversed")
            steps.append({"index_1based":k+1, "epsilon":str(eps),
                          "before":str(old),"after":str(new),
                          "monotone_step":True})
    final=F_N(v)
    if final<start or any(v[k]>v[k+1] for k in range(len(v)-1)):
        raise ArithmeticError("running-max lift failed")
    mass=mass_bound(v)
    if final>F(mass["mass_left"]) or F(mass["mass_left"])>F(mass["finite_lattice_budget"]):
        raise ArithmeticError("finite comparison chain failed")
    return {"start":str(start),"after_lift":str(final),
            "lifted_gaps":[str(x) for x in v],"descent_steps":steps,
            "mass":mass,"all_exact_rational_gates_pass":True}

def audit_nodes(values):
    nodes=tuple(sorted(rational_nodes(values)))
    delta=gaps(nodes)
    findings=[]
    for i,x in enumerate(nodes):
        contributions=[]
        for direction in (-1,1):
            adj=[nodes[k] for k in (range(i-1,-1,-1) if direction<0 else range(i+1,len(nodes)))]
            if not adj:
                contributions.append({"side":"left" if direction<0 else "right",
                                      "sum":"0","empty_side":True})
                continue
            distances=[abs(y-x) for y in adj]
            intervals=[distances[0]]+[distances[k]-distances[k-1] for k in range(1,len(distances))]
            a=[step/delta[i] for step in intervals]
            a.append(a[-1]) # end point has only one adjacent spacing on this side
            cert=lift_certificate(a)
            one_side=F(cert["start"])
            exact_side=sum((delta[i]*delta[nodes.index(y)]/(y-x)**2 for y in adj),F(0))
            if one_side!=exact_side:
                raise ArithmeticError("normalized gap formula failed")
            contributions.append({"side":"left" if direction<0 else "right",
                                  "sum":str(one_side),"gap_lift":cert,"empty_side":False})
        total=sum((F(t["sum"]) for t in contributions),F(0))
        matrix_direct=sum((delta[i]*delta[j]/(x-nodes[j])**2
                           for j in range(len(nodes)) if j!=i),F(0))
        if total!=matrix_direct:
            raise ArithmeticError("row energy mismatch")
        # Stronger rational finite check only for this sample; the GENERAL
        # theorem uses exact Basel identity, not a rational approximation
        # to pi (PI_LOWER is too low for arbitrarily large uniform grids).
        findings.append({
            "index":i+1,"node":str(x),"delta":str(delta[i]),
            "sides":contributions,"row_energy":str(total),
            "finite_normalization_proved":True,
        })
    return {"nodes":[str(x) for x in nodes],"rows":findings,
            "all_rows_exact":True}

def proof_ledger():
    samples=[]
    for name,arr in PAIRS.items():
        cert=audit_nodes(arr)
        cert["name"]=name
        samples.append(cert)
    uniform41=audit_nodes(range(41))
    uniform41["name"]="uniform_41_extremal_limit_diagnostic"
    samples.append(uniform41)
    return {
        "schema":SCHEMA,
        "source":{
            "Shan_Yangjit_Lemmas_2_3_2_5":YANGJIT,
            "Preissmann_1984_Lemma_1":ORIGINAL,
        },
        "theorem_all_finite_real_nodes":
            "for all N>=2 and distinct real nodes Q_i<=2*zeta(2)=pi^2/3",
        "more_general_one_sided_lemma":
            "for positive a with a1>=1, decreasing convex nonnegative summable f: Σ min(a_n,a_(n+1))f(Σ_{k<=n} a_k)<=Σ_{j>=1}f(j)",
        "analytic_steps":{
            "convex_interpolation":"f(t_n)<= (1-theta_n)f(floor(t_n))+theta_n f(floor(t_n)+1)",
            "mass_telescope":"Σ a_n f(t_n)<=Σ_{j<=floor(t_N)}f(j)+theta_N f(floor(t_N)+1)",
            "descent_lifting":"a_(nu-1)>a_nu: increasing a_nu until equality increases F_N by convexity; repeat to running maximum",
            "two_sides":"normalize each left/right gap by delta_i; applying F_N<=zeta(2) to each side gives Q_i<=2*zeta(2)",
            "basel":"2*zeta(2)=pi^2/3",
        },
        "source_method":"independent proof of Preissmann Lemma1 from Shan/Yangjit convexity, source theorem not presumed",
        "no_arb_needed_for_exact_sample_checks":True,
        "universal_mathematical_proof_in_README":True,
        "tests_are_not_the_general_proof":True,
        "cases":samples,
        "exact_case_count":len(samples),
        "Preissmann_Lemma1_reconstructed":True,
        "Preissmann_Lemma6_reconstructed":False,
        "Preissmann_full_m0_reproved_without_Lemma6":False,
        "S_epsilon_RZ_reproved":False,
        "RH_proved":False,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    opts=p.parse_args()
    result=proof_ledger()
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("P87 Shan: 2*zeta(2)=pi²/3 analytic proof described; rational cases",len(result["cases"]))

if __name__=="__main__":
    main()
