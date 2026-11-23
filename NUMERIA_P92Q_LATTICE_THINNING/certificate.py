"""P92-Q — variational triangular tents on a thinned integer lattice.

The theorem for all finite subsets and all positive real mesh sizes is proved
in README.md. Fraction checks here are regression, NOT the universal proof.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import json
import random

SCHEMA="NUMERIA.P92Q.PREISSMANN.LATTICE_THINNING.TENT_PACKING.1"
Z0=F(0)
ONE=F(1)


def interior_nearest_gaps(M, nodes):
    """Nearest-neighbor distance in mesh units, including endpoints 0,M."""
    if not isinstance(M,int) or M<1:
        raise ValueError("M must be a positive integer")
    nodes=tuple(sorted(nodes))
    if any(not isinstance(j,int) or not 0<j<M for j in nodes):
        raise ValueError("interior nodes must be integer mesh positions")
    if len(nodes)!=len(set(nodes)):
        raise ValueError("repeated nodes")
    centers=(0,)+nodes+(M,)
    return {j:min(j-centers[i-1],centers[i+1]-j)
            for i,j in enumerate(centers) if i>0 and i<len(centers)-1}


def tent(j,d,k):
    """Tent supported at [j-d,j+d], triangular and rational."""
    if d<=0:
        raise ValueError("positive width required")
    return max(Z0,ONE-F(abs(k-j),d))


def mesh_kernel(M,k):
    if not 0<k<M:
        raise ValueError("kernel only defined strictly between endpoints")
    return F(1,k*k*(M-k)*(M-k))


def harmonic(n,p=1):
    if n<0:
        raise ValueError("negative index")
    return sum((F(1,j**p) for j in range(1,n+1)),Z0)


def D_mesh(M):
    """Exact D(1/M) from P90, harmonic identity valid for every integer M>=1."""
    if M<1: raise ValueError("M>=1 required")
    return M*harmonic(M,2)+2*harmonic(M)-F(3,M)


def exact_case(M,nodes, details=False):
    gaps=interior_nearest_gaps(M,nodes)
    interior=sum((d*mesh_kernel(M,j) for j,d in gaps.items()),Z0)
    complete=sum((mesh_kernel(M,k) for k in range(1,M)),Z0)
    complete_via_harmonics=F(2,M*M)*harmonic(M-1,2)+F(4,M**3)*harmonic(M-1)
    if complete!=complete_via_harmonics:
        raise ArithmeticError("harmonic partial fractions disagree")
    if M**3*complete!=2*D_mesh(M):
        raise ArithmeticError("sharp uniform mesh identity failed")
    occupancy=[]
    local_jensen=[]
    for k in range(1,M):
        total=sum((tent(j,d,k) for j,d in gaps.items()),Z0)
        if not Z0<=total<=ONE:
            raise ArithmeticError("overlapping tents exceed mass one")
        occupancy.append(total)
    for j,d in gaps.items():
        mass=sum((tent(j,d,k) for k in range(1,M)),Z0)
        moment=sum((F(k)*tent(j,d,k) for k in range(1,M)),Z0)
        cost=sum((tent(j,d,k)*mesh_kernel(M,k) for k in range(1,M)),Z0)
        if mass!=d or moment!=j*d or cost<d*mesh_kernel(M,j):
            raise ArithmeticError("mass, barycentre, Jensen failed")
        local_jensen.append(cost-d*mesh_kernel(M,j))
    if interior>complete:
        raise ArithmeticError("global tent packing comparison failed")
    return {"mesh_intervals":M,"interior_nodes":list(sorted(nodes)),
            "number_of_interior_nodes":len(gaps),
            "minimum_gap_scale_units":{str(k):v for k,v in gaps.items()},
            "normalized_interior_sum_exact":str(M**3*interior),
            "normalized_complete_sum_exact":str(M**3*complete),
            "thinning_margin_exact":str(M**3*(complete-interior)),
            "twice_D_of_1_over_M_exact":str(2*D_mesh(M)),
            "min_unused_slot_mass_exact":str(min((ONE-o for o in occupancy),default=ONE)),
            "total_Jensen_slack_exact":str(sum(local_jensen,Z0))}


def atan_lower(x,n):
    """Leibniz alternating bounds for arctan on 0<x<1."""
    x=F(x)
    if not Z0<x<ONE or n<1:
        raise ValueError("atan fraction in (0,1), n>=1")
    partial=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),Z0)
    nxt=x**(2*n+1)/F(2*n+1)
    return partial if n%2==0 else partial-nxt


def atan_upper(x,n):
    x=F(x)
    partial=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),Z0)
    nxt=x**(2*n+1)/F(2*n+1)
    return partial+nxt if n%2==0 else partial


def pi_lower():
    """Machin: pi/4=4 atan(1/5)-atan(1/239), exact lower bound."""
    return 4*(4*atan_lower(F(1,5),16)-atan_upper(F(1,239),8))


def full_L6_left_right(xs,s,t):
    """Calculate both L6 sides with exact rational arithmetic, pi replaced by lower enclosure."""
    x=[F(y) for y in xs]
    if any(x[i]>=x[i+1] for i in range(len(x)-1)):
        raise ValueError("points must be strictly increasing")
    d=[min(abs(x[i]-x[j]) for j in range(len(x)) if j!=i)
       for i in range(len(x))]
    h=x[t]-x[s]
    if h<=0:raise ValueError("indices out of order")
    lhs=sum((d[k]/((x[k]-x[s])**2*(x[k]-x[t])**2)
             for k in range(len(x)) if k not in (s,t)),Z0)
    K=pi_lower()**2/3
    rhs=K/(h*h)*(1/d[s]+1/d[t])-3*(d[s]+d[t])/h**4
    return lhs,rhs,d


def certify():
    rng=random.Random(9212026)
    tests=[]
    widths=(2,3,4,5,8,11,17,26,40,60,89,120)
    for M in widths:
        allslots=list(range(1,M))
        options=[
            (),
            tuple(allslots),
            tuple(j for j in allslots if j%2==0),
            tuple(j for j in allslots if j%3==1),
            tuple(j for j in allslots if j<=min(2,M-1) or j>=max(1,M-2)),
            tuple(sorted(rng.sample(allslots,rng.randrange(0,M)))),
            tuple(sorted(rng.sample(allslots,rng.randrange(0,M))))
        ]
        for nodes in options:
            tests.append(exact_case(M,nodes))
    # A highly irregular configuration with isolated points and 3 multi-node islands
    adversarial=(1,2,5,19,20,46,73,74,99)
    tests.append(exact_case(100,adversarial))
    # Direct checks of the sharp Lema 6 with exterior points and endpoints delta<=a.
    full=[]
    for M,nodes in ((5,(2,)),(9,(1,3,6,8)),(12,(2,7,9)),
                    (34,(1,2,17,30,33)),(100,adversarial)):
        coords=[F(-1),F(0)]+[F(j) for j in nodes]+[F(M),F(M+1)]
        lhs,rhs,d=full_L6_left_right(coords,1,len(coords)-2)
        if not lhs<rhs or d[1]>1 or d[-2]>1:
            raise ArithmeticError("failed exact full L6 on sparse grid with arbitrary external tails")
        full.append({"M":M,"nodes":list(nodes),
                     "full_L6_margin_with_pi_lower_exact":str(rhs-lhs)})
    assert pi_lower()>F(3)
    return {
       "schema":SCHEMA,
       "theorem":"All finite thinned arithmetic grids, arbitrary interior holes, endpoint nearest gaps <= mesh step: full sharp Preissmann L6",
       "proof":"triangular tent capacity at each mesh site; discrete Jensen for convex two-pole kernel; P89 exterior bound and P90 D monotonicity + exact harmonic identity",
       "analytic_universal_for_stated_lattice_class":True,
       "general_Lemma5_and_6_for_arbitrary_real_positions":False,
       "general_Riemann_Hypothesis_proved":False,
       "exact_rational_regressions":len(tests),
       "full_L6_checks_with_rational_pi_lower":len(full),
       "full_checks":full,
       "cases":tests,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    x=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92Q PASS",x["exact_rational_regressions"],"Fraction tent cases and",
          x["full_L6_checks_with_rational_pi_lower"],"full L6 cases")
    print("SHA256:",hashlib.sha256(args.out.read_bytes()).hexdigest())
    print("No general Lemmas 5/6 or RH claim.")


if __name__=="__main__":
    main()
