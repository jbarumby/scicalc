"""P85: independently bracket each b4* root for printed & S2 robust B1.

Uses Arb / acb.integral at workprec=192,256,384. K1 is integrated only
ONCE per row and precision; all root tests use exact algebraic factorization
B1(b4) = D * (alpha*b4+beta)^e for K3>1. These are thresholds for the
implemented P77/P79 B1 *majorant*, not mathematical truths about zeros.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
from NUMERIA_P76_COEFFICIENTS.table import TABLE, TABLE_S2_PUBLIC_ROBUST, validate_table_s2_robust
from NUMERIA_P79_B1_PUBLIC import audit_public as p79

ROOT=Path(__file__).resolve().parents[1]
SCHEMA="NUMERIA.P85.EXACT_ARB_B4STAR_BRACKETS.1"
BITS=(192,256,384)
SCALE=10**10 # outward decimal rational grid, NOT a floating approximation

def cap(fr):
    v=F(fr)
    return arb(v.numerator)/arb(v.denominator)

def frac_from_pair(pair):
    return p79.exact_from_dyadic(pair)

def floor_fraction(fr):
    q=F(fr)
    return q.numerator//q.denominator

def ceil_fraction(fr):
    q=F(fr)
    return -((-q.numerator)//q.denominator)

def source_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def factor_row(row,cs):
    """Re-evaluate P79 upper-majorant with K3-dependent factors factored out.

    K1's upper enclosure is exactly the one selected in P79; the same
    branch 'max' structure remains for all factors independent of b4.
    """
    d,kap,delta=p79.q(row.d),p79.q(row.kappa),p79.q(row.delta)
    s1,s2=p79.rat(row.sigma_start,32),p79.rat(row.sigma_start+1,32)
    L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
    A0=(arb(10).log()*p79.q(row.log10_a0)).exp()
    if not arb(10)**10<A0<arb(p79.T0):
        raise ArithmeticError("invalid A0")
    La=A0.log()
    half=(A0/2).log()
    a1=A0/(2*pi*pi)*half**4+p79.q("34.842")*A0*half**3+arb("1.3897e16")
    r=1-arb.const_log2()/La
    a2=r**4/(2*pi*pi)+p79.q("34.842")/La*r**3+arb("9.5079e15")/(A0*La**4)
    if not (a1>0 and a2>0):
        raise ArithmeticError("invalid a1/a2")
    c=3/(4*kap)
    k1lower,intupper,tail=p79.k1_integral(c,cs["c1"],cs["c2"],L)
    coeff=p79.pw(2*a2*c.sqrt(),p79.rat(1,3))
    correction=(p79.pw(2*a1,p79.rat(1,3))*
                p79.pw(cs["c2"]+A0*cs["c1"]/p79.T0,p79.rat(2,3))/
                (p79.pw(arb(p79.T0),p79.rat(1,3))*p79.pw(L,p79.rat(4,3))))
    k1=(coeff*intupper+correction).upper()
    if not k1>0:raise ArithmeticError("invalid K1")
    two13=p79.pw(arb(2),p79.rat(1,3))
    k2=two13*c.sqrt()*pi.sqrt()/2+two13*(4*kap/3)*p79.pw(A0/arb(p79.T0),arb(3))
    coupling=pi.sqrt()/(2*(2*kap).sqrt())+pi*m0/p79.T0
    alpha=pi*m0/(2*delta)*(2*delta*arb.const_euler()/L).exp()+coupling/(5*delta*delta.exp())
    beta=coupling*p79.q(".605")*(-2*delta).exp()/L
    exp1=(2*s1-1)/(2-s1)
    exp2=(2*s2-1)/(2-s2)
    if not (alpha>0 and beta>0 and exp2>exp1>0):
        raise ArithmeticError("bad K3 affine or monotone exponent")

    l0=(kap*(1+(1/(4*arb(p79.T0))*(1+2*delta/L))**2)).exp()/(1-arb(1)/p79.T0)
    l1=max(arb(1),p79.pw(k1,3*delta/L))*(3*delta+4*delta*L.log()/L).exp()
    l2=max(arb(1),p79.pw(k2,3*delta/L))*(3*delta).exp()
    common0=max(arb(1),p79.pw(l0,-2*d/L))
    common0*=max(p79.pw(l0,2/(2-s1)),p79.pw(l0,2/(2-s2)))
    D1=l1*(d*(3+4*L.log()/L)).exp()*common0
    D1*=max(arb(1),p79.pw(k1,3*d/L))
    D1*=max(p79.pw(k1,3*(1-s1)/(2-s1)),p79.pw(k1,3*(1-s2)/(2-s2)))
    D2=l2*(3*d).exp()*common0
    D2*=max(arb(1),p79.pw(k2,3*d/L))
    D2*=max(p79.pw(k2,3*(1-s1)/(2-s1)),p79.pw(k2,3*(1-s2)/(2-s2)))
    D=(D1+D2*p79.pw(L,-6*(1-s2)/(2-s2)))/(2*pi*d)
    if not (D>0 and k1>0 and k2>0):
        raise ArithmeticError("nonpositive D, K1 or K2")
    return D,alpha,beta,exp2,{"K1_upper_dyadic":p79.dyadic(k1.upper()),
                              "K2_upper_dyadic":p79.dyadic(k2.upper()),
                              "D_lower_dyadic":p79.dyadic(D.lower()),
                              "D_upper_dyadic":p79.dyadic(D.upper()),
                              "alpha_lower_dyadic":p79.dyadic(alpha.lower()),
                              "beta_lower_dyadic":p79.dyadic(beta.lower())}

def root_bracket(D,alpha,beta,e,target):
    """Construct rational brackets, verify strict inequalities with Arb."""
    t=p79.q(target)
    exponent_inv=arb(1)/e
    root_ball=(p79.pw(t/D,exponent_inv)-beta)/alpha
    if not root_ball.is_finite() or not root_ball>0:
        raise ArithmeticError("nonpositive root or unbounded interval")
    rl=frac_from_pair(p79.dyadic(root_ball.lower()))
    rh=frac_from_pair(p79.dyadic(root_ball.upper()))
    lo=F(floor_fraction(rl*SCALE)-1,SCALE)
    hi=F(ceil_fraction(rh*SCALE)+1,SCALE)
    if not F(0)<lo<hi<F(1):
        raise ArithmeticError("root bracket outside 0..1")
    def eval_at(x):
        k3=alpha*cap(x)+beta
        if not k3>1:
            raise ArithmeticError("K3>1 branch missing at root endpoint")
        return D*p79.pw(k3,e)
    below,above=eval_at(lo),eval_at(hi)
    if not (below.upper()<t and above.lower()>t):
        raise ArithmeticError("outward root bracket failed strict interval signs")
    return {
        "threshold_B1":target,
        "root_lower_rational":str(lo),
        "root_upper_rational":str(hi),
        "root_lower_dyadic":p79.dyadic(root_ball.lower()),
        "root_upper_dyadic":p79.dyadic(root_ball.upper()),
        "bracket_width":str(hi-lo),
        "B1_at_lower_upper_dyadic":p79.dyadic(below.upper()),
        "B1_at_upper_lower_dyadic":p79.dyadic(above.lower()),
        "strict_signs_proved":True
    }

def audit(bits=256):
    if bits not in BITS:
        raise ValueError("precision must be 192/256/384")
    validate_table_s2_robust()
    source,sha=p79.load_source()
    result=[]
    with ctx.workprec(bits):
        cs=p79.constants()
        L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
        cs["c1"]=6/pi**2+p79.q("1.048")/L
        cs["c2"]=pi*m0*p79.q(".62")/L+6*m0/(pi*p79.T0)+pi*m0*p79.q("1.048")/(arb(p79.T0)*L)
        for i in range(10,17):
            original=TABLE[i-1]
            version=TABLE_S2_PUBLIC_ROBUST[i-1]
            public=F(source[i-1]["b4"])
            D,alpha,beta,e,components=factor_row(original,cs)
            old=root_bracket(D,alpha,beta,e,original.b1)
            robust=root_bracket(D,alpha,beta,e,version.b1)
            if not (F(old["root_upper_rational"])<public<
                    F(robust["root_lower_rational"])):
                raise ArithmeticError(f"b4_public not strictly between old and new roots at row {i}")
            k3=alpha*cap(public)+beta
            B=D*p79.pw(k3,e)
            if not (k3>1 and B.lower()>p79.q(original.b1) and B.upper()<p79.q(version.b1)):
                raise ArithmeticError(f"majorant at public b4 not between printed and robust: row {i}")
            derivative=D*e*alpha*p79.pw(k3,e-1)
            if not derivative>0:
                raise ArithmeticError("nonpositive certified B1 derivative")
            result.append({
                "row":i, "b4_public_exact":str(public),
                "B1_old_printed":original.b1,"B1_new_robust":version.b1,
                "e_exact":str(F(2*i,48-i)),
                "old_print_threshold":old,
                "new_robust_threshold":robust,
                "B1_public_lower_dyadic":p79.dyadic(B.lower()),
                "B1_public_upper_dyadic":p79.dyadic(B.upper()),
                "dB1_db4_lower_dyadic":p79.dyadic(derivative.lower()),
                "K3_public_lower_dyadic":p79.dyadic(k3.lower()),
                "K3_gt_1":True,
                "public_b4_between_strictly_certified_roots":True,
                "P79_factors":components
            })
    return {
        "schema":SCHEMA,
        "bits":bits,"rows":result,
        "P78_certificate_sha256":sha,
        "P79_method_sha256":source_sha(ROOT/"NUMERIA_P79_B1_PUBLIC"/"audit_public.py"),
        "P84_table_sha256":source_sha(ROOT/"NUMERIA_P76_COEFFICIENTS"/"table.py"),
        "root_meaning":"threshold of explicit P77/P79 B1 majorant, not root of RH",
        "original_S3_float_root_claims_replaced_by":"outward rational strict-sign Arb brackets",
        "proof_status":"REPO_INTERVAL_CERTIFIED_CONDITIONALLY",
        "Preissmann_weighted_Hilbert_theorem_reproved":False,
        "RZ_double_sum_reproved":False,
        "lower_X_3e10_3e12_covered":False,
        "global_RH_proved":False
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bits",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=audit(a.bits)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("P85 roots: bits",a.bits,"rows",len(result["rows"]),
          "all strict old-root < b4 public < new-root")

if __name__=="__main__":
    main()
