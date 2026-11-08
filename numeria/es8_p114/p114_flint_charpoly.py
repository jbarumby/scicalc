#!/usr/bin/env python3
"""NUMERIA P114: fast EXACT integer positivity certificate using FLINT charpoly.

Given an interval matrix W_ij = midpoint_ij +/- radius_ij, construct integer
M and E at a power-of-two scale S. Then ||W-M/S||_2 <= ||E||_F / S.
Using the exact fmpz_mat characteristic polynomial of A=M-shift*I,
check positivity of every coefficient in det(x*I+A). Because A is symmetric
all roots are real; positive coefficients imply no nonnegative root, hence
A is positive definite. This is a rigorous INTEGER test, not a numeric
eigensolver. A lower bound for the INTERVAL MATRIX requires shift>rho.
All bounds are conditional on Arb/Cauchy GL inclusions in input.
"""
import argparse,hashlib,json,math,sys,time,zipfile
from fractions import Fraction as Q
from pathlib import Path
from flint import fmpz_mat

SCALE=1<<80

def exact_dyadic(a,e):
    return Q(int(a))*Q(2)**int(e)

def round_enclosure(record):
    if len(record)!=4:raise ValueError("bad Arb dyadic record")
    a,e,r,re=record
    mid=exact_dyadic(a,e);rad=exact_dyadic(r,re)
    if rad<0:raise ValueError("negative radius")
    v=round(mid*SCALE)
    err=abs(mid-Q(v,SCALE))+rad
    i=(err.numerator*SCALE+err.denominator-1)//err.denominator
    return int(v),int(i)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--bound",required=True,type=int,choices=(500,1000))
    p.add_argument("--input-zip",required=True,type=Path)
    p.add_argument("--factor",default=None)
    p.add_argument("--out",default=None,type=Path)
    args=p.parse_args()
    if args.factor is None:
        args.factor="0.99999" if args.bound==500 else "0.90"
    with zipfile.ZipFile(args.input_zip) as z:
        matrix=json.loads(z.read(f"P109_{args.bound}_COMPACT_MATRIX.json"))
        audit=json.loads(z.read(f"P109_{args.bound}_AUDIT.json"))
    n=len(matrix["centers"])
    if n!={500:96,1000:169}[args.bound]:
        raise ValueError("incorrect dimension")
    if (matrix["max_prime"]!=args.bound or matrix["precision_bits"]!=256
        or matrix["epsilon"]!="1/100000000000"):
        raise ValueError("incorrect experiment metadata")
    if len(matrix["entries"])!=n*(n+1)//2:
        raise ValueError("incomplete interval matrix")
    M=[[0]*n for _ in range(n)]
    E=[[0]*n for _ in range(n)]
    for i,a in enumerate(matrix["centers"]):
        for j in range(i,n):
            b=matrix["centers"][j]
            m,rad=round_enclosure(matrix["entries"][f"{a},{b}"])
            M[i][j]=M[j][i]=m
            E[i][j]=E[j][i]=rad
    ss=sum(t*t for row in E for t in row)
    rho=math.isqrt(ss)
    rho+=rho*rho<ss
    guess=Q(repr(audit["full_approx_lambda_min_not_certified"]))
    factor=Q(args.factor)
    if not 0<factor<1:raise ValueError("factor must lie strictly in (0,1)")
    numerator=guess*factor*SCALE
    shift=numerator.numerator//numerator.denominator
    if shift<=0:raise ValueError("nonpositive proposed shift")
    # No floating-point operations enter the matrix or sign decisions.
    t0=time.monotonic()
    A=fmpz_mat(n,n,[int(M[i][j]-(shift if i==j else 0))
                   for i in range(n) for j in range(n)])
    polynomial=A.charpoly()
    coeff=[int(c) for c in polynomial.coeffs()]
    if len(coeff)!=n+1 or coeff[-1]!=1:
        raise ArithmeticError("wrong degree/non-monic characteristic polynomial")
    # det(xI+A) = (-1)^n det(-xI-A) = (-1)^n chi_A(-x).
    transformed=[(-1)**(n+k)*c for k,c in enumerate(coeff)]
    positive_coefficients=sum(c>0 for c in transformed)
    midpoint_shift_pd=positive_coefficients==n+1
    certified=bool(midpoint_shift_pd and shift>rho)
    result={
       "phase":"NUMERIA_ES8_P114_FLINT_CHARPOLY_EXACT",
       "source_zip":args.input_zip.name,
       "source_sha256":hashlib.sha256(args.input_zip.read_bytes()).hexdigest(),
       "dimension":n,"max_prime":args.bound,
       "integer_scale_denominator":SCALE,
       "shift_exact_numerator":shift,
       "shift_factor_exact":str(factor),
       "frob_radius_upper_exact":str(Q(rho,SCALE)),
       "shifted_midpoint_strictly_pd":midpoint_shift_pd,
       "positive_charpoly_coefficients":positive_coefficients,
       "charpoly_degree":n,
       "charpoly_sha256":hashlib.sha256(("\n".join(map(str,coeff))+"\n").encode()).hexdigest(),
       "elapsed_seconds":time.monotonic()-t0,
       "interval_matrix_positive_conditional":certified,
       "lower_lambda_rational":str(Q(shift-rho,SCALE)) if certified else None,
       "lower_lambda_display":float(Q(shift-rho,SCALE)) if certified else None,
       "status":("POSITIVE_169_OR_96_CONDITIONAL" if certified
          else ("MIDPOINT_ONLY_PD__INTERVAL_INCONCLUSIVE" if midpoint_shift_pd
                else "INCONCLUSIVE_MIDPOINT_SHIFT")),
       "warning":"Strictly positive charpoly coefficients prove shifted MIDPOINT PD; the entire interval matrix requires shift>rho. Finite-dimension conditional result, no RH."
    }
    out=args.out or Path(f"P114_{args.bound}_FLINT_EXACT.json")
    out.write_text(json.dumps(result,indent=2)+"\n")
    Path(f"P114_{args.bound}_CHARPOLY_SIGNED_COEFFS.txt").write_text(
         "\n".join(map(str,transformed))+"\n")
    print("EXACT_FLINT_CERT",json.dumps({
        "dimension":n,"status":result["status"],
        "positive_coefficients":positive_coefficients,
        "required_coefficients":n+1,"frob_radius":float(Q(rho,SCALE)),
        "shift":float(Q(shift,SCALE)),"lambda_lower":result["lower_lambda_display"],
        "elapsed_seconds":result["elapsed_seconds"]
    }),flush=True)
    if not midpoint_shift_pd:
        raise RuntimeError("FLINT exact characteristic polynomial failed midpoint positivity")
    if args.bound==500 and not certified:
        raise RuntimeError("500 x 500 certificate required from original P109 data")
if __name__=="__main__":
    sys.exit(main())
