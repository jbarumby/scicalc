#!/usr/bin/env python3
"""P122: independently reconstruct 169x169 from original P111 Arb balls;
recalculate exact Frobenius error, then prove all leading minors positive
with C-backed FLINT fraction-free Bareiss with checked exact divisions.

A finite-dimensional CONDITIONAL proof; not Riemann Hypothesis.
"""
import argparse,hashlib,json,math,time,zipfile
from fractions import Fraction as Q
from pathlib import Path

S=1<<88
MATRIX_SHA="747f6a045ff6279ac3677a6ca77c40fa3816537d9f98d3f7af3ebe4c2355e51b"
AUDIT_SHA="f2fc5f17a70832286b0b8969d9647aaa3ef366df62c075f2d180735cb54b2ae9"

def ball(v):
    if not isinstance(v,list) or len(v)!=4:raise ValueError("Invalid Arb ball")
    mid=Q(int(v[0]))*Q(2)**int(v[1])
    rad=Q(int(v[2]))*Q(2)**int(v[3])
    if rad<0:raise ValueError("Negative radius")
    im=round(mid*S)
    er=rad+abs(mid-Q(im,S))
    re=(er.numerator*S+er.denominator-1)//er.denominator
    return int(im),int(re)

def inputs(source):
    with zipfile.ZipFile(source) as z:
        if z.testzip():raise ValueError("Damaged ZIP")
        mat_bytes=z.read("P111_1000_COMPACT_MATRIX.json")
        audit_bytes=z.read("P111_1000_AUDIT.json")
    if hashlib.sha256(mat_bytes).hexdigest()!=MATRIX_SHA:
        raise ValueError("Original 169x169 matrix SHA differs")
    if hashlib.sha256(audit_bytes).hexdigest()!=AUDIT_SHA:
        raise ValueError("Original audit SHA differs")
    d=json.loads(mat_bytes);a=json.loads(audit_bytes)
    centers=d["centers"];n=len(centers)
    if (n!=169 or centers[0]!=1 or centers[-1]!=997 or
        len(d["entries"])!=14365 or d["epsilon"]!="1/10000000000000000"
        or d["precision_bits"]!=256 or a["status"]!="POSITIVE_CONDITIONAL"):
        raise ValueError("Mismatched experiment")
    M=[[0]*n for _ in range(n)]
    E=[[0]*n for _ in range(n)]
    for i,c in enumerate(centers):
        for j in range(i,n):
            v,r=ball(d["entries"][f"{c},{centers[j]}"])
            M[i][j]=M[j][i]=v
            E[i][j]=E[j][i]=r
    ss=sum(e*e for row in E for e in row)
    rho=math.isqrt(ss)
    rho+=int(rho*rho<ss)
    oldrho=Q(a["full_error_radius_upper_rational"])*S
    if oldrho.denominator!=1 or rho!=oldrho.numerator:
        raise ArithmeticError("Independently recomputed Frobenius rho mismatched")
    lb=Q(a["full_lambda_lower_rational"])*S
    if lb.denominator!=1 or lb<=0:raise ValueError("Invalid lower bound")
    shift=lb.numerator+rho
    if not shift>rho:raise ValueError("No positive margin")
    return M,rho,shift

def bareiss(A):
    from flint import fmpz,fmpz_mat
    n=len(A)
    if any(A[i][j]!=A[j][i] for i in range(n) for j in range(i+1,n)):
        raise ValueError("Not symmetric")
    T=fmpz_mat(n,n,[int(x) for row in A for x in row])
    prev=fmpz(1)
    deltas=[]
    divisions=0
    for k in range(n):
        p=T[k,k]
        if p<=0:raise ArithmeticError(f"Nonpositive leading minor at {k}: {p}")
        deltas.append(str(int(p)))
        for i in range(k+1,n):
            mik=T[i,k]
            for j in range(k+1,n):
                num=T[i,j]*p-mik*T[k,j]
                q,r=divmod(num,prev)
                if r!=0:raise ArithmeticError(f"Nonzero Bareiss division remainder at {k},{i},{j}")
                T[i,j]=q
                divisions+=1
        prev=p
        if k%30==0 or k==n-1:print("P122_MINOR",k+1,"of",n,"positive",flush=True)
    return deltas,divisions

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--prefix",type=int,default=169)
    p=ap.parse_args()
    M,rho,shift=inputs(p.input)
    n=p.prefix
    if not 1<=n<=169:raise ValueError("Bad prefix")
    A=[[M[i][j]-(shift if i==j else 0) for j in range(n)] for i in range(n)]
    start=time.monotonic()
    deltas,divisions=bareiss(A)
    elapsed=time.monotonic()-start
    lower=Q(shift-rho,S)
    result={"schema":"NUMERIA.P122.BAREISS.FMPZ.1",
       "n":n,"full_matrix_size":169,"full_input_entries":14365,
       "p111_matrix_json_sha256":MATRIX_SHA,
       "p111_audit_json_sha256":AUDIT_SHA,
       "arbitrary_precision_engine":"FLINT fmpz_mat exact integer Bareiss",
       "integer_scale":S,"shift_integer":shift,"rho_integer":rho,
       "rho_full_rational":str(Q(rho,S)),
       "lambda_lower_full_rational":str(lower),
       "positive_leading_minors":len(deltas),
       "checked_bareiss_exact_divisions":divisions,
       "leading_minors_integerized":[str(v) for v in deltas],
       "computed_seconds":elapsed,
       "status":"169X169_POSITIVE_CONDITIONAL" if n==169 else "PREFIX_PIVOT_TEST_ONLY",
       "note":"Verifies original 169-centre finite matrix conditional on Arb/GL/Cauchy inclusions; no RH proof."
    }
    p.output.write_text(json.dumps(result,indent=2)+"\n")
    print("P122_POSITIVE",json.dumps({
       "n":n,"rho":float(Q(rho,S)),"lambda_lower":float(lower),
       "positive_minors":len(deltas),"checked_divisions":divisions,
       "seconds":round(elapsed,3)}),flush=True)
if __name__=="__main__":main()
