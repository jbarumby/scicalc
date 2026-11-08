#!/usr/bin/env python3
"""ES8-P108 — fail-closed 47x47 matrix audit using exact integer Bareiss LDL.

Numerical eigvalsh values are ONLY guesses for a rational shift. Certification
uses integer leading-principal-minor signs and exact dyadic error bounds.
A green workflow does not imply positive definiteness; inspect status flags.
"""
import csv, glob, json, math, sys
from fractions import Fraction as Q
from pathlib import Path

SCALE = 1 << 80
CENTERS = tuple([1] + [p for p in range(2, 200)
                        if all(p % d for d in range(2, math.isqrt(p)+1))])

def dyadic(m, e):
    return Q(int(m)) * Q(2)**int(e)

def ceil_q(x):
    return -(-x.numerator // x.denominator)

def exact_mid_err(record):
    m, e, r, re = record
    mid=dyadic(m,e)
    rad=dyadic(r,re)
    if rad < 0:
        raise ValueError("negative radius")
    scaled_mid=round(mid*SCALE)
    distance=abs(mid-Q(scaled_mid,SCALE))
    err_int=ceil_q((rad+distance)*SCALE)
    return scaled_mid,err_int

def bareiss_strict_positive(A):
    """Fraction-free Sylvester test: every pivot is a positive leading minor."""
    n=len(A)
    a=[row[:] for row in A]
    prev=1
    for k in range(n):
        p=a[k][k]
        if p<=0: return False
        for i in range(k+1,n):
            aik=a[i][k]
            for j in range(i,n):
                num=p*a[i][j]-aik*a[k][j]
                q,r=divmod(num,prev)
                if r:
                    raise ArithmeticError("nonexact Bareiss division")
                a[i][j]=a[j][i]=q
        for i in range(k+1,n):
            a[i][k]=a[k][i]=0
        prev=p
    return True

def certified_sqrt_upper_sum(err):
    total=sum(x*x for row in err for x in row)
    r=math.isqrt(total)
    return r if r*r==total else r+1

def collect(paths):
    if not paths:
        raise ValueError("missing P108 shard outputs")
    entries={}
    row_shards=[]
    sha_set=set()
    intervals=[]
    for file in paths:
        obj=json.loads(Path(file).read_text(encoding="utf-8"))
        if obj.get("status")!="P108_ARBITRARY_MATRIX_ENTRY_INTERVALS_ANALYTICALLY_CONDITIONAL":
            raise ValueError(f"incomplete shard {file}")
        if tuple(obj.get("centers",[]))!=CENTERS or obj.get("matrix_dimension")!=47:
            raise ValueError(f"incompatible centers or matrix dimension {file}")
        if obj.get("epsilon")!="1/100000000000" or obj.get("bits")!=256:
            raise ValueError(f"incompatible eps or precision {file}")
        if obj.get("primepower_upper_integer")!=500:
            raise ValueError(f"incompatible prime-power truncation {file}")
        if obj.get("pilot_pair") is not None:
            raise ValueError("pilot result cannot be used for matrix certificate")
        first,last=obj["row_start"],obj["row_end"]
        if not (0<=first<last<=47):raise ValueError("shard outside range")
        row_shards.extend(range(first,last))
        sha_set.add(obj["source_sha256"])
        expected={f"{CENTERS[i]},{CENTERS[j]}"
                  for i in range(first,last) for j in range(i,47)}
        if set(obj["entries"])!=expected:
            raise ValueError("missing or extra entry within shard")
        d=obj["diagonal_reference"]
        intervals.append(exact_mid_err(d))
        for key,v in obj["entries"].items():
            if key in entries:raise ValueError(f"overlapping shard: {key}")
            entries[key]=v
            proof=obj["proof_records"].get(key)
            if proof is None or proof.get("entry")!=v:
                raise ValueError("entry/proof mismatch: "+key)
    if len(sha_set)!=1 or sorted(row_shards)!=list(range(47)):
        raise ValueError("missing/overlapping shard or source hash inconsistency")
    if len(entries)!=47*48//2:
        raise ValueError("missing entries in 47x47 matrix")
    # The repeated W(1,1) computations must have overlapping intervals.
    # They need not have identical endpoints, but must cover same true value.
    first_mid,first_err=intervals[0]
    low=max(Q(m-e,SCALE) for m,e in intervals)
    high=min(Q(m+e,SCALE) for m,e in intervals)
    if low>high:
        raise ValueError("inconsistent repeated diagonal references")
    M=[[0]*47 for _ in range(47)]
    E=[[0]*47 for _ in range(47)]
    for i,x in enumerate(CENTERS):
        for j in range(i,47):
            y=CENTERS[j]
            mid,err=exact_mid_err(entries[f"{x},{y}"])
            M[i][j]=M[j][i]=mid
            E[i][j]=E[j][i]=err
    return M,E,sha_set.pop(),len(paths),len(entries)

def audit(paths):
    import numpy as np
    M,E,sha,nfiles,nentries=collect(paths)
    rows=[]
    for k in range(1,48):
        m=[row[:k] for row in M[:k]]
        e=[row[:k] for row in E[:k]]
        rho=certified_sqrt_upper_sum(e)
        heuristic=float(np.linalg.eigvalsh(
            np.asarray(m,dtype="float64")/SCALE
        )[0])
        # Candidate shift chosen with floats is NOT trusted: the integer
        # Bareiss strict-PD test below confirms or rejects it exactly.
        sigma=math.floor(heuristic*0.94*SCALE)
        checks=0
        proved=False
        while sigma>rho and checks<24:
            checks+=1
            mat=[[m[i][j]-(sigma if i==j else 0)
                  for j in range(k)] for i in range(k)]
            if bareiss_strict_positive(mat):
                proved=True
                break
            sigma//=2
        lower_int=sigma-rho if proved else None
        rows.append({
            "n":k,"largest_center":CENTERS[k-1],
            "numerical_lambda_heuristic":heuristic,
            "Frobenius_error_upper_int":rho,
            "Frobenius_error_upper":float(Q(rho,SCALE)),
            "candidate_shift_int":sigma,
            "lambda_lower_dyadic_numerator":lower_int,
            "lambda_lower_dyadic_denominator":SCALE,
            "lambda_lower_for_plot":float(Q(lower_int,SCALE)) if proved else None,
            "strictly_positive_certified_conditional":proved and lower_int>0,
            "bareiss_tests":checks
        })
        print("PREFIX",k,"THROUGH",CENTERS[k-1],"LOWER",
              rows[-1]["lambda_lower_for_plot"],"CERT",
              rows[-1]["strictly_positive_certified_conditional"],flush=True)
    return {
        "status":("ALL_47_PREFIXES_POSITIVE_CONDITIONAL" if all(
             r["strictly_positive_certified_conditional"] for r in rows)
             else "SOME_PREFIXES_INCONCLUSIVE"),
        "centers":CENTERS,"epsilon":"1/100000000000",
        "precision_bits":256,
        "distinct_entries":nentries,"shard_count":nfiles,
        "source_sha256":sha,"arithmetic":"integer Bareiss LDL plus dyadic outward error radius",
        "numeric_eigenvalues_only_heuristics":True,
        "lambda_min_interlacing_note":"Nested principal submatrices have nonincreasing exact smallest eigenvalue",
        "rows":rows,
        "warning":"Finite Weil 47x47 only. All interval conclusions depend on analytic Cauchy GL bounds and FLINT/Arb correctness. No RH proof."
    }

def main():
    paths=sorted(glob.glob("P108_ROWS_*.json"))
    result=audit(paths)
    Path("P108_47X47_EXACT_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
    with open("P108_LAMBDA_TREND.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["n","largest_center","lambda_lower_for_plot",
              "numerical_lambda_heuristic","Frobenius_error_upper",
              "strictly_positive_certified_conditional"])
        w.writeheader()
        for r in result["rows"]:
            w.writerow({key:r[key] for key in w.fieldnames})
    certified=sum(bool(r["strictly_positive_certified_conditional"]) for r in result["rows"])
    last=result["rows"][-1]
    report=(
       "# NUMERIA ES8-P108 — interval matrix through prime centre 199\n\n"
       f"Centres: {len(CENTERS)}; interval entries: {result['distinct_entries']}; "
       f"epsilon=1e-11; 256-bit Arb; prime powers <=500.\n\n"
       f"Exact integer Bareiss positive prefixes: {certified}/47.\n"
       f"Full 47x47 lambda lower: {last['lambda_lower_for_plot']} "
       f"({last['strictly_positive_certified_conditional']}).\n\n"
       "Numerical eigenvalues were used only to propose candidates; "
       "positive signs are established by exact integer Sylvester tests "
       "after a rigorous Frobenius enclosure.\n\n"
       "**Limits:** finite 47x47 subspace; conditional on analytic Cauchy GL "
       "and Arb enclosures; no RH proof.\n"
    )
    Path("P108_REPORT.md").write_text(report)
    print("RESULT",result["status"],"CERTIFIED_PREFIXES",certified,"/47")
    return 0

if __name__=="__main__":
    sys.exit(main())
