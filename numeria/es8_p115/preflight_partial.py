#!/usr/bin/env python3
"""P115: fail-closed preflight of available P111 row shards, not a matrix proof.

Run repeatedly on a fixed P111 Actions run. Accept a partial set only for
progress reporting. Never claim 169x169 positivity until all 24 row shards
and the independent FLINT sign certificate are complete.
"""
import argparse,glob,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
SCALE=1<<80
NUM_SHARDS=24
N=169
def dyadic(a,e):return F(int(a))*F(2)**int(e)
def interval(v):
    if not isinstance(v,list) or len(v)!=4:raise ValueError("invalid Arb encoding")
    m,e,r,re=v
    mid=dyadic(m,e);rad=dyadic(r,re)
    if rad<0:raise ValueError("negative Arb radius")
    return mid-rad,mid+rad
def error_at_scale(v):
    lo,hi=interval(v)
    mid=(lo+hi)/2;rad=(hi-lo)/2
    vi=round(mid*SCALE)
    extra=abs(mid-F(vi,SCALE))+rad
    return (extra.numerator*SCALE+extra.denominator-1)//extra.denominator
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--shards",required=True,type=Path)
    ap.add_argument("--old-matrix",required=True,type=Path)
    ap.add_argument("--github-output",type=Path,default=None)
    ap.add_argument("--out",type=Path,default=Path("P115_PREFLIGHT.json"))
    a=ap.parse_args()
    old=json.loads(a.old_matrix.read_text(encoding="utf-8"))
    centers=[1]+[i for i in range(2,1001) if all(i%d for d in range(2,math.isqrt(i)+1))]
    if old.get("centers")!=centers or len(old["entries"])!=14365:
        raise ValueError("old comparison matrix does not match prime centres")
    files=sorted(a.shards.glob("P111_1000_ROWS_*.json"))
    if not files:raise RuntimeError("no P111 Arb data available; no certification")
    found_rows=set();allkeys=set();digest_set=set()
    max_radius=F(0);ss=0;proofs=0
    for path in files:
        d=json.loads(path.read_text(encoding="utf-8"))
        if (d.get("status")!="P109_ROW_SHARD_INTERVALS_CONDITIONAL_ON_GL_CAUCHY_ARITHMETIC"
            or d.get("bits")!=256 or d.get("epsilon")!="1/10000000000000000"
            or d.get("matrix_dimension")!=169 or d.get("max_prime")!=1000
            or d.get("primepower_upper_integer")!=2500 or d.get("shards")!=24
            or d.get("pilot_pair") is not None or d.get("centers")!=centers):
            raise ValueError("incompatible shard metadata: "+path.name)
        idx=d["shard_index"];i=idx*N//NUM_SHARDS;j=(idx+1)*N//NUM_SHARDS
        if not(0<=idx<NUM_SHARDS) or d["row_start"]!=i or d["row_end"]!=j:
            raise ValueError("row shard index/range mismatch")
        if idx in found_rows:raise ValueError("duplicate shard index")
        found_rows.add(idx); digest_set.add(d["source_sha256"])
        expected={f"{centers[k]},{centers[m]}" for k in range(i,j) for m in range(k,N)}
        if set(d["entries"])!=expected:raise ValueError("missing or extra entries")
        if d.get("diagonal_full_proof_record",{}).get("entry")!=d["diagonal_reference"]:
            raise ValueError("diagonal proof not retained")
        for key,v in d["entries"].items():
            if key in allkeys:raise ValueError("duplicate interval")
            allkeys.add(key)
            if d["proof_records"].get(key,{}).get("entry")!=v:
                raise ValueError("interval/proof discrepancy")
            lo,hi=interval(v)
            oldlo,oldhi=interval(old["entries"][key])
            if max(lo,oldlo)>min(hi,oldhi):raise ValueError("new/old interval conflict "+key)
            r=(hi-lo)/2
            max_radius=max(r,max_radius)
            err=error_at_scale(v)
            u,w=map(int,key.split(","))
            ss+=err*err*(1 if u==w else 2)
            proofs+=1
        print("P115_PREFLIGHT_SHARD",idx,"pairs",len(expected),flush=True)
    if len(digest_set)!=1:raise ValueError("different P111 source binaries")
    radius_int=math.isqrt(ss)
    radius_int+=radius_int*radius_int<ss
    complete=(found_rows==set(range(NUM_SHARDS)) and len(allkeys)==N*(N+1)//2)
    result={
       "status":"ALL_24_SHARDS_VALID__READY_FOR_EXACT_169_CERTIFICATE" if complete
                else "PARTIAL_INTERVAL_EVIDENCE__NOT_CERTIFIED",
       "complete":complete,"shard_count":len(found_rows),
       "missing_shard_indices":sorted(set(range(NUM_SHARDS))-found_rows),
       "distinct_entries_verified":len(allkeys),"expected_entries":14365,
       "exact_old_new_interval_intersections":proofs,
       "same_sha256_source_across_shards":True,
       "maximum_new_radius_rational":str(max_radius),
       "maximum_new_radius_display":float(max_radius),
       "partial_Frobenius_upper_rational":str(F(radius_int,SCALE)),
       "partial_Frobenius_upper_display":float(F(radius_int,SCALE)),
       "radius_scope":"FULL_MATRIX" if complete else "PARTIAL_POSITIONS_ONLY",
       "important":"Positivity is NEVER inferred from successful preflight; exact charpoly and shift>global rho required."
    }
    a.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    if a.github_output is not None:
        with a.github_output.open("a",encoding="utf-8") as f:
            f.write(f"ready={'true' if complete else 'false'}\n")
    print("P115_PREFLIGHT_RESULT",json.dumps({k:result[k] for k in
        ("status","shard_count","distinct_entries_verified","partial_Frobenius_upper_display")}),flush=True)
if __name__=="__main__":main()
