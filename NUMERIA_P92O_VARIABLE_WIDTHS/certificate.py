"""NUMERIA P92-O. Exact Fraction audit of the variable-width cluster theorem.

The proof for ALL real admissible parameters is in README.md.  Computer tests
check specific rational geometries, not general mathematical universality.
"""
from fractions import Fraction as F
from math import isqrt
import argparse, json, hashlib

SCHEMA="NUMERIA.P92O.DELTA.VARIABLE_WIDTHS.1"
T=F(18,5)

def root(q):
    q=F(q)
    if q<0: raise ValueError("negative root")
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    if a*a!=q.numerator or b*b!=q.denominator:
        raise ValueError("not a square Fraction")
    return F(a,b)

def clusters(specs,gaps):
    if not specs or len(gaps)!=len(specs)-1: raise ValueError("bad groups")
    out=[]; x=F(0)
    for k,(m,h) in enumerate(specs):
        m,h=int(m),F(h)
        if m<2 or h<=0: raise ValueError("nonpositive step/cardinality")
        w=(m-1)*h
        if k:
            G=F(gaps[k-1])
            if G<max(h,out[-1]["h"]):
                raise ValueError("inter-group gap below neighboring steps")
            x=out[-1]["end"]+G
        out.append({"m":m,"h":h,"w":w,"start":x,"end":x+w})
    return out

def dist(gs,i,j):
    return gs[j]["start"]-gs[i]["end"]

def comparison(gs):
    k=len(gs); a=[[F(0)]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1,k):
            x=gs[i]["w"]+gs[j]["w"]
            a[i][j]=a[j][i]=x*x/(2*dist(gs,i,j)**2)
    return a

def pairwise_gate(gs):
    for i in range(len(gs)):
        for j in range(i+1,len(gs)):
            D=dist(gs,i,j)
            w=gs[i]["w"]+gs[j]["w"]
            if D**4<(j-i)**3*w**4:
                return False
    return True

def ldlt_positive(m):
    a=[[F(z) for z in row] for row in m]
    n=len(a); least=None
    for i in range(n):
        if len(a[i])!=n: raise ValueError("not square")
        for j in range(i):
            if a[i][j]!=a[j][i]: raise ValueError("not symmetric")
    for k in range(n):
        d=a[k][k]
        if d<=0: raise AssertionError("nonpositive rational LDL pivot")
        least=d if least is None else min(d,least)
        for i in range(k+1,n):
            v=a[i][k]
            for j in range(i,n):
                a[j][i]-=v*a[j][k]/d
                a[i][j]=a[j][i]
    return least

def hs_square(gs,i,j):
    a,b=gs[i],gs[j]
    ra,rb=root(a["h"]),root(b["h"])
    beta=ra*rb*(ra-rb)**2/2
    D=dist(gs,i,j)
    return sum((beta**2/(D+u*a["h"]+v*b["h"])**4
                for u in range(a["m"]) for v in range(b["m"])),F(0))

def rational_pivot_lower(p):
    # Small, exact decimal rational <= the exact positive LDL pivot.
    scale=10**9
    value=p*scale
    return str(F(value.numerator//value.denominator,scale))

def examples():
    e={}
    def wide(specs,factor):
        widths=[(m-1)*F(h) for m,h in specs]
        return specs,[F(factor)*(widths[i]+widths[i+1])
                      for i in range(len(specs)-1)]
    e["six_irregular_not_N"]=wide([
        (2,F(100)),(2,F(1)),(2,F(1,100)),(5,F(4)),
        (10,F(1,16)),(3,F(9))],7)
    s=[(2,F(36)),(2,F(1)),(5,F(1)),(2,F(9)),(2,F(16))]
    ws=[(m-1)*F(h) for m,h in s]
    e["short_spectral_only"]=(s,[max(F(3,4)*(ws[i]+ws[i+1]),F(s[i][1]),F(s[i+1][1]))
        for i in range(len(s)-1)])
    e["anchor_plus_12_microgroups"] = (
        [(2,F(100))]+[(2,F(1,1024))]*12,
        [F(100)]+[F(1)]*11
    )
    for k in (3,4,8,14,25):
        s=[(2+(j*7)%23,F((j%4+1)**2,(j%5+2)**2))
            for j in range(1,k+1)]
        e["variable_K"+str(k)]=wide(s,k+1)
    s=[(2000+(j*137)%2701,F((j%7+1)**2,(j%13+11)**2))
       for j in range(1,31)]
    e["thirty_giant_clusters"]=wide(s,32)
    return e

def full_delta(gs,cap=35):
    n=sum(g["m"] for g in gs)
    if n>cap: raise ValueError("full matrix is intentionally capped")
    x=[]; rt=[]
    for g in gs:
        v=root(g["h"])
        for j in range(g["m"]):
            x.append(g["start"]+j*g["h"])
            rt.append(v)
    mat=[[F(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            z=rt[i]*rt[j]*(rt[i]-rt[j])**2/(2*(x[i]-x[j])**2)
            mat[i][j]=mat[j][i]=z
    return mat

def certify():
    rows=[]; n_hs=0
    for name,(specs,gaps) in examples().items():
        gs=clusters(specs,gaps); k=len(gs); A=comparison(gs)
        good=pairwise_gate(gs)
        rowsums=[sum(row,F(0)) for row in A]
        if good and max(rowsums)>=3:
            raise AssertionError("O.1 exact analytic bound failed")
        p=ldlt_positive([[T*(i==j)-A[i][j] for j in range(k)] for i in range(k)])
        if p<=0: raise AssertionError("LDL")
        checks=0
        if name!="thirty_giant_clusters":
            for i in range(k):
                for j in range(i+1,k):
                    if gs[i]["m"]*gs[j]["m"]<=5000:
                        if hs_square(gs,i,j)>A[i][j]**2:
                            raise AssertionError("HS pair bound failed")
                        checks+=1
        n_hs+=checks
        full=None
        if name=="short_spectral_only":
            D=full_delta(gs); n=len(D)
            pfull=ldlt_positive([[T*(i==j)-D[i][j] for j in range(n)]
                                  for i in range(n)])
            full=rational_pivot_lower(pfull)
        wmax=max(g["w"] for g in gs)
        old=all(F(gaps[i])>=max(F(13,16)*wmax,gs[i]["h"],gs[i+1]["h"])
                for i in range(k-1))
        rows.append({"name":name,"K":k,
                    "points":sum(g["m"] for g in gs),
                    "pairwise_3over4":good,"old_P92N":old,
                    "max_row":str(max(rowsums)),
                    "ldlt_pivot_lower":rational_pivot_lower(p),
                    "pair_HS_checks":checks,"full_delta_ldlt_lower":full})
    byname={r["name"]:r for r in rows}
    star=byname["anchor_plus_12_microgroups"]
    assert F(star["max_row"])>T and not star["pairwise_3over4"]
    assert not byname["six_irregular_not_N"]["old_P92N"]
    assert byname["six_irregular_not_N"]["pairwise_3over4"]
    assert byname["thirty_giant_clusters"]["points"]>90000
    assert F(5477,5000)**2<F(6,5)
    assert F(157,50)**2/F(3)*F(5477,5000)>T
    return {"schema":SCHEMA,
        "theorem":"ALL real widths under pairwise D_rs >= |r-s|^(3/4)*(width_r+width_s): ||Delta||_2<3",
        "all_configurations":False,"Preissmann_general":False,
        "Riemann_hypothesis":False,
        "cases":rows,"case_count":len(rows),
        "HS_checks":n_hs,"exact_comparison_LDL":len(rows),
        "full_delta_LDL":sum(r["full_delta_ldlt_lower"] is not None for r in rows)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    c=certify()
    with open(args.out,"w",encoding="utf8") as f:
        json.dump(c,f,indent=2,ensure_ascii=False)
        f.write("\n")
    print("P92O exact rational:",c["case_count"],"cases,",
          c["HS_checks"],"HS,",
          c["exact_comparison_LDL"],"comparison LDL")
    print("SHA256:",hashlib.sha256(open(args.out,"rb").read()).hexdigest())

if __name__=="__main__":
    main()
