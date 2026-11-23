"""P92-N. Exact Fraction tests; analytic proof for all real parameters in README."""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
SCHEMA="NUMERIA.P92N.MACROSPAN.LOCAL_GAPS.1"
C=F(13,16); Q=F(5,7); T=F(18,5)
PI_L=16*(F(1,5)-F(1,375))-4*F(1,239)
PI_U=16*(F(1,5)-F(1,375)+F(1,15625))-4*(F(1,239)-F(1,3*239**3))

def rational_threshold():
    if not (PI_L>F(157,50) and PI_U<F(22,7)): raise ArithmeticError("P92N exact regression check 1 failed")
    if not (F(5477,5000)**2<F(6,5) and PI_L**2/F(3)*F(5477,5000)>T): raise ArithmeticError("P92N exact regression check 2 failed")
    if not (2*Q**2>1 and 125*16**4<294*13**4): raise ArithmeticError("P92N exact regression check 3 failed")
    b=Q/C**2*F(22,7)**2/F(3)
    if not (b==F(619520,173901) and b<T): raise ArithmeticError("P92N exact regression check 4 failed")
    return b

def build(roots,sizes,materialize=False):
    roots=list(map(F,roots)); sizes=list(sizes)
    if not (len(roots)==len(sizes)>0): raise ArithmeticError("P92N exact regression check 5 failed")
    if not (all(x>0 for x in roots) and all(m>=2 for m in sizes)): raise ArithmeticError("P92N exact regression check 6 failed")
    steps=[v*v for v in roots]
    if not (all((m-1)*h<=1 for m,h in zip(sizes,steps))): raise ArithmeticError("P92N exact regression check 7 failed")
    g=[]; x=[]; at=F(0)
    for i,(v,h,m) in enumerate(zip(roots,steps,sizes)):
        if i: at+=max(C,steps[i-1],h)
        start=at; end=start+(m-1)*h
        g.append((v,h,m,start,end))
        if materialize: x += [start+j*h for j in range(m)]
        at=end
    return g,x

def exact_ldl(a):
    M=[list(map(F,row)) for row in a];n=len(M);least=None
    if not (all(len(row)==n for row in M)): raise ArithmeticError("P92N exact regression check 8 failed")
    for k in range(n):
        p=M[k][k]
        if p<=0:raise ArithmeticError("nonpositive exact LDL pivot")
        least=p if least is None else min(least,p)
        for i in range(k+1,n):
            for j in range(i,n):
                M[j][i]-=M[i][k]*M[j][k]/p
                M[i][j]=M[j][i]
    return least

def full_delta(groups):
    x=[];roots=[]
    for v,h,m,start,end in groups:
        x.extend(start+j*h for j in range(m))
        roots.extend([v]*m)
    n=len(x)
    return [[F(0) if i==j else roots[i]*roots[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
             for j in range(n)] for i in range(n)]

def audit(name, roots, sizes,full=False):
    g,_=build(roots,sizes)
    K=len(g); A=[[F(0)]*K for _ in range(K)]
    pairs=0
    for r,s in combinations(range(K),2):
        vr,hr,mr,start_r,end_r=g[r]
        vs,hs,ms,start_s,end_s=g[s]
        D=start_s-end_r
        if not (D>=C*(s-r) and D>=max(hr,hs)): raise ArithmeticError("P92N exact regression check 9 failed")
        beta=vr*vs*(vr-vs)**2/F(2)
        HS2=sum((beta**2/(D+j*hr+k*hs)**4
                 for j in range(mr) for k in range(ms)),F(0))
        if not (2*D**4*HS2<=1): raise ArithmeticError("P92N exact regression check 10 failed")
        A[r][s]=A[s][r]=Q/D**2
        if not (HS2<=A[r][s]**2): raise ArithmeticError("P92N exact regression check 11 failed")
        pairs+=1
    rowmax=max(map(sum,A))
    if not (rowmax<rational_threshold()): raise ArithmeticError("P92N exact regression check 12 failed")
    block_pivot=exact_ldl([[T*(i==j)-A[i][j] for j in range(K)] for i in range(K)])
    full_pivot=None
    if full:
        Delta=full_delta(g);n=len(Delta)
        full_pivot=exact_ldl([[T*(i==j)-Delta[i][j] for j in range(n)] for i in range(n)])
    return dict(name=name,groups=K,points=sum(sizes),distinct_scales=len(set(v*v for v in map(F,roots))),
                pair_hs_checks=pairs,max_comparison_row_sum=str(rowmax),
                comparison_ldl_min_pivot=str(block_pivot),
                full_ldl_min_pivot=None if full_pivot is None else str(full_pivot))

def large_macro():
    roots=[F(1,k) for k in range(1,41)]
    sizes=[k*k+1 for k in range(1,41)]
    g,_=build(roots,sizes)
    D=g[-1][3]-g[0][4]
    if not (D==F(559,8) and D**4<39**5 and sum(sizes)==22180): raise ArithmeticError("P92N exact regression check 13 failed")
    A=[[F(0)]*40 for _ in range(40)]
    for i,j in combinations(range(40),2):
        d=g[j][3]-g[i][4]
        A[i][j]=A[j][i]=Q/d**2
    rowmax=max(map(sum,A))
    if not (rowmax<rational_threshold()): raise ArithmeticError("P92N exact regression check 14 failed")
    return {"groups":40,"points":22180,"distinct_scales":40,
            "max_group_size":1601,"end_to_end_gap":str(D),
            "P92M_separation_violated":True,"full_matrix_constructed":False,
            "macro_comparison_row_sum_exact":str(rowmax)}

def certify():
    configs=[
        ("two_close",[F(1,2),F(1,3)],[5,10],True),
        ("three_irregular",[F(1,2),F(1,4),F(1,3)],[3,7,5],True),
        ("four_contrasts",[F(1),F(1,3),F(1,5),F(1,2)],[2,4,8,3],True),
        ("six_variable",[F(1,2),F(1,3),F(1,4),F(1,5),F(1,6),F(1,7)],[5,10,17,26,37,50],False),
        ("large_density",[F(1),F(1,7),F(1,30)],[2,50,300],False),
        ("eight_mixed",[F(1,2),F(1,8),F(1,4),F(1,6),F(1,10),F(1,3),F(1,9),F(1,7)],[4,5,3,6,7,3,5,6],False)
    ]
    cases=[audit(*config) for config in configs]
    if not (len(cases)==6 and sum(c["pair_hs_checks"] for c in cases)==56): raise ArithmeticError("P92N exact regression check 15 failed")
    if not (sum(c["full_ldl_min_pivot"] is not None for c in cases)==3): raise ArithmeticError("P92N exact regression check 16 failed")
    return dict(schema=SCHEMA,
                analytic_theorem_for_all_real_macroscale_configurations=True,
                unrestricted_point_configurations_proved=False,
                full_Preissmann_Lemma6_proved=False,
                Riemann_Hypothesis_proved=False,
                universal_rational_upper=str(rational_threshold()),
                rational_regressions=len(cases),exact_pair_hs_checks=56,
                full_exact_ldl_count=3,cases=cases,macro40=large_macro())

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    d=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("P92N PASS: 6 cases, 56 exact block HS, 3 full LDL, 40-group macro certificate")
    print("SHA256",hashlib.sha256(args.out.read_bytes()).hexdigest())
    print("Universal Preissmann lemma and RH remain OPEN")
if __name__=="__main__":main()
