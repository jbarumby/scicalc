"""P92-L exact rational checks. The analytic theorem is in README.md."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json

SCHEMA='NUMERIA.P92L.DELTA.LINEAR_SUBLINEAR_ADAPTIVE_BLOCKS.1'
ZERO=F(0)
TARGET=F(3,2)
PAIR_HS_UP=F(5,14)
ROOT_POOL=(F(1), F(2), F(1,2), F(3,2), F(1,3), F(5,2),
           F(3,5), F(7,4), F(4,3), F(2,3), F(7,5), F(1,4),
           F(5,3), F(3,4), F(8,5), F(1,5))


def positive_fraction_ldlt(mat):
    n=len(mat)
    if not n or any(len(r)!=n for r in mat): raise ValueError('square matrix required')
    a=[[F(v) for v in row] for row in mat]
    if any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('symmetric required')
    lo=[[ZERO]*n for _ in range(n)]
    piv=[ZERO]*n
    for i in range(n):
        lo[i][i]=F(1)
        for j in range(i):
            lo[i][j]=(a[i][j]-sum((lo[i][k]*piv[k]*lo[j][k] for k in range(j)),ZERO))/piv[j]
        piv[i]=a[i][i]-sum((lo[i][k]**2*piv[k] for k in range(i)),ZERO)
        if piv[i]<=0: raise ArithmeticError('nonpositive exact pivot '+str(i))
    return min(piv)


def build(roots,sizes,multiplier=F(1),gap_mode='local'):
    roots=tuple(map(F,roots)); sizes=tuple(sizes); mul=F(multiplier)
    if not roots or len(roots)!=len(sizes) or mul<1:
        raise ValueError('invalid group setup')
    if any(r<=0 for r in roots) or any(not isinstance(m,int) or m<2 for m in sizes):
        raise ValueError('root positive and group >=2')
    if gap_mode not in ('local','global'): raise ValueError('invalid gap mode')
    hs=tuple(r*r for r in roots); H=max(hs)
    xs=[]; gs=[]; last=None
    for j,(h,m) in enumerate(zip(hs,sizes)):
        first=ZERO if last is None else last+mul*(H if gap_mode=='global' else max(hs[j-1],h))
        grp=tuple(first+k*h for k in range(m))
        gs.append(tuple(range(len(xs),len(xs)+m)))
        xs.extend(grp);last=grp[-1]
    return tuple(xs),tuple(gs),hs,roots


def nearest(xs):
    if len(xs)<2 or tuple(sorted(set(xs)))!=tuple(xs):
        raise ValueError('points must be strictly increasing')
    return tuple(min(abs(xs[i]-xs[j]) for j in range(len(xs)) if j!=i) for i in range(len(xs)))


def cross_entry(xi,xj,si,sj):
    return si*sj*(si-sj)**2/(2*(xi-xj)**2)


def get_D(xs,gs,r,s):
    if not 0<=r<s<len(gs):raise ValueError('r<s required')
    return xs[gs[s][0]]-xs[gs[r][-1]]


def comparison(xs,gs,hs,roots):
    K=len(gs); A=[[ZERO]*K for _ in range(K)]; checks=[]
    for r,s in combinations(range(K),2):
        D=get_D(xs,gs,r,s)
        rho=D/max(hs[r],hs[s])
        if rho<1:raise ArithmeticError('group overlap/insufficient gap')
        if hs[r]==hs[s]:b=ZERO
        else:b=min(PAIR_HS_UP/rho,F(len(gs[r])+len(gs[s]),27)/rho**2)
        A[r][s]=A[s][r]=b
        hs_sq=sum((cross_entry(xs[i],xs[j],roots[r],roots[s])**2
                   for i in gs[r] for j in gs[s]), ZERO)
        if hs_sq>b*b:raise ArithmeticError('invalid pair envelope %d,%d'%(r,s))
        checks.append({'r':r,'s':s,'rho':str(rho),
                       'HS_squared':str(hs_sq),'rational_bound':str(b),
                       'HS_bound_integral':str(PAIR_HS_UP/rho),
                       'HS_bound_cardinality':str(F(len(gs[r])+len(gs[s]),27)/rho**2),
                       'linear_separation':rho>=s-r,
                       'three_quarters_separation':rho**4>=(s-r)**3})
    return tuple(tuple(row) for row in A),checks,max(map(len,gs))


def delta_exact(xs,gs,roots):
    node_roots=tuple(roots[r] for r,g in enumerate(gs) for _ in g)
    if tuple(v*v for v in node_roots)!=nearest(xs):
        raise ArithmeticError('nearest gaps mismatch')
    return [[ZERO if i==j else cross_entry(xs[i],xs[j],node_roots[i],node_roots[j])
             for j in range(len(xs))] for i in range(len(xs))]


def audit(roots,sizes,multiplier=F(1),full_ldlt=False,target=TARGET):
    xs,gs,hs,roots=build(roots,sizes,multiplier,gap_mode='global')
    delta=nearest(xs)
    if any(delta[i]!=hs[r] for r,g in enumerate(gs) for i in g):
        raise ArithmeticError('not nearest steps')
    A,pairs,M=comparison(xs,gs,hs,roots)
    if not all(p['linear_separation'] and p['three_quarters_separation'] for p in pairs):
        raise ArithmeticError('linear configuration invalid')
    mr=max(sum(row,ZERO) for row in A)
    if M<=14 and not F(2*M,27)**2<F(6,5):
        raise ArithmeticError('universal residual check failed')
    ldlt_A=(positive_fraction_ldlt([[F(target)*(i==j)-A[i][j] for j in range(len(A))]
                                   for i in range(len(A))])
            if mr<F(target) else None)
    ldlt_full=None
    if full_ldlt:
        delta_mat=delta_exact(xs,gs,roots)
        n=len(xs)
        ldlt_full=positive_fraction_ldlt([[F(target)*(i==j)-delta_mat[i][j]
                                          for j in range(n)] for i in range(n)])
    return {'clusters':len(gs),'points':len(xs),'max_group_size':M,
            'ratio_max_min_h':str(max(hs)/min(hs)),
            'separation_multiplier':str(F(multiplier)),
            'max_comparison_row_sum':str(mr),
            'comparison_LDL_min_pivot':None if ldlt_A is None else str(ldlt_A),
            'full_Delta_LDL_min_pivot':None if ldlt_full is None else str(ldlt_full),
            'pairs':pairs}


def star_gate(leaves=12,large=64,small_root=F(1,1000)):
    roots=(F(1),)+(F(small_root),)*leaves
    sizes=(large,)+(2,)*leaves
    xs,gs,hs,roots=build(roots,sizes,gap_mode='local')
    A,pairs,M=comparison(xs,gs,hs,roots)
    if any(A[r][s]!=ZERO for r in range(1,len(A)) for s in range(1,len(A))):
        raise ArithmeticError('comparison not a star')
    rows=[sum(row,ZERO) for row in A]
    if rows[0]<=4:raise ArithmeticError('star row did not exceed 4')
    piv=positive_fraction_ldlt([[TARGET*(i==j)-A[i][j]
                                 for j in range(len(A))] for i in range(len(A))])
    norm_sq=sum((v*v for v in A[0][1:]),ZERO)
    if norm_sq>=TARGET**2: raise ArithmeticError('star spectral gate failed')
    return {'groups':len(gs),'total_points':len(xs),'leaves':leaves,
            'max_row_sum_gt_4':True,'max_row_sum':str(rows[0]),
            'spectral_radius_squared':str(norm_sq),
            'spectral_radius_lt_3_over_2':True,
            'exact_comparison_LDL_min_pivot':str(piv),'cross_pairs':len(pairs)}


def certify():
    recipes=((2,2),(3,4,2),(2,4,6,9),(14,2,3,4,6),
             (2,3,4,6,9,14),(4,4,4,4,4,4,4,4),
             (2,3,4,6,9,14,2,3,4,6,9,14),
             (3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4))
    cases=[audit([ROOT_POOL[j%len(ROOT_POOL)] for j in range(len(sizes))],
                 sizes,full_ldlt=i<4,target=F(10))
           for i,sizes in enumerate(recipes)]
    count=sum(len(c['pairs']) for c in cases)
    if len(cases)!=8 or count!=249:raise ArithmeticError('census mismatch')
    star=star_gate()
    return {'schema':SCHEMA,
            'analytic_theorem_linear':'For any K and m_r<=M, D_rs>=|r-s| max(h_r,h_s) entails ||Delta||<=2M*pi²/81. M<=14 is below P88 residual.',
            'analytic_theorem_sublinear':'For any K, m_r<=M, D_rs>=|r-s|^(3/4) max(h_r,h_s), ||Delta||<4M/9. M<=4 is below P88 residual.',
            'analytic_theorem_adaptive':'For arbitrary finite uniform clusters with adjacent D>=max(h), rational nonnegative block A dominates Delta; exact positive LDL of T I-A proves ||Delta||<T.',
            'analytic_proof_all_real_under_hypotheses':True,
            'arbitrary_finite_cluster_count':True,
            'arbitrary_cluster_size_without_spectral_gate':False,
            'arbitrary_point_configurations_proved':False,
            'general_Preissmann_Lemma6_proved':False,
            'Riemann_Hypothesis_proved':False,
            'fraction_cases':len(cases),'pair_checks':count,
            'full_fraction_LDL_cases':sum(c['full_Delta_LDL_min_pivot'] is not None for c in cases),
            'spectral_star_counterexample_to_row_gate':star,'cases':cases}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',required=True,type=Path)
    args=ap.parse_args()
    x=certify();args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('P92L PASS:',x['fraction_cases'],'cases;',x['pair_checks'],'pairs;',
          x['full_fraction_LDL_cases'],'full LDL cases')
    print('Weighted star max row >4, exact comparison spectral norm <3/2')
    print('All unstructured configurations, general Preissmann and RH remain OPEN')
    print('SHA256:',hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=='__main__':
    main()
