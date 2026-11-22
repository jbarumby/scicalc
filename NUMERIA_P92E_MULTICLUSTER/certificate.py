"""P92-E exact rational block-HS certificate, for uniform multi-cluster arrays.

The analytical theorem for all real parameters is proved in README.md.
Fraction sample checks DO NOT replace the theorem.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

SCHEMA = 'NUMERIA.P92E.DELTA.MULTICLUSTER.BLOCK_HS.1'
FZ = F(0)
FONE = F(1)


def sqrt_rational_upper(q, digits=25):
    """Exact rigorous *upper* rational enclosure sqrt(q), q>=0."""
    q=F(q)
    if q<0 or digits<1:
        raise ValueError('invalid radicand/precision')
    if not q:
        return FZ
    n=10**digits
    floored=(q.numerator*n*n)//q.denominator
    lo=isqrt(floored)
    return F(lo+1,n)


def two_cluster_coefficient(s1,s2):
    """Exact rational kappa, where s1=sqrt(h1), s2=sqrt(h2) positive."""
    s1,s2=F(s1),F(s2)
    if s1<=0 or s2<=0:
        raise ValueError('positive square-roots required')
    return s1*s2*(s1-s2)**2/2


def pair_hs_square_exact(s1,s2,n1,n2,gap):
    """Full finite, exact, rational Hilbert-Schmidt square of cross-block Delta."""
    s1,s2,gap=F(s1),F(s2),F(gap)
    if min(s1,s2)<=0 or not all(type(n) is int and n>=2 for n in (n1,n2)):
        raise ValueError('clusters have >=2 nodes, sqrt gaps positive')
    h1,h2=s1*s1,s2*s2
    if gap<max(h1,h2):
        raise ValueError('boundary gap must be >= both internal gaps')
    k=two_cluster_coefficient(s1,s2)
    return k*k*sum((FONE/(gap+i*h1+j*h2)**4
                    for i in range(n1) for j in range(n2)),FZ)


def pair_hs_integral_ceiling_square(s1,s2,gap):
    """Exact rational upper bound from infinite double monotone sum."""
    s1,s2,gap=F(s1),F(s2),F(gap)
    if min(s1,s2)<=0 or gap<max(s1*s1,s2*s2):
        raise ValueError('gap>=max(h1,h2)>0 required')
    h1,h2=s1*s1,s2*s2
    k=two_cluster_coefficient(s1,s2)
    return k*k*(FONE/gap**4+(FONE/h1+FONE/h2)/(3*gap**3)+FONE/(6*h1*h2*gap**2))


def build_uniform_clusters(sqrt_gaps, sizes, boundary_multipliers):
    """Build exact rational nodes / deltas; gap_{r,r+1}=multiplier*max(h_r,h_{r+1})."""
    ss=tuple(F(s) for s in sqrt_gaps)
    nn=tuple(sizes)
    mult=tuple(F(t) for t in boundary_multipliers)
    k=len(ss)
    if not k or len(nn)!=k or len(mult)!=k-1:
        raise ValueError('shape mismatch')
    if any(s<=0 for s in ss) or any(type(n) is not int or n<2 for n in nn) or any(t<1 for t in mult):
        raise ValueError('positive uniform gaps, sizes>=2 and multipliers>=1 required')
    hs=tuple(s*s for s in ss)
    blocks=[]
    start=FZ
    for i in range(k):
        points=tuple(start+j*hs[i] for j in range(nn[i]))
        blocks.append(points)
        if i<k-1:
            start=points[-1]+mult[i]*max(hs[i],hs[i+1])
    if any(blocks[i][-1]>=blocks[i+1][0] for i in range(k-1)):
        raise ArithmeticError('overlap')
    return ss,nn,hs,tuple(blocks)


def block_certificate(sqrt_gaps,sizes,boundary_multipliers, digits=25):
    ss,nn,hs,blocks=build_uniform_clusters(sqrt_gaps,sizes,boundary_multipliers)
    n=len(ss)
    rows=[FZ]*n
    pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            S=blocks[j][0]-blocks[i][-1]
            if S<max(hs[i],hs[j]):
                raise ArithmeticError('separation >= both gaps did not propagate')
            actual=pair_hs_square_exact(ss[i],ss[j],nn[i],nn[j],S)
            bound=pair_hs_integral_ceiling_square(ss[i],ss[j],S)
            if not (FZ<=actual<=bound<=F(1,8)):
                raise ArithmeticError('invalid pairwise HS analytic bound')
            upper=sqrt_rational_upper(bound,digits=digits)
            if not (upper*upper>bound or bound==0):
                raise ArithmeticError('invalid sqrt upper enclosure')
            rows[i]+=upper
            rows[j]+=upper
            pairs.append({'r':i,'s':j,'boundary_distance':str(S),
                          'finite_HS_square':str(actual),
                          'integral_HS_ceiling_square':str(bound),
                          'rational_sqrt_upper':str(upper)})
    row_bound=max(rows)
    # Restricted Preissmann sufficient norm threshold: global residual bound.
    # For K<=11 the analytic bound alone suffices, using pi > 25/8.
    universal_cluster_count_gate=(n<=11)
    if universal_cluster_count_gate:
        if not (F(25,8)**4>F(375,4)):
            raise ArithmeticError('incorrect rigorous pi comparison')
    return {'cluster_count':n,'node_count':sum(nn),
            'sqrt_internal_gaps':[str(s) for s in ss],
            'nodes_per_cluster':list(nn),
            'boundary_multipliers':[str(m) for m in boundary_multipliers],
            'pairs':pairs,'row_sum_upper_exact':str(row_bound),
            'global_restricted_theorem_applies':universal_cluster_count_gate,
            'pair_count':len(pairs),
            'max_finite_HS_square':str(max((F(t['finite_HS_square']) for t in pairs),default=FZ)),
            'min_pairwise_margin_to_1_over_8_exact':str(min((F(1,8)-F(t['integral_HS_ceiling_square']) for t in pairs),default=F(1,8))),
            'row_enclosure_below_target_sufficient_by_25_over_8':row_bound*row_bound < F(2,15)*F(25,8)**4}


def certify():
    cases=[
        ('single_uniform', (F(1),), (7,), ()),
        ('two_equal_gaps', (F(1), F(1)), (7,15), (F(1),)),
        ('two_unequal_adjacent', (F(1),F(1,2)), (10,15), (F(1),)),
        ('two_widely_spaced', (F(1),F(1,10)), (4,100), (F(5),)),
        ('three_mixed', (F(1),F(1,3),F(2)), (5,7,3), (F(1),F(1))),
        ('five_mixed', (F(1),F(1,2),F(2),F(1,3),F(3,2)), (4,6,5,8,4),(F(1),)*4),
        ('eleven_close', (F(1),F(1,2),F(2),F(1,3),F(3,2),F(1),F(1,2),F(2),F(1,3),F(3,2),F(1)), (2,3,2,4,3,2,3,2,4,3,2), (F(1),)*10),
        ('fourteen_far', (F(1),F(1,2),F(2),F(1,3),F(3,2),F(1),F(1,2),F(2),F(1,3),F(3,2),F(1),F(1,2),F(2),F(1,3)),(2,)*14,(F(50),)*13),
    ]
    rows=[]
    for name,sizesqrt,sizes,mult in cases:
        result=block_certificate(sizesqrt,sizes,mult)
        if result['pair_count']!=len(sizesqrt)*(len(sizesqrt)-1)//2:
            raise ArithmeticError('wrong pair enumeration')
        if result['cluster_count']<=11 and not result['row_enclosure_below_target_sufficient_by_25_over_8']:
            raise ArithmeticError('example failed restricted-preissmann gate')
        result['name']=name
        rows.append(result)
    return {
        'schema':SCHEMA,
        'analytic_theorem':'K uniform separated clusters, arbitrary real gap sizes and lengths; ||Delta|| <= (K-1)/sqrt8',
        'pairwise_analytic_HS_norm_bound':'1/sqrt8',
        'cluster_count_unconditional_residual_gate':11,
        'pi_lower_used_for_comparison':'25/8',
        'full_Preissmann_constant_for_up_to_11_uniform_clusters_under_gap_condition':True,
        'all_arbitrary_spacings_global_Delta_bound_proved':False,
        'arbitrary_number_of_clusters_global_Delta_bound_proved':False,
        'full_Preissmann_Lemma6_proved':False,
        'Riemann_Hypothesis_proved':False,
        'case_count':len(rows),
        'finite_exact_pair_count':sum(x['pair_count'] for x in rows),
        'cases':rows,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',required=True,type=Path)
    args=ap.parse_args()
    result=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    digest=hashlib.sha256(args.out.read_bytes()).hexdigest()
    print('P92E PASS',result['case_count'],'cases',result['finite_exact_pair_count'],'rational block pairs')
    print('SHA-256',digest)
    print('Scope: <=11 separated uniform clusters; NOT general Preissmann or RH')

if __name__=='__main__':
    main()
