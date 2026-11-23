"""P92-M: rational regression for block comparison and the 1/rho obstruction.

The ALL-real theorem is proved analytically in README.md. Finite Fraction
regressions audit arithmetic and do not establish the theorem themselves.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

SCHEMA = 'NUMERIA.P92M.SUBQUADRATIC.M_FREE.SPECTRAL.1'
T = F(18, 5)


def make_configuration(roots, sizes, gap_factor=3):
    roots = tuple(F(x) for x in roots)
    sizes = tuple(sizes)
    if len(roots) != len(sizes) or not roots or any(x <= 0 for x in roots) or any(not isinstance(m,int) or m < 2 for m in sizes):
        raise ValueError('at least one uniform group; roots>0; sizes>=2')
    if F(gap_factor) < 1:
        raise ValueError('gap factor must be >=1')
    steps = [s*s for s in roots]
    max_step = max(steps)
    boundary = F(0)
    xs, ids, group_idx = [], [], []
    for g,(h,m) in enumerate(zip(steps,sizes)):
        indices = []
        for j in range(m):
            indices.append(len(xs))
            xs.append(boundary+j*h)
            ids.append(g)
        group_idx.append(indices)
        boundary=xs[-1]+F(gap_factor)*max_step
    dd = [min(abs(xs[i]-xs[j]) for j in range(len(xs)) if j!=i) for i in range(len(xs))]
    if any(dd[i] != steps[ids[i]] for i in range(len(xs))):
        raise ArithmeticError('incorrect nearest-neighbor gap')
    return xs,ids,group_idx,roots,steps


def delta_matrix(xs, ids, roots):
    n=len(xs)
    return [[F(0) if i==j else
             roots[ids[i]]*roots[ids[j]]*(roots[ids[i]]-roots[ids[j]])**2 /
             (2*(xs[i]-xs[j])**2)
             for j in range(n)] for i in range(n)]


def exact_ldl_positive(A):
    """Fraction-only unpivoted LDL test for symmetric positive definite A."""
    B=[list(map(F,row)) for row in A]
    n=len(B)
    if any(len(row)!=n for row in B) or any(B[i][j]!=B[j][i] for i in range(n) for j in range(i)):
        raise ValueError('matrix must be symmetric square')
    least=None
    for k in range(n):
        pivot=B[k][k]
        if pivot<=0:
            raise ArithmeticError(f'nonpositive exact pivot {k}')
        least=pivot if least is None else min(least,pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                B[j][i] -= B[i][k]*B[j][k]/pivot
                B[i][j]=B[j][i]
    return least


def audit_case(name, roots, sizes, gap_factor=3, full_ldl=False):
    xs, ids, grp, roots, steps=make_configuration(roots,sizes,gap_factor)
    M=delta_matrix(xs,ids,roots)
    K=len(grp)
    block_bounds = [[F(0) for _ in range(K)] for _ in range(K)]
    pairs=0
    for r in range(K):
        for s in range(r+1,K):
            a,b=steps[r],steps[s]
            D=xs[grp[s][0]]-xs[grp[r][-1]]
            H=max(a,b)
            d=s-r
            # Exact check D >= (s-r)^(5/4) H without algebraic roots.
            if D**4 < d**5 * H**4:
                raise ArithmeticError('P92-M separation condition violated')
            rho=D/H
            hs=sum((M[i][j]**2 for i,j in product(grp[r],grp[s])),F(0))
            # Pair bound imported from analytic P92-K: HS^2 <= 1/(8rho^2).
            if hs > 1/(8*rho**2):
                raise ArithmeticError('Hilbert-Schmidt pair comparison violated')
            # Fully rational upper bound 5/(14rho) > 1/(2 sqrt2 rho).
            B=F(5,14)/rho
            block_bounds[r][s]=block_bounds[s][r]=B
            if hs > B*B:
                raise ArithmeticError('fractional block enclosure violated')
            pairs+=1
    # If K=1, zero matrix; otherwise A is nonnegative symmetric block comparison.
    positive=exact_ldl_positive([[T*(r==s)-block_bounds[r][s] for s in range(K)] for r in range(K)])
    full=exact_ldl_positive([[T*(i==j)-M[i][j] for j in range(len(M))] for i in range(len(M))]) if full_ldl else None
    return {'name':name,'K':K,'N':len(xs),'distinct_delta_scales':len(set(steps)),
            'pair_interactions':pairs,'largest_comparison_row_sum':str(max(map(sum,block_bounds))),
            'comparison_ldl_min_pivot':str(positive),
            'full_ldl_min_pivot':None if full is None else str(full)}


def obstruction(n, exact_sum=False):
    """Two groups: a=1, b=n^-2, D=n, m=n, M=n^3.

    For all j,k, D+j+k/n^2<3n, so one gets a finite-dimensional
    operator-norm lower bound (1-1/n)^2/(18n) by constant unit vectors.
    """
    if not isinstance(n,int) or n < 2:
        raise ValueError('n>=2')
    b=F(1,n*n)
    beta=F(1,n)*(1-F(1,n))**2/2
    lower=(1-F(1,n))**2/(18*n)
    # A strict entrywise lower bound is obtained from D+j+k*b<3n.
    min_entry=beta/(F(n)+n-1+(n**3-1)*b)**2
    if not (min_entry > beta/(9*n*n) and lower>=F(1,72*n)):
        raise ArithmeticError('invalid obstruction bound')
    value=None
    if exact_sum:
        value = sum((beta/(F(n)+j+k*b)**2 for j in range(n) for k in range(n**3)),F(0))/n**2
        if not value>=lower:
            raise ArithmeticError('Rayleigh lower bound failed')
    return {'n':n,'rho':n,'a':'1','b':str(b),'m_anchor':n,'m_fine':n**3,
            'analytic_lower_on_operator_norm':str(lower),
            'universal_weaker_lower':str(F(1,72*n)),
            'bilinear_exact':None if value is None else str(value)}


def certify():
    cases=[]
    configurations=[
        ('two_balanced',[F(1),F(1,2)],[5,7],3,True),
        ('three_contrasts',[F(1),F(1,3),F(1,5)],[2,3,4],3,True),
        ('four_irregular',[F(2),F(1,2),F(1),F(1,4)],[3,4,2,5],3,True),
        ('five_small',[F(1),F(2),F(1,3),F(3),F(1,2)],[2,2,3,2,4],3,True),
        ('equal_scales',[F(1),F(1),F(1),F(1)],[2,3,4,5],3,True),
        ('six_mixed',[F(1),F(1,2),F(1,3),F(2),F(1,4),F(3)],[2,3,2,2,2,2],3,True),
        ('eight_small',[F(1),F(1,3),F(1,2),F(2),F(1,5),F(1),F(3),F(1,4)],[2]*8,3,True),
        ('twelve_groups',[F(1),F(1,2),F(1,3)]*4,[2]*12,2,False),
        ('two_large',[F(1),F(1,30)],[90,25],3,False),
        ('three_unbalanced',[F(1),F(1,100),F(3)],[55,2,3],3,False),
    ]
    for name, roots, sizes, gap, full_ldl in configurations:
        cases.append(audit_case(name,roots,sizes,gap,full_ldl))
    obstructions=[obstruction(n,exact_sum=n<=5) for n in (2,3,4,5,8,12)]
    if len(cases)!=10 or sum(x['full_ldl_min_pivot'] is not None for x in cases)!=7 or len(obstructions)!=6:
        raise ArithmeticError('invalid regression census')
    # Machin pi>25/8 with alternating atan: lower atan(1/5)>1/5-1/(3*125)
    # and upper atan(1/239)<1/239.
    pi_lower=16*(F(1,5)-F(1,375))-4*F(1,239)
    if not (pi_lower>F(25,8) and F(25,8)**4>F(375,4)):
        raise ArithmeticError('failed rational strict comparison 5/sqrt2 vs P88 target')
    return {'schema':SCHEMA,
        'theorem_p92m1':'For arbitrary K,m_r>=2,h_r>0,D_rs >= abs(r-s)^(5/4)*max(h_r,h_s), ||Delta||_2 < 5/sqrt(2) < pi^2/3*sqrt(6/5)',
        'theorem_p92m2':'No exponent 1+eps can replace 1/rho in a uniform two-block operator norm bound independent of block sizes and step ratios',
        'unrestricted_point_configurations_proved':False,
        'general_Preissmann_Lemma6_proved':False,
        'Riemann_Hypothesis_proved':False,
        'exact_machin_pi_lower':str(pi_lower),
        'rational_configurations':len(cases),
        'full_exact_ldl':sum(x['full_ldl_min_pivot'] is not None for x in cases),
        'rational_pair_checks':sum(x['pair_interactions'] for x in cases),
        'cases':cases,'sharp_exponent_obstructions':obstructions}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    obj=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('P92M PASS:',obj['rational_configurations'],'configurations,',obj['rational_pair_checks'],
          'pairs,',obj['full_exact_ldl'],'full LDLs and',len(obj['sharp_exponent_obstructions']),'obstructions')
    print('SHA256',hashlib.sha256(args.out.read_bytes()).hexdigest())
    print('Global unrestricted Preissmann and RH remain OPEN')

if __name__=='__main__':
    main()
