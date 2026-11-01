#!/usr/bin/env python3
"""Exact rational Gibbs Gram matrices and an obstruction for short-support Weil modes."""
from fractions import Fraction as Q
from gibbs_euler import log_interval,exact_decimal_bound


def gram(beta:int,dim:int,number_of_states:int):
    if beta<=1 or dim<1 or number_of_states<1:raise ValueError('domain')
    weights=[Q(1,n**beta) for n in range(1,number_of_states+1)]
    z=sum(weights,Q(0))
    # features: f_i(n)=(1/n)^i, i=0..dim-1; all real and exactly rational.
    return [[sum((weights[n-1]*Q(1,n**(i+j)) for n in range(1,number_of_states+1)),Q(0))/z
             for j in range(dim)] for i in range(dim)]


def ldl_symmetric(matrix):
    dim=len(matrix)
    if any(len(row)!=dim for row in matrix):raise ValueError('not square')
    if any(matrix[i][j]!=matrix[j][i] for i in range(dim) for j in range(dim)):
        raise ValueError('not symmetric')
    L=[[Q(int(i==j)) for j in range(dim)] for i in range(dim)]
    d=[]
    for i in range(dim):
        pivot=matrix[i][i]-sum((L[i][k]**2*d[k] for k in range(i)),Q(0))
        if pivot<=0:raise ValueError(f'nonpositive LDL pivot at i={i}')
        d.append(pivot)
        for j in range(i+1,dim):
            L[j][i]=(matrix[j][i]-sum((L[j][k]*L[i][k]*d[k] for k in range(i)),Q(0)))/pivot
    for i in range(dim):
        for j in range(dim):
            reconstructed=sum((L[i][k]*d[k]*L[j][k] for k in range(dim)),Q(0))
            assert reconstructed==matrix[i][j], 'LDL reconstruction failure'
    return d


def short_support_obstruction():
    L=Q(9,20)
    lower2,upper2=log_interval(2)
    lower3,upper3=log_interval(3)
    assert Q(0)<lower2<upper2<2*L<lower3<upper3
    # g_k(log 1)=g_k(0)=sin 0=0; only n=2 can be nonzero.
    # Every finite Gibbs/Gram sum equals a single nonnegative rank-1 outer product.
    return {'L':'9/20','2L':'9/10', 'log2_upper':exact_decimal_bound(upper2,True),'log3_lower':exact_decimal_bound(lower3,False),
            'allowed_state_indices':[1,2], 'nonzero_possible_states':[2],
            'rank_upper':1,'projected_rank_upper':1,
            'consequence':'Cannot prove positive definiteness on a projected dimension >=2 using this direct Gibbs Gram kernel.'}


def run():
    cases=[]
    for beta in (2,3):
        for dim,M in ((3,5),(4,8),(6,9)):
            G=gram(beta,dim,M)
            pivots=ldl_symmetric(G)
            assert all(p>0 for p in pivots)
            cases.append({'beta':beta,'dim':dim,'states':M,'ldl_pivots_lower':[exact_decimal_bound(x,False) for x in pivots],
                          'exact_reconstruction':True})
    # Deliberate negative Gibbs weight invalidates Gram positivity.
    bad=gram(3,3,5)
    bad[0][0]=Q(-1)
    try:ldl_symmetric(bad)
    except ValueError: pass
    else:raise AssertionError('negative-definite mutation was not rejected')
    return {'status':'GIBBS_GRAM_EXACT_LDL_PASS_WITH_SHORT_SUPPORT_RANK_OBSTRUCTION',
            'positive_definite_independent_cases':cases,
            'short_support_obstruction':short_support_obstruction(),
            'assumptions':['positive Gibbs weights beta>1','finite rational Gram feature family'],
            'not_equal_to_weil':True}
