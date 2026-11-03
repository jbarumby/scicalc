"""Tests distinguish algebraic/interval proofs from independent floating diagnostics."""
import sys
from pathlib import Path
from fractions import Fraction as Q
import pytest
import mpmath as mp
from es8p2_prime_cert import (L,SCALE,M,I,Iv,atan_bounds,atanh_bounds,constants,
    sincos_rational_at,trig_mode,overlap_sym,raw_prime,projection_column,
    projected_prime,audit,format_bound)

def mpq(q):return mp.mpf(q.numerator)/q.denominator

def within(v:Iv,x):return mpq(v.lo)<=x<=mpq(v.hi)

def independent_mode(k,x):
    l=mp.mpf(9)/20
    return mp.sin(k*mp.pi*(x+l)/(2*l))/mp.sqrt(l)

def independent_overlap(i,j):
    l=mp.mpf(9)/20;u=mp.log(2)
    return mp.quad(lambda x:independent_mode(i,x)*independent_mode(j,x-u),[-l+u,l])

def numeric_prime(i,j):
    u=mp.log(2)
    return -2*u/mp.sqrt(2)*(independent_overlap(i,j)+independent_overlap(j,i))/2

def numeric_projected(i,j):
    l=mp.mpf(9)/20; pi=mp.pi
    def column(k):
        p=1 if k%2 else 2
        b=-mp.mpf(k)/p*(p*p*pi*pi+l*l)/(k*k*pi*pi+l*l)
        return {p:b,k:mp.mpf(1)}
    return mp.fsum(a*b*numeric_prime(u,v) for u,a in column(i).items() for v,b in column(j).items())

def test_odd_and_even_alternating_arctan_directions():
    with mp.workdps(80):
        for terms in [44,45,46,47]:
            iv=atan_bounds(Q(1,5),terms)
            assert within(iv,mp.atan(mp.mpf(1)/5))

def test_constants_enclose_direct_mp():
    c=constants()
    with mp.workdps(80):
        for name,target in [('pi',mp.pi),('ln2',mp.log(2)),('ln3',mp.log(3)),('sqrt2',mp.sqrt(2))]:
            assert within(c[name],target),name
        assert c['ln2'].hi<Q(9,10)<c['ln3'].lo

def test_taylor_remainder_sine_cosine():
    with mp.workdps(80):
        for q in [Q(0),Q(1,3),Q(7,1),Q(24,1)]:
            s,c=sincos_rational_at(q)
            assert within(s,mp.sin(mpq(q)))
            assert within(c,mp.cos(mpq(q)))

def test_taylor_trig_all_ten_modes():
    with mp.workdps(70):
        for k in range(1,11):
            a,b=trig_mode(k); theta=mp.mpf(k)*mp.pi*mp.log(2)/(2*mpq(L))
            assert within(a,mp.sin(theta)), k
            assert within(b,mp.cos(theta)),k

def test_reflection_parity_exact_zeros():
    for i in range(1,11):
        for j in range(1,11):
            if (i-j)%2:
                assert overlap_sym(i,j)==Iv.point(0)
                assert raw_prime(i,j)==Iv.point(0)

def test_real_symmetric_raw_prime():
    for i in range(1,11):
        for j in range(1,11):assert raw_prime(i,j)==raw_prime(j,i)

def test_raw_prime_numeric_independent_checks():
    with mp.workdps(65):
        for i,j in [(1,1),(1,3),(2,4),(3,3),(3,5),(4,4)]:
            assert within(raw_prime(i,j),numeric_prime(i,j)), (i,j)

def test_sample_prime_diagonal_and_offdiagonal_signs():
    assert raw_prime(1,1).hi<0
    assert raw_prime(1,3).hi<0
    assert raw_prime(2,4).lo>0

def test_pole_projector_rational_identity_for_sample_pi_square():
    for X in [Q(9),Q(10),Q(11)]:
        for k in range(3,11):
            p=1 if k%2 else 2
            Ck=Q(k)/(k*k*X+L*L)
            Cp=Q(p)/(p*p*X+L*L)
            b=-Q(k,p)*(p*p*X+L*L)/(k*k*X+L*L)
            assert Ck+b*Cp==0

def test_projection_coefficients_contain_numeric_values():
    c=constants()
    with mp.workdps(75):
        for k in range(3,11):
            p=1 if k%2 else 2
            b=-mp.mpf(k)/p*(p*p*mp.pi**2+mpq(L)**2)/(k*k*mp.pi**2+mpq(L)**2)
            assert within(projection_column(k,10)[p],b)

def test_projected_prime_parity_and_symmetry():
    for N in [6,8,10]:
        for i in range(3,N+1):
            for j in range(3,N+1):
                assert projected_prime(i,j,N)==projected_prime(j,i,N)
                if (i-j)%2:assert projected_prime(i,j,N)==Iv.point(0)

def test_projected_prime_interval_independent_checks():
    with mp.workdps(65):
        for i,j in [(3,3),(3,5),(4,4),(6,8)]:
            assert within(projected_prime(i,j,10),numeric_projected(i,j)),(i,j)

def test_projected_diagonal_negative_example_not_full_weil():
    assert projected_prime(3,3,10).hi<0
    # This is the PRIME-ONLY matrix, and negativity says nothing about the
    # full Weil matrix (archimedean plus prime plus pole).

def test_negative_control_wrong_sign_or_factor():
    with mp.workdps(55):
        proper=numeric_prime(1,1)
        wrong=-mp.log(2)*independent_overlap(1,1)
        assert abs(proper-wrong)>mp.mpf('0.005')

def test_halfline_cross_phase_is_pure_imaginary_not_real():
    with mp.workdps(55):
        l=mp.mpf(9)/20
        def FT(k,r):return mp.quad(lambda x:independent_mode(k,x)*mp.exp(1j*r*x),[-l,l])
        p=FT(1,1)*mp.conj(FT(2,1))
        assert abs(mp.re(p))<mp.mpf('1e-45')
        assert abs(mp.im(p))>mp.mpf('0.1')

def test_audit_reproducible_and_no_rh_claim():
    a=audit();b=audit()
    assert a==b
    assert a['schema']=='NUMERIA.ES8P2.PRIME_PROJECTION_INTERVAL.1'
    assert 'NOT full Weil' in a['scientific_status']
    assert set(a['finite_prime_matrices'])=={'6','8','10'}

def test_bad_inputs():
    with pytest.raises(ValueError):trig_mode(11)
    with pytest.raises(ValueError):overlap_sym(0,2)
    with pytest.raises(ValueError):projected_prime(2,3,6)
    with pytest.raises(ZeroDivisionError):I(1)/Iv(Q(-1),Q(1))

def test_legacy_45_term_atan_direction_was_wrong_but_rounded_pi_contains_fixed():
    # Prior ES8-P's atan_bounds(1/5,45) returned [S,S+R] for odd 45,
    # instead of [S-R,S]. Check its *rounded composite pi interval* with
    # scale 10^30 still contains a rigorously enclosed corrected pi.
    terms=45
    def wrong(q):
        s=sum(((-1)**j*q**(2*j+1)/Q(2*j+1) for j in range(terms)),Q())
        r=q**(2*terms+1)/Q(2*terms+1)
        return s,s+r
    a,b=wrong(Q(1,5));c,d=wrong(Q(1,239))
    S=10**30
    oldlo=Q((16*a-4*d).numerator*S//(16*a-4*d).denominator,S)
    x=16*b-4*c
    oldhi=Q(-((-x.numerator*S)//x.denominator),S)
    fixed=constants()['pi']
    assert oldlo<=fixed.lo<=fixed.hi<=oldhi
