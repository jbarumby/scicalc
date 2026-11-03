"""Tests of ES7 interval arithmetic + independent mpmath diagnostics."""
from fractions import Fraction as F
import mpmath as mp
import pytest
from es7_archimedean import (S,up,enc,mul,pi_interval,ln2_interval,lnpi_interval,
    gamma_interval,omega0_interval,omega_interval,exp_negative_interval,
    gaussian_tail,arch_certificate)


def within(b,x):
    return mp.mpf(b[0])/S <= x <= mp.mpf(b[1])/S


def test_01_arithmetic_directed():
    assert up(3,2)==2 and up(-3,2)==-1
    assert enc(F(1,3))==(S//3,up(S,3))
    assert mul(enc(F(-2,3)),enc(F(5,7)))[0] <= -10*S//21


def test_02_pi_interval():
    with mp.workdps(60): assert within(pi_interval(),mp.pi)


def test_03_logs():
    with mp.workdps(60):
        assert within(ln2_interval(),mp.log(2))
        assert within(lnpi_interval(pi_interval()),mp.log(mp.pi))


def test_04_gamma_and_special_value():
    with mp.workdps(60):
        assert within(gamma_interval(),mp.euler)
        assert within(omega0_interval(),mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi))


def test_05_digamma_series_intervals():
    base=omega0_interval()
    with mp.workdps(55):
        for t in (F(0),F(1,10),F(1),F(8),F(24),F(48)):
            val=mp.re(mp.digamma(mp.mpf(1)/4+mp.mpf(t.numerator)/t.denominator*1j/2))-mp.log(mp.pi)
            assert within(omega_interval(t,192,base),val)


def test_06_exp_scaling():
    with mp.workdps(60):
        for x in (F(0),F(1,20),F(1),F(3),F(20),F(48)):
            assert within(exp_negative_interval(x),mp.exp(-mp.mpf(x.numerator)/x.denominator))


def test_07_monotone_arch_convergence():
    slow=arch_certificate(F(1,50),M=256,terms=64,T=48)
    finer=arch_certificate(F(1,50),M=512,terms=128,T=48)
    a=list(map(int,slow['interval_scaled_exact']))
    b=list(map(int,finer['interval_scaled_exact']))
    assert max(a[0],b[0]) <= min(a[1],b[1])
    assert b[1]-b[0] < a[1]-a[0]


def test_08_contains_arch_gaussian_and_differential_independently():
    x=arch_certificate(F(1,50),M=1024,terms=128,T=48)
    with mp.workdps(45):
        a=mp.mpf(1)/50
        def omega(t):return mp.re(mp.digamma(mp.mpf(1)/4+1j*t/2))-mp.log(mp.pi)
        knots=[0,1,3,7,15,30,60]
        arch=mp.quad(lambda t:mp.exp(-a*t*t)*omega(t)/mp.pi,knots)
        darch=mp.quad(lambda t:-t*t*mp.exp(-a*t*t)*omega(t)/mp.pi,knots)
        assert within(tuple(map(int,x['interval_scaled_exact'])),arch)
        assert within(tuple(map(int,x['derivative_interval_scaled_exact'])),darch)


def test_09_parameter_dependence():
    for a in (F(1,50),F(1,25)):
        x=arch_certificate(a,M=512,terms=128,T=48)
        lo,hi=map(int,x['interval_scaled_exact'])
        dl,dh=map(int,x['derivative_interval_scaled_exact'])
        assert lo<hi and dl<dh
        assert x['status'].startswith('CERTIFIED')


def test_10_tail_is_positive_tiny():
    t1=gaussian_tail(F(1,50),48,pi_interval()[0])
    t2=gaussian_tail(F(1,25),48,pi_interval()[0])
    assert 0<t2<=t1 < S//10**11


def test_11_invalid_inputs():
    for m,t,T in [(4,128,48),(128,12,48),(128,100,8)]:
        with pytest.raises(ValueError):arch_certificate(F(1,50),M=m,terms=t,T=T)
    with pytest.raises(ValueError):arch_certificate(F(-1,50),M=128)


def test_12_scientific_scope_not_global():
    x=arch_certificate(F(1,50),M=128,terms=32,T=48)
    assert 'ONLY' in x['status']
    assert 'ONLY' in x['derivative_status']
