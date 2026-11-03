import math
import mpmath as mp
import pytest
from es4_spectral_count import (smooth_count,energy,energy_derivative,model_count,
                                diagonal_truncation,gaussian_trace,li_like_partial,benchmark,shifted_model_count)


def test_inverse_formula_first_30():
    with mp.workdps(65):
        for n in range(1,31):
            assert abs(smooth_count(energy(n))-n) < mp.mpf('1e-59')


def test_positive_and_strictly_increasing_levels():
    with mp.workdps(60):
        v=diagonal_truncation(150)
        assert all(e>0 for e in v)
        assert all(b>a for a,b in zip(v,v[1:]))


def test_gaps_shrink_and_derivative():
    with mp.workdps(55):
        gaps=[energy(n+1)-energy(n) for n in range(1,40)]
        assert all(gaps[i+1]<gaps[i] for i in range(len(gaps)-1))
        assert energy_derivative(100000)<energy_derivative(100)<energy_derivative(1)
        for n in [1,5,20]:
            h=mp.mpf('1e-9')
            x=mp.mpf(n)-mp.mpf(7)/8
            def f(z):return 2*mp.pi*(z-mp.mpf(7)/8)/mp.lambertw((z-mp.mpf(7)/8)/mp.e)
            derivative=(f(n+h)-f(n-h))/(2*h)
            assert abs(derivative-energy_derivative(n))<mp.mpf('1e-16')


def test_exact_counting_by_direct_enumeration():
    with mp.workdps(50):
        for t in [1,10,18,20,50,100,300,500]:
            actual=sum(energy(n)<=t for n in range(1,500))
            assert model_count(t)==actual


def test_counting_gap_below_one():
    with mp.workdps(60):
        for n in [1,2,10,100,1000]:
            t=(energy(n)+energy(n+1))/2
            assert model_count(t)==n
            assert 0<=smooth_count(t)-n<1


def test_diagonal_selfadjoint_finite_truncations():
    with mp.workdps(55):
        for n in [1,10,100]:
            levels=diagonal_truncation(n)
            assert len(levels)==n
            assert all(mp.im(e)==0 for e in levels)


def test_zero_reference_not_input_to_model():
    with mp.workdps(50):
        z=mp.im(mp.zetazero(1))
        assert abs(energy(1)-z)>1


def test_finite_gaussian_not_zero_trace_identity():
    with mp.workdps(60):
        mod=diagonal_truncation(10)
        zs=[mp.im(mp.zetazero(n)) for n in range(1,11)]
        assert gaussian_trace(mod,30)>0
        assert abs(gaussian_trace(mod,30)-gaussian_trace(zs,30))>mp.mpf('1e-4')


def test_li_like_positivity_is_tautology():
    with mp.workdps(55):
        for k in [1,2,5,10,20]:
            assert li_like_partial(diagonal_truncation(15),k)>=0
            assert li_like_partial([3,7,11,15],k)>=0


def test_invalid_inputs():
    for bad in [0,-1,2.5,'1']:
        with pytest.raises(ValueError):energy(bad)
    with pytest.raises(ValueError):smooth_count(0)
    with pytest.raises(ValueError):diagonal_truncation(0)
    with pytest.raises(ValueError):gaussian_trace([1,2],0)
    with pytest.raises(ValueError):li_like_partial([1],0)
    with pytest.raises(ValueError):benchmark(3,60)
    with pytest.raises(ValueError):benchmark(10,10)


def test_benchmark_schema_and_limits():
    r=benchmark(8,45)
    assert 'not RH proof' in r['status']
    assert r['j']==8
    assert len(r['samples'])>8
    assert len(r['finite_gaussian_traces'])==3
    assert abs(float(r['rmse_level_error']))>0
    assert len(r['li_like_partial'])==4


def test_shifted_operator_same_asymptotic_different_levels():
    with mp.workdps(60):
        for t in [20,50,100,500,1000]:
            assert shifted_model_count(t,1)==sum(energy(n)+1<=t for n in range(1,1100))
        assert abs((energy(3)+1)-energy(3))==1
        with pytest.raises(ValueError):shifted_model_count(100,-1)
