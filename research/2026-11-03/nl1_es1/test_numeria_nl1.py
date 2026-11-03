"""pytest -q test_numeria_nl1.py"""
from fractions import Fraction as Q
import json
import pytest
from numeria_nl1 import Interval,krawczyk_2d,numerical_spectral_newton,zeta_newton_2d,v24_checks,perturbation_gate_exact


def test_interval_products_all_signs():
    assert Interval(Q(-2),Q(3))*Interval(Q(-4),Q(5)) == Interval(Q(-12),Q(15))


def test_exact_krawczyk_first():
    p=krawczyk_2d([Q(2),Q(1,2)],[Q(1,1000),Q(1,1000)])
    assert p['certified_unique_root_in_box']
    assert p['strict_inclusion']


def test_exact_krawczyk_second():
    p=krawczyk_2d([Q(1,2),Q(2)],[Q(1,1000),Q(1,1000)])
    assert p['certified_unique_root_in_box']
    assert Q(p['contraction_infinity_exact'])<1


def test_box_not_containing_zero_rejected():
    p=krawczyk_2d([Q(3),Q(1,2)],[Q(1,1000),Q(1,1000)])
    assert not p['certified_unique_root_in_box']


def test_singular_center_rejected():
    with pytest.raises(ValueError):
        krawczyk_2d([Q(1),Q(1)],[Q(1,100),Q(1,100)])


def test_bad_width_rejected():
    with pytest.raises(ValueError):
        krawczyk_2d([Q(2),Q(1,2)],[Q(0),Q(1,100)])


def test_newton_separates_two_roots():
    a,b,it,res=numerical_spectral_newton(1.8,.7)
    assert abs(a-2)<1e-12 and abs(b-.5)<1e-12
    assert res<1e-12
    a,b,it,res=numerical_spectral_newton(.7,1.8)
    assert abs(a-.5)<1e-12 and abs(b-2)<1e-12


def test_zeta_first_root_numerical_only():
    result=zeta_newton_2d(.62,14.10,digits=50)
    assert result['converged_numerically']
    assert not result['certified_by_interval']
    assert abs(float(result['sigma'])-.5)<1e-25
    assert abs(float(result['t'])-14.134725141734693)<1e-12


def test_v24_envelope_if_fixture_present():
    from pathlib import Path
    p=Path(__file__).parent/'EXTRATO_V24_COTAS.json'
    if not p.exists():
        pytest.skip('optional V24 fixture absent')
    r=v24_checks(p)
    assert r['quadratic_scaling_consistent_to_1e-25']
    assert not r['observations'][2]['under_illustrative_margin_0p001']
    assert r['observations'][3]['under_illustrative_margin_0p001']


def test_exact_weyl_gate_hypothetical_only():
    assert perturbation_gate_exact(Q(1,1000),Q(1,10000),Q(1))['positive_definiteness_preserved_conditional']
    assert not perturbation_gate_exact(Q(1,1000),Q(1,100),Q(1))['positive_definiteness_preserved_conditional']
    with pytest.raises(ValueError):
        perturbation_gate_exact(Q(1,1000),Q(-1),Q(1))
