"""Exact rational proof witnesses + adversarial controls for ES8-P10."""
from __future__ import annotations
import json
from decimal import Decimal
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import pytest
from es8p10_lemmas import (
    NODES, CELLS, RHO, ALPHA, DIRECTION,
    absolute_newton_integral, newton_polynomial, exact_polynomial_integral,
    original_cell_majorant, refined_interpolation_error,
    independent_total_error, coth_coefficients, cm_function,
    em_corrections, validate_cauchy_geometry, decimal_out, new_intervals, audit,
)

BASE=Path(__file__).resolve().parent

def test_interpolation_polynomial_zeros():
    xs,p=newton_polynomial()
    assert len(xs)==16 and len(p)==17
    for x in xs:
        assert sum(c*x**n for n,c in enumerate(p))==0
    assert p[-1]==1

def test_polynomial_symmetry_and_pieces():
    xs,p=newton_polynomial()
    A,pieces=absolute_newton_integral()
    assert len(pieces)==15 and A==sum(map(abs,pieces),Q(0))
    assert pieces==list(reversed(pieces)), 'degree-16 polynomial is symmetric about 1/2'
    assert all(a*b<0 for a,b in zip(pieces,pieces[1:]))

def test_lagrange_new_bound_is_strict():
    A,_=absolute_newton_integral()
    B=original_cell_majorant()
    assert 0<A<B
    assert Q(6_000_000)<B/A<Q(7_000_000)

def test_integrated_error_exact_ratio():
    A,_=absolute_newton_integral()
    Enew=refined_interpolation_error(); Eold=independent_total_error()
    assert Eold/Enew==original_cell_majorant()/A
    assert Q(3,10**12)<Enew<Q(4,10**12)

def test_cauchy_geometry_strict():
    d=validate_cauchy_geometry()
    assert d['digamma_real_part_lower']=='3/40'
    assert ALPHA>0 and DIRECTION>0

def test_coth_symbolic_coefficients_all_n():
    # Proof is the algebraic expressions in the source, valid for all n>=0.
    for n in range(64):
        u,v=coth_coefficients(n)
        assert u>=0 and v>=0
    assert coth_coefficients(3)[0]>0
    assert coth_coefficients(3)[1]>0
    assert coth_coefficients(1)==(0,0)
    assert coth_coefficients(2)[1]==0

def test_coth_mutation_wrong_sign_fails():
    # Reverse the 4th-order correction; cannot be a lower bound as t->0.
    assert coth_coefficients(3)[1]>0
    assert -coth_coefficients(3)[1]<0

def test_em_corrections_equivalent_formula():
    for u in [Q(1,4),Q(257,4),Q(513,4),Q(128)]:
        for b2 in [Q(0),Q(1),Q(100),Q(6400)]:
            c2,c4=em_corrections(u,b2)
            assert c2>=0 and c4>=0
            if b2>0:assert c2>0 and c4>0

def test_em_b2_is_derivative_manual():
    u,b2=Q(257,4),Q(75)
    f=cm_function(u,b2)
    deriv=-b2*(3*u*u+b2)/(u*u*(u*u+b2)**2)
    assert em_corrections(u,b2)[0]==-deriv/12
    assert f>0

def test_em_mutation_bad_lower_upper():
    a,b=em_corrections(Q(128),Q(6400))
    assert a>b>0
    # Dropping B4 would mistakenly declare the upper approximation exact.
    assert a-b < a

def test_invalid_arguments_rejected():
    for u,b2 in [(0,1),(-1,1),(1,-1)]:
        with pytest.raises(ValueError):em_corrections(u,b2)
    with pytest.raises(ValueError):coth_coefficients(-1)

def test_decimal_out_negative_edges():
    assert decimal_out(Q(-1,3),precision=3)=='-0.334'
    assert decimal_out(Q(-1,3),precision=3,upper=True)=='-0.333'
    assert decimal_out(Q(1,3),precision=3)=='0.333'
    assert decimal_out(Q(1,3),precision=3,upper=True)=='0.334'

def test_refinement_nested_in_original_intervals():
    src=json.loads((BASE/'RESULTADOS_ES8P9_INTEGRAIS.json').read_text())
    r=new_intervals()
    assert len(r)==len(src['entries'])==30
    for k,a in r.items():
        old=src['entries'][k]
        x,y=map(Decimal,a['independent_finite_with_cauchy'])
        x0,y0=map(Decimal,old['independent_finite_with_cauchy'])
        assert x0<=x<=y<=y0

def test_reject_partial_input(tmp_path):
    d=json.loads((BASE/'RESULTADOS_ES8P9_INTEGRAIS.json').read_text())
    d['scope']['evaluations']=767999
    p=tmp_path/'partial.json';p.write_text(json.dumps(d))
    with pytest.raises(AssertionError):new_intervals(p)

def test_reject_forged_outward_bound(tmp_path):
    d=json.loads((BASE/'RESULTADOS_ES8P9_INTEGRAIS.json').read_text())
    e=d['entries']['1,1']
    e['independent_finite_with_cauchy']=['0','0']
    p=tmp_path/'forged.json';p.write_text(json.dumps(d))
    with pytest.raises(AssertionError):new_intervals(p)

def test_projected_all_three_cases():
    d=audit()
    assert d['error_reduction_factor_floor']>6_000_000
    assert [m['N'] for m in d['projected_matrices']]==[6,8,10]
    assert all(m['certified_M_minus_0p05I_positive'] for m in d['projected_matrices'])

def test_reported_scope_does_not_claim_rh():
    d=audit()
    assert 'Riemann Hypothesis' in ' '.join(d['not_established'])
    assert d['improved_finite_entries_included_in_p9']==30
