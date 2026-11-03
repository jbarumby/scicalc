import math
import numpy as np
import pytest
import sympy as sp
from es3_spin_zeta import (spin, spin_generator, local_spin_audit, kron_sum_spectrum,
                           harmonic_finite_matrix, oscillator_exact, zeta_zeros,
                           li_partial, fit_affine_overclass, analyze)


def test_all_eight_spin_states_exact():
    assert len(local_spin_audit()) == 8
    for nz in (-1, 1):
        for nx in (-1, 1):
            for ny in (-1, 1):
                s = spin(nz,nx,ny)
                assert s*s == s
                assert sp.trace(s) == 1
                assert set((s+s.H).eigenvals()) == {1-sp.sqrt(2),1+sp.sqrt(2)}
                assert set((s.H*s).eigenvals()) == {0,2}


def test_tensor_spin_spectrum_exact():
    for m in range(1,5):
        result=kron_sum_spectrum(m)
        assert sum(mult for eigen,mult in result['eigenvalue_multiplicities'])==2**m
        assert result['zero_eigenvalue'] == (m%2==0)


def test_nonhermitian_example_and_automatic_positive():
    s=spin(1,1,1)
    assert s!=s.H
    assert (s.H*s).det()==0
    assert set((s.H*s).eigenvals())=={0,2}


def test_harmonic_matrix_symmetric_and_convergence():
    alpha,beta,gamma=0.4,3.,0.5
    h=harmonic_finite_matrix(48,alpha,beta,gamma)
    assert np.allclose(h,h.T)
    spec=np.linalg.eigvalsh(h)[:6]
    exact=[oscillator_exact(n,alpha,beta,gamma) for n in range(6)]
    assert np.max(np.abs(spec-np.array(exact))) < 1e-8


def test_no_go_constant_gap():
    a,b,g=0.4,3.0,0.5
    energies=[oscillator_exact(n,a,b,g) for n in range(6)]
    differences=np.diff(energies)
    assert np.allclose(differences,math.sqrt(b*b-4*g*g))
    with pytest.raises(ValueError): oscillator_exact(0,0.,1.,1.)


def test_li_partial_double_conjugate_and_nonnegativity():
    zeros=zeta_zeros(5)
    p=li_partial(zeros,indices=(1,2))
    assert float(p['1']) > 0
    assert float(p['2']) > 0
    assert float(p['1']) < 1.0


def test_held_out_test_distinct_from_training():
    zeros=zeta_zeros(10)
    fit=fit_affine_overclass(zeros,training=6,holdout=4)
    assert fit['gap']>0
    assert fit['holdout_rmse']>0
    assert len(fit['holdout_errors'])==4


def test_invalid_inputs():
    with pytest.raises(AssertionError): spin(0,1,1)
    with pytest.raises(AssertionError): kron_sum_spectrum(5)
    with pytest.raises(ValueError): analyze(J=7,training=5,holdout=5)
