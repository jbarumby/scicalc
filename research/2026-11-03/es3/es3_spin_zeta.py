"""NUMERIA ES3: falsifiable spin/zeta toy-model audit (not RH proof).

CLI outputs reproducible JSON for symbolic spin identities, harmonic-oscillator
truncations, an intentionally over-flexible affine-spectrum benchmark, and
partial Li sums from known critical-line zeros. No interval certification.
"""
import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

I2 = sp.eye(2)
SX = sp.Matrix([[0, 1], [1, 0]])
SY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
SZ = sp.Matrix([[1, 0], [0, -1]])


def spin(nz=1, nx=1, ny=1):
    """Sanctuary idempotent; ± orientations only."""
    assert all(v in (-1, 1) for v in (nz, nx, ny))
    return (I2 + nz*SZ + nx*SX + sp.I*ny*SY)/2


def spin_generator(nz=1, nx=1, ny=1):
    return 2*spin(nz, nx, ny)-I2


def local_spin_audit():
    out = []
    for nz in (-1, 1):
        for nx in (-1, 1):
            for ny in (-1, 1):
                s = spin(nz, nx, ny)
                a = 2*s-I2
                assert s*s == s
                assert sp.trace(s) == 1 and s.det() == 0
                assert a*a == I2
                assert sorted(s.eigenvals(), key=float) == [0, 1]
                assert sorted((s+s.H).eigenvals(), key=float) == [1-sp.sqrt(2), 1+sp.sqrt(2)]
                assert sorted((s.H*s).eigenvals(), key=float) == [0, 2]
                out.append({'nz': nz, 'nx': nx, 'ny': ny, 's_idempotent': True,
                            'spec_s': ['0', '1'], 'spec_s_plus_adjoint': ['1-sqrt(2)', '1+sqrt(2)'],
                            'spec_s_adjoint_s': ['0', '2'], 'spec_a': ['-1', '1']})
    return out


def kron_sum_spectrum(m):
    """Exact finite-site tensor sum, not a construction on an adele space."""
    assert 1 <= m <= 4
    a = spin_generator()
    total = sp.zeros(2**m)
    for site in range(m):
        term = sp.ones(1, 1)
        for k in range(m):
            term = sp.kronecker_product(term, a if k == site else I2)
        total += term
    spec = sorted([(int(v), int(n)) for v,n in total.eigenvals().items()])
    predicted = sorted([(-m+2*k, math.comb(m,k)) for k in range(m+1)])
    assert spec == predicted
    return {'sites':m,'dimension':2**m,'eigenvalue_multiplicities':spec,
            'zero_eigenvalue':any(v==0 for v,n in spec)}


def harmonic_finite_matrix(M, alpha, beta, gamma):
    """H=alpha(a+a*) + beta a*a + gamma(a*a*+aa); Hermite basis 0..M-1."""
    assert M>=3
    h=np.diag(beta*np.arange(M,dtype=float))
    for j in range(M-1):
        off=alpha*math.sqrt(j+1)
        h[j,j+1]=h[j+1,j]=off
    for j in range(M-2):
        off=gamma*math.sqrt((j+1)*(j+2))
        h[j,j+2]=h[j+2,j]=off
    return h


def oscillator_exact(n, alpha,beta,gamma):
    """Stable confining sector beta>2 abs(gamma)."""
    if not beta>2*abs(gamma):
        raise ValueError('Confining spectrum requires beta>2*abs(gamma)')
    omega=math.sqrt(beta*beta-4*gamma*gamma)
    e0=(omega-beta)/2-alpha*alpha/(beta+2*gamma)
    return e0+omega*n


def zeta_zeros(J,dps=65):
    with mp.workdps(dps):
        return [mp.im(mp.zetazero(k)) for k in range(1,J+1)]


def li_partial(zeros,indices=(1,2,5,10)):
    """Sums BOTH conjugate zeros: 2 sum_{gamma>0}(1-cos(n theta)).
    They are NOT the complete Li coefficients and NOT evidence of RH."""
    a={}
    for n in indices:
        terms=[]
        for t in zeros:
            theta=mp.pi-2*mp.atan(2*t)
            terms.append(2*(1-mp.cos(n*theta)))
        a[str(n)]=mp.nstr(mp.fsum(terms),45)
    return a


def fit_affine_overclass(zeros,training=8,holdout=6):
    """A+omega*j, with arbitrary A: a superset of stable quadratic spectra.
    Fit gamma_j (not theta_j), and score disjoint holdout data."""
    if training<3 or training+holdout>len(zeros):
        raise ValueError('Need at least 3 training zeros and enough holdout zeros')
    y=np.asarray([float(x) for x in zeros],dtype=float)
    x=np.arange(len(zeros),dtype=float)
    omega=float(np.sum((x[:training]-np.mean(x[:training]))*(y[:training]-np.mean(y[:training])))/np.sum((x[:training]-np.mean(x[:training]))**2))
    offset=float(y[:training].mean()-omega*x[:training].mean())
    pred=offset+omega*x
    return {'training':training,'holdout':holdout,'offset':offset,'gap':omega,
            'train_rmse':float(np.sqrt(np.mean((y[:training]-pred[:training])**2))),
            'holdout_rmse':float(np.sqrt(np.mean((y[training:training+holdout]-pred[training:training+holdout])**2))),
            'holdout_errors':[float(e) for e in (y[training:training+holdout]-pred[training:training+holdout])],
            'model':'unconstrained affine spectrum; more flexible than given stable quadratic Hamiltonian'}


def rv_m_main(t):
    return t/(2*math.pi)*(math.log(t/(2*math.pi))-1)+7/8


def analyze(J=20,training=8,holdout=6):
    if J<training+holdout: raise ValueError('J must >= training+holdout')
    spin_checks=local_spin_audit()
    site_models=[kron_sum_spectrum(m) for m in (1,2,3,4)]
    alpha,beta,gamma=0.4,3.0,0.5
    oscs=[]
    for M in (12,24,48):
        calculated=np.linalg.eigvalsh(harmonic_finite_matrix(M,alpha,beta,gamma))[:6]
        predicted=np.array([oscillator_exact(n,alpha,beta,gamma) for n in range(6)])
        oscs.append({'truncation_size':M,'eigenvalues_first_six':calculated.tolist(),
                     'max_abs_difference_from_exact':float(np.max(np.abs(calculated-predicted)))})
    with mp.workdps(65):
        zeros=zeta_zeros(J)
        partial=li_partial(zeros)
        lambda1=1+mp.euler/2-mp.log(4*mp.pi)/2
        li_1_gap=lambda1-mp.mpf(partial['1'])
        zeros_out=[mp.nstr(x,55) for x in zeros]
        li_one=mp.nstr(lambda1,55)
        li_gap=mp.nstr(li_1_gap,40)
    affine=fit_affine_overclass(zeros,training,holdout)
    # Mean spacing of zeros goes to zero; stable harmonic gap is constant.
    # R-v-M formula represents only the main asymptotic count, not certified counts.
    counts=[{'T':T,'RvM_main_approx':rv_m_main(T),
             'harmonic_count_with_exact_parameters':max(0,math.floor((T-oscillator_exact(0,alpha,beta,gamma))/math.sqrt(beta**2-4*gamma**2))+1)}
            for T in (100,300,1000,3000,10000)]
    return {'schema':'NUMERIA.ES3.SPIN_ZETA_AUDIT.1',
            'scientific_status':'toy-model no-go and numerical exploration; not RH proof; not interval certified',
            'sources':{'main_input':'Qwen_latex_20261103_m86z47xnk.txt','reference':'Sanctuary arXiv:0908.3219'},
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'mpmath':mp.__version__,'sympy':sp.__version__},
            'spin':{'eight_orientations':spin_checks,'finite_site_tensor_sum':site_models,
                    'conclusion':'finite tensor sum eigenvalues are -m,-m+2,...,m; NOT ±sqrt(2**m)'},
            'oscillator':{'parameters':{'alpha':alpha,'beta':beta,'gamma':gamma},
                          'omega':math.sqrt(beta**2-4*gamma**2),
                          'exact_E0':oscillator_exact(0,alpha,beta,gamma),'convergence':oscs,
                          'no_go':'stable quadratic oscillator has affine eigenvalue sequence; R-v-M count grows as T log T'},
            'zeta':{'positive_ordinates_first_J':zeros_out,'precision_digits':65,
                    'nature':'mpmath numerically computed known critical-line zeros; no interval or completeness certificate'},
            'affine_fit':affine,
            'li':{'partial_paired_critical_line_zeros':partial,'exact_lambda_1_formula':li_one,
                  'lambda_1_minus_partial_lambda_1':li_gap,
                  'warning':'finite partial Li sums automatically nonnegative, never complete Li coefficients'},
            'spectral_counts':counts,
            'limitations':['No adelic Hilbert space/operator constructed',
                           'No analytic identification of zeta zeros as an operator spectrum',
                           'Finite-spin spectral degeneracy is not the conjectured 2**(3*m)',
                           'Positive S* S is tautological and does not imply RH',
                           'Toy-model stable quadratic Hamiltonian cannot realize zeta zero counting law',
                           'No Krawczyk/Arb certification of any zeta zero or matrix of Weil',
                           'No inference about off-critical zeros']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='RESULTADOS_ES3.json')
    parser.add_argument('--zeros',type=int,default=20)
    parser.add_argument('--train',type=int,default=8)
    parser.add_argument('--holdout',type=int,default=6)
    args=parser.parse_args()
    data=analyze(args.zeros,args.train,args.holdout)
    p=Path(args.out); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(f'ES3 written {p} | spin=8 cases | sites=4 | zeros={args.zeros} | status=NO_RH_PROOF')

if __name__=='__main__':main()
