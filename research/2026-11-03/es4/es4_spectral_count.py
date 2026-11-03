"""NUMERIA ES4: explicit non-quadratic diagonal operator with R-v-M *smooth counting*.

This is a deliberately constructed control operator, NOT a Hilbert-Polya
operator: it uses no primes or zeta zeros and makes no RH claim.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path

import mpmath as mp

OFFSET = mp.mpf(7) / 8


def smooth_count(t):
    """N_0(t)=(t/(2*pi))*(log(t/(2*pi))-1)+7/8; not N_zeta(t)."""
    t = mp.mpf(t)
    if t <= 0:
        raise ValueError('t must be positive')
    u = t/(2*mp.pi)
    return u*(mp.log(u)-1)+OFFSET


def energy(n):
    """Positive inverse branch of N_0(E_n)=n, n>=1.

    E_n=2*pi*(n-7/8)/W_0((n-7/8)/e). No zeros are inputs.
    """
    if not isinstance(n,int) or n<1:
        raise ValueError('n must be a positive integer')
    v=mp.mpf(n)-OFFSET
    return 2*mp.pi*v/mp.lambertw(v/mp.e)


def energy_derivative(n):
    """Derivative of continuous interpolation E(x), x>=1."""
    x=mp.mpf(n)
    if x<1:
        raise ValueError('n must be >= 1')
    return 2*mp.pi/(1+mp.lambertw((x-OFFSET)/mp.e))


def model_count(t):
    """Exact floor rule above the first level; 0 below it.

    Formula assumes t>=2*pi, so only used in that monotonic region.
    """
    t=mp.mpf(t)
    if t<=0:
        return 0
    if t<energy(1):
        return 0
    return max(0,int(mp.floor(smooth_count(t))))



def shifted_model_count(t, shift=1):
    """Count of H+shift*I; same leading R-v-M asymptotic, different spectrum."""
    shift=mp.mpf(shift)
    if shift<0:
        raise ValueError('shift must be nonnegative for this diagnostic')
    return model_count(mp.mpf(t)-shift)


def diagonal_truncation(m):
    if not isinstance(m,int) or m<1:
        raise ValueError('m must be positive')
    return [energy(n) for n in range(1,m+1)]


def gaussian_trace(levels,scale):
    """Finite positive sum, not a trace formula for zeta or Weil."""
    scale=mp.mpf(scale)
    if scale<=0:
        raise ValueError('scale must be positive')
    return mp.fsum([mp.exp(-(mp.mpf(e)/scale)**2) for e in levels])


def li_like_partial(levels,k):
    """Positive tautology for arbitrary real positive levels, NOT full Li coeffs."""
    if not isinstance(k,int) or k<1:
        raise ValueError('k must be a positive integer')
    return 2*mp.fsum([1-mp.cos(2*k*mp.atan(1/(2*mp.mpf(e)))) for e in levels])


def benchmark(j=25,dps=60):
    if not isinstance(j,int) or j<5:
        raise ValueError('j must be at least 5')
    if not isinstance(dps,int) or dps<30:
        raise ValueError('dps must be >= 30')
    with mp.workdps(dps):
        es=diagonal_truncation(j)
        # Independent reference only: zeros are not used to construct H.
        zeros=[mp.im(mp.zetazero(k)) for k in range(1,j+1)]
        gaps=[es[k]-es[k-1] for k in range(1,j)]
        rmse=mp.sqrt(mp.fsum([(a-b)**2 for a,b in zip(es,zeros)])/j)
        sample=[]
        for k in [1,2,3,5,10,20,25,50,100,1000,10000]:
            if k<=0:continue
            e=energy(k)
            sample.append({'n':k,'energy':mp.nstr(e,45),
                           'smooth_count_at_energy':mp.nstr(smooth_count(e),45),
                           'continuous_gap_derivative':mp.nstr(energy_derivative(k),45)})
        counts=[]
        for t in [20,50,100,500,1000,10000]:
            counts.append({'T':t,'N_model':model_count(t),
                           'N_smooth':mp.nstr(smooth_count(t),45)})
        shifted_counts=[{'T':t,'N_base':model_count(t),'N_shift_by_one':shifted_model_count(t,1),'difference':model_count(t)-shifted_model_count(t,1)} for t in [20,50,100,500,1000,10000]]
        trace=[]
        for scale in [20,40,80]:
            t_model=gaussian_trace(es,scale)
            t_zeros=gaussian_trace(zeros,scale)
            trace.append({'scale':scale,'finite_model':mp.nstr(t_model,42),
                          'finite_known_zeros':mp.nstr(t_zeros,42),
                          'finite_difference':mp.nstr(t_model-t_zeros,42),
                          'warning':'both sums truncated at J; NOT global trace identity'})
        return {
          'schema':'NUMERIA.ES4.SPECTRAL_COUNT_CONTROL.1',
          'status':'explicit self-adjoint diagonal control; counting law only; not RH proof',
          'math':{
            'hilbert_space':'ell^2(N_{>=1})',
            'operator':'H e_n=E_n e_n, E_n=2*pi*(n-7/8)/W_0((n-7/8)/e)',
            'domain':'D(H)={c in ell^2: sum_{n>=1} E_n^2 |c_n|^2 < infinity}',
            'spectral_count':'N_H(T)=max(0,floor(N_0(T))) for T >= E_1; 0 below E_1',
            'smooth_count':'N_0(T)=T/(2*pi)*(log(T/(2*pi))-1)+7/8',
            'counting_asymptotic':'N_H(T)=T/(2*pi)*log(T/(2*pi))-T/(2*pi)+O(1)',
            'limitations':['No prime terms or arithmetic trace formula',
                           'No identification of eigenvalues with zeta zeros',
                           'Finite comparisons to known critical-line zeros not interval-certified',
                           'Model matches only Weyl/smooth counting, not fluctuations',
                           'Li-like sums positive for ANY real energies by construction']},
          'j':j,'working_precision_dps':dps,'environment':{'python':platform.python_version(),'mpmath':mp.__version__},
          'model_levels_first_J':[mp.nstr(x,45) for x in es],
          'zeta_known_critical_line_ordinates_first_J':[mp.nstr(x,45) for x in zeros],
          'zeta_reference_status':'mpmath zetazero; exploratory numerical values only, no certification of all zeros',
          'max_abs_level_error':mp.nstr(max(abs(a-b) for a,b in zip(es,zeros)),45),
          'rmse_level_error':mp.nstr(rmse,45),
          'min_gap_model_first_J':mp.nstr(min(gaps),45),
          'max_gap_model_first_J':mp.nstr(max(gaps),45),
          'samples':sample,'counts':counts,'shifted_operator_adversarial_counts':shifted_counts,
          'shifted_operator_warning':'H+I has a different spectrum and the same R-v-M leading count, hence asymptotic counting does not identify the zeta zeros',
          'finite_gaussian_traces':trace,
          'li_like_partial':{str(k):mp.nstr(li_like_partial(es,k),45) for k in [1,2,5,10]},
          'li_warning':'These sums are nonnegative even when E_n differ from zeta zeros; not Li coefficients',
          'next_gate':'Arithmetic prime/archimedean trace identity with justified domains, normalizations and global tails.'
        }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='RESULTADOS_ES4.json')
    ap.add_argument('--zeros',type=int,default=25)
    ap.add_argument('--dps',type=int,default=60)
    a=ap.parse_args()
    data=benchmark(a.zeros,a.dps)
    p=Path(a.out)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('ES4 result:',p,'J=',a.zeros,'rmse=',data['rmse_level_error'],'no RH proof')

if __name__=='__main__':
    main()
