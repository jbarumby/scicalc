#!/usr/bin/env python3
"""ES8-P6: independent rational audit of Phase-4 Cauchy interpolation and E7 tail.

Does NOT rerun all nodal evaluations and does NOT prove Weil/RH.
Analytic lemmas and their assumptions are stated in RELATORIO_ES8P6.md.
Only uses Python stdlib, not original source functions.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
from fractions import Fraction as Q
import hashlib
import io
import json
from math import factorial
from pathlib import Path
import zipfile
from functools import lru_cache

L=Q(9,20); H=Q(1,10); RHO=Q(3,10); NODES=16; CELLS=1600; R=160
D=RHO+H/2
ALPHA=Q(1,4)-D/2
DIRECTION=1-2*L*D
EXPECTED_PAIRS={(i,j) for i in range(1,11) for j in range(i,11) if (i-j)%2==0}


def lagrange_weights(m=NODES):
    """Integrate Lagrange basis polynomials DIRECTLY; no Vandermonde inversion."""
    x=[Q(j,m-1) for j in range(m)]
    out=[]
    for j in range(m):
        coefficients=[Q(1)]
        divisor=Q(1)
        for k in range(m):
            if k==j: continue
            new=[Q(0)]*(len(coefficients)+1)
            for power,coef in enumerate(coefficients):
                new[power]-=coef*x[k]
                new[power+1]+=coef
            coefficients=new
            divisor*=x[j]-x[k]
        out.append(sum((coef/Q(power+1) for power,coef in enumerate(coefficients)), Q())/divisor)
    return out


def interpolation_product(m=NODES):
    val=Q(1)
    for j in range(m):
        x=Q(j,m-1)
        val*=max(x,1-x)
    return val


def independent_cauchy_bound(center):
    """Q-bound for |Re(Fi Fbarj) Omega/pi| on a common complex region.

    Justification:
      Fi ≤ 2 sqrt(L) exp(LD) ≤ 2 sqrt(L)/(1-LD)
      |Fi Fj| ≤4L/(1-2LD)
      w=1/4 ± iz/2, Re(w)≥alpha>0
      |psi(w)| ≤1+1/alpha+2|w| from DLMF 5.7.6 and zeta(2)<2
      |Omega(z)|≤7/2+1/alpha+|center|+D
      ≤ 3+(1+1/alpha)*(1+(|center|+D)/2) since alpha≤1
      pi≥3.
    """
    center=Q(center)
    if not (0<L and 0<D<Q(1,2) and 0<ALPHA<=1 and DIRECTION>0):
        raise AssertionError('Cauchy strip conditions do not hold')
    b=3+(1+1/ALPHA)*(1+(abs(center)+D)/2)
    stronger=Q(7,2)+1/ALPHA+abs(center)+D
    assert b>=stronger
    return 4*L*b/(3*DIRECTION)


def independent_total_error():
    """O(1) closed sum: M(center) is affine for these nonnegative cell centers."""
    pref=4*L/(3*DIRECTION)
    c0=pref*(3+(1+1/ALPHA)*(1+D/2))
    c1=pref*(1+1/ALPHA)/2
    sum_midpoints=H*Q(CELLS*CELLS,2)
    return H*interpolation_product()*(H/RHO)**NODES*(CELLS*c0+c1*sum_midpoints)


def exact_error_by_cells():
    return H*interpolation_product()*(H/RHO)**NODES*sum(
        (independent_cauchy_bound((Q(k)+Q(1,2))*H) for k in range(CELLS)), Q())


def atan_inv_interval(q, n):
    """Alternating arctan(1/q) enclosure, independent of legacy atan function."""
    assert q>=2 and n>=1
    x=Q(1,q)
    s=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)),Q())
    next_term=x**(2*n+1)/Q(2*n+1)
    return (s,s+next_term) if n%2==0 else (s-next_term,s)


@lru_cache(maxsize=None)
def pi_bounds():
    a=atan_inv_interval(5,40)
    b=atan_inv_interval(239,25)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]


@lru_cache(maxsize=None)
def log_positive_bounds(x,terms=80):
    """atanh series, |(x-1)/(x+1)|<1."""
    x=Q(x)
    if x<=0:raise ValueError
    q=(x-1)/(x+1)
    s=sum((2*q**(2*k+1)/Q(2*k+1) for k in range(terms)),Q())
    rem=2*abs(q)**(2*terms+1)/(Q(2*terms+1)*(1-q*q))
    return s-rem,s+rem


@lru_cache(maxsize=None)
def ln480_bounds():
    """ln 480 = 9 ln 2 - ln(16/15), accelerated via |atanh|<=1/3."""
    a,b=log_positive_bounds(2,80)
    c,d=log_positive_bounds(Q(16,15),80)
    return 9*a-d,9*b-c


def baseline_tail_upper(i,j):
    """Independent conservative Prop.-5 bound for R=160, with interval constants."""
    assert 1<=i<=10 and 1<=j<=10
    if (i-j)%2: return Q(0)
    pl,ph=pi_bounds()
    assert 3<pl<=ph<4
    assert R*L>=10*ph             # R>=2*w_max=N*pi/L
    ln480hi=ln480_bounds()[1]
    return 4*i*j*ph/L**3*(8+ln480hi+Q(1,3))/(3*R**3)



@lru_cache(maxsize=None)
def sharp_omega_constant_witness():
    """All rational comparisons used to prove |Omega(r)|<=ln(r)+157/50 for r>=1.

    Duplication/reflection ψ identities and positivity/monotonicity of the
    S(r) series are external classical facts. No floating-point inputs.
    """
    n=256
    h_n=sum((Q(1,k) for k in range(1,n+1)),Q())
    ln2lo,ln2hi=log_positive_bounds(2,80)
    gamma_lo=h_n-8*ln2hi-Q(1,2*n)+Q(1,12*n*n)-Q(1,120*n**4)
    gamma_hi=h_n-8*ln2lo-Q(1,2*n)+Q(1,12*n*n)
    pi_lo,pi_hi=pi_bounds()
    assert pi_lo>Q(157,50) and pi_hi<Q(22,7)
    lnpi_lo=log_positive_bounds(pi_lo,100)[0]
    lnpi_hi=log_positive_bounds(pi_hi,100)[1]
    C_lo=gamma_lo+3*ln2lo+lnpi_lo
    C_hi=gamma_hi+3*ln2hi+lnpi_hi
    assert C_lo>Q(38013,10000) and C_hi<Q(38014,10000)
    ln5lo,ln5hi=log_positive_bounds(5,120)
    assert ln5lo/2>Q(8047,10000) and ln5hi/2<Q(8048,10000)
    coshpi_lo=sum((Q(157,50)**(2*k)/Q(factorial(2*k)) for k in range(9)),Q())
    correction_hi=Q(22,7)/(2*coshpi_lo)
    assert correction_hi<Q(14,100)
    omega_lower_constant=-C_hi+ln5lo/2-correction_hi
    omega_upper_offset=-C_lo+2+ln5hi/2
    assert omega_lower_constant>-Q(157,50)
    assert omega_upper_offset<-Q(99,100)
    # Thus Omega(r)>= -3.14 and Omega(r)<=ln(r)-0.99, for r>=1.
    return {'gamma_bounds':[gamma_lo,gamma_hi],
      'C_bounds':[C_lo,C_hi],'ln5_bounds':[ln5lo,ln5hi],
      'coshpi_lower':coshpi_lo,'correction_upper':correction_hi,
      'omega_lower':omega_lower_constant,'omega_upper_offset':omega_upper_offset}


def sharp_tail_upper(i,j):
    """Independent rational Proposition-6-type tail enclosure for R=160."""
    assert 1<=i<=10 and 1<=j<=10
    if (i-j)%2:return Q(0)
    sharp_omega_constant_witness()
    pi_hi=pi_bounds()[1]
    wi_hi=Q(i)*pi_hi/(2*L)
    wj_hi=Q(j)*pi_hi/(2*L)
    assert wi_hi*2<=R and wj_hi*2<=R
    deni=1-(wi_hi/R)**2
    denj=1-(wj_hi/R)**2
    assert deni>0 and denj>0
    ln2hi=log_positive_bounds(2,80)[1]
    ln5hi=log_positive_bounds(5,120)[1]
    ln160hi=5*ln2hi+ln5hi
    # 4 wi wj/(pi L) = ij*pi/L^3, with positive rational bounds.
    return Q(i*j)*pi_hi/L**3/(deni*denj)*(ln160hi+Q(157,50)+Q(1,3))/(3*R**3)


def e7_original_tail_check(e7_zip:Path):
    """Two historical E7 entries, validated by original zip digest.

    The expected E7 bound and our bound should agree at high precision.
    No claim that this independently audits *all* original nodal sums.
    """
    EXPECTED='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
    assert hashlib.sha256(e7_zip.read_bytes()).hexdigest()==EXPECTED
    with zipfile.ZipFile(e7_zip) as z:
        assert z.testzip() is None
        results=json.loads(z.read('RESULTADO_INTEGRAIS_E7.json'))['results']
        out=[]
        for i,entry in [(1,'G11'),(2,'G22')]:
            r=next(x for x in results if x['entry']==entry)
            orig=Q(r['tail_proposition6'])
            new=sharp_tail_upper(i,i)
            assert abs(orig-new)<Q(1,10**14),(entry,'tail mismatch',float(orig-new))
            out.append({'entry':entry,'historical':float(orig),'independent':float(new),'difference_less_than_1e-14':True})
        return out

def stored_hash_check(z:zipfile.ZipFile,path:str,expected_hash:str):
    blob=z.read(path)
    assert hashlib.sha256(blob).hexdigest()==expected_hash,('hash mismatch',path)


def historical_phase4_audit(phase4_zip:Path):
    """Does not execute legacy code; checks raw input files and exact arithmetic."""
    error=independent_total_error()
    assert error==exact_error_by_cells()
    pi_lo,pi_hi=pi_bounds()
    assert 3<pi_lo<pi_hi<4
    with zipfile.ZipFile(phase4_zip) as z:
        assert z.testzip() is None
        summary=json.loads(z.read('RESUMO_FASE4_1600CELULAS_D40.json'))
        assert summary['status']=='CONCLUIDO' and summary['num_pairs']==30
        assert len(summary['records'])==30
        for name,digest in summary['source_hashes'].items():
            stored_hash_check(z,name,digest)
        pairs=set()
        nodal_checks=0
        tail_checks=0
        range_checks=0
        for rec in summary['records']:
            i,j=rec['i'],rec['j'];pair=(i,j)
            assert pair in EXPECTED_PAIRS and pair not in pairs
            pairs.add(pair)
            name=f'FASE4_G{i}{j}_1600CELULAS_D40.json'
            stored_hash_check(z,name,rec['sha256'])
            d=json.loads(z.read(name))
            assert d['source_hashes']==summary['source_hashes']
            assert (d['i'],d['j'],d['digits'],d['series_terms'],d['cells'],d['nodes_per_cell'])==(i,j,40,64,1600,16)
            assert d['L']=='9/20' and d['end']=='160' and d['cauchy_radius']=='3/10'
            assert Q(d['interpolation_error_exact'])==error
            lo=Q(Decimal(d['integral_lower'])); hi=Q(Decimal(d['integral_upper']))
            assert lo<hi
            a=Q(Decimal(d['nodal_integral_lower'])); b=Q(Decimal(d['nodal_integral_upper']))
            assert a<=b
            # These outward 24-decimal displays each add up to one ulp.
            tol=Q(2,10**24)
            assert abs(lo-(a-error))<=tol
            assert abs(hi-(b+error))<=tol
            assert d['status']=='rigorous_finite_segment_given_phase3_analytic_majorant'
            assert d['cell_width']=='1/10'
            assert (d['integral_lower'],d['integral_upper'])==(rec['lower'],rec['upper'])
            assert baseline_tail_upper(i,j)>sharp_tail_upper(i,j)>0
            nodal_checks+=1;tail_checks+=1;range_checks+=1
        assert pairs==EXPECTED_PAIRS
        return {
            'status':'ANALYTIC_CAUCHY_INTERPOLATION_AND_TAIL_REDERIVED;_NODES_NOT_RECOMPUTED',
            'input_sha256':hashlib.sha256(phase4_zip.read_bytes()).hexdigest(),
            'pairs':len(pairs),'source_hashes':len(summary['source_hashes']),
            'stored_numeric_intervals_checked':range_checks,
            'exact_interpolation_errors_checked':nodal_checks,
            'baseline_tail_entries_checked':tail_checks,
            'sharp_tail_entries_checked':tail_checks,
            'sharp_omega_rational_constant_verified':True,
            'parameters':{'L':str(L),'h':str(H),'rho':str(RHO),'d':str(D),
                          'alpha':str(ALPHA),'one_minus_2Ld':str(DIRECTION),
                          'R':R,'nodes_per_cell':NODES,'cells':CELLS},
            'interpolation_error_exact':str(error),
            'interpolation_error_decimal':f'{float(error):.18g}',
            'lagrange_weight_count':len(lagrange_weights()),
            'lagrange_exact_moment_checks':NODES,
            'pi_in_3_4': True,
            'independent_analytic_proofs':[
               'holomorphic digamma symmetrization in |Im z|<1/2',
               'absolute F_k bound from compact [-L,L] support and exp(Ld)',
               'DLMF 5.7.6 digamma series yielding |psi(w)|<=1+1/alpha+2|w|',
               'Cauchy error for real degree-15 polynomial interpolation',
               'exact rational closed sum of 1600 cell error majorants',
               'conservative E7 tail via classical psi reflection/duplication and monotone series',
               'sharp E7 Proposition 6 tail using rational gamma, ln, cosh estimates'],
            'not_proved':[
               'independent re-evaluation of all 25600 Phase4 nodal function values per entry',
               'a fully formal proof in Lean of special-function identities',
               'ES8-P3 Euler-Maclaurin implementation, apart from the E7 sharp tail rederived here',
               'equivalence to universal Weil positivity or Riemann Hypothesis']
        }


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--phase4',type=Path,required=True)
    a.add_argument('--out',type=Path,required=True)
    a.add_argument('--e7',type=Path,default=None)
    args=a.parse_args()
    output=historical_phase4_audit(args.phase4)
    if args.e7 is not None:output['historical_e7_tail_comparisons']=e7_original_tail_check(args.e7)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('ES8-P6 independent analytic audit', output['pairs'],'entries', output['interpolation_error_decimal'])
    print('NOT CERTIFIED: universal Weil positivity / Riemann Hypothesis')

if __name__=='__main__': main()
