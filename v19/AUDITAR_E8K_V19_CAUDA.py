#!/usr/bin/env python3
"""V19: independently replay the analytic archimedean tail assumptions/bounds.

All comparisons and bounds use fractions.Fraction, and elementary *stated*
analytical inequalities. No float decision. Not a proof of global Weil / RH.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import zipfile

E7_SHA = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
E7_TAIL_SHA = '5d3b0b199a1fd257def3aa07b7b49bca4665f078200966456a6bd715e9402e06'
SUBMITTED_SHA = 'f96300fc5dc0e5e88c15ae5928bdfec1eb34f0730768beedb6b3f7f9b2208cc6'
L = F(9,20)
R = F(160)


def digest(b):
    return hashlib.sha256(b).hexdigest()


def outward_decimal(value, upward=True, digits=40):
    """Integer-only directed quantization; decimal represents an exact bound."""
    value=F(value)
    scale=10**digits
    q=(value.numerator*scale)//value.denominator
    if upward and F(q,scale)<value:
        q+=1
    sign='-' if q<0 else ''
    q=abs(q)
    return f'{sign}{q//scale}.{q%scale:0{digits}d}'


def arctan_positive_recip(d, n=16):
    """Alternating series arctan(1/d); bound is next term (valid for d>1)."""
    assert d > 1 and n > 0
    s = sum(((-1)**j * F(1, (2*j+1)*d**(2*j+1)) for j in range(n)), F(0))
    nxt = F(1, (2*n+1)*d**(2*n+1))
    return (s, s+nxt) if n%2 == 0 else (s-nxt, s)


def pi_bracket():
    a0,a1=arctan_positive_recip(5)
    b0,b1=arctan_positive_recip(239)
    return 16*a0-4*b1,16*a1-4*b0  # Machin pi/4=4atan(1/5)-atan(1/239)


def log_interval(x, terms=40):
    """Exact rational interval using 2atanh q with a positive geometric tail.

    For x>0 reduce x=2**k*y, 1<=y<2. Then q=(y-1)/(y+1)<=1/3.
    """
    x=F(x)
    if x<=0:raise ValueError('log requires positive rational')
    if x==1:return F(0),F(0)
    k=x.numerator.bit_length()-x.denominator.bit_length()
    two=F(2)**k
    while two>x:k-=1;two/=2
    while 2*two<=x:k+=1;two*=2
    y=x/two
    assert F(1)<=y<F(2)
    q=(y-1)/(y+1)
    def series(z):
        s=sum((2*z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
        t=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return s,s+t
    ylo,yhi=series(q)
    llo,lhi=series(F(1,3))
    return (k*llo+ylo,k*lhi+yhi) if k>=0 else (k*lhi+ylo,k*llo+yhi)


def gamma_bracket(n=1200):
    """Euler gamma: H_n-log(n)-1/(2n) < gamma < H_n-log(n)-1/(2(n+1)).

    The classical harmonic/log inequalities are analytic inputs, NOT proven
    by this executable. Their intervals are recomputed exactly here.
    """
    h=sum((F(1,j) for j in range(1,n+1)),F(0))
    ll,lh=log_interval(F(n))
    return h-lh-F(1,2*n),h-ll-F(1,2*(n+1))


def analytic_constant_witnesses():
    pi0,pi1=pi_bracket()
    assert F(3)<pi0<pi1<F(22,7)
    ln2=log_interval(2)
    half_ln5=tuple(x/2 for x in log_interval(5))
    assert F(8047,10000)<half_ln5[0] and half_ln5[1]<F(8048,10000)
    gamma=gamma_bracket()
    log_pi=(log_interval(pi0)[0],log_interval(pi1)[1])
    c0=gamma[0]+3*ln2[0]+log_pi[0]
    c1=gamma[1]+3*ln2[1]+log_pi[1]
    assert F(38013,10000)<c0<c1<F(38014,10000), 'Gamma+3ln2+ln pi enclosure failed'
    cosh_lower=sum((pi0**(2*n)/factorial(2*n) for n in range(6)),F(0))
    cosh_term_upper=pi1/(2*cosh_lower)
    assert cosh_term_upper < F(1356,10000), 'pi/(2cosh pi) enclosure failed'
    # From monotone positive digamma series and Binet/duplication/reflection.
    assert -F(38014,10000)+F(8047,10000)-F(1356,10000)>-F(314,100)
    assert -F(38013,10000)+2+F(8048,10000)<-F(99,100)
    return {'pi_lower':outward_decimal(pi0,False),'pi_upper':outward_decimal(pi1,True),
            'gamma_lower':outward_decimal(gamma[0],False),'gamma_upper':outward_decimal(gamma[1],True),
            'log2_lower':outward_decimal(ln2[0],False),'log2_upper':outward_decimal(ln2[1],True),
            'constant_lower':outward_decimal(c0,False),'constant_upper':outward_decimal(c1,True),
            'cosh_correction_upper':outward_decimal(cosh_term_upper,True),
            'proven_inputs':['Machin arctan alternating remainder','atanh positive geometric remainder'],
            'assumed_analytic_lemmas':['Euler constant H_n-ln n enclosure','digamma series, duplication, reflection',
                                      'monotonicity/integral comparison for S(r)','cosh power-series positivity']},(pi0,pi1)


def radii(i,j,rr,pi_bounds):
    """Outward *upper* rational bounds for published P5 and E7 P6."""
    low,high=pi_bounds
    q=F(rr)
    assert i>=1 and j>=1 and q>=1
    # Require R>=N*pi/L, for both analytic estimates.
    if q <= max(i,j)*high/L:
        raise ValueError('R too small: R >= max(i,j)*pi/L not certified')
    if (i-j)%2:return {'p5':F(0),'p6':F(0),'tight_p5':F(0),'parity_zero':True}
    ln3r=log_interval(3*q)[1]
    lnr=log_interval(q)[1]
    base=F(i*j)*high/L**3/(3*q**3)
    p5=4*base*(8+ln3r+F(1,3))
    tighter=4*base*(F(314,100)+ln3r+F(1,3))
    epsi=(F(i)*high/(2*L*q))**2
    epsj=(F(j)*high/(2*L*q))**2
    assert epsi <=F(1,4) and epsj<=F(1,4)
    p6=base*(lnr+F(314,100)+F(1,3))/((1-epsi)*(1-epsj))
    return {'p5':p5,'p6':p6,'tight_p5':tighter,'parity_zero':False}


def structure_and_replay(e7_path,submitted_path):
    e7=e7_path.read_bytes()
    assert digest(e7)==E7_SHA, 'E7 SHA mismatch'
    assert digest(submitted_path.read_bytes())==SUBMITTED_SHA, 'Submitted TAIL_PROOF mismatch'
    with zipfile.ZipFile(e7_path) as archive:
        internal=archive.read('fontes/TAIL_PROOF.md')
        assert digest(internal)==E7_TAIL_SHA, 'E7 tail proof provenance changed'
        txt=internal.decode('utf-8')
        assert '## 8. Proposition 6' in txt and '1-\\varepsilon_i' in txt, 'E7 P6 missing'
        assert submitted_path.read_text(encoding='utf-8').count('## 8. Proposition 6')==0
        entry=json.loads(archive.read('RESULTADO_INTEGRAIS_E7.json'))
        prop6=json.loads(archive.read('RESULTADO_E7.json'))['prop6']
        for n,total in [(6,18),(8,32),(10,50)]:
            assert next(x['entries_same_parity'] for x in prop6 if x['N']==n)==total
        for n in (6,8,10):
            j=json.loads(archive.read(f'referencias_proposicao6/weil-interval-audit-v2-N{n}.json'))
            assert j['parameters']['N']==n and j['parameters']['tail_bound']=='lemma5'
            assert j['projected_dimension']==n-2
        return entry


def audit(e7_path,submitted_path):
    source=structure_and_replay(e7_path,submitted_path)
    witnesses,pi_bounds=analytic_constant_witnesses()
    assert R>10*pi_bounds[1]/L
    records=[]
    old_by_name={x['entry']:x for x in source['results']}
    checked=0
    for i in range(1,11):
        for j in range(i,11):
            v=radii(i,j,R,pi_bounds)
            if not v['parity_zero']:
                assert v['p6'] < v['p5']
                assert v['tight_p5'] <= v['p5']
                checked+=1
            rec={'i':i,'j':j,'same_parity':not v['parity_zero'],
                 'p5_upper_exact':outward_decimal(v['p5']), 'p6_upper_exact':outward_decimal(v['p6']),
                 'tight_p5_upper_exact':outward_decimal(v['tight_p5'])}
            if i==j and i in (1,2):
                previous=F(old_by_name[f'G{i}{i}']['tail_proposition6'])
                # Independent P6 bounds from differing rational pi/log algorithms
                # are comparable but not necessarily nested bitwise.
                discrepancy=previous-v['p6']
                rec['e7_p6_previous_exact']=outward_decimal(previous)
                rec['reconstruction_difference_exact']=outward_decimal(discrepancy,True)
                rec['relative_difference_approx']=format(float(discrepancy/previous),'.4e')
                assert abs(discrepancy) < F(1,10**22), 'Historical P6 bound disagrees materially'
            records.append(rec)
    assert checked==30 and len(records)==55
    # Mutations: same-parity pair must not certify with subcritical R; also N=30 must fail at R160.
    failed_cases=[]
    for n,r in [(10,20),(20,160),(30,160)]:
        try:
            radii(n,n,r,pi_bounds)
        except ValueError:
            failed_cases.append(f'N={n},R={r}')
    assert failed_cases==['N=10,R=20','N=30,R=160']
    # N=20 feasible at 160 and monotonicity with R increase for two configurations.
    assert radii(20,20,160,pi_bounds)['p6'] > 0
    assert radii(10,10,160,pi_bounds)['p6'] > radii(10,10,256,pi_bounds)['p6']
    return {'schema':'NUMERIA.E8K.V19.ANALYTIC_TAIL_REPLAY.1',
            'status':'P5_P6_CONSTANTS_AND_TAIL_ENCLOSURES_REPLAYED_CONDITIONALLY',
            'parameters':{'L':'9/20','R':'160','N_list':[6,8,10], 'upper_mode': 'Fraction exact'},
            'source_hashes':{'e7_zip':E7_SHA,'e7_tail_proof':E7_TAIL_SHA,'submitted_tail_proof':SUBMITTED_SHA},
            'witnesses':witnesses,
            'same_parity_independent_entries':checked,'parity_zero_pairs':25,
            'historical_comparison':'G11 and G22 P6 tails agree to <1e-22 (not byte equality)',
            'records':records,
            'negative_checks':{'subcritical_N10_R20_rejected':True,
                               'N30_R160_rejected':True,'N20_R160_accepted':True,
                               'R256_p6_tighter_than_R160':True},
            'not_certified':['V19 independently checks exact arithmetic and P5/P6 analytic consequences; not a mechanized proof of all external lemmas',
                             'Phase4 1600-cell derivative/interpolation bounds NOT independently recomputed here',
                             'Fase4 complete matrix LDL and all prime/pole normalization not rerun here',
                             'Applicability of Guinand-Weil as a global identity is external',
                             'No positivity for all N or a proof of RH']}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--e7',required=True,type=Path)
    p.add_argument('--submitted-proof',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    try:
        result=audit(args.e7,args.submitted_proof)
    except (AssertionError, ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        print('FAIL:',exc)
        return 2
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: 30 pares de mesma paridade; 25 zeros de paridade; P5/P6; constantes analiticas; 4 controles condicionais')
    print('STATUS:',result['status'])
    return 0

if __name__=='__main__':
    raise SystemExit(main())
