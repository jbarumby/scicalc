#!/usr/bin/env python3
"""P237 conditional 377x377 REAL matrix interval LDL using Arb.

Requires FULL 1008-panel C for BOTH sectors. Adds analytical Gauss GL60
remainder ONCE, then beta and rank-one poles, then interval LDL. The GL60
radius and entrywise positivity (if any) are conditional on imported classic
analytic lemmas; a finite-dimensional LDL is NOT a Weil/RH certificate.

Requires python-flint 0.9.0 with Arb ball arithmetic.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,math,time
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,arb_mat,ctx
from p237_core import N,UPPER,upper_indices,as_endpoint_interval
from p237_beta_rational import beta_bounds
from p237_beta_tight import beta_tight,decimal_enclosure
from p237_gate_preflight import check_shadow


def ball(x):
    m,e=x.mid().man_exp();r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]


def arb_dyadic(m,e):
    return arb(m)*(arb(2)**e)


def enclosing_ball_wire(wire):
    lo,hi,e=as_endpoint_interval(wire)
    a=arb_dyadic(lo,e).lower();b=arb_dyadic(hi,e).upper()
    return (a+b)/2+arb(0,((b-a)/2).upper())


def bases():
    df=1;res=[];L=Q(161,200)
    for n in range(754):
        df*=2*n+1
        B=min(Q(2),Q(2*406**n,df))
        target=2*L*(2*n+1)
        b=math.isqrt((target.numerator+target.denominator-1)//target.denominator)
        while b*b*target.denominator<target.numerator:b+=1
        res.append(B*b)
    assert max(res)<=84
    return res


def robust_radius(q):
    if q<0:raise ValueError('negative error')
    return (arb(q.numerator)/arb(q.denominator)).upper()


def validate_full(d,sector):
    if d.get('schema')!='NUMERIA.ES8P237.EXACT_DYADIC_QUADRATURE_OUTER.1':
        raise ValueError('unsupported complete quadrature schema')
    checks={'sector':sector,'L':'161/200','Tsharp':504,
            'n_modes':N,'max_degree':752 if sector=='even' else 753,
            'panels_covered':[0,1008],'panel_count':1008,'GL60_nodes_covered':60480,
            'upper_triangle_count':UPPER,'all_1008_panels_covered':True,
            'analytic_GL60_error_added':False,'beta_and_poles_added':False,
            'LDL_tested':False,'analytic_bridge_verified':False,'RH_proved':False}
    for k,v in checks.items():
        if d.get(k)!=v:raise ValueError(f'full C not authorized: {k}')
    chunks=d.get('chunks')
    if not isinstance(chunks,list) or len(chunks)!=41:
        raise ValueError('full matrix must have 41 chunks including P236')
    if chunks[0].get('panel_range')!=[0,25] or chunks[-1].get('panel_range')!=[1000,1008]:
        raise ValueError('bad first/last chunk')
    prev=0
    for z in chunks:
        a,b=z['panel_range']
        if a!=prev or not a<b<=1008:raise ValueError('gap or overlapping chunk')
        if len(z.get('sha256',''))!=64:raise ValueError('missing source sha256')
        prev=b
    if prev!=1008:raise ValueError('incomplete source coverage')
    upper=d.get('C_quadrature_outer_upper_triangle')
    if not isinstance(upper,list) or len(upper)!=UPPER:
        raise ValueError('bad matrix entries')
    return upper


def form_matrix(d,sector,bits):
    ctx.prec=bits
    entries=validate_full(d,sector)
    L=arb(161)/200;T=arb(504);two=arb(2);pi=arb.pi()
    l2=two.log();l3=arb(3).log();l5=arb(5).log()
    AL=two*l2/two.sqrt()+two*l3/arb(3).sqrt()+l2+two*l5/arb(5).sqrt()
    beta=(T/(two*pi)).log()-arb(1)/504-AL
    # Independent rational upper/lower beta enclosure: use for cross-validation.
    beta_lo,beta_hi,_,_=beta_bounds()
    from p237_core import as_ball_interval
    from fractions import Fraction as F
    def dyadq(m,e):return F(m*(1<<e),1) if e>=0 else F(m,1<<-e)
    mb=ball(beta);bl,bh,be=as_ball_interval(mb)
    if not F(beta_lo)<=dyadq(bl,be)<=dyadq(bh,be)<=F(beta_hi):
        raise ArithmeticError('Arb beta outside independent rational certificate')
    # A second, narrower independent exact certificate uses Machin's pi.
    tight_lo,tight_hi,_,_=beta_tight()
    if not tight_lo<=dyadq(bl,be)<=dyadq(bh,be)<=tight_hi:
        raise ArithmeticError('Arb beta not enclosed by second rational Machin certificate')
    if beta.lower()<=0 or beta.upper()>=2:
        raise ArithmeticError('beta gate inconclusive')
    z0=L/2
    pole=[two*((2*n+1)*L/2).sqrt()*(pi/(two*z0)).sqrt()*
          z0.bessel_i(arb(n)+arb(1)/2)
          for n in range(int(sector=='odd'),754,2)]
    if len(pole)!=N:raise AssertionError('pole rank wrong')
    abases=bases()
    aa=[abases[2*i+int(sector=='odd')] for i in range(N)]
    fac=Q(1008*50,3**120)
    # Tight enough to retain small high-mode GL60 errors, with directed radii.
    aa_rad=[robust_radius(a) for a in aa]
    fac_rad=robust_radius(fac)
    R=[[None]*N for _ in range(N)]
    witnesses=[];max_error_bound=Q(0)
    for (i,j),row in zip(upper_indices(),entries):
        if type(row)!=list or len(row)!=3 or (row[0],row[1])!=(i,j):
            raise ValueError('missing/duplicate/reordered full matrix upper entry')
        C=enclosing_ball_wire(row[2])
        # Outward error bound for the GLOBAL 1008-panel sum, added once.
        E=fac_rad*aa_rad[i]*aa_rad[j]
        e_exact=fac*aa[i]*aa[j]
        if E.lower() < 0 or not E.lower()>=0:
            raise ArithmeticError('nonpositive budget')
        if max_error_bound<e_exact:max_error_bound=e_exact
        # E.upper() is guaranteed >= exact degreewise budget.
        val=C+arb(0,E.upper())
        if i==j:val+=beta
        pp=two*pole[i]*pole[j]
        val=val+pp if sector=='even' else val-pp
        if not val.is_finite():raise ArithmeticError('nonfinite assembled matrix')
        R[i][j]=R[j][i]=val
        if (i,j) in {(0,0),(0,16),(0,376),(16,16),(100,100),(199,199),(376,376)}:
            witnesses.append([i,j,ball(C),ball(E),ball(val)])
    return R,{'beta_ball':ball(beta),
              'beta_independent_machin_decimal_bounds':decimal_enclosure(tight_lo,tight_hi),
              'AL_ball':ball(AL),
              'pole_vector_balls':[ball(x) for x in pole],
              'GL60_global_error_row_analytic_fraction_coefficient':str(fac),
              'E_max_entry_exact':str(max_error_bound),
              'witnesses_i_j_C_E_R':witnesses}


def interval_ldl(R):
    # Mathematically sound but may be inconclusive due to dependency blowup.
    # Every positive interval pivot implies strict positivity for the FULL
    # symmetric interval matrix defined by R. No test of reduced Weil done here.
    ll=[[arb(0) for _ in range(N)] for _ in range(N)]
    dd=[];pivots=[]
    for k in range(N):
        p=R[k][k]
        for a in range(k):p-=ll[k][a]*ll[k][a]*dd[a]
        if p.lower()<=0 or not p.is_finite():
            return False,k,pivots,ball(p) if p.is_finite() else None
        dd.append(p);pivots.append(ball(p))
        for i in range(k+1,N):
            t=R[i][k]
            for a in range(k):t-=ll[i][a]*ll[k][a]*dd[a]
            ll[i][k]=t/p
        if (k+1)%25==0:print('LDL positive pivots',k+1,'/',N,flush=True)
    return True,N,pivots,None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--folder',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--bits',type=int,default=1152)
    ap.add_argument('--shadow-audit',type=Path,required=True,help='independent exact summation audit, bound to full C SHA256')
    args=ap.parse_args()
    if args.bits<900:raise ValueError('precision < 900 forbidden')
    # Fail closed: no matrix is formed unless BOTH full sectors have an
    # authenticated independent exact re-aggregation audit of their bytes.
    verified_digests=check_shadow(args.shadow_audit,args.folder)
    # Fail closed: no matrix is formed unless BOTH full sectors have been checked.
    source={}
    for sec in ('even','odd'):
        p=args.folder/f'P237_FULL_C_{sec}.json.gz'
        raw=p.read_bytes();d=json.loads(gzip.decompress(raw));validate_full(d,sec)
        if hashlib.sha256(raw).hexdigest()!=verified_digests[sec]:
            raise ArithmeticError('shadow digest not bound to LDL input')
        source[sec]={'data':d,'sha256':hashlib.sha256(raw).hexdigest()}
    results={}
    for sec in ('even','odd'):
        t=time.monotonic()
        R,meta=form_matrix(source[sec]['data'],sec,args.bits)
        success,pivotcount,pivots,failed=interval_ldl(R)
        results[sec]={'source_full_C_sha256':source[sec]['sha256'],
          'beta_poles_and_GL60_added_once':True,
          'interval_LDL_positive_pivot_count':pivotcount,
          'interval_LDL_complete':success,
          'first_inconclusive_pivot_index':None if success else pivotcount,
          'first_inconclusive_pivot_ball':failed,
          'pivots':pivots,
          'assembled_full_R_upper_triangle':[[i,j,ball(R[i][j])]
             for i,j in upper_indices()],
          'assembly_metadata':meta,
          'elapsed_seconds_diagnostic':round(time.monotonic()-t,2)}
    report={'schema':'NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1',
       'L':'161/200','Tsharp':504,'n_modes_per_sector':N,
       'bits':args.bits,'python_flint_version':__import__('flint').__version__,
       'both_full_C_1008_panels_required':True,
       'full_quadrature_coverage_verified':True,
       'independent_shadow_exact_audit_sha256':hashlib.sha256(args.shadow_audit.read_bytes()).hexdigest(),
       'all_142506_interval_sums_shadow_audited':True,
       'global_GL60_error_added_once_per_matrix_entry':True,
       'analytical_Binet_Bessel_Bernstein_lemmas_imported':True,
       'independent_Weil_Zhu_analytic_bridge_formalized':False,
       '377_mode_finite_positivity_certified':all(r['interval_LDL_complete'] for r in results.values()),
       'infinite_dimensional_Weil_positivity_proved':False,
       'RH_proved':False,'sectors':results}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_bytes(gzip.compress(json.dumps(report,separators=(',',':')).encode(),9,mtime=0))
    print('P237 FINITE LDL',report['377_mode_finite_positivity_certified'],
          [(sec,results[sec]['interval_LDL_positive_pivot_count']) for sec in results],
          'SHA256',hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=='__main__':main()
