#!/usr/bin/env python3
"""P236 real Arb 377-mode ONE-PANEL stress test (not integral or LDL certificate).

Checks high-order spherical Bessel with native Arb special functions and
forms a 377x377 interval quadrature matrix over a single GL60 panel.
Exports selected upper-triangle positions, 16x16 corners, and cross block.
No analytic Gauss remainder, beta or pole has been added to the output C.
"""
from __future__ import annotations
import argparse,hashlib,json,os,time
from pathlib import Path
from flint import arb, acb, arb_mat, ctx

def ball(x):
    m,e=x.mid().man_exp();r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]

def run(sector,panel,bits):
    if sector not in ('even','odd') or panel not in (0,503,1007) or bits<900:
        raise ValueError('unsupported gate parameters')
    ctx.prec=bits;t0=time.monotonic()
    N=377;T=504;L=arb(161)/200;two=arb(2);pi=arb.pi()
    l2=two.log();l3=arb(3).log();l5=arb(5).log()
    prime=[(two*l2/two.sqrt(),l2),
           (two*l3/arb(3).sqrt(),l3),
           (l2,two*l2),
           (two*l5/arb(5).sqrt(),l5)]
    A=sum((a for a,v in prime),arb(0))
    beta=(arb(T)/(two*pi)).log()-arb(1)/T-A
    if not 0<beta.lower():raise ArithmeticError('beta is not positive')
    modes=[2*i+int(sector=='odd') for i in range(N)]
    norm=[(-1)**i*(two*L*(2*n+1)).sqrt() for i,n in enumerate(modes)]
    roots=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    if not all(roots[k][0].lower()>roots[k+1][0].upper() for k in range(59)):
        raise ArithmeticError('roots not separated')
    if not sum((w for x,w in roots),arb(0)).contains(2):
        raise ArithmeticError('bad GL60 weights')
    ctr=arb(2*panel+1)/4;quarter=arb(1)/4
    V=[];WV=[];bsamples=[]
    observed_max_rad=arb(0);ordered_radii=[]
    for index,(root,wgl) in enumerate(roots):
        t=ctr+quarter*root;x=L*t
        if x.lower()<=0 or x.upper()>=arb(406):
            raise ArithmeticError('Bessel argument outside domain')
        pre=(pi/(two*x)).sqrt()
        vals=[pre*x.bessel_j(arb(n)+arb(1)/2) for n in modes]
        if index in (0,30,59):
            j0=vals[0] if sector=='even' else pre*x.bessel_j(arb(1)/2)
            if not j0.overlaps(x.sin()/x):
                raise ArithmeticError('native j0 mismatch')
            bsamples.append({'node_index':index,'t':ball(t),
               'degrees':[[n,ball(vals[k])] for k,n in enumerate(modes)
                          if k in (0,15,31,63,127,199,299,375,376)]})
        # Probe high-mode rounding precision without dictating a tolerance
        # that might invalidate genuine zeros or extremely small values.
        high=vals[-1]
        if index in (0,30,59):
            ordered_radii.append([index,ball(high)])
        dig=acb(arb(1)/4,t/2).digamma().real
        psi=dig-pi.log()-sum((a*(t*ln).cos() for a,ln in prime),arb(0))
        weight=quarter*wgl*(psi-beta)/pi
        u=[norm[k]*vals[k] for k in range(N)]
        V.extend(u);WV.extend([weight*v for v in u])
        if index in (0,10,20,30,40,50,59):
            print('sector',sector,'panel',panel,'node',index,
                  'elapsed',round(time.monotonic()-t0,1),'s',flush=True)
    J=arb_mat(60,N,V);WJ=arb_mat(60,N,WV)
    C=J.transpose()*WJ
    samples=[]
    indices=set()
    for i in range(16):
        for j in range(i,16):indices.add((i,j))
    for i in range(N-16,N):
        for j in range(i,N):indices.add((i,j))
    for i in range(16):
        for j in range(N-16,N):indices.add((i,j))
    for i,j in [(0,0),(0,16),(0,376),(16,16),(32,32),
                (100,100),(199,199),(299,299),(375,375),(376,376)]:
        indices.add((i,j))
    for i,j in sorted(indices):
        z=(C[i,j]+C[j,i])/2
        if not z.is_finite():raise ArithmeticError(f'nonfinite interval {i},{j}')
        samples.append([i,j,ball(z)])
    return {
      'schema':'NUMERIA.ES8P236.ARB_377_SINGLE_GL60_PANEL.1',
      'precision_bits':bits,'flint_version':__import__('flint').__version__,
      'sector':sector,'L':'161/200','Tsharp':T,'mode_count':N,
      'Legendre_degree_max':modes[-1],'panel_index':panel,
      'panel_range_t_exact':[f'{panel}/2',f'{panel+1}/2'],
      'root_count':60,'matrix_dimension':[N,N],
      'all_C_entries_computed':True,'sample_entries_exported':len(samples),
      'GL60_integral_analytic_remainder_added':False,
      'beta_and_poles_added':False,
      'LDL_or_positive_definiteness_tested':False,
      'gauss_interval_matrix_sample':samples,
      'Bessel_probes':bsamples,
      'high_degree_probes':ordered_radii,
      'github_sha':os.environ.get('GITHUB_SHA','local'),
      'elapsed_seconds_diagnostic':round(time.monotonic()-t0,2),
      'Riemann_hypothesis_proved':False}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--sector',choices=['even','odd'],required=True)
    ap.add_argument('--panel',type=int,required=True)
    ap.add_argument('--bits',type=int,default=900)
    ap.add_argument('--out',type=Path,required=True)
    x=ap.parse_args()
    d=run(x.sector,x.panel,x.bits)
    x.out.write_text(json.dumps(d,separators=(',',':'))+'\n')
    print('PASS native Arb GL60 full 377-mode single-panel',d['sector'],
          d['panel_index'],'elapsed',d['elapsed_seconds_diagnostic'],
          'SHA256',hashlib.sha256(x.out.read_bytes()).hexdigest(),flush=True)
