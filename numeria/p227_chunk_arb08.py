#!/usr/bin/env python3
"""P227 Arb 0.8 independent native Bessel/digamma 50-panel chunks.
Encloses quadrature sums only; no quadrature theorem, beta or poles added.
"""
import argparse,gzip,hashlib,json,os,time
from pathlib import Path
from flint import arb,acb,arb_mat,ctx
from p227_replay_arb08 import sphere,ball
SCHEMA='NUMERIA.ES8P227.ARB08_NATIVE_50PANELS.1'
def compute(sector,start,stop,bits=1024):
    if sector not in ('even','odd') or bits<1024 or start%50!=0 or stop-start!=50:
        raise ValueError('P227 geometry/precision')
    T=200 if sector=='even' else 150
    if not 0<=start<stop<=2*T:raise ValueError('outside sector')
    ctx.prec=bits;p=arb.pi();L=arb(4)/5;two=arb(2)
    l2=two.log();l3=arb(3).log();l4=2*l2
    a2=2*l2/two.sqrt();a3=2*l3/arb(3).sqrt();a4=l2
    beta=(arb(T)/(2*p)).log()-arb(1)/T-a2-a3-a4
    modes=[2*k+int(sector=='odd') for k in range(200)]
    norms=[(-1)**(n//2)*(2*L*(2*n+1)).sqrt() for n in modes]
    gl=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    assert sum((w for x,w in gl),arb(0)).contains(2)
    assert all(gl[i][0].lower()>gl[i+1][0].upper() for i in range(59))
    quarter=arb(1)/4;result=arb_mat(200,200)
    samples=[];t0=time.monotonic()
    for panel in range(start,stop):
        center=arb(2*panel+1)/4;V=[];WV=[]
        for root,w in gl:
            t=center+quarter*root;j=sphere(L*t)
            u=[norms[i]*j[n] for i,n in enumerate(modes)]
            dig=acb(arb(1)/4,t/2).digamma().real
            psi=dig-p.log()-a2*(t*l2).cos()-a3*(t*l3).cos()-a4*(t*l4).cos()
            wt=quarter*w*(psi-beta)/p
            V.extend(u);WV.extend([wt*a for a in u])
        mat=arb_mat(60,200,V).transpose()*arb_mat(60,200,WV)
        if panel in (start,stop-1):
            samples.append({'panel':panel,'values':[[i,j,ball(mat[i,j])] for i,j in [(0,0),(0,1),(1,1),(0,199),(10,120),(199,199)]]})
        result+=mat
        print(f'{sector} Arb08 panel {panel+1}/{stop} {time.monotonic()-t0:.1f}s',flush=True)
    ent=[[i,j,ball(result[i,j])] for i in range(200) for j in range(i,200)]
    if len(ent)!=20100:raise AssertionError('entry count')
    return {'schema':SCHEMA,'backend':'python-flint/Arb 0.8.0','version':__import__('flint').__version__,
       'bits':bits,'sector':sector,'Tsharp':T,'L':'4/5','panel_range':[start,stop],
       'modes':200,'nodes':3000,'gauss_degree':60,
       'independent_methods':'native_Bessel_J_half_order_native_complex_digamma',
       'entries_count':len(ent),'quadrature_error_added':False,
       'beta_and_poles_added':False,'status':'QUADRATURE_SUM_ENCLOSURES_NOT_INTEGRAL',
       'RH_proved':False,'samples':samples,'entries':ent,
       'github_sha':os.environ.get('GITHUB_SHA','local'),
       'elapsed_seconds_diagnostic':round(time.monotonic()-t0,1)}
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--sector',choices=['even','odd'],required=True)
    parser.add_argument('--start',type=int,required=True)
    parser.add_argument('--stop',type=int,required=True)
    parser.add_argument('--bits',type=int,default=1024)
    parser.add_argument('--out',required=True)
    a=parser.parse_args();d=compute(a.sector,a.start,a.stop,a.bits)
    path=Path(a.out)
    path.write_bytes(gzip.compress(json.dumps(d,separators=(',',':')).encode(),compresslevel=9,mtime=0))
    print('PASS Arb08',d['sector'],d['panel_range'],
          'SHA256',hashlib.sha256(path.read_bytes()).hexdigest(),flush=True)
