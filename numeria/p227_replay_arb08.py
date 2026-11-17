#!/usr/bin/env python3
"""P227: alternate Arb-0.8.0 native Bessel J and native digamma panel replay.
The result contains intervals for quadrature sums, NOT certified integrals.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
from flint import arb,acb,arb_mat,ctx
SCHEMA='NUMERIA.ES8P227.ARB08_NATIVE_SPECIAL_PANEL.1'
def ball(x):
    m,e=x.mid().man_exp();r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]
def sphere(x):
    if not (x.lower()>0 and x.upper()<161):raise ValueError('Bessel domain')
    c=(arb.pi()/(2*x)).sqrt(); j=[None]*402
    j[400]=c*x.bessel_j(arb(400)+arb(1)/2)
    j[401]=c*x.bessel_j(arb(401)+arb(1)/2)
    for n in range(400,0,-1):j[n-1]=(2*n+1)*j[n]/x-j[n+1]
    assert j[0].overlaps(x.sin()/x)
    assert j[1].overlaps(x.sin()/(x*x)-x.cos()/x)
    if j[0].rad().upper()>=arb(10)**-46 or j[1].rad().upper()>=arb(10)**-46:
        raise ArithmeticError('Bessel radii too wide')
    return j
def calculate(sector,panel,bits):
    T=200 if sector=='even' else 150
    assert sector in ('even','odd') and 0<=panel<2*T and bits>=800
    ctx.prec=bits;pi=arb.pi();L=arb(4)/5;two=arb(2)
    l2=two.log();l3=arb(3).log();l4=2*l2
    a2=2*l2/two.sqrt();a3=2*l3/arb(3).sqrt();a4=l2
    beta=(arb(T)/(2*pi)).log()-arb(1)/T-a2-a3-a4
    modes=[2*k+int(sector=='odd') for k in range(200)]
    norms=[(-1)**(n//2)*(2*L*(2*n+1)).sqrt() for n in modes]
    quarter=arb(1)/4;center=arb(2*panel+1)/4
    v=[];wv=[]
    roots=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    for root,w in roots:
        t=center+quarter*root;j=sphere(L*t)
        u=[norms[i]*j[n] for i,n in enumerate(modes)]
        dig=acb(arb(1)/4,t/2).digamma().real
        psi=dig-pi.log()-a2*(t*l2).cos()-a3*(t*l3).cos()-a4*(t*l4).cos()
        fac=quarter*w*(psi-beta)/pi
        v.extend(u);wv.extend([fac*a for a in u])
    C=arb_mat(60,200,v).transpose()*arb_mat(60,200,wv)
    inds=[(0,0),(0,1),(1,1),(0,199),(10,120),(199,199)]
    for i,j in inds:
        if not C[i,j].is_finite() or C[i,j].rad().upper()>=arb(10)**-44:
            raise ArithmeticError(f'wide/nonfinite ({i},{j})')
    return {'schema':SCHEMA,'version':__import__('flint').__version__,
        'backend':'Arb native J_(n+1/2) + native digamma','bits':bits,
        'sector':sector,'Tsharp':T,'L':'4/5','panel':panel,
        'nodes':60,'modes':200,'gauss_remainder_added':False,
        'beta_and_poles_added':False,
        'samples':[[i,j,ball(C[i,j])] for i,j in inds],
        'status':'INTERVAL_QUADRATURE_ONLY_NOT_GLOBAL_CERTIFICATE',
        'RH_proved':False,'commit':os.environ.get('GITHUB_SHA','local')}
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sector',choices=['even','odd'],required=True)
    p.add_argument('--panel',type=int,required=True)
    p.add_argument('--bits',type=int,default=1024)
    p.add_argument('--output',required=True)
    a=p.parse_args(); t=time.monotonic()
    result=calculate(a.sector,a.panel,a.bits)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',a.sector,a.panel,'secs',round(time.monotonic()-t,2),
          'version',result['version'],flush=True)
    print('SHA256',hashlib.sha256(Path(a.output).read_bytes()).hexdigest(),flush=True)
