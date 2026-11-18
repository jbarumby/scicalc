#!/usr/bin/env python3
"""P236: actual 377x377 real Arb GL60 QUADRATURE CHUNK for first 25 panels.
Export ALL 71,253 upper entries as exact dyadic balls. No poles, beta,
analytical Gauss remainder, LDL, or RH claim. These are only PARTIAL sums.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,os,time
from pathlib import Path
from flint import arb,acb,arb_mat,ctx
from p236_377mode_panel import ball

N=377
def prefix_positions():
    s=set()
    for i in range(16):
        for j in range(i,16):s.add((i,j))
    for i in range(N-16,N):
        for j in range(i,N):s.add((i,j))
    for i in range(16):
        for j in range(N-16,N):s.add((i,j))
    for i,j in [(0,0),(0,16),(0,376),(16,16),(32,32),
                (100,100),(199,199),(299,299),(375,375),(376,376)]:
        s.add((i,j))
    if len(s)!=534:raise RuntimeError('sample indices wrong')
    return sorted(s)

def compute(sector,start=0,stop=25,bits=900):
    if sector not in ('even','odd') or (start,stop)!=(0,25) or bits!=900:
        raise ValueError('only first 25-panel chunk at 900 bits supported')
    ctx.prec=bits;start_time=time.monotonic()
    T=504;L=arb(161)/200;two=arb(2);pi=arb.pi()
    l2=two.log();l3=arb(3).log();l5=arb(5).log()
    prime=[(two*l2/two.sqrt(),l2),(two*l3/arb(3).sqrt(),l3),
           (l2,two*l2),(two*l5/arb(5).sqrt(),l5)]
    A=sum((a for a,v in prime),arb(0))
    beta=(arb(T)/(two*pi)).log()-arb(1)/T-A
    if beta.lower()<=0:raise ArithmeticError('beta positivity failure')
    modes=[2*i+int(sector=='odd') for i in range(N)]
    norm=[(-1)**i*(two*L*(2*n+1)).sqrt() for i,n in enumerate(modes)]
    roots=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    if not all(roots[i][0].lower()>roots[i+1][0].upper() for i in range(59)):
        raise ArithmeticError('GL60 roots order failure')
    if not sum((w for x,w in roots),arb(0)).contains(2):
        raise ArithmeticError('GL60 weight sum failure')
    accum=arb_mat(N,N);first_panel_sample=None
    quarter=arb(1)/4
    for panel in range(start,stop):
        center=arb(2*panel+1)/4;V=[];WV=[]
        for root,wgl in roots:
            t=center+quarter*root
            x=L*t
            if x.lower()<=0 or x.upper()>=arb(406):
                raise ArithmeticError('Bessel domain')
            pref=(pi/(two*x)).sqrt()
            vals=[pref*x.bessel_j(arb(n)+arb(1)/2) for n in modes]
            dig=acb(arb(1)/4,t/2).digamma().real
            psi=dig-pi.log()-sum((a*(t*ln).cos() for a,ln in prime),arb(0))
            w=quarter*wgl*(psi-beta)/pi
            u=[norm[k]*vals[k] for k in range(N)]
            V.extend(u);WV.extend([w*v for v in u])
        J=arb_mat(60,N,V);WJ=arb_mat(60,N,WV)
        P=J.transpose()*WJ
        if panel==start:
            first_panel_sample=[[i,j,ball((P[i,j]+P[j,i])/2)]
                 for i,j in prefix_positions()]
        accum+=P
        print('P236 377-mode partial chunk',sector,'panel',panel+1,'of',stop,
              'elapsed',round(time.monotonic()-start_time,1),'s',flush=True)
    entries=[[i,j,ball((accum[i,j]+accum[j,i])/2)]
          for i in range(N) for j in range(i,N)]
    if len(entries)!=71253:raise AssertionError('wrong upper triangle dimension')
    return {'schema':'NUMERIA.ES8P236.ARB_377_25PANEL_CHUNK.1',
       'sector':sector,'Tsharp':T,'L':'161/200','n_modes':N,
       'degree_max':modes[-1],'bits':ctx.prec,
       'flint_version':__import__('flint').__version__,
       'panel_range':[start,stop],'GL60_node_evaluations':60*(stop-start),
       'upper_triangle_count':len(entries),
       'source_method':'direct native bessel_J_half and digamma, GL60 roots independent',
       'C_partial_upper_triangle':entries,
       'first_panel_534_sample':first_panel_sample,
       'whole_integral_computed':False,'analytic_GL60_error_added':False,
       'beta_and_poles_added':False,'LDL_tested':False,
       'RH_proved':False,'github_sha':os.environ.get('GITHUB_SHA','unknown'),
       'elapsed_seconds_diagnostic':round(time.monotonic()-start_time,2)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--sector',choices=['even','odd'],required=True)
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--stop',type=int,default=25)
    p.add_argument('--bits',type=int,default=900)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();d=compute(a.sector,a.start,a.stop,a.bits)
    raw=json.dumps(d,separators=(',',':')).encode()
    a.out.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
    print('PASS 71253 intervals for 25 panels; SHA256',
       hashlib.sha256(a.out.read_bytes()).hexdigest(),flush=True)

if __name__=='__main__':main()
