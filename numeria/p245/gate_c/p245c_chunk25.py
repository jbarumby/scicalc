#!/usr/bin/env python3
"""P245 Gate C: independent rational interval reconstruction of 25 GL60 panels.

Computes selected entries for EVERY panel 0..24 and verifies whole-chunk sums
against independently generated P236 Arb intervals. The certified items are
GL60 quadrature sums, NOT the exact integral or complete 377x377 matrix.
"""
from __future__ import annotations
import argparse, hashlib, json, gzip, sys
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parent))
from p245b_certified_gl60 import (I, iv, S, PI, LN2, PREC, log_interval,
    certified_nodes, re_digamma, sin_cos, spherical_j_small,
    double_factorial_odd, source_ball_to_fraction, ceildiv)

N=377
LQ=Q(161,200)
PAIRS=((0,0),(0,1),(0,2),(0,3),(1,1),(1,2),(2,2),(3,3),
       (0,16),(16,16),(0,376),(1,376),(376,376))
INDICES=sorted(set(j for p in PAIRS for j in p))
PANEL_START=0
PANEL_STOP=25
HASHES={
 'even':'6cfc2d3080389751b08d62dbdcdabdbd16a8be64ea5b23993dac0f27f61383f5',
 'odd':'76a7c147bf805fca4b67f65cbfd01464cf95134b30e990d3b5f24226a4de1f3d'}

def bessel_series(x,n,terms=70):
    """Rigorous all-real j_n for 0<x<=11 using the entire Bessel series.

    |next terms| is bounded using x²<=121 and geometric tail ratio after K.
    """
    x=iv(x)
    if x.l<=0 or x.u>11*S:raise ValueError('Bessel series domain 0<x<=11')
    y=x*x
    v=iv(1);s=v
    q=Q(1)
    for k in range(1,terms):
        v=-v*y/Q(2*k*(2*n+2*k+1));s=s+v
        q*=Q(121,2*k*(2*n+2*k+1))
    k=terms
    q*=Q(121,2*k*(2*n+2*k+1))
    ratio=Q(121,2*(k+1)*(2*n+2*k+3))
    if ratio>=1:raise ArithmeticError('geometric tail not decaying')
    tail=q/(1-ratio)
    # add analytic bound as a rational interval, with outward rounding
    s=s+I.endpoints(-tail,tail)
    return x**n*s/double_factorial_odd(n)

def bound_bessel_high(panel,parity):
    n=2*376+parity
    x_max=Q(161*(panel+1),400)
    # j_n(x) <= x^n/(2n+1)!!; factor 50 is a universal norm bound
    return Q(50)*x_max**n/double_factorial_odd(n)

def triangular_pos(i,j):
    if i>j:i,j=j,i
    return i*N-i*(i-1)//2+(j-i)

def get_source(chunk_dir,sector):
    path=chunk_dir/f'P236_CHUNK25_{sector}.json.gz'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=HASHES[sector]:
        raise ArithmeticError(f'P236 source changed: {sector}')
    obj=json.loads(gzip.decompress(path.read_bytes()))
    if obj['schema']!='NUMERIA.ES8P236.ARB_377_25PANEL_CHUNK.1' or obj['panel_range']!=[0,25] or obj['sector']!=sector or len(obj['C_partial_upper_triangle'])!=N*(N+1)//2:
        raise ArithmeticError('P236 chunk metadata inconsistent')
    return obj,path

def evaluate_chunk(chunks:Path,roots_file:Path):
    roots=certified_nodes(roots_file)
    log_pi=log_interval(PI)
    logs=[LN2,log_interval(3),2*LN2,log_interval(5)]
    coeff=[2*LN2/iv(2).sqrt(),2*logs[1]/iv(3).sqrt(),LN2,2*logs[3]/iv(5).sqrt()]
    beta=log_interval(iv(504)/(2*PI))-iv(Q(1,504))-sum(coeff,iv(0))
    L=iv(LQ)
    norm={par:{i:(2*L*(2*(2*i+par)+1)).sqrt() for i in INDICES if i!=376}for par in (0,1)}
    sums={par:{pair:iv(0) for pair in PAIRS} for par in (0,1)}
    panel_results=[]
    for panel in range(PANEL_START,PANEL_STOP):
        vals={par:{pair:iv(0) for pair in PAIRS} for par in (0,1)}
        high={par:bound_bessel_high(panel,par) for par in (0,1)}
        for root,w,_ in roots:
            t=(2*panel+1+root)/4
            x=L*t
            dig=re_digamma(t/2)
            ps=dig-log_pi-sum((coef*sin_cos(t*ln)[1] for coef,ln in zip(coeff,logs)),iv(0))
            z=w*(ps-beta)/(4*PI)
            js=spherical_j_small(x)
            for par in (0,1):
                u={}
                for i in INDICES:
                    if i==376:
                        u[i]=I.endpoints(-high[par],high[par]);continue
                    deg=2*i+par
                    v=js[deg] if deg<=3 else bessel_series(x,deg)
                    if i%2:v=-v
                    u[i]=norm[par][i]*v
                for i,j in PAIRS:
                    vals[par][i,j]=vals[par][i,j]+z*u[i]*u[j]
        for par in (0,1):
            for pair in PAIRS:
                sums[par][pair]=sums[par][pair]+vals[par][pair]
        panel_results.append({'panel':panel,'entries_recomputed':len(PAIRS)*2})
        if panel%5==4:
            print(f'interval GL60 panels {panel-4}..{panel} completed',flush=True)
    records=[]
    for par in (0,1):
        sector='even' if par==0 else 'odd'
        src,path=get_source(chunks,sector)
        rec={'sector':sector,'source_hash':HASHES[sector], 'interval_pairs':[]}
        for pair in PAIRS:
            i,j=pair
            ii,jj,b=src['C_partial_upper_triangle'][triangular_pos(i,j)]
            if [ii,jj]!=[i,j]:raise ArithmeticError('triangle order mismatch')
            lo,hi=source_ball_to_fraction(b)
            result=sums[par][pair]
            inside=result.contain_fraction(lo) and result.contain_fraction(hi)
            if not inside:
                raise ArithmeticError(f'P245C mismatch: {sector} pair={pair} independent={result} Arb span={float(lo)}, {float(hi)}')
            rec['interval_pairs'].append({'pair':[i,j],'Legendre_degrees':[2*i+par,2*j+par],
              'interval':result.desc(),'source_interval_contained':inside})
        records.append(rec)
    return {'schema':'NUMERIA.ES8P245C.FIRST_25_PANEL_RATIONAL_REPLAY.1',
      'parameters':{'L':'161/200','Tsharp':504,'panels':[0,25],
      'nodes_per_panel':60,'precision_bits':PREC,'modes_per_parity':377},
      'source_root_brackets_sha256':hashlib.sha256(roots_file.read_bytes()).hexdigest(),
      'method':'fixed point outward 256 bits, classical digamma-Binet and Bessel series',
      'source_chunks':HASHES,'selected_pairs_per_sector':len(PAIRS),
      'all_panel_pair_evaluations':(PANEL_STOP-PANEL_START)*len(PAIRS)*2,
      'verified_chunk_sum_enclosures':len(PAIRS)*2,
      'chunk_sums':records,'panel_records':panel_results,
      'actual_integrals_fully_certified':False,'GL60_analytic_remainder_verified_from_first_principles':False,
      'all_142506_entries_certified_independently':False,'RH_proved':False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--chunks',type=Path,required=True)
    ap.add_argument('--roots',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    report=evaluate_chunk(args.chunks,args.roots)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print('P245C: certified',report['verified_chunk_sum_enclosures'],'chunk intervals and',report['all_panel_pair_evaluations'],'panel entries')
if __name__=='__main__': main()
