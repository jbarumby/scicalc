#!/usr/bin/env python3
"""P245C stratified independent interval GL60 panel sums (not integral errors).

Panels spread across the full frequency interval. Where P236 individual
samples exist (503, 1007), independently computed boxes contain the Arb boxes.
The other panel enclosures have no external Arb reference in this package.
"""
from __future__ import annotations
import json, hashlib, argparse
from pathlib import Path
from fractions import Fraction as Q
from p245b_certified_gl60 import (I,iv,PI,LN2,PREC,certified_nodes,log_interval,
    re_digamma,sin_cos,spherical_j_small,source_ball_to_fraction,PILOT_HASHES,PAIRS as PILOT_PAIRS)
from p245c_chunk25 import LQ,PAIRS,INDICES,bound_bessel_high
PANELS=(100,250,503,750,875,1007)

def forward_all(x,kmax):
    s=spherical_j_small(x)
    a=[s[0],s[1]]
    for n in range(1,kmax):a.append((2*n+1)*a[n]/x-a[n-1])
    return a

def run(roots_file:Path,pilot_dir:Path):
    roots=certified_nodes(roots_file)
    log_pi=log_interval(PI)
    logs=[LN2,log_interval(3),2*LN2,log_interval(5)]
    coeff=[2*LN2/iv(2).sqrt(),2*logs[1]/iv(3).sqrt(),LN2,2*logs[3]/iv(5).sqrt()]
    beta=log_interval(iv(504)/(2*PI))-iv(Q(1,504))-sum(coeff,iv(0))
    L=iv(LQ)
    norms={par:{i:(2*L*(2*(2*i+par)+1)).sqrt() for i in INDICES if i!=376} for par in (0,1)}
    report={'schema':'NUMERIA.ES8P245C.STRATIFIED_PANEL_RATIONAL_ENCLOSURES.1',
      'parameters':{'L':'161/200','Tsharp':504,'panel_indices':list(PANELS),'bits':PREC},
      'panel_entries_reconstructed':len(PANELS)*len(PAIRS)*2,
      'externally_compared_panel_entries':0,'panels':[],
      'all_GL60_panels_independently_certified':False,'RH_proved':False}
    for panel in PANELS:
        vals={par:{pair:iv(0) for pair in PAIRS} for par in (0,1)}
        highs={par:bound_bessel_high(panel,par) for par in (0,1)}
        for root,w,_ in roots:
            t=(2*panel+1+root)/4
            x=L*t
            if x.l<30*(1<<PREC):raise ArithmeticError('forward Bessel requires x>=30')
            dig=re_digamma(t/2)
            ps=dig-log_pi-sum((coef*sin_cos(t*ln)[1] for coef,ln in zip(coeff,logs)),iv(0))
            z=w*(ps-beta)/(4*PI)
            js=forward_all(x,33)
            for par in (0,1):
                u={}
                for i in INDICES:
                    if i==376:
                        u[i]=I.endpoints(-highs[par],highs[par]);continue
                    n=2*i+par
                    v=js[n]
                    if i%2:v=-v
                    u[i]=norms[par][i]*v
                for i,j in PAIRS:
                    vals[par][i,j]=vals[par][i,j]+z*u[i]*u[j]
        for par in (0,1):
            sec='odd' if par else 'even'
            bounds=[]
            source_map={}
            source_name=None
            if panel in (503,1007):
                source_name=f'P236_{sec}_{panel}_{900 if panel==503 else 1152}.json'
                p=pilot_dir/source_name
                if hashlib.sha256(p.read_bytes()).hexdigest()!=PILOT_HASHES[source_name]:
                    raise ArithmeticError('source SHA mismatch '+source_name)
                d=json.loads(p.read_text())
                source_map={(i,j):b for i,j,b in d['gauss_interval_matrix_sample']}
            for pair in PAIRS:
                box=vals[par][pair]
                confirmed=False
                if pair in source_map:
                    lo,hi=source_ball_to_fraction(source_map[pair])
                    confirmed=box.contain_fraction(lo) and box.contain_fraction(hi)
                    if not confirmed:raise ArithmeticError(f'{sec}/{panel}/{pair} fails Arb inclusion')
                    report['externally_compared_panel_entries']+=1
                bounds.append({'pair':list(pair),'interval':box.desc(),
                               'external_Arb_interval_compared':pair in source_map,
                               'external_Arb_interval_contained':confirmed})
            report['panels'].append({'panel':panel,'sector':sec,'P236_source':source_name,
                                    'entries':bounds})
        print('stratified GL60 panel',panel,'both parities done',flush=True)
    return report

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--roots',type=Path,required=True)
    ap.add_argument('--pilots',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    r=run(a.roots,a.pilots);a.out.write_text(json.dumps(r,indent=2)+'\n')
    print('P245C stratified panel entries',r['panel_entries_reconstructed'],'external comparisons',r['externally_compared_panel_entries'])
if __name__=='__main__':main()
