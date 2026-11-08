#!/usr/bin/env python3
"""Fail-closed gate: all five P110 pilots must overlap the original P109 Arb intervals.

The original and tight approximations use different support truncations;
overlap checks expose gross inconsistencies before a full 169x169 run.
Overlap alone is not an independent proof of analytic inclusion.
"""
import json, hashlib, math
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parent
PAIRS={'diagonal':'1,1','prime-neighbors':'2,3','near-contact':'197,433',
       'high-centres':'499,997','largest-separation':'1,997'}
def endpts(v):
    assert len(v)==4
    a,e,r,f=v
    m=F(int(a))*F(2)**int(e)
    rad=F(int(r))*F(2)**int(f)
    if rad<0:raise ValueError('negative Arb radius')
    return m-rad,m+rad

def gate():
    old=json.loads((ROOT/'old'/'P109_1000_COMPACT_MATRIX.json').read_text())
    if old['precision_bits']!=256 or len(old['entries'])!=14365:
        raise ValueError('bad original P109 169-centre input')
    src=set()
    records=[]
    for label,pair in PAIRS.items():
        fn=ROOT/'pilots'/f'P110_{label}.json'
        d=json.loads(fn.read_text())
        if (d['status']!='P109_ROW_SHARD_INTERVALS_CONDITIONAL_ON_GL_CAUCHY_ARITHMETIC'
            or d['bits']!=256 or d['epsilon']!='1/10000000000000000'
            or d['max_prime']!=1000 or d['pilot_pair']!=pair or len(d['entries'])!=1):
            raise ValueError('invalid tight pilot '+label)
        t=d.get('test_profile',{})
        if (t.get('z0')!='127/128' or t.get('inner_GL_error_per_panel')!='1e-22'
            or t.get('outer_GL_error_per_panel')!='1e-20'):
            raise ValueError('wrong GL tight configuration')
        source=d['source_sha256']
        if len(source)!=64:raise ValueError('source SHA invalid')
        src.add(source)
        a,b=endpts(d['entries'][pair])
        c,e=endpts(old['entries'][pair])
        if max(a,c)>min(b,e):
            raise ArithmeticError('NEW OLD INTERVALS DISJOINT: '+pair)
        tight_rad=(b-a)/2; old_rad=(e-c)/2
        if tight_rad >= F(1,10**16):
            raise ArithmeticError('TIGHT PILOT ABOVE PRECISION TARGET: '+pair)
        if not d['stats']['outer_complex_verified'] or not d['stats']['inner_complex_verified']:
            raise ArithmeticError('missing geometric ellipse records')
        records.append({
            'label':label,'pair':pair,'old_radius_rational':str(old_rad),
            'new_radius_rational':str(tight_rad),
            'improvement_factor_for_display':float(old_rad/tight_rad),
            'overlap_rational':True,
            'source_sha256':source,
            'verified_outer_ellipses':d['stats']['outer_complex_verified'],
            'verified_inner_ellipses':d['stats']['inner_complex_verified'],
            'artifact_sha256':hashlib.sha256(fn.read_bytes()).hexdigest()})
    if len(src)!=1:raise ValueError('mixed pilot sources')
    out={'status':'ALL_FIVE_TIGHT_PILOTS_COMPLETE_AND_EXACTLY_OVERLAP_P109',
         'source_sha256':src.pop(),'pilots':records,
         'target':'full 169x169 finite Weil matrix; conditional proof only',
         'warning':'No universal Weil positivity or RH proof.'}
    (ROOT/'P111_GATE_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print('P111_GATE_SUCCESS',len(records),'pilot enclosures; max new radius',
          max(float(F(p['new_radius_rational'])) for p in records),flush=True)

if __name__=='__main__':gate()
