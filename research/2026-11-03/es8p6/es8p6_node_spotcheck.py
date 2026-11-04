#!/usr/bin/env python3
"""Exact-Fraction sample cross-check of the two historical nodal evaluators.

Both share several analytic identities, so intersection is a regression check,
not an independent full proof of the special-function inequalities.
"""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys

SNAPSHOT=Path(__file__).resolve().parent/'source_snapshot'
sys.path.insert(0,str(SNAPSHOT))
from fase4_fixedpoint_sem_arb import FixedRationalNodes
from fase4_nos_sem_arb import omega_interval, amplitude_interval

EXPECTED={
 'fase4_fixedpoint_sem_arb.py':'0935bb988f7a7bd4c83b36a4a873a71b1612eda00c77ac9bd6ba7daeef600bda',
 'fase4_nos_sem_arb.py':'151ec9a0a6c249f443bcbac278e8f004e02d268ef7590f81e78e5280872cbe49'
}

def intersection(a,b): return max(a[0],b[0])<=min(a[1],b[1])

def spotcheck():
    for name,sha in EXPECTED.items():
        assert hashlib.sha256((SNAPSHOT/name).read_bytes()).hexdigest()==sha
    obj=FixedRationalNodes(digits=40,N=64)
    events=[]
    for r in (Q(0),Q(1),Q(4),Q(10),Q(160)):
        z=obj.omega(r)
        a=(Q(z[0],obj.S),Q(z[1],obj.S))
        b=omega_interval(r,N=256)
        assert intersection(a,b),(r,'omega')
        events.append({'category':'omega','r':str(r),'overlap':True,
            'source_width':float(a[1]-a[0]),'independent_width':float(b[1]-b[0])})
    for k,r in ((1,Q(0)),(1,Q(4)),(2,Q(4)),(3,Q(10))):
        z=obj.amp(r,k)
        a=(Q(z[0],obj.S),Q(z[1],obj.S))
        b=amplitude_interval(r,k)
        assert intersection(a,b),(r,k,'amplitude')
        events.append({'category':'amplitude','mode':k,'r':str(r),'overlap':True,
            'source_width':float(a[1]-a[0]),'independent_width':float(b[1]-b[0])})
    return {'status':'9_EXACT_RATIONAL_INTERVAL_INTERSECTIONS',
            'n_cases':len(events),'checks':events,'source_sha256':EXPECTED,
            'limitation':'Both codes use shared classical identities; not a 30x1600x16 nodal replay'}

if __name__=='__main__':
    p=Path(__file__).resolve().parent/'RESULTADOS_NOS_ES8P6.json'
    p.write_text(json.dumps(spotcheck(),indent=2,ensure_ascii=False)+'\n')
    print('ES8-P6 nodes',9,'intersection checks passed')
