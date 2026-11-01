#!/usr/bin/env python3
"""Deterministic and negative checks for V20 Cheng q-digamma extension."""
from fractions import Fraction as F
import json,subprocess,sys,tempfile
from pathlib import Path
import q_digamma_bose as m
for q in (F(1,2),F(3,4),F(9,10)):
    assert m.overlaps(m.lambert(q,1,128),m.divisor_enclosure(q,128))
    for x in (1,2,3):
        a,b=m.lambert(q,x,128)
        assert 0<a<b and b-a<F(1,10**8)
        assert m.psi_enclosure(q,x,128)[0]<m.psi_enclosure(q,x,128)[1]
for args in ((0,1,128),(1,1,128),(2,1,128),(F(1,2),0,128)):
    try:m.lambert(*args)
    except ValueError:pass
    else:raise AssertionError('invalid accepted')
q=F(1,2)
assert sum((q**n for n in range(1,129)),F(0))<m.lambert(q,1,128)[0]
with tempfile.TemporaryDirectory() as t:
    p1,p2=(Path(t)/'a.json',Path(t)/'b.json')
    for p in (p1,p2):
        r=subprocess.run([sys.executable,str(Path(__file__).with_name('q_digamma_bose.py')),'--out',str(p)],capture_output=True,text=True,timeout=60)
        assert r.returncode==0,r.stderr
    assert p1.read_bytes()==p2.read_bytes()
    report=json.loads(p1.read_text())
    assert report['checks']=={'evaluations':9,'divisor_crosschecks':3,'recurrences':6}
    assert report['not_certified']
print('PASS: independent q-series/recurrences and negative cases; deterministic evidence')
