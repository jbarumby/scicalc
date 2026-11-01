#!/usr/bin/env python3
"""Negative/integrity/replay regressions for V19 using stdlib only."""
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('v19', ROOT / 'v19/AUDITAR_E8K_V19_CAUDA.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
E7 = ROOT / 'NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip'
PROOF = ROOT / 'referencias/TAIL_PROOF_ENVIADO_2026-11-01.md'
assert mod.digest(E7.read_bytes()) == mod.E7_SHA
assert mod.digest(PROOF.read_bytes()) == mod.SUBMITTED_SHA

pi = mod.pi_bracket()
assert abs(float((pi[0]+pi[1])/2)-math.pi)<1e-14  # float sanity only
for n in [2, 5, 160, 480]:
    lo, hi = mod.log_interval(n)
    assert abs(float((lo+hi)/2)-math.log(n))<1e-12  # float sanity only
witnesses, _ = mod.analytic_constant_witnesses()
assert F(witnesses['constant_lower']) < F('3.8013872') < F(witnesses['constant_upper'])
for N, R in [(10, 20), (30, 160)]:
    try:
        mod.radii(N, N, R, pi)
        raise AssertionError('rejected configuration was accepted')
    except ValueError:
        pass
assert mod.radii(20, 20, 160, pi)['p6'] > 0
assert mod.radii(10, 10, 256, pi)['p6'] < mod.radii(10, 10, 160, pi)['p6']
assert mod.radii(1, 2, 160, pi)['p6'] == 0

with tempfile.TemporaryDirectory() as tmp:
    temp = Path(tmp)
    def replay(e7, proof, out):
        return subprocess.run([
            sys.executable, str(ROOT / 'v19/AUDITAR_E8K_V19_CAUDA.py'),
            '--e7', str(e7), '--submitted-proof', str(proof), '--out', str(out)
        ], capture_output=True, text=True, timeout=90)
    a, b = temp/'a.json', temp/'b.json'
    assert replay(E7, PROOF, a).returncode == 0
    assert replay(E7, PROOF, b).returncode == 0
    assert a.read_bytes() == b.read_bytes(), 'non-deterministic V19'
    doc = json.loads(a.read_text())
    assert len(doc['records']) == 55
    assert doc['same_parity_independent_entries'] == 30
    assert doc['parity_zero_pairs'] == 25
    assert doc['not_certified']
    bad_e7 = temp/'invalid_E7.zip'
    bad_e7.write_bytes(E7.read_bytes()[:100])
    assert replay(bad_e7, PROOF, temp/'bad1.json').returncode == 2
    bad_proof = temp/'bad_tail.md'
    bad_proof.write_bytes(PROOF.read_bytes()+b'\nunauthorized revision\n')
    assert replay(E7, bad_proof, temp/'bad2.json').returncode == 2
print('PASS: V19 witness, 55 independent entries, 2 replay runs, 4 negative conditions, source integrity')
