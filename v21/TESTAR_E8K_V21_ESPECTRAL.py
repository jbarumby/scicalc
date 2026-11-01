#!/usr/bin/env python3
"""V21 adversarial and exact-arithmetic regression tests (Python standard library)."""
import json
import subprocess
import sys
import tempfile
from fractions import Fraction as Q
from pathlib import Path

import AUDITAR_E8K_V21_ESPECTRAL as v


def rejects(call):
    try:call()
    except (ValueError,AssertionError):return True
    return False


assert rejects(lambda:v.exp_neg_enclosure(Q(-1)))
assert rejects(lambda:v.exp_neg_enclosure(Q(1),degree=16))
assert rejects(lambda:v.heat_trace_enclosure(Q(0),8))
assert rejects(lambda:v.spectral_zeta_enclosure(0,8))
assert rejects(lambda:v.spectral_zeta_enclosure(1,1))
assert rejects(lambda:v.normalized_mellin_integral(0,2))
assert rejects(lambda:v.reciprocal_power_series([Q(0),Q(1)],1))

# A finite rational oscillator/Dirichlet test cannot contain the zero mode:
assert v.spectral_zeta_enclosure(1,16)[0]>Q(1)
assert v.spectral_zeta_enclosure(1,16)[1]<Q(2)
assert v.normalized_mellin_integral(2,3)==Q(1,81)
assert v.normalized_mellin_integral(3,3)==Q(1,729)
assert v.normalized_mellin_integral(1,3)!=Q(1,3)  # n^-2 vs mutated n^-1

# Altering the Casimir coefficient sign or removing counterterm must fail.
coeff=v.abel_casimir_certificate()
assert coeff['abel_finite_part']=='-1/12'
assert coeff['zeta_minus_one_from_bernoulli']=='-1/12'
assert coeff['normalized_casimir_energy_E_over_pi_hbar_c_div_a']=='-1/24'
assert coeff['abel_finite_part']!='1/12'
assert coeff['laurent_coefficients_x_minus2_minus1_zero'][1]=='0'

for x in (Q(0),Q(1,10),Q(1,2),Q(4),Q(20)):
    a,b=v.exp_neg_enclosure(x)
    assert 0<a<=b<=1
for t in (Q(1,4),Q(1,2),Q(1)):
    a,b=v.heat_trace_enclosure(t,6)
    c,d=v.heat_trace_enclosure(t,10)
    assert a<=d and c<=b
    assert d-c<b-a
    assert c>0

with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)
    reports=[]
    for i in (1,2):
        path=p/f'v21-{i}.json'
        run=subprocess.run([sys.executable,str(Path(__file__).with_name('AUDITAR_E8K_V21_ESPECTRAL.py')),
            '--out',str(path)],capture_output=True,text=True,timeout=120)
        assert run.returncode==0,run.stderr+run.stdout
        reports.append(path.read_bytes())
    assert reports[0]==reports[1],'non deterministic V21 report'
    report=json.loads(reports[0])
    assert report['schema']=='NUMERIA.E8K.V21.SPECTRAL_ZETA_MELLIN_CASIMIR.1'
    assert report['checks']['mellin_single_eigenmode_identities']==24
    assert report['checks']['heat_trace_t_cases']==3
    assert len(report['not_certified'])>=5
    if (p0:=Path(__file__).with_name('RESULTADO_E8K_V21.json')).is_file():
        assert p0.read_bytes()==reports[0],'versioned report not reproduced'
print('PASS: exact spectral zeta, Mellin, heat tail and Casimir witness; 7 domain failures; 2 identical replays')
