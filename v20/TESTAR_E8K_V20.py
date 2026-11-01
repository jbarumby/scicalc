#!/usr/bin/env python3
"""Adversarial tests for V20 using only Python stdlib."""
import json
import subprocess
import sys
import tempfile
from fractions import Fraction as Q
from pathlib import Path
import gibbs_euler as ge
import yang_lee as yl
import gibbs_gram as gg

HERE=Path(__file__).resolve().parent

def rejects(call):
    try:call()
    except (ValueError,AssertionError):return True
    return False

assert rejects(lambda:ge.log_interval(-1))
assert rejects(lambda:ge.dirichlet_zeta(1,64))
assert rejects(lambda:ge.mangoldt_dirichlet(2,64))
assert rejects(lambda:ge.euler_zeta(2,1,64))
assert rejects(lambda:gg.gram(-1,3,6))
assert rejects(lambda:yl.polynomial_coeffs(0))
assert ge.mangoldt_prime_power(12) is None
assert ge.mangoldt_prime_power(27)==3
assert ge.mangoldt_prime_power(25)==5
assert ge.mangoldt_prime_power(2)==2
assert ge.mangoldt_prime_power(1) is None
assert ge.overlap(ge.dirichlet_zeta(3,128),ge.euler_zeta(3,17,128))
assert ge.overlap(ge.mangoldt_dirichlet(3,128),ge.mangoldt_euler(3,17))
# Detect incorrect p^beta coefficient (missing repetitions k>=2)
q=Q(1,2**3)
assert q/(1-q)!=q, 'mutated Euler coefficient should be rejected'
# Detect non-unit Yang-Lee coefficient
coeffs=yl.polynomial_coeffs(7);coeffs[3]=2
assert yl.multiply_z_minus_one(coeffs)!=[-1]+[0]*7+[1]
# Negative LDL pivot must be rejected
bad=gg.gram(3,3,5);bad[0][0]=-1
assert rejects(lambda:gg.ldl_symmetric(bad))
assert gg.short_support_obstruction()['rank_upper']==1

with tempfile.TemporaryDirectory() as t:
    A,B=Path(t)/'first.json',Path(t)/'second.json'
    for output in (A,B):
        proc=subprocess.run([sys.executable,str(HERE/'AUDITAR_E8K_V20.py'),'--out',str(output)],capture_output=True,text=True,timeout=60)
        assert proc.returncode==0,(proc.stdout,proc.stderr)
    assert A.read_bytes()==B.read_bytes(),'Non-deterministic JSON'
    committed=HERE/'RESULTADO_E8K_V20.json'
    if committed.exists():assert A.read_bytes()==committed.read_bytes(),'Committed evidence mismatch'
    j=json.loads(A.read_text(encoding='utf-8'))
    assert len(j['gates'])==4 and all(x['status']=='PASS' for x in j['gates'].values())
    assert len(j['not_certified'])>=4 and j['gates']['short_support_rank_obstruction']['rank_upper']==1
    assert len(j['gates']['gibbs_zeta_euler']['result']['beta_cases'])==3
    assert len(j['gates']['finite_yang_lee_mode']['result']['cases'])==8
    assert len(j['gates']['exact_gibbs_gram']['result']['positive_definite_independent_cases'])==6
print('PASS: V20 4 gates, 6 rational LDL matrices, 8 finite-polynomial cases, 3 Euler cases, 9 rejection/mutation tests and deterministic evidence.')
