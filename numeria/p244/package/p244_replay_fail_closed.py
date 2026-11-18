#!/usr/bin/env python3
"""P244 strict replay: regenerate P242 MPFR traces from source and compare bytes.

Unlike p242_verify_exact.py alone, tampered inverse bounds or pivots cannot be
accepted: the proof trace is regenerated from the hashed original P237 matrix.
Exact byte equality is intentionally strict: a different MPFR version may fail
closed even if both certificates are mathematically correct.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from p242_verify_exact import audit

ROOT = Path(__file__).resolve().parent
EXPECTED_CODE_HASHES = {
    'p242_mpfr_interval_ldl.c':'3fe68cac0ba12d72b1036d42b4b3c30933c3c2db2f3f520edbf846071d217009',
    'p242_verify_exact.py':'afec8dd53f4960f74f0ea12a48fbc6589f56f9fded98325fbc8506b61abc1652',
}
EXPECTED_TRACES={
   'even':'0080384607c00fd63d3b4c481c17555118791b36b19372b91ea47cb40984fd2b',
   'odd':'d61f312cd24c3f2b518692798022ea90dfcc81b0f4d42d9fad321d16dd1ba5a5',
}

def digest(p:Path)->str: return hashlib.sha256(p.read_bytes()).hexdigest()

def validate_replays(work:Path, traces:Path, p242_cert:Path)->dict:
    prior=json.loads(p242_cert.read_text())
    result={}
    for sec in ('even','odd'):
        filename=f'P242_PIVOTS_{sec}_1024.txt'
        rebuilt=work/filename
        original=traces/filename
        e=EXPECTED_TRACES[sec]
        if digest(original)!=e or prior['sectors'][sec]['independent_trace_sha256']!=e:
            raise ArithmeticError(f'{sec} original trace SHA mismatch')
        if digest(rebuilt)!=e:
            raise ArithmeticError(f'{sec} regenerated trace SHA mismatch (refuse unbound data)')
        if rebuilt.read_bytes()!=original.read_bytes():
            raise ArithmeticError(f'{sec} replay byte mismatch')
        result[sec]={'trace_sha256':e,'identical_reexecution':True}
    return result

def run(data:Path,traces:Path,cert:Path,work:Path,recompute:bool)->dict:
    for name,expected in EXPECTED_CODE_HASHES.items():
        if digest(ROOT/name)!=expected:
            raise ArithmeticError(f'P244 refuses modified source: {name}')
    work.mkdir(parents=True,exist_ok=True)
    if recompute:
        cmd=['python3',str(ROOT/'export_centers.py'),'--source',str(data/'P237_FINITE_LDL_2048.json.gz'),'--out-dir',str(work)]
        subprocess.run(cmd,check=True)
        subprocess.run(['gcc','-O2','-Wall','-Wextra',str(ROOT/'p242_mpfr_interval_ldl.c'),'-o',str(work/'mpfr_ldl'),'-l:libmpfr.so.6','-lgmp'],check=True)
        for sec in ('even','odd'):
            subprocess.run([str(work/'mpfr_ldl'),str(work/f'P242_CENTERS_{sec}.txt'),str(work/f'P242_PIVOTS_{sec}_1024.txt'),'1024','377'],check=True)
    linked=validate_replays(work,traces,cert)
    full=audit(data,work)
    if full['total_independent_positive_pivots']!=754:
        raise ArithmeticError('P242 independent trace incomplete')
    return {
      'schema':'NUMERIA.ES8P244.STRICT_REPLAY.1',
      'P237_source_sha256':digest(data/'P237_FINITE_LDL_2048.json.gz'),
      'original_p242_certificate_sha256':digest(cert),
      'recomputed_using_compiled_MPFR_C_in_this_run':bool(recompute),
      'MPFR_traces':linked,
      'positive_pivots':754,
      'margin_lower_exact':'1/10^18',
      'regenerated_trace_matches_original':True,
      'original_P242_Python_verifier_accepts_forged_inverse_norm':True,
      'all_original_P237_matrix_assembly_integrals_independently_redone':False,
      'Weil_RH_proved':False}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--traces',type=Path,required=True)
    p.add_argument('--cert',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--check-existing',action='store_true',help='verify already recomputed traces; default executes MPFR')
    a=p.parse_args()
    x=run(a.data,a.traces,a.cert,a.work,not a.check_existing)
    a.out.write_text(json.dumps(x,indent=2)+'\n')
    print('P244 strict source-bound replay: PASS, 754 pivots, inverse norms, and Weyl bound')
