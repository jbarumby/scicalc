#!/usr/bin/env python3
"""Gate para V11: distingue testes aprovados de SKIP esperados; nunca os reclassifica."""
import argparse, json, subprocess, sys
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['independent', 'sourced'])
    p.add_argument('--script', type=Path, required=True)
    p.add_argument('--zip', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.mode == 'sourced' and a.zip is None:
        p.error('--zip exigido em sourced')
    a.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable,str(a.script),'--out',str(a.out)]
    if a.zip is not None: cmd += ['--zip',str(a.zip)]
    proc = subprocess.run(cmd,check=False)
    if not a.out.is_file():
        raise SystemExit('V11 não produziu JSON; código de saída '+str(proc.returncode))
    report = json.loads(a.out.read_text(encoding='utf-8'))
    gates=report.get('approval_gates',{})
    expected={
        'independent': {
            'independent_transcendental_cotejo':'PASS',
            'ast_methods':'SKIP','mutations':'SKIP',
            'edge_cases':'SKIP','rational_tests':'SKIP',
            'b5_pi_bounds_overlap':'SKIP',
        },
        'sourced': {
            'independent_transcendental_cotejo':'PASS',
            'ast_methods':'PASS','mutations':'PASS',
            'edge_cases':'SKIP','rational_tests':'SKIP',
            'b5_pi_bounds_overlap':'SKIP',
        },
    }[a.mode]
    observed={key:gates.get(key,{}).get('status') for key in expected}
    if any(value=='FAIL' for value in (v.get('status') for v in gates.values())):
        raise SystemExit('GATE V11 FALHOU; consulte JSON')
    if observed != expected or len(gates)!=len(expected):
        raise SystemExit(f'Padrão de gates V11 diferente do conhecido: {observed}')
    if proc.returncode != 2:
        raise SystemExit(f'V11 mudou o código de saída esperado (2): {proc.returncode}')
    print(f'V11 {a.mode}: comportamento conhecido confirmado; '
          f'{sum(v=="PASS" for v in observed.values())} PASS, '
          f'{sum(v=="SKIP" for v in observed.values())} SKIP (NÃO aprovados).')
    print('Gate de regressão aceito; NÃO é aprovação de suíte completa.')

if __name__=='__main__': main()
