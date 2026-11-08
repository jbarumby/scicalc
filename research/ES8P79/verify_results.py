#!/usr/bin/env python3
"""Fail-closed verification of six actual Arb outputs, no synthetic success fallback."""
import argparse,hashlib,json
from pathlib import Path
from replay import DEGREES,compare,interval,sample_pairs

HERE=Path(__file__).resolve().parent

def verify(folder):
    folder=Path(folder)
    b=(HERE/'sample_reference.json').read_bytes()
    fixture=json.loads(b)
    digest=hashlib.sha256(b).hexdigest()
    errors=[];per_window={};relations={}
    total=0
    for label in DEGREES:
        seen=set()
        expected=set(fixture['entries'][label])
        for shard in (0,1):
            expected_shard=set(sample_pairs(fixture,label,shard,2))
            path=folder/f'{label}_sample_{shard}.json'
            if not path.is_file():
                errors.append('MISSING:'+path.name);continue
            try:
                run=json.loads(path.read_text())
                required={'window':label,'mode':'sample','shard':shard,
                          'shards':2,'bits':256,'status':'ARB_REPLAY_COMPLETE',
                          'fixture_sha256':digest}
                for key,value in required.items():
                    if run.get(key)!=value:errors.append('BAD_META:'+path.name+':'+key)
                keys=set()
                for entry in run.get('records',[]):
                    key=entry['pair']
                    if key not in expected_shard or key in keys or key in seen:
                        errors.append('BAD_PAIR:'+label+':'+key)
                        continue
                    keys.add(key);seen.add(key)
                    recomputed=compare(entry['replay'],fixture['entries'][label][key])
                    if recomputed=='DISJOINT' or entry['relation']!=recomputed:
                        errors.append('BAD_INTERVAL:'+label+':'+key)
                    relations[recomputed]=relations.get(recomputed,0)+1
                    total+=1
                if keys!=expected_shard:errors.append('INCOMPLETE_SHARD:'+path.name)
            except (ValueError,KeyError,TypeError,ArithmeticError) as exc:
                errors.append('INVALID:'+path.name+':'+str(exc))
        if seen!=expected:errors.append('INCOMPLETE_WINDOW:'+label)
        per_window[label]={'expected':len(expected),'seen':len(seen)}
    return {'status':'VERIFIED_SAMPLE_ARTIFACTS' if not errors else 'FAIL_CLOSED',
            'total_checked':total,'windows':per_window,
            'relations':relations,'errors':errors,
            'scientific_limit':'checks runner artifacts, not analytic proof of Weil or RH'}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('folder')
    p.add_argument('--out',default='ES8P79_AUDIT.json')
    a=p.parse_args()
    report=verify(a.folder)
    Path(a.out).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if report['errors']:raise SystemExit(1)
