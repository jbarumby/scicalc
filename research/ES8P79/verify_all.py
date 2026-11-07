#!/usr/bin/env python3
"""ES8-P81: fail-closed independent structure audit for 1328 Arb entries.
References 30 historical samples only. Full historical comparisons need
the external historical reference files from the ES8-P78 package.
"""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction
from replay import DEGREES,interval,compare

COUNTS={'a030':156,'a034':182,'a040':990}

def check(folder):
    folder=Path(folder)
    base=Path(__file__).resolve().parent
    fixture_bytes=(base/'sample_reference.json').read_bytes()
    sample=json.loads(fixture_bytes)
    fixture_hash=hashlib.sha256(fixture_bytes).hexdigest()
    errors=[]
    summaries={}
    sha={}
    for win,N in [(s,sum(v)) for s,v in DEGREES.items()]:
        all_keys={f'{i},{j}' for i in range(N) for j in range(i,N)}
        if len(all_keys)!=COUNTS[win]:raise AssertionError('wrong dimensions')
        seen=set();checked_sample=set()
        for sh in range(8):
            p=folder/f'{win}_all_{sh}.json'
            if not p.is_file():
                errors.append('MISSING:'+p.name);continue
            raw=p.read_bytes()
            sha[p.name]=hashlib.sha256(raw).hexdigest()
            try:
                d=json.loads(raw)
                req={'status':'ARB_REPLAY_COMPLETE','window':win,
                     'mode':'all','bits':256,'shard':sh,'shards':8,
                     'fixture_sha256':fixture_hash,'flint_version':'0.9.0'}
                for k,v in req.items():
                    if d.get(k)!=v:errors.append(f'META:{p.name}:{k}')
                expected={f'{i},{j}' for i in range(N) for j in range(i,N)}
                expected=sorted(expected,key=lambda s:tuple(map(int,s.split(','))))[sh::8]
                keys=set()
                for row in d.get('records',[]):
                    k=row.get('pair')
                    if k not in all_keys or k in seen or k in keys:errors.append(f'DUPLICATE:{p.name}:{k}');continue
                    keys.add(k);seen.add(k)
                    interval(row['replay'])
                    if k in sample['entries'][win]:
                        if compare(row['replay'],sample['entries'][win][k])!='REPLAY_SUBSET_REFERENCE':
                            errors.append(f'BAD_HISTORICAL_SAMPLE:{win}:{k}')
                        elif row.get('relation')!='REPLAY_SUBSET_REFERENCE':
                            errors.append(f'BAD_RELATION:{win}:{k}')
                        checked_sample.add(k)
                    elif row.get('relation')!='NOT_COMPARED':
                        errors.append(f'UNEXPECTED_RELATION:{win}:{k}')
                if set(keys)!=set(expected):errors.append('SHARD_INCOMPLETE:'+p.name)
            except (ValueError,KeyError,TypeError,ArithmeticError) as ex:
                errors.append(f'INVALID:{p.name}:{ex}')
        if seen!=all_keys:errors.append('WINDOW_INCOMPLETE:'+win)
        if checked_sample!=set(sample['entries'][win]):errors.append('SAMPLE_INCOMPLETE:'+win)
        summaries[win]={'all_expected':len(all_keys),'seen':len(seen),
                        'historical_samples_checked':len(checked_sample),
                        'full_historical_comparison':'NOT_DONE'}
    return {'status':'RECORDED_1328_INTERVALS' if not errors else 'FAIL_CLOSED',
            'windows':summaries,'artifacts_sha256':sha,'errors':errors,
            'scientific_limit':'Full historical 1328-entry comparisons and mathematical proof of analytic error bounds remain to be audited; not RH.'}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('folder')
    p.add_argument('--out',default='ES8P81_AUDIT.json')
    a=p.parse_args()
    res=check(a.folder)
    Path(a.out).write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps(res,indent=2))
    if res['errors']:raise SystemExit(1)
