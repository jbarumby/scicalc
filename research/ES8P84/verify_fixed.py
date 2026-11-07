#!/usr/bin/env python3
"""ES8-P84 fail-closed verification of six fixed-110-node Arb jobs."""
import json,hashlib,argparse
from pathlib import Path
from fixed_order_replay import sample_pairs,compare,interval,DEGREES,HERE
def run(folder):
    root=Path(folder)
    raw=(HERE/'sample_reference.json').read_bytes()
    base=json.loads(raw)
    digest=hashlib.sha256(raw).hexdigest()
    errors=[]; counts={}
    for win in DEGREES:
        seen=set(); historical=set(base['entries'][win]); counts[win]=0
        for sh in range(2):
            file=root/f'{win}_fixed_{sh}.json'
            if not file.exists():errors.append(f'MISSING:{file.name}');continue
            try:
                d=json.loads(file.read_bytes())
                for name,value in {'window':win,'mode':'sample','status':'ARB_REPLAY_COMPLETE',
                      'shard':sh,'shards':2,'bits':256,'flint_version':'0.9.0',
                      'fixture_sha256':digest,'quadrature':'fixed-110-nodes-global-M600',
                      'n_per_segment':110,'uniform_M':600}.items():
                    if d.get(name)!=value:errors.append(f'METADATA:{file.name}:{name}')
                expected=set(sample_pairs(base,win,sh,2));shard_seen=set()
                for r in d.get('records',[]):
                    p=r.get('pair')
                    if p not in expected or p in seen or p in shard_seen:
                        errors.append(f'BAD_PAIR:{file.name}:{p}');continue
                    seen.add(p);shard_seen.add(p)
                    rr=compare(r['replay'],base['entries'][win][p])
                    if rr=='DISJOINT' or rr!=r.get('relation'):
                        errors.append(f'BAD_INTERVAL:{file.name}:{p}')
                    counts[win]+=1
                if shard_seen!=expected:errors.append(f'INCOMPLETE_SHARD:{file.name}')
            except (ValueError,KeyError,TypeError,ZeroDivisionError) as e:
                errors.append(f'INVALID:{file.name}:{e}')
        if seen!=historical:errors.append(f'INCOMPLETE_WINDOW:{win}')
    return {'status':'VERIFIED_FIXED_110_ARB_SAMPLE' if not errors else 'FAIL_CLOSED',
            'checks':counts,'total':sum(counts.values()),
            'errors':errors,
            'scope':'30 fixed-order intervals compared exactly to historical boxes; adaptive replay comparison and domain/Schur require separate review'}
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('folder');p.add_argument('--out',default='ES8P84_VERDICT.json')
    a=p.parse_args();v=run(a.folder)
    Path(a.out).write_text(json.dumps(v,indent=2)+'\n')
    print(json.dumps(v,indent=2))
    if v['errors']:raise SystemExit(1)
