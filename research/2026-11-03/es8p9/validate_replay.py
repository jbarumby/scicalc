#!/usr/bin/env python3
"""Validate recomputed 32 weighted-sum hashes and 30 output intervals."""
import json,sys
from pathlib import Path

def verify(original,again):
    old=json.loads(Path(original).read_text());new=json.loads(Path(again).read_text())
    assert old['status']==new['status']=='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED'
    assert old['scope']==new['scope']
    assert old['quadrature']==new['quadrature']
    assert old['compatibility']==new['compatibility']
    assert old['entries']==new['entries']
    assert old['batches']==new['batches']
    return {'status':'REPLAY_EXACT_MATCH','batches':len(new['batches']),'entries':len(new['entries']),
            'evaluations':new['scope']['evaluations']}

if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit('usage: validate_replay.py expected.json replay.json')
    print(verify(sys.argv[1],sys.argv[2]))
