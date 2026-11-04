#!/usr/bin/env python3
"""Strictly compare a fresh ES8-P8 replay to the committed 32-block certificate."""
import argparse, json
from pathlib import Path

def validate(a,b):
    assert a == b, 'replayed interval ledger differs from versioned certificate'
    assert a['status']=='FULL_GRID_PAIRWISE_INTERVAL_INTERSECTIONS_CHECKED'
    assert a['count_positions']==25600 and a['count_compared']==768000
    assert a['failures_total']==0 and not a['failure_examples']
    assert len(a['chunks'])==32
    for n,x in enumerate(a['chunks']):
        assert x['start_cell']==n*50 and x['end_cell_exclusive']==(n+1)*50
        assert x['comparisons']==24000 and x['incompatible']==0
        assert len(x['sha256_pairwise_intervals'])==64
    assert sum(x['comparisons'] for x in a['chunks'])==768000
    return True

def main():
    p=argparse.ArgumentParser()
    p.add_argument('expected',type=Path)
    p.add_argument('replay',type=Path)
    x=p.parse_args()
    a=json.loads(x.expected.read_text()); b=json.loads(x.replay.read_text())
    validate(a,b)
    print('PASS: 768000 intervals, 32 matching SHA-256 blocks, exact replay equality')
if __name__=='__main__': main()
