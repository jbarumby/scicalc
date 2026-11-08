#!/usr/bin/env python3
"""Fail-closed exact-rational check for 8 Arb pole/prime shards; no archimedean claim."""
import json,hashlib,argparse
from fractions import Fraction
from pathlib import Path
def endpoints(x):
    m,p,r,q=map(int,x)
    if r<0:raise ValueError('negative radius')
    c=Fraction(m)*Fraction(2)**p
    rad=Fraction(r)*Fraction(2)**q
    return c-rad,c+rad
def check(folder):
    root=Path(folder);seen=set();errors=[];hashes={}
    expected=[(i,j) for i in range(8) for j in range(i,8)]
    widths={'poles':Fraction(0),'primes':Fraction(0)}
    for shard in range(8):
        f=root/f'shard_{shard}.json'
        if not f.is_file():errors.append('MISSING '+f.name);continue
        raw=f.read_bytes();hashes[f.name]=hashlib.sha256(raw).hexdigest()
        try:
            data=json.loads(raw)
            metadata={'status':'CERTIFIED_POLE_PRIME_COMPONENTS_ONLY',
                      'shard':shard,'shards':8,'precision_bits':256,
                      'python_flint':'0.9.0',
                      'center_integers':[1,2,3,5,7,11,13,17],
                      'h':'2/5'}
            for k,v in metadata.items():
                if data.get(k)!=v:errors.append('BAD_META '+f.name+' '+k)
            expected_shard=set(expected[shard::8]);subseen=set()
            for entry in data['records']:
                ij=entry['i'],entry['j']
                if ij not in expected_shard or ij in seen:
                    errors.append('BAD_PAIR '+str(ij));continue
                seen.add(ij);subseen.add(ij)
                for name in ('poles','primes'):
                    lo,hi=endpoints(entry[name]);width=hi-lo
                    widths[name]=max(widths[name],width)
                    if width>Fraction(1,10**8):
                        errors.append('TOO_WIDE '+name+' '+str(ij))
            if subseen!=expected_shard:errors.append('INCOMPLETE_SHARD '+str(shard))
        except (ValueError,KeyError,TypeError,OverflowError,ZeroDivisionError) as exc:
            errors.append('INVALID '+str(shard)+' '+str(exc))
    if seen!=set(expected):errors.append('INCOMPLETE_36_PAIRS')
    return {'status':'VALIDATED_POLE_PRIME_COMPONENTS' if not errors else 'FAIL_CLOSED',
            'pairs':len(seen),'errors':errors,'artifact_sha256':hashes,
            'width_bounds':{k:str(v) for k,v in widths.items()},
            'archimedean_status':'NOT_CERTIFIED',
            'weils_matrix_positive_definite':'NOT_CERTIFIED'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('directory')
    ap.add_argument('--out',default='ES8P89_VERDICT.json')
    args=ap.parse_args();result=check(args.directory)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if result['errors']:raise SystemExit(1)
