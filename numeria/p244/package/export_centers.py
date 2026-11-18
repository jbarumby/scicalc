#!/usr/bin/env python3
"""Export P237 dyadic center entries for the independently coded C/MPFR LDL.

A center (midpoint) is stored as ((mantissa, exponent), (radius_mantissa, radius_exponent)).
For the SPD center replay we use only the first dyadic pair.
"""
import argparse, gzip, hashlib, json
from pathlib import Path

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args()
    raw=a.source.read_bytes();src=json.loads(gzip.decompress(raw))
    expected='2695768a3cca11237dd63cb40c0b6c34427b390329e5ff2ec089e31111c7898a'
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('P237 matrix SHA mismatch')
    if src.get('schema')!='NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1':
        raise ValueError('source matrix schema mismatch')
    a.out_dir.mkdir(parents=True,exist_ok=True)
    N=377
    for sec in ('even','odd'):
        rows=src['sectors'][sec]['assembled_full_R_upper_triangle']
        if len(rows)!=N*(N+1)//2: raise ValueError('matrix not 377x377')
        output=a.out_dir/f'P242_CENTERS_{sec}.txt'
        with output.open('w') as f:
            ix=0
            for j in range(N):
                for i in range(j+1):
                    # Stored data are upper triangular, in row-major (i,j) with i<=j.
                    # Index into the original row-major triangular representation.
                    k=i*N-i*(i-1)//2+(j-i)
                    ri,rj,b=rows[k]
                    if (ri,rj)!=(i,j):raise ValueError(f'wrong upper triangular row at {k}')
                    if not isinstance(b,list) or len(b)!=2:raise ValueError('invalid ball')
                    f.write(f'{i} {j} {b[0][0]} {b[0][1]}\n')
                    ix+=1
        print(sec,ix,'source',hashlib.sha256(raw).hexdigest(),'export_sha',hashlib.sha256(output.read_bytes()).hexdigest(),flush=True)

if __name__=='__main__':main()
