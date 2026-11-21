#!/usr/bin/env python3
"""Seeded, NON-PROOF discovery experiments for the P73 Suzuki spectral window.

All scoring uses IEEE-754 floats. The final proof is a separate program:
certify_adaptive.py, which uses exact fractions and a full cell covering.
"""
import argparse
from fractions import Fraction as F
import hashlib
import heapq
import json
import math
from pathlib import Path
import random

H=3*10**12
HI_TAIL=5/H+38/H**2

def ordinates(data,n):
    vals=[]
    for z in data['zeros'][:n]:
        def dyad(x):return float(int(x[0])*F(2)**int(x[1]))
        vals.append((dyad(z['lo_dyad'])+dyad(z['hi_dyad']))/2)
    if len(vals)!=n:raise ValueError('Not enough ordinates')
    return vals

def score(t, gammas):
    truncated=sum(2*(1-math.cos(g*t))/(g*g) for g in gammas)
    estimated_tail=(3+math.exp(t/2))*HI_TAIL
    return truncated-estimated_tail

def discover(approach, gammas, left=44.35, right=47.48, budget=512, seed=73):
    if budget<3 or not left<right:raise ValueError('Invalid experiment parameters')
    visited={}
    def measure(x):
        if x not in visited:visited[x]=score(x,gammas)
        return visited[x]
    measure(left); measure(right)
    if approach=='uniform':
        for i in range(1,budget-1):measure(left+(right-left)*i/(budget-1))
    elif approach=='random':
        rng=random.Random(seed)
        while len(visited)<budget:measure(left+(right-left)*rng.random())
    elif approach=='adaptive':
        # Optimistic priority from endpoints plus a curvature uncertainty penalty.
        # This is a discovery heuristic, NOT a certificate for the interval.
        def priority(a,b):
            w=b-a
            return -w*(1+2/(1+max(0,min(measure(a),measure(b)))*100))
        heap=[(priority(left,right),left,right)]
        while len(visited)<budget:
            _,a,b=heapq.heappop(heap)
            mid=(a+b)/2
            if mid in visited:
                # Underflow of float midpoint may occur only for extreme budgets.
                raise ArithmeticError('Adaptive floating resolution exhausted')
            measure(mid)
            heapq.heappush(heap,(priority(a,mid),a,mid))
            heapq.heappush(heap,(priority(mid,b),mid,b))
    else:raise ValueError('Unknown strategy')
    worst_t=min(visited,key=lambda x:(visited[x],x))
    canonical=''.join(f'{t:.17g},{s:.17g}\n' for t,s in sorted(visited.items()))
    return {'strategy':approach,'seed':seed if approach=='random' else None,
            'evaluations':len(visited),'best_t':worst_t,'minimum_heuristic_score':visited[worst_t],
            'negative_samples':sum(v<0 for v in visited.values()),
            'first_negative_by_t':next((t for t in sorted(visited) if visited[t]<0),None),
            'trace_sha256':hashlib.sha256(canonical.encode()).hexdigest()}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--zeros',type=Path,required=True)
    p.add_argument('--n',type=int,default=20)
    p.add_argument('--from-t',type=float,default=44.35)
    p.add_argument('--to',type=float,default=47.48)
    p.add_argument('--budget',type=int,default=512)
    p.add_argument('--seed',type=int,default=73)
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    data=json.loads(args.zeros.read_text())
    G=ordinates(data,args.n)
    r=[discover(x,G,args.from_t,args.to,args.budget,args.seed) for x in ('uniform','random','adaptive')]
    out={'schema':'NUMERIA.P73.DISCOVERY_ONLY.1','zero_count':args.n,
         'start_t':args.from_t,'end_t':args.to,'evaluations_per_strategy':args.budget,
         'strategies':r,'disclaimer':'Floating-point exploration only; no randomized result certifies a continuous interval.'}
    content=json.dumps(out,indent=2,ensure_ascii=False)+'\n'
    if args.out:args.out.write_text(content,encoding='utf-8')
    print(content)
if __name__=='__main__':main()
