#!/usr/bin/env python3
"""Testes racionais independentes para OS ENUNCIADOS Lean V13; nao compilam Lean.
Confere largura fina por extremos, clausulas de inclusao, reciprocal e
tick do B5, usando Fraction e arredondamento dirigido. Deterministico.
"""
import argparse
import json
import random
from fractions import Fraction as F
from pathlib import Path
S=10**42

def bounds(a,b,c,d):
    x=(a*c,a*d,b*c,b*d)
    return min(x),max(x)

def test(trials,seed):
    rng=random.Random(seed)
    count={'mul_contains':0,'width_mul_le':0,'width_mul_abs_bound':0,
           'directed_singleton':0,'directed_monotonicity':0,
           'recip_pos':0,'recip_neg':0,'zero_guard':0}
    for i in range(trials):
        vals=[F(rng.randint(-20000,20000),rng.randint(1,20000)) for _ in range(4)]
        a,b=sorted(vals[:2]);c,d=sorted(vals[2:]); v,w=bounds(a,b,c,d)
        x=(a+b)/2;y=(c+d)/2
        assert v<=x*y<=w
        count['mul_contains']+=1
        maxa=max(abs(a),abs(b));maxb=max(abs(c),abs(d))
        assert w-v<=(b-a)*maxb+(d-c)*maxa
        count['width_mul_le']+=1
        assert w-v<=2*maxa*maxb
        count['width_mul_abs_bound']+=1
        for z in (a,b,c,d,x,y):
            dn=F((z*S).numerator//(z*S).denominator,S)
            up=F(-((-z*S).numerator//(-z*S).denominator),S)
            assert dn<=z<=up
            assert up-dn<=F(1,S)
            count['directed_singleton']+=1
        z1,z2=sorted([x,y]);
        dn=lambda z: F((z*S).numerator//(z*S).denominator,S)
        up=lambda z: F(-((-z*S).numerator//(-z*S).denominator),S)
        assert dn(z1)<=dn(z2) and up(z1)<=up(z2)
        count['directed_monotonicity']+=1
        if a>0:
            assert F(1,b)<=F(1,x)<=F(1,a);count['recip_pos']+=1
        if b<0:
            assert F(1,b)<=F(1,x)<=F(1,a);count['recip_neg']+=1
        if a<=0<=b: count['zero_guard']+=1
    return {'schema':'NUMERIA.E8K.V13.LEAN_STATEMENT_SANITY.1',
            'trials':trials,'seed':seed,'results':count,
            'status':'RATIONAL_SANITY_PASSED_NOT_LEAN_COMPILED',
            'not_certified':['Lean 4 typechecking e semantica exata B5 pendentes']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=5000)
    p.add_argument('--seed',type=int,default=20261028);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args(); assert args.trials>0
    result=test(args.trials,args.seed)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('OK:',result['results'])
if __name__=='__main__':main()
