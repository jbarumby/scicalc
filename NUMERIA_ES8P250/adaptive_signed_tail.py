#!/usr/bin/env python3
"""P250: adaptive fixed-point certificate for the signed-prime tail on L=161/200.

Only integer-directed arithmetic and declared classical digamma/Binet identities.
Independent of P249 Arb 24960-cell search; uses pinned P245B library.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / 'numeria/p245/package/p245b_certified_gl60.py'
EXPECTED_LIBRARY_SHA256 = '29a6ad24fe6ce073531eccd033a8ac9690010f5d76cc304377c3f6bbc5e544b8'
L=F(161,200)
BETA_CUTOFF=504
CERT_CUTOFF=309
NEGATIVE_WITNESS=F(6169,20)


def load_p245b():
    if hashlib.sha256(LIBRARY.read_bytes()).hexdigest() != EXPECTED_LIBRARY_SHA256:
        raise ArithmeticError('P245B fixedpoint source not pinned')
    spec = importlib.util.spec_from_file_location('p250_p245b_pinned',LIBRARY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.PREC != 256:
        raise ArithmeticError('unsupported interval precision')
    return module


class SignedSymbol:
    def __init__(self, arith):
        g=arith
        self.g=g
        self.logs=[g.LN2,g.log_interval(3),2*g.LN2,g.log_interval(5)]
        self.coeff=[2*g.LN2/g.iv(2).sqrt(),2*self.logs[1]/g.iv(3).sqrt(),g.LN2,2*self.logs[3]/g.iv(5).sqrt()]
        self.A=sum(self.coeff,g.iv(0))
        self.M2=sum((a*omega*omega for a,omega in zip(self.coeff,self.logs)),g.iv(0))
        self.beta=g.log_interval(g.iv(BETA_CUTOFF)/(2*g.PI))-g.iv(F(1,BETA_CUTOFF))-self.A
        self.log_pi=g.log_interval(g.PI)

    def symbol_and_derivative(self,t):
        g=self.g
        t=g.iv(t)
        sym=g.iv(0)
        der=g.iv(0)
        for a,w in zip(self.coeff,self.logs):
            sin,cos=g.sin_cos(t*w)
            sym=sym+a*cos
            der=der-a*w*sin
        return sym,der

    def gap_lower(self,left,right):
        """Enclose B(t)-beta-S(t) uniformly on [left,right]."""
        g=self.g
        left,right=F(left),F(right)
        if not 0<left<right<=BETA_CUTOFF:
            raise ValueError('invalid interval for the P250 analytic lower bound')
        c=(left+right)/2
        radius=(right-left)/2
        s,d=self.symbol_and_derivative(c)
        # B(t) = log(t/(2pi)) - 1/t is strictly increasing for t>0.
        Bleft=g.log_interval(g.iv(left)/(2*g.PI))-g.iv(1)/g.iv(left)
        # S(t) <= S(c)+r|S'(c)|+r^2/2 * sum(a_q omega_q^2).
        err=(g.iv(radius)*g.iv(F(d.abs_upper(),g.S)) +
             g.iv(radius*radius)*g.iv(F(self.M2.u,g.S))/2)
        return Bleft-self.beta-s-err

    def true_symbol_at(self,t):
        """Digamma-based Phi-beta, with P245B explicit Binet truncation bounds."""
        g=self.g
        t=F(t)
        s,_=self.symbol_and_derivative(t)
        return g.re_digamma(g.iv(t)/2)-self.log_pi-s-self.beta


def adaptive_search(symbol,start=CERT_CUTOFF,end=BETA_CUTOFF,max_depth=18):
    g=symbol.g
    leaves=[]
    attempts=0
    max_depth_observed=0
    minimum_lower=None
    def recurse(a,b,depth):
        nonlocal attempts,max_depth_observed,minimum_lower
        attempts+=1
        bound=symbol.gap_lower(a,b)
        if bound.l>0:
            leaves.append((a,b,bound.l))
            max_depth_observed=max(max_depth_observed,depth)
            minimum_lower=bound.l if minimum_lower is None else min(minimum_lower,bound.l)
            return
        if depth>=max_depth:
            raise ArithmeticError(f'cannot prove strict symbol positivity on [{a},{b}] at depth {depth}')
        middle=(a+b)/2
        recurse(a,middle,depth+1)
        recurse(middle,b,depth+1)
    for k in range(start,end):
        recurse(F(k),F(k+1),0)
    if not leaves or leaves[0][0]!=start or leaves[-1][1]!=end:
        raise ArithmeticError('coverage endpoints missing')
    for prior,current in zip(leaves,leaves[1:]):
        if prior[1]!=current[0]:raise ArithmeticError('coverage gap')
    return {'cells':len(leaves),'attempts':attempts,
            'max_depth':max_depth_observed,
            'min_lower_dyadic_numerator':str(minimum_lower),
            'interval_scale_power_of_two':g.PREC,
            'leaves':[[str(a),str(b),str(lo)] for a,b,lo in leaves]}


def run():
    g=load_p245b()
    s=SignedSymbol(g)
    witness=s.true_symbol_at(NEGATIVE_WITNESS)
    if witness.u>=0:
        raise ArithmeticError('negative symbol witness failed')
    edge=s.true_symbol_at(F(CERT_CUTOFF))
    if edge.l<=0:
        raise ArithmeticError('boundary point positivity failed')
    result=adaptive_search(s)
    # Finite interval is covered; for t>504, B(t)>B(504) and S(t)<=A.
    return {'schema':'NUMERIA.ES8P250.ADAPTIVE_SIGNED_TAIL.1',
            'L':'161/200','beta_frozen_at':BETA_CUTOFF,
            'pin_p245b_sha256':EXPECTED_LIBRARY_SHA256,
            'arithmetic':'outward integer fixedpoint 256-bit, exact dyadic partitions',
            'phase_formula':'S(t)=sum(2 log p / sqrt(q))*cos(t log q), q=2,3,4,5',
            'analytic_input':'Binet lower bound Re digamma(1/4+it/2)-log(pi)>=B(t) for t>=309',
            'certified_positive_frequency_range':[CERT_CUTOFF,'infinity'],
            'negative_exact_rational_frequency':str(NEGATIVE_WITNESS),
            'negative_actual_symbol_interval':witness.desc(),
            'positive_at_309_interval':edge.desc(),
            'smallest_integer_tail_cutoff_for_this_fixed_beta':CERT_CUTOFF,
            'count_cells':result['cells'],'count_attempts':result['attempts'],
            'max_depth':result['max_depth'],
            'min_lower_dyadic_numerator':result['min_lower_dyadic_numerator'],
            'interval_scale_power_of_two':g.PREC,
            'leaves':result['leaves'],
            'new_finite_block_certified':False,
            'RH_proved':False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    obj=run()
    a.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('P250 PASS; cells',obj['count_cells'],'attempts',obj['count_attempts'],'max_depth',obj['max_depth'])

if __name__=='__main__':main()
