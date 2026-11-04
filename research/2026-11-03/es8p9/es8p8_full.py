#!/usr/bin/env python3
"""ES8-P8: independently enclosed entire historical 1600x16x30 nodal grid.

Certifies pairwise interval intersections, NOT equality and NOT the true-value
inclusion of either evaluator independently of their stated analytic lemmas.

Independent evaluator uses 10^-35 fixed-point outward rounding, 128 positive
terms of digamma and the elementary integral/trapezoid tail. The historical
source is retained verbatim and pinned by SHA-256, with 10^-40 precision and 64
terms plus Euler correction. No float is used to compute interval bounds.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys, time
from fractions import Fraction as Q
from pathlib import Path
from math import factorial, isqrt
from functools import lru_cache
from es8p7_independent import pi as independent_pi, omega0 as independent_omega0
from fase4_fixedpoint_sem_arb import FixedRationalNodes

L=Q(9,20)
PAIRS=tuple((i,j) for i in range(1,11) for j in range(i,11) if (j-i)%2==0)
HIST_SHA={
 'fase4_fixedpoint_sem_arb.py':'0935bb988f7a7bd4c83b36a4a873a71b1612eda00c77ac9bd6ba7daeef600bda',
 'fase4_nos_sem_arb.py':'151ec9a0a6c249f443bcbac278e8f004e02d268ef7590f81e78e5280872cbe49'
}
assert len(PAIRS)==30

class IndependentFixed:
    def __init__(self,digits=35,N=128):
        if N<2:raise ValueError('N>=2')
        self.S=10**digits;self.N=N;self.one=(self.S,self.S)
        self.zero=(0,0)
        self.pi=self.enclose(*independent_pi())
        self.omega_at_zero=self.enclose(*independent_omega0())
        self.rootL=self.sqrtq(L)
        self.w={k:self.div(self.scale(self.pi,k),self.num(2*L)) for k in range(1,11)}
        self.cache={};self.omega_cache={}
        self.sin_err=(4**(2*32+1)*self.S+factorial(2*32+1)-1)//factorial(2*32+1)
        self.sinc_err=(self.S+factorial(2*28+1)-1)//factorial(2*28+1)
    def enclose(self,a,b):
        a=Q(a);b=Q(b);S=self.S
        return a.numerator*S//a.denominator,-((-b.numerator*S)//b.denominator)
    def num(self,a):return self.enclose(a,a)
    def add(self,a,b):return a[0]+b[0],a[1]+b[1]
    def neg(self,a):return -a[1],-a[0]
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        S=self.S;p=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]);lo=min(p);hi=max(p)
        return lo//S,-((-hi)//S)
    def scale(self,a,b):return self.mul(a,self.num(b))
    def invpos(self,a):
        if a[0]<=0:raise ValueError('nonpositive interval')
        T=self.S*self.S
        return T//a[1],-((-T)//a[0])
    def div(self,a,b):
        if b[0]>0:return self.mul(a,self.invpos(b))
        if b[1]<0:return self.neg(self.div(a,self.neg(b)))
        raise ValueError('zero divisor')
    def sqrtq(self,x):
        x=Q(x);i=isqrt((x.numerator*self.S*self.S)//x.denominator)
        return i,i+1
    def sin_center(self,y):
        """|y|<=4, polynomial degree 63, remainder <=4^65/65!."""
        if abs(y)>4*self.S:raise ValueError('reduced argument outside [-4,4]')
        v=(y,y);sq=self.mul(v,v);term=v;acc=term
        for k in range(1,32):
            term=self.neg(self.div(self.mul(term,sq),self.num((2*k)*(2*k+1))))
            acc=self.add(acc,term)
        return acc[0]-self.sin_err,acc[1]+self.sin_err
    def sin(self,x):
        mid=(x[0]+x[1])//2
        pi_mid=(self.pi[0]+self.pi[1])//2
        two_pi_mid=2*pi_mid
        # Closest multiple; exact integer choice, with verification of the resulting enclosure.
        k=(mid+two_pi_mid//2)//two_pi_mid
        reduced=self.sub(x,self.scale(self.pi,2*k))
        c=(reduced[0]+reduced[1])//2
        rad=max(c-reduced[0],reduced[1]-c)
        val=self.sin_center(c)
        return val[0]-rad,val[1]+rad
    def sinc_center(self,y):
        if abs(y)>self.S:raise ValueError('sinc small argument outside [-1,1]')
        yy=self.mul((y,y),(y,y));term=self.one;acc=term
        for k in range(1,28):
            term=self.neg(self.div(self.mul(term,yy),self.num((2*k)*(2*k+1))))
            acc=self.add(acc,term)
        return acc[0]-self.sinc_err,acc[1]+self.sinc_err
    def sinc(self,x):
        if x[0]>=-self.S and x[1]<=self.S:
            c=(x[0]+x[1])//2;rad=max(c-x[0],x[1]-c)
            v=self.sinc_center(c)
            return v[0]-((rad+1)//2),v[1]+((rad+1)//2)
        if x[0]<=0<=x[1]:raise ValueError('large argument straddles zero')
        return self.div(self.sin(x),x)
    def amp(self,r,k):
        key=(r,k)
        if key in self.cache:return self.cache[key]
        rint=self.num(r)
        left=self.sinc(self.scale(self.sub(rint,self.w[k]),L))
        right=self.sinc(self.scale(self.add(rint,self.w[k]),L))
        comb=self.add(left,right) if k%2 else self.sub(left,right)
        out=self.mul(self.rootL,comb)
        self.cache[key]=out
        return out
    def _log_one_plus(self,b2,u):
        """ln(1+b²/u²)=2atanh(t), t=b²/(2u²+b²), t>=0. No library log."""
        x=b2/(2*u*u+b2)
        tx=self.enclose(x,x)
        tx2=self.mul(tx,tx)
        power=tx;acc=self.zero
        for j in range(45):
            acc=self.add(acc,self.div(power,self.num(2*j+1)))
            power=self.mul(power,tx2)
        value=self.scale(acc,2)
        remainder=self.div(self.scale(power,2),self.mul(self.num(91),self.sub(self.one,tx2)))
        return value[0],value[1]+remainder[1]
    def omega(self,r):
        if r in self.omega_cache:return self.omega_cache[r]
        r=Q(r);p=r.numerator;q=r.denominator;S=self.S
        if p==0:
            self.omega_cache[r]=self.omega_at_zero
            return self.omega_at_zero
        # Re psi(1/4+ir/2)-psi(1/4) = sum 16p²/[v(v²q²+4p²)], v=4n+1.
        s0=s1=0;num=16*p*p*S;pp=4*p*p;qq=q*q
        for n in range(self.N):
            v=4*n+1;den=v*(v*v*qq+pp)
            c,rem=divmod(num,den)
            s0+=c;s1+=c+bool(rem)
        b2=r*r/4;u=Q(self.N)+Q(1,4)
        integral=self.div(self._log_one_plus(b2,u),self.num(2))
        g=b2/(u*(u*u+b2))
        correction=b2*(3*u*u+b2)/(12*u*u*(u*u+b2)**2)
        # Integral + g/2 <= sum tail <= integral+g/2-corrected g'/12
        tail=self.add(integral,self.enclose(g/2,g/2+correction))
        out=self.add(self.omega_at_zero,self.add((s0,s1),tail))
        self.omega_cache[r]=out
        return out
    def qij(self,r,i,j):
        if (j-i)%2:return self.zero
        sign=-1 if ((j-i)//2)%2 else 1
        v=self.div(self.mul(self.mul(self.amp(r,i),self.amp(r,j)),self.omega(r)),self.pi)
        return v if sign==1 else self.neg(v)

def full_grid():
    for cell in range(1600):
        for j in range(16):
            yield cell,j,Q(cell,10)+Q(j,150)

def digest_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def scan(start=0,end=1600,chunk=50,terms=128):
    if not 0<=start<end<=1600:raise ValueError('bad cell bounds')
    if chunk<=0:raise ValueError('chunk must positive')
    base=Path(__file__).resolve().parent
    hashes={name:digest_file(base/name) for name in HIST_SHA}
    if hashes!=HIST_SHA:raise ValueError('historical source not SHA-pinned')
    chunks=[];failures=[];compared=0;max_width_hist=0;max_width_ind=0;total_nodes=0
    started=time.monotonic()
    for first in range(start,end,chunk):
        last=min(end,first+chunk)
        reader=FixedRationalNodes(digits=40,N=64)
        other=IndependentFixed(digits=35,N=terms)
        md=hashlib.sha256();c_checks=0;bad=0;h_width=i_width=0
        for cell in range(first,last):
            for j in range(16):
                r=Q(cell,10)+Q(j,150)
                total_nodes+=1
                for i,k in PAIRS:
                    a=reader.qij(r,i,k)
                    b=other.qij(r,i,k)
                    # Fast exact comparison without making Fractions
                    overlap=(a[0]*other.S <= b[1]*reader.S and b[0]*reader.S <= a[1]*other.S)
                    c_checks+=1
                    if not overlap:
                        bad+=1
                        if len(failures)<30:
                            failures.append({'cell':cell,'node':j,'r':str(r),'mode_i':i,'mode_j':k,
                                             'historical':[str(x) for x in a],
                                             'independent':[str(x) for x in b]})
                    h_width=max(h_width,a[1]-a[0]);i_width=max(i_width,b[1]-b[0])
                    # Canonical chain hash of both entire compared endpoint pairs.
                    md.update(f'{cell},{j},{i},{k}:{a[0]},{a[1]}:{b[0]},{b[1]}\n'.encode())
        compared+=c_checks;max_width_hist=max(max_width_hist,h_width);max_width_ind=max(max_width_ind,i_width)
        chunks.append({'start_cell':first,'end_cell_exclusive':last,'nodes':(last-first)*16,
                       'comparisons':c_checks,'incompatible':bad,'sha256_pairwise_intervals':md.hexdigest(),
                       'max_width_historical_grid_units_1e_minus_40':str(h_width),
                       'max_width_independent_grid_units_1e_minus_35':str(i_width)})
        print(f'cells [{first},{last}), {c_checks} intervals, failures={bad}, time={time.monotonic()-started:.1f}s',flush=True)
    result={'schema':'NUMERIA.ES8P8.EXHAUSTIVE_NODE_CROSSCHECK.1',
            'status':'FULL_GRID_PAIRWISE_INTERVAL_INTERSECTIONS_CHECKED' if not failures and start==0 and end==1600 else 'PARTIAL_OR_MISMATCH',
            'scope':{'start_cell':start,'end_cell_exclusive':end,'cells_total':1600,'nodes_per_cell':16,'mode_pairs':30,
                     'full_expected_comparisons':768000,'historical_digits':40,'independent_digits':35,
                     'independent_digamma_terms':terms},
            'historical_source_sha256':hashes,'count_positions':total_nodes,
            'count_compared':compared,'failures_total':sum(x['incompatible'] for x in chunks),
            'failure_examples':failures,'chunks':chunks,
            'max_width_historical':str(Q(max_width_hist,10**40)),
            'max_width_independent':str(Q(max_width_ind,10**35)),
            'limitations':[
               'Interval intersection is compatibility, not identity of computed enclosures',
               'This audit assumes analytic series and tail lemmas; the code is not a Lean proof',
               'It checks pointwise integrands, not quadrature integration or universal Weil positivity',
               'Representations of Boyadzhiev and Kanemitsu not implemented as third evaluator']}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--end',type=int,default=1600)
    p.add_argument('--chunk',type=int,default=50)
    p.add_argument('--terms',type=int,default=128)
    p.add_argument('--out',type=Path,default=Path('RESULTADOS_ES8P8.json'))
    a=p.parse_args();out=scan(a.start,a.end,a.chunk,a.terms)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f"STATUS {out['status']}: {out['count_compared']} checks, {out['failures_total']} incompatibilities")
    return 0 if out['failures_total']==0 else 1
if __name__=='__main__':raise SystemExit(main())
