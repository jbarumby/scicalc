#!/usr/bin/env python3
"""Enclosures fixed-point directed, no Arb, for real q_ij(r).

All intervals are pairs of integers representing multiples of 10**-digits.
Outward rounding at every multiplication/division preserves inclusion.
The special-function identities and phase convention are those of Fases 3/4.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import isqrt, factorial
from functools import lru_cache
from fase4_nos_sem_arb import pi_interval, psi_quarter_interval

class FixedRationalNodes:
    def __init__(self,digits=40,N=64,L=Q(9,20)):
        if not 25 <= digits <= 90: raise ValueError('digits 25..90')
        if N < 2: raise ValueError('N >= 2')
        self.S=10**digits; self.digits=digits; self.N=N; self.L=Q(L)
        self.one=(self.S,self.S)
        self.zero=(0,0)
        self.pi=self.q(pi_interval()[0],pi_interval()[1])
        self.rootL=self.sqrt_q(self.L)
        self.lnpi=self.ln_pos(self.pi,120)
        self.psi0=self.q(*psi_quarter_interval())
        self.w={k:self.div(self.scale(self.pi,k),self.q(2*self.L)) for k in range(1,11)}
        self.cache={}
        self.sinc_cache={}
        self.omega_cache={}
    def floor(self,n,d): return n//d
    def ceil(self,n,d): return -((-n)//d)
    def q(self,a,b=None):
        a=Q(a);b=a if b is None else Q(b)
        if a>b:raise ValueError('interval inverted')
        s=self.S
        return (a.numerator*s//a.denominator, self.ceil(b.numerator*s,b.denominator))
    def toQ(self,x):return Q(x,self.S)
    def add(self,a,b):return a[0]+b[0],a[1]+b[1]
    def neg(self,a):return -a[1],-a[0]
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        p=[u*v for u in a for v in b];s=self.S
        return min(p)//s,self.ceil(max(p),s)
    def scale(self,a,k):return self.mul(a,self.q(k))
    def inv_pos(self,b):
        if b[0] <= 0:raise ZeroDivisionError('positive interval required')
        s2=self.S*self.S
        return s2//b[1],self.ceil(s2,b[0])
    def div(self,a,b):
        if b[0]>0:return self.mul(a,self.inv_pos(b))
        if b[1]<0:return self.neg(self.div(a,self.neg(b)))
        raise ZeroDivisionError('interval crosses zero')
    def sqrt_q(self,v):
        v=Q(v)
        if v<=0:raise ValueError
        n=v.numerator*self.S*self.S//v.denominator
        root=isqrt(n)
        return root,root+1
    def ln_pos(self, x, terms=100):
        """ln(x)=2(sum_{k<terms} t^(2k+1)/(2k+1)+positive remainder).
        For x>=1, t=(x-1)/(x+1) lies in [0,1).
        """
        if x[0]<self.S:raise ValueError('ln_pos needs x >=1')
        t=self.div(self.sub(x,self.one),self.add(x,self.one))
        if t[0] < 0 or t[1]>=self.S:raise ValueError('bad atanh input')
        square=self.mul(t,t)
        power=t
        ans=self.zero
        for k in range(terms):
            ans=self.add(ans,self.div(power,self.q(2*k+1)))
            power=self.mul(power,square)
        # remaining positive series bounded by 2*t^(2n+1)/((2n+1)*(1-t²))
        res=self.div(self.scale(power,2),self.mul(self.q(2*terms+1),self.sub(self.one,square)))
        value=self.scale(ans,2)
        return value[0],value[1]+res[1]
    def sin_poly(self,y):
        """Sin(mid) via 38-term Taylor, with one-ulp remainder for |y|<=4.
        Bounds include directed rounding of each recurrence.
        """
        if abs(y) > 4*self.S:raise ValueError('Taylor sin requires |y|<=4')
        term=(y,y); total=term
        yy=self.mul(term,term)
        for k in range(1,38):
            term=self.neg(self.div(self.mul(term,yy),self.q((2*k)*(2*k+1))))
            total=self.add(total,term)
        # Taylor-Lagrange at degree 76 (degree 76 coefficient is 0):
        # |R| <= 4^77/77!. Include rounding at the fixed-point grid.
        residual=self.ceil(4**77*self.S,factorial(77))
        return total[0]-residual,total[1]+residual
    def sinc_poly(self,y):
        if abs(y)>self.S:raise ValueError('sinc Taylor requires |y|<=1')
        term=self.one;total=term; yy=self.mul((y,y),(y,y))
        for k in range(1,30):
            term=self.neg(self.div(self.mul(term,yy),self.q((2*k)*(2*k+1))))
            total=self.add(total,term)
        # Alternating sinc Taylor tail for |y|<=1: <=1/61!.
        residual=self.ceil(self.S,factorial(61))
        return total[0]-residual,total[1]+residual
    def sinc(self,x):
        if min(x)>=-self.S and max(x)<=self.S:
            mid=(x[0]+x[1])//2; rad=max(mid-x[0],x[1]-mid)
            val=self.sinc_poly(mid)
            half=self.ceil(rad,2)
            return val[0]-half,val[1]+half
        if x[0]<=0<=x[1]: raise ValueError('crosses zero in high-range sinc')
        mid=(x[0]+x[1])//2
        two_pi=self.scale(self.pi,2)
        pi_mid=(two_pi[0]+two_pi[1])//2
        k=(2*mid+pi_mid)//(2*pi_mid)
        residual=self.sub(x,self.scale(two_pi,k))
        if max(abs(residual[0]),abs(residual[1]))>4*self.S:
            raise ValueError('bad pi reduction')
        y=(residual[0]+residual[1])//2
        rad=max(y-residual[0],residual[1]-y)
        sin=self.sin_poly(y)
        sin=(sin[0]-rad,sin[1]+rad)  # |sin'|<=1
        return self.div(sin,x)
    def amp(self,r,k):
        key=(r,k)
        if key in self.cache:return self.cache[key]
        rv=self.q(r)
        L=self.q(self.L)
        sm=self.sinc(self.mul(L,self.sub(rv,self.w[k])))
        sp=self.sinc(self.mul(L,self.add(rv,self.w[k])))
        combo=self.add(sm,sp) if k%2 else self.sub(sm,sp)
        v=self.mul(self.rootL,combo)
        self.cache[key]=v
        return v
    def omega(self,r):
        if r in self.omega_cache:return self.omega_cache[r]
        r=Q(r)
        bn=r.numerator**2;bd=4*r.denominator**2
        # exact rational positive series using the identity in the source paper
        s0=0;s1=0;S=self.S
        for n in range(self.N):
            v=4*n+1
            num=64*bn
            den=v*(v*v*bd+16*bn)
            t=num*S//den
            s0+=t;s1+=t+(1 if num*S%den else 0)
        sums=(s0,s1)
        b2=self.q(Q(bn,bd))
        u=self.q(Q(self.N)+Q(1,4)); u2=self.mul(u,u)
        du=self.add(u2,b2)
        logarg=self.add(self.one,self.div(b2,u2))
        ilog=self.div(self.ln_pos(logarg,90),self.q(2))
        g=self.div(b2,self.mul(u,du))
        B2=self.div(self.mul(b2,self.add(self.scale(u2,3),b2)),self.mul(self.q(12),self.mul(u2,self.mul(du,du))))
        u4=self.mul(u2,u2); b4=self.mul(b2,b2)
        numerator=self.add(self.sub(u4,self.scale(self.mul(u2,b2),6)),b4)
        quotient=self.div(numerator,self.mul(self.mul(du,du),self.mul(du,du)))
        B4=self.div(self.sub(self.div(self.one,u4),quotient),self.q(120))
        # Analytically B4 >= 0 from phase 3; choose nonnegative upper bound.
        B4hi=max(0,B4[1])
        two=(self.add(ilog,self.div(g,self.q(2))))
        tail=self.add(two,B2)
        tail=(tail[0]-B4hi,tail[1])
        out=self.sub(self.add(self.psi0,self.add(sums,tail)),self.lnpi)
        self.omega_cache[r]=out
        return out
    def qij(self,r,i,j):
        if (i-j)%2:return self.zero
        phase=-1 if (abs(i-j)//2)%2 else 1
        v=self.mul(self.mul(self.amp(r,i),self.amp(r,j)),self.omega(r))
        v=self.div(v,self.pi)
        return v if phase==1 else self.neg(v)

