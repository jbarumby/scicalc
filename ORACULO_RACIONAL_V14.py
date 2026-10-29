#!/usr/bin/env python3
"""NUMERIA E8-K V14 — oráculo racional auxiliar para a cota global da cauda da digama e oraculo racional.

Requer o E7 ZIP fixado por SHA-256. Executa a fonte B5 localmente (NÃO e sandbox
 de seguranca). O oraculo de omega é independente das rotinas B5: usa Machin
com mais termos, gamma EM em N=128, log via atanh para constantes e serie
alternada de log(1+t) para o termo integral da cauda.

Alcance: prova analitica da desigualdade da cauda sob teoremas declarados,
comparacoes racionais exatas para nos amostrados e conformidade dos registros
B5 selecionados. NAO prova todos os caminhos nem a formula global de Weil.
"""
from __future__ import annotations
import argparse, ast, hashlib, json, math, platform
from fractions import Fraction as F
from pathlib import Path
from zipfile import ZipFile

ZIP_SHA='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
SRC_SHA='a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
SOURCE='fontes/AUDITOR_PONTOS_B5.py'
S=10**42
T=10**65

def sha(x: bytes): return hashlib.sha256(x).hexdigest()
def lower(x,scale=T):return F(x.numerator*scale//x.denominator,scale)
def upper(x,scale=T):return F(-((-x.numerator*scale)//x.denominator),scale)
def bounds(iv): return F(iv.a,S),F(iv.b,S)
def contains(outer,inner): return outer[0]<=inner[0] and inner[1]<=outer[1]
def display(x):
    # Approximate display only; comparisons use exact Fraction.
    return f'{float(x):.12e}'

def atan_bounds_pos(x, n):
    assert 0<x<1
    p=x; z=F(0)
    for j in range(n):
        z += (p if j%2==0 else -p)/F(2*j+1)
        p *= x*x
    rem=p/F(2*n+1)
    return (z,z+rem) if n%2==0 else (z-rem,z)

def pi_reference():
    a=atan_bounds_pos(F(1,5),110)
    b=atan_bounds_pos(F(1,239),55)
    lo=16*a[0]-4*b[1];hi=16*a[1]-4*b[0]
    return lower(lo),upper(hi)

def atanh_log1p_reference(t,n=120,coarse=True):
    """Bounds for log(1+t), t>=0. Series in z=t/(2+t)."""
    assert t>=0
    z=t/(2+t)
    p=z; total=F(0)
    for j in range(n):
        total+=2*p/F(2*j+1)
        p*=z*z
    rem=2*p/F((2*n+1)) /(1-z*z)
    return (lower(total),upper(total+rem)) if coarse else (total,total+rem)

def log1p_alternating_reference(t, n=160):
    """Independent alternating log series for 0<=t<1; rational enclosure."""
    assert F(0)<=t<F(1)
    p=t; z=F(0)
    for j in range(1,n+1):
        z+= (p if j%2 else -p)/F(j)
        p*=t
    rem=p/F(n+1)
    lo,hi=(z,z+rem) if n%2==0 else (z-rem,z)
    return lower(lo),upper(hi)

def ln2_reference():
    return atanh_log1p_reference(F(1),130,coarse=False)

def logpi_reference(pi,ln2):
    # pi in (2,4): log(pi)=log2 + log(pi/2) using monotonicity.
    assert 2<pi[0]<pi[1]<4
    alo=atanh_log1p_reference(pi[0]/2-1,120)
    ahi=atanh_log1p_reference(pi[1]/2-1,120)
    return lower(alo[0]+ln2[0]),upper(ahi[1]+ln2[1])

_B=[F(1,6),-F(1,30),F(1,42),-F(1,30),F(5,66),-F(691,2730),F(7,6),-F(3617,510)]
def gamma_reference(ln2, N=128):
    H=sum((F(1,k) for k in range(1,N+1)),F(0))
    corr=[_B[j]/F(2*(j+1)*N**(2*j+2)) for j in range(8)]
    upper_g=H-F(1,2*N)+sum(corr[:7])  # no log yet
    lower_g=upper_g+corr[7]
    # log N = 7 log 2 for N=128
    assert N==128
    return lower(lower_g-7*ln2[1]),upper(upper_g-7*ln2[0])

def omega_zero_reference(gamma,pi,ln2,logpi):
    return (lower(-gamma[1]-pi[1]/2-3*ln2[1]-logpi[1]),
            upper(-gamma[0]-pi[0]/2-3*ln2[0]-logpi[0]))

def omega_reference(r, omega0, M=128):
    """Exact finite digamma terms + exact analytic tail bracket at rational r."""
    assert r>=0
    nr,dr=r.numerator,r.denominator
    total=F(0)
    for j in range(M):
        u=4*j+1
        total+=F(16*nr*nr,u*(u*u*dr*dr+4*nr*nr))
    u=F(4*M+1,4)
    t=r*r/(4*u*u)
    # all selected points satisfy t<1 (r<=160)
    loglo,loghi=log1p_alternating_reference(t)
    extra=r*r/(4*u*(u*u+r*r/4))
    return lower(omega0[0]+total+loglo/2),upper(omega0[1]+total+loghi/2+extra)

def multiply_intervals(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def divide_pos(a,b):
    assert b[0]>0
    return multiply_intervals(a,(1/b[1],1/b[0]))


# Motor aritmetico independente do B5: escala, implementacao e ordens distintas.
# Todas as operacoes usam inteiros e arredondamento dirigido; sem mpmath/float.
U = 10**62
class Iv:
    __slots__=('lo','hi')
    def __init__(self,lo,hi):
        if lo>hi: raise ValueError('intervalo invertido')
        self.lo,self.hi=lo,hi
    @staticmethod
    def q(x):
        x=F(x)
        return Iv((x.numerator*U)//x.denominator, -((-x.numerator*U)//x.denominator))
    def bounds(self): return F(self.lo,U),F(self.hi,U)
    def __add__(self,b):
        b=co(b); return Iv(self.lo+b.lo,self.hi+b.hi)
    def __neg__(self): return Iv(-self.hi,-self.lo)
    def __sub__(self,b):return self+-co(b)
    def __mul__(self,b):
        b=co(b); vals=(self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi)
        return Iv(min(vals)//U,-((-max(vals))//U))
    def sq(self):
        if self.lo<=0<=self.hi: return Iv(0,-((-max(self.lo**2,self.hi**2))//U))
        return Iv(min(self.lo**2,self.hi**2)//U,-((-max(self.lo**2,self.hi**2))//U))
    def recip(self):
        if self.lo<=0<=self.hi: raise ZeroDivisionError('intervalo contem zero')
        if self.hi<0:return -(-self).recip()
        return Iv((U*U)//self.hi,-((-U*U)//self.lo))
    def __truediv__(self,b): return self*co(b).recip()
    def absmax(self):return F(max(abs(self.lo),abs(self.hi)),U)
def co(x):return x if isinstance(x,Iv) else Iv.q(x)
def enclosing(interval): return Iv.q(interval[0]) + Iv(0,0) if interval[0]==interval[1] else Iv(Iv.q(interval[0]).lo,Iv.q(interval[1]).hi)
def expanded(v,eps):
    return Iv(v.lo-Iv.q(eps).hi,v.hi+Iv.q(eps).hi)

def trig_reference(x,kind,pi,terms=43):
    if kind not in ('sin','cos'):raise ValueError(kind)
    p0,p1=pi.bounds()
    mid=(p0+p1)/2
    k=round(float(x.bounds()[0]/(2*mid)))
    a=x-Iv.q(2*k)*pi
    if not a.absmax()<4:raise AssertionError('reducao nao comprovada')
    xx=-a.sq()
    v=Iv.q(1)
    if kind=='cos':
        for j in range(terms,0,-1):v=Iv.q(1)+xx*v/Iv.q((2*j)*(2*j-1))
        rem=F(4**(2*terms+2), math.factorial(2*terms+2))
    else:
        for j in range(terms,0,-1):v=Iv.q(1)+xx*v/Iv.q((2*j)*(2*j+1))
        v=a*v
        rem=F(4**(2*terms+3),math.factorial(2*terms+3))
    return expanded(v,rem)

def sinc_reference(x,terms=43):
    if not x.absmax()<2: raise ValueError('sinc requer |x|<2')
    xx=-x.sq();v=Iv.q(1)
    for j in range(terms,0,-1):v=Iv.q(1)+xx*v/Iv.q((2*j)*(2*j+1))
    rem=F(3**(2*terms+2),math.factorial(2*terms+3))
    return expanded(v,rem)

def amp_reference(r,k,pi):
    L=F(9,20);w=Iv.q(k)*pi/Iv.q(2*L)
    x=Iv.q(L*r)-Iv.q(k)*pi/Iv.q(2)
    if x.absmax()<2:
        return Iv.q(4*L)*w.sq()*sinc_reference(x).sq()/(Iv.q(r)+w).sq()
    g=trig_reference(Iv.q(L*r),'cos' if k%2 else 'sin',pi)
    d=Iv.q(r*r)-w.sq()
    if d.lo<=0<=d.hi: raise AssertionError('polo aparente nao removido')
    return Iv.q(k*k)*pi.sq()*g.sq()/(Iv.q(L**3)*d.sq())

def make_constants():
    pi=pi_reference(); ln2=ln2_reference(); lp=logpi_reference(pi,ln2)
    ga=gamma_reference(ln2);o0=omega_zero_reference(ga,pi,ln2,lp)
    return {'pi':pi,'ln2':ln2,'logpi':lp,'gamma':ga,'omega0':o0,'pi_iv':enclosing(pi)}

def evaluate_independent(r,k,const):
    rr=F(r)
    om=omega_reference(rr,const['omega0'])
    ai=amp_reference(rr,k,const['pi_iv'])
    result=ai*enclosing(om)/const['pi_iv']
    return result.bounds(),om
