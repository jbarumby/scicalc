#!/usr/bin/env python3
"""NUMERIA fase 4: enclosures RACIONAIS, SEM flint/Arb/mpmath.

Aplica-se a r REAL RACIONAL em poucos nós de prova, não à integração inteira.
Usa: arctan de Machin, ln via arctanh, Euler-Maclaurin para gamma,
série positiva e cauda integral para Re(psi(1/4+ir/2)),
sinc por Taylor + limite global |sinc'|<=1/2.
Todas as decisões e intervalos são Fraction exatos.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import factorial, isqrt
import argparse, json, math
from functools import lru_cache

I = tuple[Q, Q]

def interval(lo:Q, hi:Q)->I:
    if lo>hi: raise ValueError('limites invertidos')
    return (lo,hi)

def add(a:I,b:I)->I: return a[0]+b[0], a[1]+b[1]
def neg(a:I)->I: return -a[1],-a[0]
def sub(a:I,b:I)->I: return add(a,neg(b))
def mul(a:I,b:I)->I:
    terms=[x*y for x in a for y in b]
    return min(terms),max(terms)
def div_positive(a:I,b:I)->I:
    if b[0]<=0: raise ValueError('divisor intervalar deve ser positivo')
    return mul(a,(1/b[1],1/b[0]))
def scal(a:I,k:Q)->I: return mul(a,(k,k))
def sq(a:I)->I:
    if a[0]<=0<=a[1]: return Q(0),max(a[0]**2,a[1]**2)
    return min(a[0]**2,a[1]**2),max(a[0]**2,a[1]**2)

def atan_inv_interval(q:int, terms:int=32)->I:
    assert q>=2 and terms>=1
    x=Q(1,q)
    s=sum(((-1)**j * x**(2*j+1)/Q(2*j+1) for j in range(terms)),Q(0))
    nxt=x**(2*terms+1)/Q(2*terms+1)
    return (s,s+nxt) if terms%2==0 else (s-nxt,s)

@lru_cache(maxsize=None)
def pi_interval()->I:
    return sub(scal(atan_inv_interval(5),Q(16)),scal(atan_inv_interval(239),Q(4)))

def ln_interval(x:Q,terms:int=80)->I:
    """Intervalo exato para log(x), x>0, via arctanh de argumento racional."""
    x=Q(x)
    if x<=0:raise ValueError('log de x<=0')
    t=(x-1)/(x+1)
    s=sum((2*t**(2*k+1)/Q(2*k+1) for k in range(terms)),Q(0))
    rem=Q(2)*abs(t)**(2*terms+1)/(Q(2*terms+1)*(1-t*t))
    return s-rem,s+rem

@lru_cache(maxsize=None)
def ln2_interval()->I: return ln_interval(Q(2),50)

@lru_cache(maxsize=None)
def euler_gamma_interval(n:int=256)->I:
    """H_n-log(n)-1/(2n)+1/(12n^2)-1/(120n^4) <= gamma <= sem ultimo termo.
    Ver desigualdade de Euler-Maclaurin para H_n-log(n)-gamma.
    n DEVE ser potencia de 2 para usar log(n)=k log2 com série rápida.
    """
    if n<=0 or n&(n-1): raise ValueError('n deve ser potencia de 2')
    h=sum((Q(1,k) for k in range(1,n+1)),Q(0))
    ln_n=scal(ln2_interval(),Q(n.bit_length()-1))
    approx=sub((h,h),ln_n)
    approx=add(approx,( -Q(1,2*n)+Q(1,12*n*n),)*2)
    return approx[0]-Q(1,120*n**4),approx[1]

@lru_cache(maxsize=None)
def psi_quarter_interval()->I:
    """Identidade psi(1/4)=-gamma - pi/2 -3 log 2."""
    return sub(sub(neg(euler_gamma_interval()),scal(pi_interval(),Q(1,2))),scal(ln2_interval(),Q(3)))

@lru_cache(maxsize=None)
def delta_psi_quarter(r:Q,N:int=1024)->I:
    """Re psi(1/4+ir/2) - psi(1/4), com cauda acelerada rigorosa.

    g(u)=b²/[u(u²+b²)] é completamente monotona, pois
    g(u)=integral_0^inf exp(-ut)*(1-cos(bt)) dt.
    Para u=N+1/4 e I=integral_u^inf g(v)dv=log(1+b²/u²)/2,
    I+g(u)/2 <= sum_{n=N}^inf g(n+1/4)
                 <= I+g(u)/2 -g'(u)/12.
    Justificacao: 1/t+1/2 <= 1/(1-e^-t) <= 1/t+1/2+t/12.
    """
    b=Q(r,2); b2=b*b; a=Q(1,4)
    if N<1: raise ValueError('N deve ser >=1')
    if b2==0: return Q(0),Q(0)
    s=Q(0)
    for n in range(N):
        y=Q(n)+a
        s+=b2/(y*(y*y+b2))
    u=Q(N)+a
    integ=scal(ln_interval(1+b2/(u*u),terms=35),Q(1,2))
    g=b2/(u*(u*u+b2))
    deriv_upper=b2*(3*u*u+b2)/(12*u*u*(u*u+b2)**2)
    return (s+integ[0]+g/2, s+integ[1]+g/2+deriv_upper)

@lru_cache(maxsize=None)
def omega_interval(r:Q,N:int=1024)->I:
    lp=interval(*pi_interval())
    # log(pi) is monotonic, use independent interval arithmetic on endpoints.
    lnpi=(ln_interval(lp[0])[0],ln_interval(lp[1])[1])
    return sub(add(psi_quarter_interval(),delta_psi_quarter(r,N)),lnpi)

@lru_cache(maxsize=None)
def sqrt_rational_interval(q:Q,digits:int=48)->I:
    if q<0:raise ValueError('raiz de numero negativo')
    if q==0:return Q(0),Q(0)
    scale=10**digits
    n=q.numerator*scale*scale//q.denominator
    t=isqrt(n)
    lo=Q(t,scale); hi=Q(t+1,scale)
    return lo,hi

def sinc_rational_point(x:Q, terms:int=45)->I:
    """sinc(x)=sum_{k<terms}(-1)^k x^(2k)/(2k+1)! + erro <= |x|^(2terms)/(2terms+1)! ."""
    s=sum((Q((-1)**k)*x**(2*k)/factorial(2*k+1) for k in range(terms)),Q(0))
    err=abs(x)**(2*terms)/factorial(2*terms+1)
    return s-err,s+err

def sinc_interval(x:I,terms:int=45)->I:
    """Para x real, sinc(x)=∫_0^1 cos(tx)dt: |sinc'(x)|<=1/2."""
    mid=(x[0]+x[1])/2; rad=(x[1]-x[0])/2
    s=sinc_rational_point(mid,terms)
    return s[0]-rad/2,s[1]+rad/2

@lru_cache(maxsize=None)
def amplitude_interval(r:Q,k:int=1,L:Q=Q(9,20))->I:
    """A_k(r)=sqrt(L)[sinc(L(r-w_k))-(-1)^k sinc(L(r+w_k))]."""
    w=div_positive(scal(pi_interval(),Q(k)),scal((L,L),Q(2)))
    rm=(r,r)
    left=sinc_interval(scal(sub(rm,w),L))
    right=sinc_interval(scal(add(rm,w),L))
    combo=add(left,right) if k%2 else sub(left,right)
    return mul(sqrt_rational_interval(L),combo)

def integrand_11_interval(r:Q,N:int=1024,L:Q=Q(9,20))->I:
    """q_11(r)=A_1(r)^2 Omega(r)/pi; integrando, nao integral."""
    return div_positive(mul(sq(amplitude_interval(r,1,L)),omega_interval(r,N)),pi_interval())

def float_str(q:Q, digits:int=19)->str:
    from decimal import Decimal, localcontext
    with localcontext() as c:
        c.prec=digits
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def summary(a:I, digits:int=26)->dict:
    # Exportacao decimal dirigida: piso da ponta inferior e teto da superior
    # inteiramente por divisoes de inteiros, sem Decimal/float.
    scale=10**digits
    lo=(a[0].numerator*scale)//a[0].denominator
    hi=-((-a[1].numerator*scale)//a[1].denominator)
    outward=(Q(lo,scale),Q(hi,scale))
    assert outward[0] <= a[0] <= a[1] <= outward[1]
    return {'lower':f'{float_str(outward[0],digits+5)}',
            'upper':f'{float_str(outward[1],digits+5)}',
            'width_upper':float_str(outward[1]-outward[0],18),
            'decimal_places':digits,'outward_inclusion_checked':True}


def main()->int:
    p=argparse.ArgumentParser(description='Fase 4: sinc, digamma e integrando G11 sem Arb')
    p.add_argument('--r',nargs='*',default=['0','1','4','10'])
    p.add_argument('--series',type=int,default=1024,help='Numero de termos positivos da digamma')
    p.add_argument('--out',default='FASE4_NOS_INDEPENDENTES.json')
    args=p.parse_args()
    res={'metodo':'arithmetic_fraction_no_arb','status':'intervals_constructed_not_full_integral',
         'series_terms':args.series,'nodes':[]}
    for rv in args.r:
        r=Q(rv); omega=omega_interval(r,args.series); a=amplitude_interval(r)
        g=integrand_11_interval(r,args.series)
        rec={'r':str(r),'omega':summary(omega),'A_1':summary(a),'q_11':summary(g)}
        print(f'r={str(r):<5} Omega=({float_str(omega[0],13)}, {float_str(omega[1],13)}) '
              f'width={float_str(omega[1]-omega[0],5)}')
        res['nodes'].append(rec)
    from pathlib import Path
    Path(args.out).write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('ARQUIVO:',args.out)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
