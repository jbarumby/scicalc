"""P92-R: rational verification of the proposed independent Preissmann L5/L6 proofs.

Analytic arguments are in README.md. Finite Fraction checks are regression only.
"""
from fractions import Fraction as F
from random import Random
from pathlib import Path
import json, argparse, hashlib

SCHEMA = "NUMERIA.P92R.LEMMA5.BELLMAN.1"


def harmonic(n):
    return sum((F(1, k) for k in range(1, n+1)), F(0))


def nearest(gaps):
    g = tuple(F(x) for x in gaps)
    if not g or any(x <= 0 for x in g):
        raise ValueError("all gaps must be positive")
    x = []
    total = F(0)
    for d in g:
        total += d
        x.append(total)
    delta = [min(g[i], g[i+1]) if i+1<len(g) else g[i] for i in range(len(g))]
    return x, delta


def lemma5(gaps):
    x, delta = nearest(gaps)
    a, X = x[0], x[-1]
    k = int(X//a)
    S = sum((d/z for d,z in zip(delta,x)), F(0))
    mass = sum(delta, F(0))
    bound = harmonic(k) + (X-k*a)/X
    tighter = harmonic(k) + (mass-k*a)/X
    if not (S <= tighter <= bound and mass <= X):
        raise ArithmeticError("quantitative Lemma 5 regression failed")
    return {"sum":str(S), "bound":str(bound), "margin":str(bound-S),
            "number_of_points":len(x), "harmonic_index":k}


def ramp_U0(L,h):
    L,h=F(L),F(h)
    if L<0 or h<=0: raise ValueError("L >=0, h>0")
    if L==0:return F(0)
    n=int(L//h)+int(L%h>0)
    return h*n*L-h*h*n*(n-1)/2


def ramp_bellman(L,h,b):
    L,h,b=F(L),F(h),F(b)
    if L<0 or h<=0 or b<=0:raise ValueError()
    return ramp_U0(L,h)-(min(h,b)*L+ramp_U0(max(F(0),L-b),b))


def pi_lower():
    """Exact rational lower enclosure using Machin and alternating atan."""
    def lo(x,n):
        S=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
        return S if n%2==0 else S-x**(2*n+1)/F(2*n+1)
    def hi(x,n):
        S=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
        return S+x**(2*n+1)/F(2*n+1) if n%2==0 else S
    return 4*(4*lo(F(1,5),16)-hi(F(1,239),8))


def l6(points,s,t):
    z=[F(x) for x in points]
    if not 0<=s<t<len(z) or any(x>=y for x,y in zip(z,z[1:])):
        raise ValueError("need sorted points, 0 <= s < t < N")
    d=[min(abs(z[i]-z[j]) for j in range(len(z)) if i!=j)
       for i in range(len(z))]
    h=z[t]-z[s]
    left=sum((d[j]/((z[j]-z[s])**2*(z[j]-z[t])**2)
              for j in range(len(z)) if j!=s and j!=t),F(0))
    right=pi_lower()**2/(3*h*h)*(1/d[s]+1/d[t])-3*(d[s]+d[t])/h**4
    if left>right:raise ArithmeticError("L6 rational regression failed")
    return {"lhs":str(left), "rhs_pi_lower":str(right),"margin":str(right-left)}


def finite_compensator(gaps,A,L,R):
    g=[F(x) for x in gaps];L,R=F(L),F(R)
    if not g or A<1 or not 0<L<=g[0] or not 0<R<=g[-1]:
        raise ValueError("invalid mixed-index truncation")
    h=sum(g,F(0));Z=h+A*L;J=int(Z//R);Zp=h+J*R
    x=[-j*L for j in range(A,0,-1)]+[F(0)]
    p=F(0)
    for d in g:
        p+=d
        x.append(p)
    x += [h+j*R for j in range(1,J+1)]
    t,s=A,A+len(g)
    d=[min(x[i]-x[i-1] if i else x[1]-x[0],
           x[i+1]-x[i] if i+1<len(x) else x[-1]-x[-2])
       for i in range(len(x))]
    mixed=sum((d[j]*(1/(x[j]-x[t])-1/(x[j]-x[s]))
               for j in range(len(x)) if j not in (s,t)),F(0))
    remainder=g[-1]/Z+g[0]/Zp+2*h/(A*L)
    margin=-(L+R)/h+remainder-mixed
    if margin<0:raise ArithmeticError("mixed-index compensated bound failed")
    return {"A":A,"right_cutoff":J,"mixed":str(mixed),
            "remainder":str(remainder),"margin":str(margin)}


def verify():
    r=Random(101738)
    l5=[]
    for n in (1,2,3,4,5,7,12,25):
        for _ in range(26):
            g=[F(1)]+[F(r.randint(1,24),r.randint(1,16)) for j in range(n-1)]
            l5.append(lemma5(g))
    for n in range(1,15):
        l5.append(lemma5([F(1)]*n))
    bellman=[]
    for _ in range(400):
        margin=ramp_bellman(F(r.randint(0,100),r.randint(1,20)),
              F(r.randint(1,35),r.randint(1,15)),
              F(r.randint(1,35),r.randint(1,15)))
        if margin<0:raise ArithmeticError("Bellman ramp failure")
        bellman.append(str(margin))
    mix=[]
    for g in ([F(1)],[F(1),F(2)],[F(1),F(4),F(3),F(1)],
              [F(3,2),F(1,4),F(2),F(7,3),F(1)]):
        for A in (3,10,40):
            mix.append(finite_compensator(g,A,min(g[0],1),min(g[-1],1)))
    checks=[]
    for n in (3,5,9,16,30,45):
        for _ in range(10):
            z=[F(0)]
            for j in range(n):
                z.append(z[-1]+F(r.randint(1,30),r.randint(1,20)))
            checks.append(l6(z,0,len(z)-1))
    return {"schema":SCHEMA,"l5":l5,"bellman":bellman,"mixed":mix,"l6":checks,
            "counts":{"L5":len(l5),"ramp":len(bellman),"mixed":len(mix),"L6":len(checks)},
            "analytic_independent_L5_proof_in_README":True,
            "analytic_general_L6_reconstruction_proposed_for_audit":True,
            "independent_peer_review_completed":False,
            "Riemann_Hypothesis_proved":False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    x=verify()
    args.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92R PASSED:",x["counts"],"SHA256",hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__=="__main__":
    main()
