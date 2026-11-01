#!/usr/bin/env python3
"""Rational finite enclosures for q-digamma and ideal Bose excited occupation."""
import json
from pathlib import Path
from fractions import Fraction as F
from gibbs_euler import log_interval,exact_decimal_bound

def lambert(q,x,N):
    q=F(q)
    if not 0<q<1 or x<1 or N<2:raise ValueError('domain')
    s=sum((q**n/(1-q**n) for n in range(x,x+N)),F(0))
    t=q**(x+N)/(1-q)
    return s+t,s+t/(1-q**(x+N))

def divisor_enclosure(q,N):
    q=F(q)
    d=[0]*(N+1)
    for j in range(1,N+1):
        for n in range(j,N+1,j):d[n]+=1
    s=sum((d[n]*q**n for n in range(1,N+1)),F(0))
    return s,s+q**(N+1)*((N+1)-N*q)/(1-q)**2

def psi_enclosure(q,x,N):
    a,b=lambert(q,x,N)
    la,lb=log_interval(q)
    ca,cb=log_interval(1-q)
    assert lb<0
    return -cb+la*b,-ca+lb*a

def overlaps(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def fmt(p):return [exact_decimal_bound(p[0],False,30),exact_decimal_bound(p[1],True,30)]

def run():
    cases=[]
    for q in (F(1,2),F(3,4),F(9,10)):
        d=divisor_enclosure(q,128)
        assert overlaps(d,lambert(q,1,128))
        for x in (1,2,3):
            s=lambert(q,x,128);psi=psi_enclosure(q,x,128)
            assert 0<s[0]<s[1] and s[1]-s[0]<F(1,10**8)
            assert psi[0]<psi[1]
            if x<=2:
                nxt=psi_enclosure(q,x+1,128)
                lq=log_interval(q)
                factor=q**x/(1-q**x)
                jump=(-lq[1]*factor,-lq[0]*factor)
                assert overlaps((nxt[0]-psi[1],nxt[1]-psi[0]),jump)
            r={'q':str(q),'x':x,'S_q':fmt(s),'psi_q':fmt(psi)}
            if x==1:r['independent_divisor_sum']=fmt(d)
            cases.append(r)
    return {'schema':'NUMERIA.V20.QDIGAMMA_BOSE.1',
            'status':'RATIONAL_ENCLOSURES_AND_QDIGAMMA_RECURRENCES_PASS',
            'checks':{'evaluations':9,'divisor_crosschecks':3,'recurrences':6},
            'cases':cases,
            'not_certified':['Full text unavailable for paper equation audit',
                'Real q-digamma is not complex Weil digamma',
                'No Gibbs-Weil identity or RH proof']}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();o=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: 9 q-digamma enclosures, 3 independent divisor crosschecks, 6 recurrences')
