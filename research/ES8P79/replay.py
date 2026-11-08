#!/usr/bin/env python3
"""ES8-P79: independent 256-bit Gauss-Legendre/Bernstein replay, 30 selected pairs.
No acb.integral; historical values are only used for exact interval comparison.
Does not certify the underlying Weil identification or prove RH.
"""
import argparse, hashlib, json, platform, sys, time
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
DEGREES={'a030':(24,), 'a034':(26,), 'a040':(11,22,11)}
RHO_N,RHO_D=34,25
PARTS=12

def interval(v):
    m,p,r,q=map(int,v)
    if r<0: raise ValueError('negative radius')
    return F(m)*F(2)**p-F(r)*F(2)**q,F(m)*F(2)**p+F(r)*F(2)**q

def compare(a,b):
    x,y=interval(a);z,w=interval(b)
    if y<z or w<x:return 'DISJOINT'
    if z<=x and y<=w:return 'REPLAY_SUBSET_REFERENCE'
    if x<=z and w<=y:return 'REFERENCE_SUBSET_REPLAY'
    return 'OVERLAP'

def sample_pairs(data,label,shard,shards=2):
    if not 0<=shard<shards:raise ValueError('invalid shard')
    if label not in DEGREES or tuple(data['degrees'][label])!=DEGREES[label]:
        raise ValueError('unexpected basis')
    d=data['entries'][label]
    pairs=sorted(d,key=lambda s:tuple(map(int,s.split(','))))
    for s in pairs:
        i,j=map(int,s.split(','))
        if not 0<=i<=j<sum(DEGREES[label]):raise ValueError(s)
        interval(d[s])
    return pairs[shard::shards]

def replay(label,shard=0,shards=2,bits=256,target_exp=23,out=None,mode='sample'):
    if bits!=256 or not 15<=target_exp<=50:
        raise ValueError('audit policy requires 256 bits and 15<=target_exp<=50')
    import flint
    from flint import arb,acb,ctx
    ctx.prec=bits
    if ctx.prec!=bits:raise AssertionError('failed to set precision')
    data_bytes=(HERE/'sample_reference.json').read_bytes()
    data=json.loads(data_bytes)
    if shards<1 or not 0<=shard<shards: raise ValueError('invalid shard configuration')
    if mode=='sample':
        selected=sample_pairs(data,label,shard,shards)
    elif mode=='all':
        N=sum(DEGREES[label])
        selected=[f'{i},{j}' for i in range(N) for j in range(i,N)][shard::shards]
    else:
        raise ValueError('unknown replay mode')
    if out is not None and Path(out).exists():raise FileExistsError(str(out))
    begin=time.monotonic()
    a={'a030':arb(3)/10,'a034':arb(17)/50,'a040':arb(2)/5}[label]
    t=arb(2).log()
    pieces=([(-t/2,a-t/2,11),(arb(0),t-a,22),(t/2,a-t/2,11)]
            if label=='a040' else [(arb(0),a,sum(DEGREES[label]))])
    idx=[(i,p) for i,(_,_,n) in enumerate(pieces) for p in range(n)]
    assert len(idx)==sum(DEGREES[label])
    pi=arb.pi();s0=acb(arb(1)/4)
    psiR=acb(arb(1)/4,arb(15)).digamma().real
    mR=psiR-pi.log()
    def psi_symmetric(z):
        shift=acb(0,1)*z/2
        return ((s0+shift).digamma()+(s0-shift).digamma())/2
    def odd_factorial(p):
        r=arb(1)
        for k in range(1,2*p+2,2):r*=k
        return r
    def j(p,z):
        return z**p/odd_factorial(p)*(-z*z/4).hypgeom_0f1(arb(p)+arb(3)/2)
    def jj(p,z):
        return z**p/odd_factorial(p)*(z*z/4).hypgeom_0f1(arb(p)+arb(3)/2)
    uv=[]
    for I,p in idx:
        c,h,_=pieces[I]
        factor=((arb(2*p+1)*h/2).sqrt())*2*jj(p,h/2)
        uv.append(((c/2).exp()*factor,(-c/2).exp()*factor*(-1 if p%2 else 1)))
    rho=arb(RHO_N)/RHO_D
    L=arb(5)/4
    ellipse_A=L*(rho+1/rho)/2
    ellipse_B=L*(rho-1/rho)/2
    sigma=arb(1)/4-ellipse_B/2
    if not bool(sigma.lower()>0):raise ArithmeticError('ellipse hits digamma poles')
    psiR_upper=arb(103) # proved loose bound from convergent digamma series
    target=arb(10)**(-target_exp)
    root_cache={}
    counts={'quadrature_nodes':0,'max_order':0,'structural_zeros':0,'integrated_pairs':0}
    def majorant(p,h,q,k,delta,seg):
        rmax=arb(5)*seg/2+L+ellipse_A
        psi_bound=arb(1)+(rmax/2+arb(5)/4)*(1/sigma+2)
        def j_bound(deg,half):
            return (half*rmax)**deg/odd_factorial(deg)*(half*ellipse_B).cosh()
        return (j_bound(p,h)*j_bound(q,k)*
                (psiR_upper+psi_bound)*(abs(delta)*ellipse_B).cosh()).upper()
    def integrate(p,h,q,k,delta,even):
        total=arb(0);error=arb(0)
        for seg in range(PARTS):
            M=majorant(p,h,q,k,delta,seg)
            n=16
            err=8*L*M*rho**(1-2*n)/(rho-1)
            while not bool(err.upper()<target):
                n+=2
                if n>512:raise ArithmeticError('Gauss order exceeded 512')
                err=8*L*M*rho**(1-2*n)/(rho-1)
            if n not in root_cache:
                root_cache[n]=[arb.legendre_p_root(n,i,weight=True) for i in range(n)]
            center=arb(5)*seg/2+L
            s=acb(0)
            for x,w in root_cache[n]:
                r=center+L*x
                tr=(r*delta).cos() if even else (r*delta).sin()
                s+=w*j(p,h*r)*j(q,k*r)*(psiR-psi_symmetric(r))*tr
            term=L*s
            if not term.is_finite():raise ArithmeticError('nonfinite quadrature')
            total+=term.real
            error+=err.upper()
            counts['quadrature_nodes']+=n
            counts['max_order']=max(n,counts['max_order'])
        return total+arb(0,error.upper())
    def encode(x):
        m,e=x.mid().man_exp();r,f=x.rad().man_exp()
        return [str(m),str(e),str(r),str(f)]
    records=[]
    for key in selected:
        i,j_ix=map(int,key.split(','))
        (I,p),(J,q)=idx[i],idx[j_ix]
        c,h,_=pieces[I];d,k,_=pieces[J]
        delta=c-d;even=(p+q)%2==0
        if I==J and not even:
            arch=arb(0);counts['structural_zeros']+=1
        else:
            integral=integrate(p,h,q,k,delta,even)
            if label=='a040':
                phase=(-1)**((p-q)//2 if even else (p-q+1)//2)
                fac=arb(2)/pi*(h*k*arb((2*p+1)*(2*q+1))).sqrt()
            else:
                phase=(-1)**((q-p)//2)
                fac=2*h/pi*arb((2*p+1)*(2*q+1)).sqrt()
            arch=phase*fac*integral
            counts['integrated_pairs']+=1
        pole=uv[i][0]*uv[j_ix][1]+uv[i][1]*uv[j_ix][0]
        prime=t/arb(2).sqrt() if label=='a040' and {I,J}=={0,2} and p==q else arb(0)
        value=(mR if i==j_ix else arb(0))+pole-prime-arch
        ball=encode(value)
        rel=compare(ball,data['entries'][label][key]) if key in data['entries'][label] else 'NOT_COMPARED'
        records.append({'pair':key,'replay':ball,'relation':rel})
        print(label,shard,key,rel,flush=True)
        if rel=='DISJOINT':raise ArithmeticError('disjoint result '+label+' '+key)
        if int(ctx.prec)!=bits:raise AssertionError('precision mutated')
    result={'status':'ARB_REPLAY_COMPLETE','window':label,'mode':mode,
            'shard':shard,'shards':shards,'bits':int(ctx.prec),
            'flint_version':flint.__version__,
            'fixture_sha256':hashlib.sha256(data_bytes).hexdigest(),
            'records':records,'counts':counts,'elapsed_seconds':round(time.monotonic()-begin,2),
            'platform':platform.platform(),'python':sys.version,
            'disclaimer':'interval comparison does not prove analytic majorants or RH'}
    if out:Path(out).write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--window',choices=sorted(DEGREES),required=True)
    ap.add_argument('--shard',type=int,required=True)
    ap.add_argument('--shards',type=int,default=2)
    ap.add_argument('--mode',choices=['sample','all'],default='sample')
    ap.add_argument('--out',required=True)
    ap.add_argument('--bits',type=int,default=256)
    ap.add_argument('--target-exp',type=int,default=23)
    args=ap.parse_args()
    replay(args.window,args.shard,args.shards,args.bits,args.target_exp,args.out,args.mode)
