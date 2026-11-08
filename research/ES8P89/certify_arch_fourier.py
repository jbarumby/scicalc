#!/usr/bin/env python3
"""ES8-P89: rigorous finite archimedean integral via FLINT/Arb.
Each shard computes 30 of 240 segments of [0,600], integrating the Fourier
transform of the actual smooth compact-support bump by certified complex
integration, using degree-14 Taylor shifts on 95 certified GL nodes per segment.
The infinite Fourier tail is NOT included in the shard outputs; see
certify_derivative.py for a separate rigorous r>600 tail.
"""
import argparse,json,math
from pathlib import Path
from flint import acb,arb,ctx
import flint
CENTERS=(1,2,3,5,7,11,13,17)
PARTS=240
NGL=95
K=14
def run(shard,shards=8,pilot=0):
    if not 0<=shard<shards:raise ValueError('bad shard')
    ctx.prec=256
    h=arb(2)/5;L=arb(5)/4
    rho=arb(34)/25
    ellipseB=L*(rho-1/rho)/2
    ellipseA=L*(rho+1/rho)/2
    sigma=arb(1)/4-ellipseB/2
    if not bool(sigma.lower()>0):raise ArithmeticError('ellipse intersects digamma poles')
    rmax=arb(600)+ellipseA
    # Digamma bound from convergent series; b-hat bound from 0 <= b <=1.
    psi_upper=arb(1)+(rmax/2+arb(5)/4)*(1/sigma+2)+2
    b_upper=2*h*(h*ellipseB).cosh()
    majorant=psi_upper*b_upper**2*(arb(17).log()*ellipseB).cosh()
    if not bool(majorant.upper()<arb(10000)):
        raise ArithmeticError('Bernstein majorant M=10000 fails')
    errseg=8*L*arb(10000)*rho**(1-2*NGL)/(rho-1)
    pi=arb.pi()
    if not bool((arb(PARTS)*errseg/pi).upper()<arb('1e-9')):
        raise ArithmeticError('GL analytic truncation too large')
    taylor=2*h*(h*L)**(K+1)/math.factorial(K+1)
    if not bool(taylor.upper()<arb('3e-17')):
        raise ArithmeticError('Taylor remainder too large')
    def b(z,analytic):
        re=z.real
        if bool(re.lower()>-h) and bool(re.upper()<h):
            return (1-acb(1)/(1-(z/h)**2)).exp()
        if bool(re.upper() < -h) or bool(re.lower()>h):
            return acb(0)
        if analytic:return acb(arb('nan'))
        if not bool(z.imag==0):return acb(arb('nan'))
        return acb(arb('0.5 +/- 0.5'))
    def fmt(v):return [str(y) for x in (v.mid(),v.rad()) for y in x.man_exp()]
    roots=[arb.legendre_p_root(NGL,i,weight=True) for i in range(NGL)]
    centers=[arb(n).log() for n in CENTERS]
    pairs=[(i,j) for i in range(8) for j in range(i,8)]
    vals={ij:arb(0) for ij in pairs}
    segs=0
    for seg in range(shard,PARTS,shards):
        if pilot and segs>=pilot:break
        c=arb(5)*seg/2+L
        mus=[]
        for k in range(K+1):
            def integrand(z,analytic):
                return b(z,analytic)*z**k*(acb(0,1)*c*z).exp()
            mu=acb.integral(integrand,-h,h,abs_tol=arb('1e-13'),
                            rel_tol=arb('1e-24'),deg_limit=180,
                            eval_limit=300000,depth_limit=90)
            if not mu.is_finite():raise ArithmeticError(f'nonfinite seg={seg} k={k}')
            if not bool(mu.real.rad().upper()<arb('2e-12') and mu.imag.rad().upper()<arb('2e-12')):
                raise ArithmeticError(f'moment radius seg={seg} k={k}: {mu}')
            mus.append(mu)
        for x,w in roots:
            r=c+L*x
            delta=r-c
            F=acb(0);term=acb(1)
            for k in range(K+1):
                F+=term*mus[k]
                term*=acb(0,delta)/(k+1)
            B=F.real+arb(0,taylor.upper())
            z=acb(0,r/2)
            t=acb(arb(1)/4)
            m=((t+z).digamma()+(t-z).digamma()).real/2-pi.log()
            common=w*L/pi*m*B*B
            for i,j in pairs:
                vals[(i,j)]+=common*((centers[j]-centers[i])*r).cos()
        segs+=1
        print('SEGMENT',seg,'PROVED finite-interval quadrature',flush=True)
    for ij in pairs:
        vals[ij]+=arb(0,arb(segs)*errseg/pi)
    return {'status':'PILOT_PARTIAL_SEGMENTS_ONLY' if pilot else 'RIGOROUS_FOURIER_0_TO_600',
            'mode':'8-center-bump','bits':ctx.prec,'python_flint':flint.__version__,
            'shard':shard,'shards':shards,'segments':segs,'expected_segments':PARTS,
            'degree':K,'gl_nodes':NGL,
            'analytic_bound_M':fmt(majorant),'taylor_error':fmt(taylor),
            'single_segment_quad_error':fmt(errseg),
            'records':[{'i':i,'j':j,'arch_finite':fmt(vals[(i,j)])} for i,j in pairs],
            'scientific_limit':'No Fourier tail included in this file, no full Weil matrix certificate'}
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--shard',type=int,required=True)
    p.add_argument('--shards',type=int,default=8)
    p.add_argument('--pilot-segments',type=int,default=0)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    out=Path(a.out)
    if out.exists():raise FileExistsError(out)
    result=run(a.shard,a.shards,a.pilot_segments)
    out.write_text(json.dumps(result,indent=2)+'\n')
