#!/usr/bin/env python3
"""ES8-P89: rigorous 256-bit upper bound ||b''''_(h=2/5)||_1 < 70000."""
import json,argparse
from pathlib import Path
from flint import arb,ctx
import flint
def run():
    ctx.prec=256
    h=arb(2)/5;u0=arb(995)/1000;n=4096
    step=u0/n;acc=arb(0)
    for k in range(n):
        u=(u0*k/n).union(u0*(k+1)/n)
        v=1-u*u
        if not bool(v.lower()>0):raise ArithmeticError('bad partition')
        poly=4*(30*u**10+45*u**8-132*u**6+58*u**4+6*u**2-3)
        acc+=step*((1-1/v).exp()*poly/v**8).abs_upper()
    v0=1-u0*u0
    if not bool(v0.upper()<arb(1)/8):raise AssertionError('invalid tail monotonicity')
    # For u0<=u<1, 0<v<=v0<1/8, exp(1-1/v)*v^-8 increasing in v
    # and |4*(30u^10+...-3)| <= 1096; endpoint bound includes full interval.
    tail=(1-u0)*arb(1096)*(1-1/v0).exp()/v0**8
    whole=2*(acc+tail)/h**3
    if not bool(whole.upper()<arb(70000)):
        raise ArithmeticError('L1 D4 bound not proved: '+str(whole))
    def bin(x):return [str(v) for y in (x.mid(),x.rad()) for v in y.man_exp()]
    return {'status':'PROVED_D4_L1_LESS_THAN_70000','bits':ctx.prec,
            'python_flint':flint.__version__,'segments':n,
            'u0':'995/1000','D4_L1_upper':bin(whole),'claimed_upper':70000,
            'tail_600_upper_formula':'14*70000^2/(6*pi*600^6)',
            'scope':'Certified fourth derivative norm and absolute Fourier tail only'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    a=p.parse_args();r=run()
    Path(a.out).write_text(json.dumps(r,indent=2)+'\n')
    print('CERTIFIED_D4_L1',r['D4_L1_upper'])
