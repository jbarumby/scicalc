#!/usr/bin/env python3
"""Generate exact DYADIC enclosures of the first N zeta zeros with FLINT/Arb.

Requires: python-flint >=0.9.0. Run e.g.:
  python -m pip install python-flint==0.9.0
  python export_zeros_flint.py --n 20 --out arb_zeros_20.json
  python verify_spectral_rational.py --zeros arb_zeros_20.json --n 20 --to 44 --out arb_cert_44.json

We do not round decimal display strings to construct the enclosure; all
endpoints are exact dyadic rationals from rigorous Arb balls.
"""
import argparse,json
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=20)
    p.add_argument('--precision',type=int,default=192)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not (1<=a.n<=50000 and a.precision>=96):raise ValueError('invalid n/precision')
    from flint import acb, arb, ctx
    rows=[]
    with ctx.workprec(a.precision):
        zeros=acb.zeta_zeros(1,a.n)
        for index,z in enumerate(zeros,1):
            if not (z.real == arb(1)/2):
                raise ArithmeticError('unexpected off-line value from FLINT')
            im=z.imag
            lo=im.lower().man_exp(); hi=im.upper().man_exp()
            if lo==hi:
                raise ArithmeticError('unexpectedly exact zero ordinate')
            rows.append({'index':index,'lo_dyad':[str(lo[0]),int(lo[1])],
                         'hi_dyad':[str(hi[0]),int(hi[1])]})
    payload={ 'schema':'NUMERIA.P72.ZEROS.1',
              'zero_certification':'arb_enclosures',
              'generator':'python-flint acb.zeta_zeros(1,n); rigorous isolation and Turing counting',
              'precision_bits':a.precision,'python_flint_version':getattr(__import__('flint'),'__version__','unknown'), 'zeros':rows}
    a.out.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print('written',a.out,'count',len(rows))

if __name__=='__main__':main()
