"""P75 full-cell directed-Arb Suzuki positivity experiment, not global RH.

The proof gate uses Chourasiya--Simonic v2 PREPRINT coefficients, a published
finite-RH verification and Bellotti--Wong zero counting. A positive verdict is
conditional on independent mathematical review of all cited inputs.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx
from NUMERIA_ES8P259.suzuki import _dyadic
from NUMERIA_P74_ZERO_SIDE.zeros import verified_ordinates
from NUMERIA_P74_ZERO_SIDE.certify import Spectrum, rational, ball
from .cs_v2 import ChourasiyaV2Tail
from .local import local_curvature_upper, lower_secant_local

SCHEMA='NUMERIA.P75.CS_V2.LOCAL_CURVATURE.1'
MAX_DEPTH=12


def cover(start,end,ordinates,tails,base_width,depth_limit=MAX_DEPTH):
    spectrum=Spectrum(ordinates)
    # Value enclosure is inherited from P74; no random sampling in the proof.
    count=len(ordinates)
    evaluated=0; accepted=0; attempts=0; deepest=0
    cached={}; sha=sha256(); best=None; failure=None; cursor=start

    def value(t):
        nonlocal evaluated
        if t not in cached:
            cached[t]=spectrum.value(t)
            evaluated+=1
        return cached[t]

    def attempt(a,b,depth):
        nonlocal attempts,accepted,deepest,best,cursor
        attempts+=1
        mid=ball((a+b)/2)
        width=ball(b-a)
        ceiling=local_curvature_upper(ordinates,mid,width)
        lower=lower_secant_local(value(a),value(b),width,ceiling)
        margin=(lower - tails.bound(ball(b))).lower()
        if margin>0:
            accepted+=1; deepest=max(deepest,depth); cursor=b
            if best is None or margin<best[0]:
                best=(margin,[str(a),str(b)])
            record=[str(a),str(b),_dyadic(margin),_dyadic(ceiling)]
            sha.update((json.dumps(record,separators=(',',':'))+'\n').encode('ascii'))
            return None
        if depth>=depth_limit:
            return {'cell':[str(a),str(b)],'lower_dyadic':_dyadic(margin),'depth':depth}
        m=(a+b)/2
        return attempt(a,m,depth+1) or attempt(m,b,depth+1)

    a=start
    while a<end:
        b=min(end,a+base_width)
        failure=attempt(a,b,0)
        if failure is not None:
            break
        a=b
    return {
        'status':'conditional_finite_positive' if failure is None else 'inconclusive',
        'positive_under_imported_lemmas':failure is None,
        'certified_through_exact':str(cursor) if accepted else None,
        'cell_count':accepted,
        'attempted_cells':attempts,
        'endpoint_evaluations':evaluated,
        'max_depth':deepest,
        'minimum':{'lower_dyadic':_dyadic(best[0]),'cell':best[1]} if best else None,
        'first_unresolved_cell':failure,
        'ordered_cells_sha256':sha.hexdigest(),
    }


def run(end='79',count=500,prec=256,segments=128,mesh='1/64',start='18'):
    start=rational(start); end=rational(end); mesh=rational(mesh)
    if not (0<start<end and 0<mesh<=1):
        raise ValueError('invalid interval/mesh')
    if type(count) is not int or count<1:
        raise ValueError('invalid count')
    if type(prec) is not int or prec<128:
        raise ValueError('precision too low')
    with ctx.workprec(prec):
        zeros,zero_meta=verified_ordinates(count)
        tails=ChourasiyaV2Tail(segments)
        outcome=cover(start,end,zeros,tails,mesh)
        return {
          'schema':SCHEMA,'begin_t':str(start),'end_t':str(end),
          'zero_count':count,'precision_bits':prec,'bands':segments,
          'mesh_exact':str(mesh),'zero_provenance':zero_meta,
          'tails':tails.metadata(),
          'overlaps_p71':bool(ball(start).exp()<10**8),
          'publication_status':'DRAFT; source v2 is a preprint and must be audited independently',
          **outcome,
        }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--end',default='79');p.add_argument('--start',default='18')
    p.add_argument('--zeros',type=int,default=500);p.add_argument('--precision',type=int,default=256)
    p.add_argument('--segments',type=int,default=128);p.add_argument('--mesh',default='1/64')
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    data=run(args.end,args.zeros,args.precision,args.segments,args.mesh,args.start)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:data[k] for k in ('status','certified_through_exact','cell_count','minimum','first_unresolved_cell')},indent=2))
    # Research reports are valid even for inconclusive results. CI gates separately.

if __name__=='__main__':main()
