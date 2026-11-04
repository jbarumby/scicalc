#!/usr/bin/env python3
"""NUMERIA ES8-P9: independent weighted integration of all Phase-4 nodal enclosures.

No floating-point operations in the enclosing algorithm. Exact rational
Newton-Cotes weights; outward-fixed 1e-35 function evaluations from ES8-P8;
independent Cauchy interpolation and E7 tail certificates from ES8-P6.

All certificates are conditional on explicitly listed analytic lemmas.
Full 1600-cell execution required to assert the 30-entry integrated crosscheck.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
from fractions import Fraction as Q
import hashlib
import json
from math import lcm
from pathlib import Path
import time

from es8p8_full import IndependentFixed, PAIRS
from es8p6_independent_audit import independent_total_error, sharp_tail_upper, baseline_tail_upper, lagrange_weights

BASE=Path(__file__).resolve().parent
REFERENCE=BASE/'HISTORICAL_PHASE4_REFERENCE.json'
NODE_COUNT=16
CELLS=1600
CELL_H=Q(1,10)


def vandermonde_weights(n=16):
    """Independently solve exact integration moment equations via Gauss-Jordan."""
    t=[Q(j,n-1) for j in range(n)]
    a=[[t[j]**k for j in range(n)]+[Q(1,k+1)] for k in range(n)]
    for k in range(n):
        p=next(i for i in range(k,n) if a[i][k])
        a[k],a[p]=a[p],a[k]
        v=a[k][k]
        a[k]=[u/v for u in a[k]]
        for i in range(n):
            if i!=k:
                factor=a[i][k]
                if factor:
                    a[i]=[u-factor*v for u,v in zip(a[i],a[k])]
    w=[row[-1] for row in a]
    assert all(sum(w[j]*t[j]**k for j in range(n))==Q(1,k+1) for k in range(n))
    assert w==lagrange_weights(n), 'independent weight derivations differ'
    return w


def reference_inputs():
    d=json.loads(REFERENCE.read_text(encoding='utf-8'))
    assert d['source_phase4_sha256']=='b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe'
    assert len(d['historical_integrals'])==30
    assert set(d['historical_integrals'])=={f'{i},{j}' for i,j in PAIRS}
    for e in d['historical_integrals'].values():
        assert e['cells']==CELLS and e['nodes_per_cell']==NODE_COUNT
        assert e['interpolation_error_exact']==str(independent_total_error())
    return d


def outward_decimal(x: Q, digits=24, upper=False):
    s=10**digits
    k=-( (-x.numerator*s)//x.denominator) if upper else x.numerator*s//x.denominator
    a=abs(k)
    return ('-' if k<0 else '')+str(a//s)+'.'+str(a%s).zfill(digits)


def total_interval(raw, shared_denominator):
    return Q(raw[0], shared_denominator),Q(raw[1], shared_denominator)


def reassemble_integrals(*,start=0,end=1600,batch=50,terms=128,digits=35):
    if not 0<=start<end<=CELLS: raise ValueError('invalid subset')
    if not 0<batch<=CELLS: raise ValueError('invalid batch')
    ref=reference_inputs()
    weights=vandermonde_weights()
    den=lcm(*(x.denominator for x in weights))
    ints=[int(w*den) for w in weights]
    assert sum(ints)==den and any(x<0 for x in ints)
    # Correct sign selection is essential for Newton-Cotes (negative weights).
    running={pair:[0,0] for pair in PAIRS}
    batches=[];timing=time.monotonic()
    scale=10**digits
    for first in range(start,end,batch):
        last=min(end,first+batch)
        E=IndependentFixed(digits=digits,N=terms)
        partial={pair:[0,0] for pair in PAIRS}
        for cell in range(first,last):
            for node,w in enumerate(ints):
                r=Q(cell,10)+Q(node,150)
                for pair in PAIRS:
                    a,b=E.qij(r,*pair)
                    if w>=0: partial[pair][0]+=w*a;partial[pair][1]+=w*b
                    else: partial[pair][0]+=w*b;partial[pair][1]+=w*a
        h=hashlib.sha256()
        for i,j in PAIRS:
            lo,hi=partial[(i,j)]
            if lo>hi:raise ArithmeticError('invalid integral bounds')
            running[(i,j)][0]+=lo;running[(i,j)][1]+=hi
            h.update(f'{i},{j}:{lo},{hi}\n'.encode('ascii'))
        batches.append({'start_cell':first,'end_cell_exclusive':last,'evaluations':(last-first)*NODE_COUNT*30,
                        'sha256_weighted_pair_sums':h.hexdigest()})
        print(f'ES8-P9 cells {first}-{last} ({last-start}/{end-start}), elapsed={time.monotonic()-timing:.2f}s',flush=True)
    common=10*den*scale
    finite={}; compatible=0; interlo=Q(0); interhi=Q(0)
    err=independent_total_error()
    for i,j in PAIRS:
        pair=(i,j)
        key=f'{i},{j}'
        low,high=total_interval(running[pair],common)
        base=ref['historical_integrals'][key]
        old_lo=Q(Decimal(base['nodal_integral_lower']))
        old_hi=Q(Decimal(base['nodal_integral_upper']))
        nodal_intersection=(low<=old_hi and old_lo<=high)
        interval=(low-err,high+err)
        phase4_lo=Q(Decimal(base['integral_lower']))
        phase4_hi=Q(Decimal(base['integral_upper']))
        finite_intersection=interval[0]<=phase4_hi and phase4_lo<=interval[1]
        if end-start==CELLS:
            if not nodal_intersection or not finite_intersection:
                raise ArithmeticError(f'Historical independent sum mismatch {(i,j)}: nodal={nodal_intersection}, finite={finite_intersection}')
            compatible+=1
        t=sharp_tail_upper(i,j)
        assert 0<t<baseline_tail_upper(i,j)
        finite[key]={'i':i,'j':j,
            'independent_nodal_0_160':[outward_decimal(low),outward_decimal(high,upper=True)],
            'independent_finite_with_cauchy':[outward_decimal(interval[0]),outward_decimal(interval[1],upper=True)],
            'historical_nodal_0_160':[base['nodal_integral_lower'],base['nodal_integral_upper']],
            'historical_finite_with_cauchy':[base['integral_lower'],base['integral_upper']],
            'nodal_intersects_historical':nodal_intersection,
            'finite_intersects_historical':finite_intersection,
            'independent_global_with_tail':[outward_decimal(interval[0]-t),outward_decimal(interval[1]+t,upper=True)],
            'sharp_tail_upper_decimal_outward':outward_decimal(t,upper=True),
            'max_nodal_interval_width':outward_decimal(high-low,upper=True),
        }
    status='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED' if start==0 and end==CELLS and compatible==30 else 'PARTIAL_OR_INCOMPATIBLE'
    return {'schema':'NUMERIA.ES8P9.EXACT_WEIGHTED_INTEGRALS.1',
        'status':status,
        'scope':{'start_cell':start,'end_cell_exclusive':end,'cells':CELLS,'nodes_per_cell':NODE_COUNT,
                 'pairs':30,'evaluations':(end-start)*NODE_COUNT*30,'batch_size':batch,
                 'digits':digits,'digamma_terms':terms},
        'provenance':{'historical_source_zip_sha256':ref['source_phase4_sha256'],
                      'historical_source_entry_sha256':ref['source_entry_sha256']},
        'quadrature':{'degree_of_exactness_at_least':NODE_COUNT-1,
                      'weights_denominator':str(den),
                      'signed_integer_weights':[str(x) for x in ints],
                      'positive_weights':sum(x>0 for x in ints),'negative_weights':sum(x<0 for x in ints),
                      'interpolation_error_fraction':str(err)},
        'compatibility':{'nodal_integrals_intersecting':sum(x['nodal_intersects_historical'] for x in finite.values()),
                         'finite_integrals_intersecting':sum(x['finite_intersects_historical'] for x in finite.values())},
        'batches':batches,
        'entries':finite,
        'limitations':[
            'conditional on analytic enclosure lemmas and error bounds in ES8-P6 and ES8-P8',
            'independent sums use same exact Newton-Cotes quadrature rule as historical Phase4',
            'historical numerical interval intersections do not independently prove common special-function identities',
            'full Weil universal positivity and RH not implied by finite matrix tests',
        ]}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--end',type=int,default=1600)
    p.add_argument('--batch',type=int,default=50)
    p.add_argument('--terms',type=int,default=128)
    p.add_argument('--digits',type=int,default=35)
    p.add_argument('--out',type=Path,default=BASE/'RESULTADOS_ES8P9_INTEGRAIS.json')
    a=p.parse_args()
    data=reassemble_integrals(start=a.start,end=a.end,batch=a.batch,terms=a.terms,digits=a.digits)
    a.out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print('STATUS:',data['status'],data['compatibility'],flush=True)
    return 0 if data['status']!='PARTIAL_OR_INCOMPATIBLE' or a.end-a.start<CELLS else 1

if __name__=='__main__':raise SystemExit(main())
