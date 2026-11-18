"""Pure-Python exact dyadic endpoints; no flint, floats or unproved linear algebra.

All enclosure decisions use integer arithmetic. Interval addition is outward-
rounded to MAX_BITS at the *endpoint level* to bound output growth; unlike
ball-to-ball summation it never silently drops an error term.
"""
from __future__ import annotations

MAX_BITS = 1200
N = 377
UPPER = N*(N+1)//2
T = 504
TOTAL_PANELS = 2*T


def _term(x):
    if not isinstance(x,list) or len(x)!=2 or any(type(v)!=str for v in x):
        raise ValueError('dyadic must be [integer-as-string, exponent-as-string]')
    if not x[0].lstrip('-').isdigit() or not x[1].lstrip('-').isdigit():
        raise ValueError('noninteger dyadic')
    m,e=int(x[0]),int(x[1])
    if m.bit_length()>5000 or not -32768<=e<=32768:
        raise ValueError('dyadic exceeds resource limits')
    return m,e


def as_ball_interval(ball):
    """P236/P237 ball [[mid,exp],[nonnegative_radius,exp]] -> (lo, hi, exp)."""
    if not isinstance(ball,list) or len(ball)!=2:raise ValueError('bad ball')
    m,e=_term(ball[0]);r,k=_term(ball[1])
    if r<0:raise ValueError('negative radius')
    a=min(e,k);M=m<<(e-a);R=r<<(k-a)
    return M-R,M+R,a


def as_endpoint_interval(endpoints):
    if not isinstance(endpoints,list) or len(endpoints)!=2:raise ValueError('bad endpoints')
    lm,le=_term(endpoints[0]);um,ue=_term(endpoints[1]);z=min(le,ue)
    lo=lm<<(le-z);hi=um<<(ue-z)
    if lo>hi:raise ValueError('inverted endpoints')
    return lo,hi,z


def outbound_reduce(l,u,e,max_bits=MAX_BITS):
    if l>u:raise ValueError('inverted endpoints')
    shift=max(0,max(abs(l).bit_length(),abs(u).bit_length())-max_bits)
    if shift:
        d=1<<shift
        l=l//d
        u=-((-u)//d)
        e+=shift
    if l==0 and u==0:return 0,0,0
    # Divide by all common powers of 2, permitting negative endpoints.
    while l%2==0 and u%2==0:
        l//=2;u//=2;e+=1
    return l,u,e


def plus(a,b,max_bits=MAX_BITS):
    lo,hi,e=a;bl,bh,be=b
    z=min(e,be)
    return outbound_reduce((lo<<(e-z))+(bl<<(be-z)),
                           (hi<<(e-z))+(bh<<(be-z)),z,max_bits)


def to_wire(z):
    l,u,e=z
    if l>u:raise ValueError('bad wire interval')
    return [[str(l),str(e)],[str(u),str(e)]]


def upper_indices():
    for i in range(N):
        for j in range(i,N):yield i,j


def expected_ranges():
    return [(0,25)]+[(i,min(i+25,TOTAL_PANELS)) for i in range(25,TOTAL_PANELS,25)]


def batch_name(sector,start,stop):
    return f'P237_CHUNK_{sector}_{start:04d}_{stop:04d}.json.gz'


def validate_headers(d,sector,start,stop,baseline=False):
    s=('NUMERIA.ES8P236.ARB_377_25PANEL_CHUNK.1' if baseline
       else 'NUMERIA.ES8P237.ARB_377_QUADRATURE_CHUNK.1')
    if d.get('schema')!=s or d.get('sector')!=sector:
        raise ValueError('wrong schema or sector')
    checks={'L':'161/200','Tsharp':504,'n_modes':N,
            'degree_max':752 if sector=='even' else 753,
            'bits':900,'flint_version':'0.9.0',
            'panel_range':[start,stop],
            'GL60_node_evaluations':60*(stop-start),
            'upper_triangle_count':UPPER}
    for k,v in checks.items():
        if d.get(k)!=v:raise ValueError(f'wrong {k}: {d.get(k)} expected {v}')
    for key in ('whole_integral_computed','analytic_GL60_error_added',
                'beta_and_poles_added','LDL_tested','RH_proved'):
        if d.get(key) is not False:raise ValueError('unjustified status '+key)
    entries=d.get('C_partial_upper_triangle')
    if not isinstance(entries,list) or len(entries)!=UPPER:
        raise ValueError('wrong triangular length')
    return entries
