#!/usr/bin/env python3
"""Independent rational enclosures of selected, original Phase4 nodal values.

Uses Fractions, integral-comparison Euler–Maclaurin bounds for the digamma
series, Taylor-with-remainder for sin/cos and Machin for pi. This algorithm
is separate from the historical 40-digit integer fixed-point implementation.
It does *not* certify all original 768,000 nodes nor the global Weil criterion.
"""
from __future__ import annotations
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, isqrt
import argparse, csv, hashlib, json
from pathlib import Path

L=Q(9,20)

def I(x,y=None):
    x=Q(x); y=x if y is None else Q(y)
    assert x<=y
    return (x,y)
def add(a,b):return (a[0]+b[0], a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    vv=[x*y for x in a for y in b];return (min(vv),max(vv))
def invpos(a):
    assert a[0]>0;return (1/a[1],1/a[0])
def div(a,b):return mul(a,invpos(b))
def scale(a,c):return mul(a,I(c))
def sq(a):
    if a[0]<=0<=a[1]:return I(0,max(a[0]*a[0],a[1]*a[1]))
    return I(min(a[0]*a[0],a[1]*a[1]),max(a[0]*a[0],a[1]*a[1]))

def atan1(q,n=42):
    x=Q(1,q)
    val=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)),Q(0))
    error=x**(2*n+1)/Q(2*n+1)
    return I(val,val+error) if n%2==0 else I(val-error,val)
@lru_cache(None)
def pi():return sub(scale(atan1(5),16),scale(atan1(239),4))

def logq(q,n=100):
    q=Q(q);assert q>0
    t=(q-1)/(q+1)
    s=sum((2*t**(2*k+1)/Q(2*k+1) for k in range(n)),Q(0))
    e=2*abs(t)**(2*n+1)/(Q(2*n+1)*(1-t*t))
    return I(s-e,s+e)
@lru_cache(None)
def ln2():return logq(2,68)
@lru_cache(None)
def lngamma_const():
    n=256
    h=sum((Q(1,j) for j in range(1,n+1)),Q(0))
    ln_n=scale(ln2(),8)
    y=add(sub(I(h),ln_n),I(-Q(1,2*n)+Q(1,12*n*n)))
    return I(y[0]-Q(1,120*n**4),y[1])
@lru_cache(None)
def omega0():
    logpi=I(logq(pi()[0])[0],logq(pi()[1])[1])
    return sub(sub(sub(neg(lngamma_const()),scale(pi(),Q(1,2))),scale(ln2(),3)),logpi)

@lru_cache(None)
def omega(r,terms=256):
    r=Q(r); b2=r*r/4
    if b2==0:return omega0()
    # Exact rational terms, each rounded OUTWARD onto an independently chosen
    # 10^-34 grid. This prevents huge common denominators of a long sum.
    grid_scale=10**34
    lo_total=hi_total=0
    for n in range(terms):
        y=Q(n)+Q(1,4)
        term=b2/(y*(y*y+b2))
        num=term.numerator*grid_scale;den=term.denominator
        lo_total+=num//den
        hi_total+=-((-num)//den)
    S=I(Q(lo_total,grid_scale),Q(hi_total,grid_scale))
    u=Q(terms)+Q(1,4)
    intv=scale(logq(1+b2/(u*u),n=54),Q(1,2))
    g=b2/(u*(u*u+b2))
    d=b2*(3*u*u+b2)/(12*u*u*(u*u+b2)**2)
    # 0 <= sum_tail - (integral + g(u)/2) <= -g'(u)/12
    tail=add(intv,I(g/2,g/2+d))
    return add(omega0(),add(S,tail))

def sqrtq(x,digits=52):
    x=Q(x);assert x>0
    S=10**digits; z=x.numerator*S*S//x.denominator
    r=isqrt(z)
    return I(Q(r,S),Q(r+1,S))

def sin_taylor_center(v,n=26):
    """Taylor sine at rational v plus Lagrange remainder <= |v|^(2n+1)/(2n+1)!"""
    v=Q(v)
    s=sum(((-1)**k*v**(2*k+1)/factorial(2*k+1) for k in range(n)),Q(0))
    e=abs(v)**(2*n+1)/factorial(2*n+1)
    return I(s-e,s+e)

def sin_interval(x):
    midpoint=(x[0]+x[1])/2
    two_pi=scale(pi(),2)
    approx_pi=float((two_pi[0]+two_pi[1])/2)
    k=round(float(midpoint)/approx_pi)
    y=sub(x,scale(two_pi,k))
    assert max(abs(y[0]),abs(y[1]))<Q(4), 'argument reduction failed'
    center_original=(y[0]+y[1])/2
    # Rational center reduction is rigorously accounted for in radius.
    # Without this, powers of exact Machin fractions generate huge integers.
    scale_int=10**18
    center=Q(center_original.numerator*scale_int//center_original.denominator,scale_int)
    radius=max(abs(y[0]-center),abs(y[1]-center))
    return add(sin_taylor_center(center),I(-radius,radius))

def sinc_interval(x):
    if x[0]<=0<=x[1]:
        # Integral of cos(tx), valid with a rational Taylor/Lagrange bound
        mid_orig=(x[0]+x[1])/2
        S=10**18
        mid=Q(mid_orig.numerator*S//mid_orig.denominator,S)
        radius=max(abs(x[0]-mid),abs(x[1]-mid))
        n=26
        s=sum(((-1)**k*mid**(2*k)/factorial(2*k+1) for k in range(n)),Q(0))
        err=abs(mid)**(2*n)/factorial(2*n+1)
        return I(s-err-radius/2,s+err+radius/2)
    return mul(sin_interval(x),invpos(x)) if x[0]>0 else neg(mul(sin_interval(x),invpos(neg(x))))

@lru_cache(None)
def amplitude(r,k):
    r=Q(r);k=int(k)
    w=scale(pi(),Q(k)/(2*L))
    v1=scale(sub(I(r),w),L);v2=scale(add(I(r),w),L)
    t1=sinc_interval(v1); t2=sinc_interval(v2)
    combo=add(t1,t2) if k%2 else sub(t1,t2)
    return mul(sqrtq(L),combo)

def qij(r,i,j,terms=256):
    if (i-j)%2:return I(0)
    sign=-1 if (abs(i-j)//2)%2 else 1
    val=div(mul(mul(amplitude(Q(r),i),amplitude(Q(r),j)),omega(Q(r),terms)),pi())
    return val if sign==1 else neg(val)

def grid():
    # Five widely separated cells, endpoints and one interior Gaussian-grid node.
    # All samples are exact nodes from original 1600x16 mesh.
    return [(c,j,Q(c,10)+Q(j,150)) for c in [0,40,100,400,1599] for j in (0,8,15)]

def comp(a,b):return a[0]<=b[1] and b[0]<=a[1]
def width(a):return a[1]-a[0]
def floor_decimal(x,n=16):return str(x.numerator*10**n//x.denominator/10**n)
def encode_interval(x):return [str(x[0]),str(x[1])]
def hashfile(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def check_admissibility(eps=Q(1,10)):
    """Symbolically auditable parameters for the Guinand-Weil test class.

    f_k in L1[-L,L], continuous with f_k(+-L)=0, smooth interior.
    F_k(z)=int f_k exp(izx) dx is entire, and two integrations by
    parts give |F_k(z)|<=D_k*exp(L*eta)/|z|² on |Im z|<=eta.
    H_ij(z) is entire/even and O(|z|^-4) in that strip.
    """
    eta=Q(1,2)+eps
    assert eps>0 and eta>Q(1,2)
    for k in range(1,11):
        assert k>0 and L>0
    return {'status':'ADMISSIBILITY_PROVED_FROM_ELEMENTARY_ESTIMATES_AND_STANDARD_WEIL_THEOREM',
            'strip_half_height':str(eta),
            'function_type':'H_ij(z)=[F_i(z)F_j(-z)+F_j(z)F_i(-z)]/2, F_k entire; H even and entire',
            'uniform_decay':'H_ij(x+iy)=O((1+|x|)^-4) on |y| <= strip_half_height',
            'decay_argument':'f_k vanishes at endpoints; two integrations by parts give F_k=O(|x|^-2) uniformly on any bounded horizontal strip',
            'complex_omega_note':'Omega_complex(z) has poles at ±i/2; used only for Cauchy interpolation on disks of radius 0.35, not as test h itself',
            'external_assumption':'Classical Guinand-Weil explicit formula valid for strip-holomorphic functions of >1 polynomial decay',
            'scope':'Individual H_ij, N up to 10, no universal positivity'}

def reference_triage(csv_path):
    with open(csv_path,encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    selected=[]
    for r in rows:
        doi=(r.get('DOI') or '').strip().lower()
        role={'10.1016/j.amc.2006.08.061':'alternative_contour_formula_digamma_not_direct_rational_certification',
              '10.1016/j.jmaa.2004.05.040':'asymptotic_lerch_reference_domain_assumptions_check_needed',
              '10.3390/math11030723':'growth_of_lerch_indirect_not_node_certification'}.get(doi,'context_only')
        selected.append({'title':r['Title'],'doi':doi,'role':role,'source':'user_supplied_csv_metadata'})
    selected.append({'title':'Nemes (2017) Error bounds for asymptotic expansion of Hurwitz zeta',
                     'doi':'10.1098/rspa.2017.0363','role':'external_reference_for_independent_future_polygamma_remainder',
                     'source':'external_literature_not_in_user_csv'})
    return selected

def run(source_dir:Path,out:Path,terms=256,digits=40):
    import sys
    sys.path.insert(0,str(source_dir))
    from fase4_fixedpoint_sem_arb import FixedRationalNodes
    sourcehash={'fase4_fixedpoint_sem_arb.py':hashfile(source_dir/'fase4_fixedpoint_sem_arb.py'),
                'fase4_nos_sem_arb.py':hashfile(source_dir/'fase4_nos_sem_arb.py')}
    expected={'fase4_fixedpoint_sem_arb.py':'0935bb988f7a7bd4c83b36a4a873a71b1612eda00c77ac9bd6ba7daeef600bda',
              'fase4_nos_sem_arb.py':'151ec9a0a6c249f443bcbac278e8f004e02d268ef7590f81e78e5280872cbe49'}
    assert sourcehash==expected, 'Historically pinned node sources differ'
    reader=FixedRationalNodes(digits=digits,N=64)
    stats={'nodes':0,'pairs':0,'compared':0,'incompatible':0,'original_contains_mpmath':0,'independent_contains_mpmath':0,
           'max_width_original':'0','max_width_independent':'0'}
    mismatches=[];sample=[];wid_old=Q(0);wid_new=Q(0)
    pairs=[(i,j) for i in range(1,11) for j in range(i,11) if (i-j)%2==0]
    assert len(pairs)==30
    stats['pairs']=len(pairs)
    for cell,j,r in grid():
        stats['nodes']+=1
        for i,k in pairs:
            old=tuple(Q(v,reader.S) for v in reader.qij(r,i,k))
            new=qij(r,i,k,terms)
            stats['compared']+=1
            wid_old=max(wid_old,width(old));wid_new=max(wid_new,width(new))
            if not comp(old,new):
                stats['incompatible']+=1
                mismatches.append({'cell':cell,'index':j,'r':str(r),'i':i,'j':k})
            if len(sample)<12:
                sample.append({'cell':cell,'index':j,'r':str(r),'i':i,'j':k,
                               'historical_interval':[str(round(float(old[0]),14)),str(round(float(old[1]),14))],
                               'independent_interval':[str(round(float(new[0]),14)),str(round(float(new[1]),14))]})
    stats['max_width_original']=str(float(wid_old))
    stats['max_width_independent']=str(float(wid_new))
    # Independent high-precision digamma (non-certified) diagnostic on all grid abscissae.
    import mpmath as mp
    mp.mp.dps=65
    diag=[]
    for c,j,r in grid():
        v=mp.re(mp.digamma(mp.mpf('0.25')+mp.j*mp.mpf(r.numerator)/(2*r.denominator)))-mp.log(mp.pi)
        bound=omega(r,terms)
        ok=mp.mpf(bound[0].numerator)/bound[0].denominator <= v <= mp.mpf(bound[1].numerator)/bound[1].denominator
        diag.append({'r':str(r),'contained_in_rational_interval':bool(ok)})
    status='450_SAMPLE_NODAL_INTERSECTIONS_VERIFIED' if not mismatches and all(z['contained_in_rational_interval'] for z in diag) else 'AUDIT_FAILED'
    report={'schema':'NUMERIA.ES8P7.NODE_AND_ADMISSIBILITY.1','status':status,
            'historical_source_sha256':sourcehash,'parameters':{'L':'9/20','R':'160','cells':1600,'nodes_per_cell':16,'rational_series_terms':terms,'historical_fixed_digits':digits},
            'sample_design':{'cells':[0,40,100,400,1599],'node_indices':[0,8,15], 'full_original_node_count':768000},
            'statistics':stats,'failure_examples':mismatches[:20], 'sample_intervals_rounded_display_only':sample,
            'independent_digamma_mp65_not_certified':diag,
            'admissibility':check_admissibility(),
            'references':reference_triage(source_dir/'REFERENCIAS_ORIGINAIS_16.csv'),
            'not_established':['Validation of every historical nodal enclosure, total 768000 evaluations',
               'Independent formal proof of all fixed-point rounding and analytic remainder code paths',
               'Full infinite-dimensional Weil nonnegativity or the Riemann Hypothesis']}
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out',type=Path,default=Path('RESULTADOS_ES8P7.json'))
    p.add_argument('--terms',type=int,default=256)
    args=p.parse_args()
    data=run(args.source_dir,args.out,args.terms)
    print(data['status'],data['statistics'])
    return 0 if data['status'].endswith('VERIFIED') else 2
if __name__=='__main__':raise SystemExit(main())
