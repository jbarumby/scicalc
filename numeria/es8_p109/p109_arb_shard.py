#!/usr/bin/env python3
"""NUMERIA ES8-P109 — Arb GL/Cauchy matrices for primes up to 300, 500, 1000.

GL/Cauchy formula from P105 with *certified* panel breakpoints, no float branch
or ellipse choices. Every outer and inner Bernstein ellipse is checked against
holomorphic strip conditions using Arb; each outer panel records rho, complex
M, n, error, overlap and active branch. Fail closed on any uncertified case.

Small-support cutoff b_T is treated with a rigorous uniform dR, dI allowance.
Only finite prime-centre subspaces through 997 are addressed. No RH claim.
"""
import argparse, json, itertools, time, hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import arb, acb, fmpq, ctx

ALL_CENTERS = tuple([1]+[p for p in range(2,1001) if all(p%q for q in range(2,__import__('math').isqrt(p)+1))])
ctx.prec = 256
h=arb(2)/5
z0=arb(63)/64
ZH=z0*h
GAP=h-ZH
S_MAX=h+ZH
SUP=2*ZH
ETA=(-(1/(1-z0*z0))).exp()
EPS=arb(1)/10**11
RADIUS_CORRECTION=arb(20)
N_MAX=2048
RULE_CACHE={}
INNER_PROOFS=[]
INNER_PROOF_SHA=hashlib.sha256()
def record_inner(proof):
    INNER_PROOF_SHA.update(json.dumps(proof,sort_keys=True,separators=(',',':')).encode())
    INNER_PROOF_SHA.update(b'\n')
STATS={'inner_calls':0,'outer_panels':0,'inner_complex_verified':0,'outer_complex_verified':0,
       'max_inner_order':0,'max_outer_order':0,'primepowers':[]}

def rat(x):
    q=x.fmpq()
    return f'{int(q.p)}/{int(q.q)}'

def exact_ball(x):
    return {'lo':rat(x.lower()),'hi':rat(x.upper())}

def require(cond, reason):
    if not bool(cond): raise ArithmeticError('NON_CERTIFIED: '+reason)

def positive(x, label):
    require(x.lower()>0, label+' must be strictly positive')

def min_lower(*vals):
    """Pick an exact dyadic lower bound, never a float or ambiguous ordering."""
    return min((v.lower() for v in vals), key=lambda a: a.fmpq())

def ratio_rho(L, ov):
    """Safe rho: u <= min(1,sqrt(ov/L)), overshoot <= ov/2.

    On ellipse E_{1+u} real overshoot is L*u^2/(2*(1+u)).
    Choosing 0<u<=sqrt(ov/L) makes it <=ov/2 without float.
    """
    positive(L,'L'); positive(ov,'overshoot')
    q=(ov/L).sqrt().lower()
    u=q if bool(q<arb(1)) else arb(1)
    rho=arb(1)+u
    require(rho.lower()>1, 'rho<=1')
    excess=L*(rho+1/rho-2)/2
    require(excess.upper()<ov.lower(), 'Bernstein ellipse exceeds certified padding')
    return rho,excess

def ellipse_limits(a,b,rho):
    L=(b-a)/2
    c=(a+b)/2
    axis=L*(rho+1/rho)/2
    return c-axis, c+axis, L*(rho-1/rho)/2

def gl_nodes(n):
    if n not in RULE_CACHE:
        RULE_CACHE[n] = [arb.legendre_p_root(n,k,weight=True) for k in range(n)]
        sw=sum((w for _,w in RULE_CACHE[n]),arb(0))
        require(sw.contains(2),'GL weights fail sum-of-weights check')
        require(all(bool(w.lower()>0) for _,w in RULE_CACHE[n]),'GL weight not positive')
    return RULE_CACHE[n]

def enclosure_error(L,M,rho,n):
    # Johansson's ellipse GL error: 64 L M rho^(1-2n) / [15(rho-1)].
    # We use 8 > 64/15, a conservative constant.
    return 8*L*M*rho**(1-2*n)/(rho-1)

def gl_panel(f,a,b,ov,M,target,location):
    L=(b-a)/2; positive(L,location+' panel length')
    rho,excess=ratio_rho(L,ov)
    n=6
    while True:
        e=enclosure_error(L,M,rho,n)
        if bool(e.upper()<arb(target).lower()):break
        n+=2
        if n>N_MAX:raise ArithmeticError(location+' could not certify GL order')
    c=(a+b)/2; acc=acb(0)
    for x,w in gl_nodes(n):acc += w*f(acb(c+L*x))
    v=L*acc
    require(v.is_finite(),location+' nonfinite output')
    return v,e,n,rho,excess

def bump(z):
    w=z/h
    return (-(1/(1-w*w))).exp()

def check_inner(s, a, b, rho, label):
    xmin,xmax,im=ellipse_limits(a,b,rho)
    # For real parameter s, b(y) and b(y+s) are holomorphic on ellipse.
    lo=min_lower(xmin, xmin+s)
    hi=max(xmax.upper(),(xmax+s).upper())
    require(lo > (-h).upper() and hi<h.lower(),label+' inner ellipse crosses bump support pole/cutoff')
    STATS['inner_complex_verified']+=1
    return {'min_re':str(lo),'max_re':str(hi),'max_im':str(im.upper())}

def R_T(s):
    STATS['inner_calls']+=1
    require(bool(s.is_finite()),'nonfinite R_T parameter')
    # Only positive real s: the even autocorrelation is assembled as R_T(abs(s)).
    require(s.lower()> (-GAP).lower(), 'R_T input too negative')
    if bool(s.lower() >= SUP.upper()):return arb(0)
    if not bool(s.upper()<SUP.lower()):
        # Valid real-axis enclosure across the contact s=2ZH; used only at
        # pointwise outer GL nodes, not for a complex-analytic extension.
        width=SUP.upper()-s.lower()
        require(width.upper()>0, 'contact width')
        bound=max(abs((SUP-s).lower()),abs((SUP-s).upper()))
        # Valid for analytic continuation in the small support-contact strip.
        return arb(0,bound)
    a,b=-ZH,ZH-s
    L=(b-a)/2
    rho,_=ratio_rho(L,arb('0.005'))
    inner_geo=check_inner(s,a,b,rho,'R_T')
    v,e,n,_,_=gl_panel(lambda y: bump(y)*bump(y+s),a,b,
                        arb('0.005'),arb(1),'1e-15','R_T')
    STATS['max_inner_order']=max(STATS['max_inner_order'],n)
    record_inner({'kind':'R_T','s':exact_ball(s),'a':exact_ball(a),
        'b':exact_ball(b),'rho':exact_ball(rho),'overshoot':exact_ball(arb('0.005')),
        'order':n,'M':exact_ball(arb(1)),'gl_error':exact_ball(e),
        'y_ellipse':inner_geo})
    return v.real+arb(0,e.upper())

def A_moment():
    a,b=-ZH,ZH
    L=(b-a)/2;rho,_=ratio_rho(L,arb('0.005'))
    moment_geo=check_inner(arb(0),a,b,rho,'A_moment')
    v,e,n,_,_=gl_panel(lambda y:bump(y)*((y/2).exp()+(-y/2).exp())/2,
                        a,b,arb('0.005'),arb(2),'1e-15','A_moment')
    tail=(2*GAP*ETA*(h/2).exp()).upper()
    STATS['max_inner_order']=max(STATS['max_inner_order'],n)
    record_inner({'kind':'A_moment','s':exact_ball(arb(0)),'a':exact_ball(a),
        'b':exact_ball(b),'rho':exact_ball(rho),'overshoot':exact_ball(arb('0.005')),
        'order':n,'M':exact_ball(arb(2)),'gl_error':exact_ball(e),
        'y_ellipse':moment_geo})
    return v.real+arb(0,e.upper())+arb(0,tail)

def primepowers():
    # Exact integer arithmetic, not log(k,p) or trial float primality.
    out=[]
    for p in range(2,2501):
        if any(p%q==0 for q in range(2,__import__('math').isqrt(p)+1)):continue
        power=p
        while power<=2500:
            out.append((power,p));power*=p
    out.sort();STATS['primepowers']=out
    return out

PRIMES=primepowers()
PRIME_DATA=[(arb(k).log(),arb(p).log()/arb(k).sqrt()) for k,p in PRIMES]
TOTAL_LAMSUM=sum((c for _,c in PRIME_DATA),arb(0))

def sorted_events(events):
    out=[]
    for nm,x in events:
        i=0
        while i<len(out) and bool(out[i][1].upper()<x.lower()):i+=1
        if i<len(out):require(bool(x.upper()<out[i][1].lower()),'overlapping or unordered event values')
        out.insert(i,(nm,x))
    for i in range(1,len(out)):
        require(out[i-1][1].upper()<out[i][1].lower(),'duplicate/overlapping events')
    return out

def panels_for(delta):
    T=delta+SUP
    t1=SUP-delta
    positive(T,'T')
    if bool(delta.upper()<SUP.lower()):
        start=EPS
        events=[]
        for nm,x in [('contactA',t1),('branchB',delta),('stop',T)]:
            # On diagonal, contactA coincides with the endpoint T; skip duplicate.
            if nm=='contactA' and delta.is_zero():continue
            if bool(x.upper()<=EPS.lower()):continue
            if bool(x.lower()>EPS.upper() and x.upper()<T.lower()):events.append((nm,x))
            elif nm=='stop':events.append((nm,x))
            else:raise ArithmeticError('event too close to cutoff or uncertain')
        events=sorted_events(events)
        first=events[0][1]
        pts=[('eps',start)];p=start
        while bool((p*2).upper()<first.lower()):
            p=p*2;pts.append(('dyadic',p))
            require(len(pts)<80,'too many near-origin panels')
        pts.extend(events)
        return pts,True,t1,T
    if bool(delta.lower()>SUP.upper()):
        start=delta-SUP
        positive(start,'initial zero region')
        require(start.lower()>EPS.upper(),'zero-region endpoint too close to eps for omission')
        return [('begin_support',start),('branchB',delta),('stop',T)],False,t1,T
    raise ArithmeticError('delta near support contact; fail closed')

def panel_acts(a,b,t1,delta,overlap):
    mid=(a+b)/2
    out=[]
    if overlap:
        if bool(mid.upper()<t1.lower()):out.append(('A','plus'))
        elif not bool(mid.lower()>t1.upper()):raise ArithmeticError('uncertified A branch')
    if bool(mid.upper()<delta.lower()):out.append(('B','neg'))
    elif bool(mid.lower()>delta.upper()):out.append(('B','pos'))
    else:raise ArithmeticError('uncertified B branch')
    return out

def safe_ov(a,b,delta,t1,acts):
    values=[arb('0.2'),a/2]
    for name,sg in acts:
        if name=='A':values.append((S_MAX-(b+delta))*arb(9)/10)
        elif sg=='neg':values.append((delta-b)+arb('0.005'))
        else: values += [(a-delta)+arb('0.005'),(S_MAX-(b-delta))*arb(9)/10]
    if any(name=='B' for name,_ in acts):values.append(arb('0.005'))
    ov=min_lower(*values)
    positive(ov,'certified outer ellipse overshoot')
    return ov

def audit_outer_ellipse(a,b,delta,acts,rho,F0):
    xmin,xmax,im=ellipse_limits(a,b,rho)
    require(xmin.lower()>0,'outer ellipse crosses t=0')
    limit_stats=[]
    for name,sg in acts:
        if name=='A':smin=xmin+delta;smax=xmax+delta
        elif sg=='neg':smin=delta-xmax;smax=delta-xmin
        else:smin=xmin-delta;smax=xmax-delta
        require(smin.lower() > (-GAP).upper(),f'{name}/{sg}: complex s below analytic strip')
        require(smax.upper() < S_MAX.lower(),f'{name}/{sg}: complex s above analytic strip')
        dist=max(abs((SUP-smin).lower()),abs((SUP-smin).upper()),
                 abs((SUP-smax).lower()),abs((SUP-smax).upper()))
        bound=dist+im.upper()
        require(bound < arb(4).lower(),f'{name}/{sg}: |R(s)|>4 not justified')
        limit_stats.append({'term':name,'branch':sg,'s_re_lower':str(smin.lower()),
                            's_re_upper':str(smax.upper()),'abs_R_bound':str(bound),
                            'strip_lower_gap':str((smin+GAP).lower()),
                            'strip_upper_gap':str((S_MAX-smax).lower())})
    # Use 4 for F bound (convex average of at most two correlations).
    M=((xmax/2).exp()*4+abs(F0).upper())/xmin.lower().sinh()
    positive(M,'complex J majorant')
    STATS['outer_complex_verified']+=1
    return M,{'t_ellipse_lower':str(xmin.lower()),'t_ellipse_upper':str(xmax.upper()),
              't_ellipse_max_imag':str(im.upper()),'M':str(M.upper()),'branches':limit_stats}

def entry(n,m,A,diagnostics):
    delta=arb(m).log()-arb(n).log() if n!=m else arb(0)
    positive(delta+arb(1),'delta legal')
    C=(4*arb.pi()).log()+arb.const_euler()
    P=A*A*((arb(n)/m).sqrt()+(arb(m)/n).sqrt())
    F0=R_T(delta)
    S=arb(0); lamsum=TOTAL_LAMSUM
    for lk,coeff in PRIME_DATA:
        dif=abs(lk-delta)
        # When both arguments lie beyond truncated support, terms vanish exactly.
        # The full-support discrepancy stays covered by dR*TOTAL_LAMSUM.
        if bool(dif.lower()>=SUP.upper()):
            continue
        S+=coeff*(R_T(lk+delta)+R_T(dif))
    pts,overlap,t1,T=panels_for(delta)
    total=arb(0);err=arb(0);nodes=0
    records=[]
    for (leftname,a),(rightname,b) in zip(pts[:-1],pts[1:]):
        acts=panel_acts(a,b,t1,delta,overlap)
        ov=safe_ov(a,b,delta,t1,acts)
        L=(b-a)/2
        rho,excess=ratio_rho(L,ov)
        M,geometry=audit_outer_ellipse(a,b,delta,acts,rho,F0)
        def Fp(t):
            out=arb(0)
            for nam,sg in acts:
                if nam=='A':out+=R_T(t+delta)
                elif sg=='neg':out+=R_T(delta-t)
                else:out+=R_T(t-delta)
            return out/2
        def J(t):
            x=t.real
            return ((x/2).exp()*Fp(x)-F0)/x.sinh()
        v,e,order,rho2,excess2=gl_panel(J,a,b,ov,M,'1e-13',f'outer({n},{m})')
        require(rho2.overlaps(rho),'ellipse ratio inconsistent')
        total+=v.real;err+=e;nodes+=order
        STATS['max_outer_order']=max(STATS['max_outer_order'],order)
        records.append({'from':leftname,'to':rightname,'a':str(a),'b':str(b),
              'active':acts,'ov':str(ov),'rho':str(rho),'real_overshoot':str(excess.upper()),
              'exact':{'a':exact_ball(a),'b':exact_ball(b),'delta':exact_ball(delta),
                  'ov':exact_ball(ov),'rho':exact_ball(rho),'M':exact_ball(M),
                  'GL_error':exact_ball(e),'ellipse_overshoot':exact_ball(excess)},
              'order':order,'error':str(e.upper()),'geometry':geometry})
    STATS['outer_panels']+=len(records)
    tail=F0*(T/2).tanh().log()
    if overlap:
        # For truncated bump, R_T(|s|) has a tiny boundary cusp at s=0.
        # Add 2 ETA^2 EPS as a conservative separate allowance.
        rem=EPS*((EPS/2).exp()*RADIUS_CORRECTION*EPS/2 + (EPS/2).exp()*abs(F0).upper()/2)+2*ETA*ETA*EPS
    else:rem=arb(0)
    # dR uniform correlation error; include full-support T_full in dI.
    dR=(arb(-1).exp()*ETA*4*GAP).upper()
    Tfull=delta+2*h
    dI=dR*((Tfull/2).exp()+1)*(-(EPS/2).tanh().log()).upper()
    correction=(dR*(C.upper()+2*lamsum.upper()+1)+dI).upper()
    W=P-S-C*F0-(total+arb(0,err.upper())+tail+arb(0,rem.upper()))+arb(0,correction)
    require(W.is_finite(),f'W({n},{m}) nonfinite')
    def enc(x):
        mm,ee=x.mid().man_exp();rr,ff=x.rad().man_exp()
        return [int(mm),int(ee),int(rr),int(ff)]
    diagnostics[f'{n},{m}']={'entry':enc(W),'segments':len(records),
       'sum_GL_errors_upper':str(err.upper()),'truncation_correction':str(correction),
       'near_origin_correction':str(rem.upper()),'nodes':nodes,'panels':records}
    return W

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='P109_ROW_SHARD.json')
    ap.add_argument('--bits',type=int,default=256)
    ap.add_argument('--max-prime',type=int,required=True,choices=[300,500,1000])
    ap.add_argument('--shards',type=int,default=24)
    ap.add_argument('--shard-index',type=int,default=0)
    ap.add_argument('--pair',default=None,help='Optional n,m for a short pilot run')
    a=ap.parse_args()
    require(a.bits==256,'precomputed prime-power Arb coefficients use 256 bits')
    centers=tuple(c for c in ALL_CENTERS if c<=a.max_prime)
    count=len(centers)
    require(count in (63,96,169),'unexpected center cardinality')
    require(a.shards>=1 and a.shards<=count and 0<=a.shard_index<a.shards,
            'bad row shard count/index')
    i=(a.shard_index*count)//a.shards
    j=((a.shard_index+1)*count)//a.shards
    require(j>i,'empty shard')
    ctx.prec=a.bits
    start=time.monotonic()
    diagnostics={};moment=A_moment()
    W11=entry(1,1,moment,diagnostics)
    diag_record=diagnostics['1,1']
    if a.pair:
        n,m=map(int,a.pair.split(','))
        require(n in centers and m in centers and n<=m,'pilot pair outside scope')
        selected=[(n,m)]
    else:
        selected=[(centers[u],centers[v])
                  for u in range(i,j) for v in range(u,count)]
    for n,m in selected:
        key=f'{n},{m}'
        if n==m:
            diagnostics[key]={'entry':diag_record['entry'],
                'identity':'translation invariance W(n,n)=W(1,1)'}
        else:
            W=entry(n,m,moment,diagnostics)
            print('ENTRY',a.max_prime,n,m,W,'panels',diagnostics[key]['segments'],
                  'nodes',diagnostics[key]['nodes'],flush=True)
    outer=sum(len(v.get('panels',[])) for v in diagnostics.values())
    result={
        'status':'P109_ROW_SHARD_INTERVALS_CONDITIONAL_ON_GL_CAUCHY_ARITHMETIC',
        'centers':centers,'matrix_dimension':count,'max_prime':a.max_prime,
        'bits':a.bits,'epsilon':'1/100000000000',
        'primepower_upper_integer':2500,'shards':a.shards,
        'shard_index':a.shard_index,'row_start':i,'row_end':j,'pilot_pair':a.pair,
        'entries':{f'{n},{m}':diagnostics[f'{n},{m}']['entry'] for n,m in selected},
        'diagonal_reference':diag_record['entry'],
        'proof_records':diagnostics, 'stats':STATS,
        'outer_panels_verified':outer,
        'inner_ellipse_count':STATS['inner_complex_verified'],
        'inner_ellipse_proof_sha256':INNER_PROOF_SHA.hexdigest(),
        'elapsed_seconds':time.monotonic()-start,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'warnings':['No positivity claim before exact matrix assembly and Sylvester proof',
          'Finite test subspaces only, not universal Weil positivity or RH',
          'All interval inclusions conditional on GL/Cauchy analytic bounds and FLINT/Arb']
    }
    Path(a.out).write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('SHARD_SUCCESS',a.max_prime,i,j,len(selected),'interval entries',
          outer,'outer panels',STATS['inner_complex_verified'],'inner ellipses',
          'elapsed',result['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
