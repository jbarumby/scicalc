"""P92: exact four-gap sharp Preissmann Lemma 6 certificate.

No symbolic CAS, floats, or numerical tolerances: all halfspaces, vertices,
facet angular ordering, tetra volumes, and coefficients use Fraction.
The full mathematical reduction is in README.md.

This certifies the P89/P90 sufficient interior budget for m=4 gaps only;
it does NOT prove the general Preissmann L5/L6 or Riemann hypothesis.
"""
import argparse
import hashlib
import json
from functools import cmp_to_key
from itertools import combinations, product
from fractions import Fraction as Q
from pathlib import Path

SCHEMA = "NUMERIA.P92.PREISSMANN.FOUR_GAPS.EXACT_POLYHEDRA.1"
ZERO = Q(0)
ONE = Q(1)
ENDPOINT_A = (Q(0), Q(1, 4), Q(1, 3), Q(1, 2))
ENDPOINT_D = (Q(0), Q(1, 4), Q(1, 3), Q(1, 2), Q(1))
ENVELOPES = (
    (Q(3, 2), Q(139, 36)),
    (Q(19, 9), Q(17, 12)),
    (Q(9, 4), Q(1)),
    (Q(5, 2), Q(1, 2)),
)

def cross(x, y):
    return (x[1]*y[2]-x[2]*y[1],
            x[2]*y[0]-x[0]*y[2],
            x[0]*y[1]-x[1]*y[0])

def sub(x, y):
    return tuple(a-b for a,b in zip(x,y))

def dot(plane, x):
    return sum((plane[i]*x[i] for i in range(3)), plane[3])

def det3(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])

def intersection(p,q,r):
    a,b,c=p[:3],q[:3],r[:3]
    determinant=det3(a,b,c)
    if not determinant:
        return None
    rhs=(-p[3],-q[3],-r[3])
    bc,ca,ab=cross(b,c),cross(c,a),cross(a,b)
    return tuple((rhs[0]*bc[i]+rhs[1]*ca[i]+rhs[2]*ab[i])/determinant
                 for i in range(3))

def region_planes(bits, ia, jd):
    """Halfspaces A*a+B*b+D*d+E>=0, c=1-a-b-d."""
    al,au=ENDPOINT_A[ia:ia+2]
    dl,du=ENDPOINT_D[jd:jd+2]
    rows=[
        (1,0,0,0), (0,1,0,0), (0,0,1,0),
        (-1,-1,-1,1),       # c>=0
        (-1,0,1,0),         # a<=d reflection
        (1,0,0,-al),(-1,0,0,au),
        (0,0,1,-dl),(0,0,-1,du),
        (-1,1,0,0) if bits[0]==0 else (1,-1,0,0),
        (-1,-2,-1,1) if bits[1]==0 else (1,2,1,-1),
        (1,1,2,-1) if bits[2]==0 else (-1,-1,-2,1),
    ]
    return tuple(tuple(Q(y) for y in p) for p in rows)

def region_vertices(planes):
    vertices=set()
    for triple in combinations(planes,3):
        x=intersection(*triple)
        if x is not None and all(dot(p,x)>=0 for p in planes):
            vertices.add(x)
    return tuple(sorted(vertices))

def facet_vertices(planes, vertices):
    """Exact circular ordering of convex face vertices, NO atan2/float."""
    seen=set()
    for p in planes:
        pts=tuple(v for v in vertices if dot(p,v)==0)
        if len(pts)<3 or frozenset(pts) in seen:
            continue
        seen.add(frozenset(pts))
        omit=max(range(3),key=lambda i:abs(p[i]))
        axes=tuple(i for i in range(3) if i!=omit)
        mean=tuple(sum((v[i] for v in pts),ZERO)/len(pts)
                   for i in range(3))
        def compare(x,y):
            xx,xy=(x[i]-mean[i] for i in axes)
            yx,yy=(y[i]-mean[i] for i in axes)
            hx=0 if xy>0 or (xy==0 and xx>=0) else 1
            hy=0 if yy>0 or (yy==0 and yx>=0) else 1
            if hx!=hy:
                return -1 if hx<hy else 1
            z=xx*yy-xy*yx
            if z:
                return -1 if z>0 else 1
            r1=xx*xx+xy*xy
            r2=yx*yx+yy*yy
            return (r1>r2)-(r1<r2)
        poly=tuple(sorted(pts,key=cmp_to_key(compare)))
        corners=[]
        for i in range(len(poly)):
            u=sub(poly[(i+1)%len(poly)],poly[i])
            v=sub(poly[(i+2)%len(poly)],poly[(i+1)%len(poly)])
            corners.append(cross(u,v)[omit])
        if not (any(x!=0 for x in corners) and
                (all(x>=0 for x in corners) or all(x<=0 for x in corners))):
            raise ArithmeticError("invalid rational polygon facet ordering")
        yield poly

def tetrahedra(planes, vertices):
    if len(vertices)<4:
        return
    center=tuple(sum((v[i] for v in vertices),ZERO)/len(vertices)
                 for i in range(3))
    if not all(dot(p,center)>0 or all(dot(p,v)>0 for v in vertices)
               for p in planes):
        # Some redundant inequalities are strict on the whole polyhedron.
        raise ArithmeticError("barycenter not in full-dimensional region")
    for face in facet_vertices(planes,vertices):
        for i in range(1,len(face)-1):
            tet=(center,face[0],face[i],face[i+1])
            volume=abs(det3(sub(tet[1],center),sub(tet[2],center),sub(tet[3],center)))/6
            if volume>0:
                yield tet,volume

# Sparse 4-variable polynomial in barycentric U,V,W,X.
UNIT={(0,0,0,0): ONE}

def padd(*args):
    out={}
    for p in args:
        for exponent,coeff in p.items():
            out[exponent]=out.get(exponent,ZERO)+coeff
    return {k:v for k,v in out.items() if v}

def pscale(p,n):
    return {k:v*n for k,v in p.items() if v*n}

def pmul(a,b):
    out={}
    for k,v in a.items():
        for l,w in b.items():
            e=tuple(k[i]+l[i] for i in range(4))
            out[e]=out.get(e,ZERO)+v*w
    return {k:v for k,v in out.items() if v}

def pprod(*args):
    out=UNIT
    for p in args:
        out=pmul(out,p)
    return out

def ppow(p,n):
    if type(n) is not int or n<0:
        raise ValueError("nonnegative exponent required")
    out=UNIT
    for _ in range(n):
        out=pmul(out,p)
    return out

def plinear(vals):
    return {tuple(int(i==j) for i in range(4)): Q(v)
            for j,v in enumerate(vals) if v}

def peval(p, vals):
    return sum((c*Q(vals[0])**e[0]*Q(vals[1])**e[1]
                *Q(vals[2])**e[2]*Q(vals[3])**e[3]
                for e,c in p.items()),ZERO)

def numerator(bits, ia, jd, tet):
    """H^13 times a^2*d^2*(1-a)^2*(1-d)^2*(a+b)^2*(c+d)^2
    times (d4(a)+d4(d)-M4), in rational barycentric coordinates.
    """
    H=plinear((1,1,1,1))
    A=plinear(v[0] for v in tet)
    B=plinear(v[1] for v in tet)
    D=plinear(v[2] for v in tet)
    C=padd(H,pscale(A,-1),pscale(B,-1),pscale(D,-1))
    L=padd(H,pscale(A,-1))
    R=padd(H,pscale(D,-1))
    T=padd(A,B)
    V=padd(C,D)
    m1=A if bits[0]==0 else B
    m2=B if bits[1]==0 else C
    m3=C if bits[2]==0 else D
    alpha,beta=ENVELOPES[ia]
    gamma,eta=ENVELOPES[jd]
    a2,d2,l2,r2,t2,v2=(ppow(p,2) for p in (A,D,L,R,T,V))
    h2,h4=ppow(H,2),ppow(H,4)
    common=pprod(a2,d2,l2,r2,t2,v2)
    poly=padd(
        pscale(pprod(h2,A,d2,l2,r2,t2,v2),alpha),
        pscale(pprod(h2,a2,D,l2,r2,t2,v2),gamma),
        pscale(pprod(H,common),beta+eta),
        pscale(pprod(padd(A,D),common),-3),
        pscale(pprod(h4,m1,d2,r2,t2,v2),-1),
        pscale(pprod(h4,m2,a2,d2,l2,r2),-1),
        pscale(pprod(h4,m3,a2,l2,t2,v2),-1),
    )
    if any(sum(e)!=13 for e in poly):
        raise ArithmeticError("nonhomogeneous degree-13 polynomial")
    return poly

def d4(s):
    s=Q(s)
    if not ZERO<s<=ONE:
        raise ValueError("0<s<=1 required")
    j=0 if s<=Q(1,4) else 1 if s<=Q(1,3) else 2 if s<=Q(1,2) else 3
    alpha,beta=ENVELOPES[j]
    return alpha/s+beta-3*s

def interior(g):
    a,b,c,d=map(Q,g)
    if not all(x>0 for x in (a,b,c,d)) or a+b+c+d!=ONE:
        raise ValueError("four positive gaps summing to one required")
    return (min(a,b)/(a*a*(1-a)**2)
            +min(b,c)/((a+b)**2*(c+d)**2)
            +min(c,d)/(d*d*(1-d)**2))

def check_barycentric_identity(poly, tet):
    weights=(1,2,3,4)
    h=sum(weights)
    a=sum((weights[i]*tet[i][0] for i in range(4)),ZERO)/h
    b=sum((weights[i]*tet[i][1] for i in range(4)),ZERO)/h
    d=sum((weights[i]*tet[i][2] for i in range(4)),ZERO)/h
    c=1-a-b-d
    factor=a*a*d*d*(1-a)**2*(1-d)**2*(a+b)**2*(c+d)**2
    expected=h**13*factor*(d4(a)+d4(d)-interior((a,b,c,d)))
    if peval(poly,weights)!=expected or expected<0:
        raise ArithmeticError("barycentric numerator not equal to original budget")

def poly_digest(poly):
    data="\n".join(f"{','.join(map(str,e))}:{c}"
                   for e,c in sorted(poly.items()))
    return hashlib.sha256(data.encode("ascii")).hexdigest()

def certify(full=False):
    regions=0
    total_tets=0
    total_positive=0
    total_vol=ZERO
    entries=[]
    for bits in product((0,1),repeat=3):
        for ia in range(3):
            for jd in range(4):
                planes=region_planes(bits,ia,jd)
                vertices=region_vertices(planes)
                tets=list(tetrahedra(planes,vertices))
                if not tets:
                    continue
                regions+=1
                records=[]
                cell_volume=ZERO
                for tet,volume in tets:
                    if volume<=0 or not all(
                        all(dot(p,v)>=0 for p in planes) for v in tet):
                        raise ArithmeticError("tetra outside its declared halfspaces")
                    poly=numerator(bits,ia,jd,tet)
                    if not poly or any(v<=0 for v in poly.values()):
                        raise ArithmeticError("nonpositive polynomial coefficient")
                    check_barycentric_identity(poly,tet)
                    cell_volume+=volume
                    total_tets+=1
                    total_positive+=len(poly)
                    item={"vertices_a_b_d":[[str(x) for x in v] for v in tet],
                          "volume":str(volume),"positive_monomials":len(poly),
                          "sha256_coefficients":poly_digest(poly)}
                    if full:
                        item["coefficients"]={
                            ",".join(map(str,e)):str(c)
                            for e,c in sorted(poly.items())}
                    records.append(item)
                total_vol+=cell_volume
                entries.append({"min_branches":"".join(map(str,bits)),
                                "a_interval":ia,"d_interval":jd,
                                "region_vertices":[[str(x) for x in v]
                                                   for v in vertices],
                                "volume":str(cell_volume),
                                "tetrahedra":records})
    if (regions,total_tets,total_positive,total_vol)!=(35,246,133363,Q(1,12)):
        raise ArithmeticError("incomplete half-simplex or changed polynomial census")
    if d4(Q(1,4))*2!=interior((Q(1,4),)*4):
        raise ArithmeticError("uniform four-gap equality lost")
    return {"schema":SCHEMA,"region_count":regions,"tetra_count":total_tets,
            "positive_exact_coefficients":total_positive,
            "negative_coefficients":0,"half_simplex_volume":str(total_vol),
            "all_four_gap_budget_proved":True,
            "all_m_ge_5_budget_proved":False,
            "full_Preissmann_Lemma5_proved":False,
            "full_Preissmann_Lemma6_proved":False,
            "Riemann_Hypothesis_proved":False,
            "uniform_gaps":"1/4,1/4,1/4,1/4",
            "uniform_M_and_two_D":str(interior((Q(1,4),)*4)),
            "regions":entries}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--full",action="store_true",
                    help="include every exact rational coefficient")
    args=ap.parse_args()
    result=certify(full=args.full)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print("P92 PASS: 35 regions, 246 tetrahedra, 133363 positive rational coefficients")
    print("Scope: sharp L6 for index separation <=4; general L6 and RH NOT proved.")

if __name__=="__main__":
    main()
