"""P91: exact, dependency-free rational polynomial certificate for 3 gaps.

Theorem (read README): for any a,b,c>0, a+b+c=1, the P90
interior budget M(a,b,c) <= D(a)+D(c). This only covers m=3 gaps.
No floating point or sympy is used in the certificate.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

SCHEMA = "NUMERIA.P91.PREISSMANN.THREE_GAPS.EXACT_COEFFICIENTS.1"
ZERO = F(0)
ONE = F(1)

# Each entry is (order, left envelope, right envelope, three vertices).
# Vertices are (a,c); b=1-a-c. The nine triangles partition
# {0 <= a <= c, a+c <= 1}, which has area 1/4.
TRIANGLES = (
    ("abc", "d1", "d2", (("1/3", "1/3"), ("1/4", "1/2"), ("0", "1/2"))),
    ("abc", "d1", "d3", (("0", "1/2"), ("1/4", "1/2"), ("0", "1"))),
    ("acb", "d1", "d1", (("0", "0"), ("1/3", "1/3"), ("0", "1/3"))),
    ("acb", "d1", "d2", (("0", "1/3"), ("1/3", "1/3"), ("0", "1/2"))),
    ("bac", "d1", "d2", (("1/3", "1/3"), ("1/3", "1/2"), ("1/4", "1/2"))),
    ("bac", "d1", "d3", (("1/4", "1/2"), ("1/3", "1/2"), ("1/3", "2/3"))),
    ("bac", "d1", "d3", (("1/4", "1/2"), ("1/3", "2/3"), ("0", "1"))),
    ("bac", "d2", "d2", (("1/3", "1/3"), ("1/2", "1/2"), ("1/3", "1/2"))),
    ("bac", "d2", "d3", (("1/3", "1/2"), ("1/2", "1/2"), ("1/3", "2/3"))),
)

# D(s) >= alpha/s + beta - 3s (0<s<=1) on its indicated interval.
# Derived from convex chords of S(s) at 0,1/3,1/2,1 and Z>3/2.
ENVELOPES = {"d1": (F(3, 2), F(13, 4)),
             "d2": (F(9, 4), F(1)),
             "d3": (F(5, 2), F(1, 2))}
RANGES = {"d1": (F(0), F(1, 3)),
          "d2": (F(1, 3), F(1, 2)),
          "d3": (F(1, 2), F(1))}

# Sparse polynomials in three barycentric variables U,V,W.
# Dict maps exponents (i,j,k) to Fraction coefficients.
def const(q):
    q = F(q)
    return {(0, 0, 0): q} if q else {}

def add(*args):
    ans = {}
    for p in args:
        for power, coeff in p.items():
            ans[power] = ans.get(power, ZERO) + coeff
            if not ans[power]:
                del ans[power]
    return ans

def scale(p, q):
    q = F(q)
    return {k: q*v for k, v in p.items() if q*v}

def mul(*args):
    ans = const(1)
    for p in args:
        nxt = {}
        for (i, j, k), a in ans.items():
            for (l, m, n), b in p.items():
                power = (i+l, j+m, k+n)
                nxt[power] = nxt.get(power, ZERO) + a*b
        ans = {k: v for k, v in nxt.items() if v}
    return ans

def power(p, n):
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative integer exponent required")
    ans = const(1)
    for _ in range(n):
        ans = mul(ans, p)
    return ans

def linear(values):
    return {exp: F(v) for exp, v in
            zip(((1, 0, 0), (0, 1, 0), (0, 0, 1)), values) if F(v)}

def evaluate(p, values):
    u, v, w = map(F, values)
    return sum((coeff*u**i*v**j*w**k for (i,j,k),coeff in p.items()), ZERO)

def triangle_area(vertices):
    (a,c),(d,e),(f,g) = vertices
    return abs((d-a)*(g-c)-(f-a)*(e-c))/2

def in_cell(a, c, order, left, right):
    b = 1-a-c
    if not(0 <= a <= c and b >= 0):
        return False
    conditions = {"abc": a <= b <= c, "acb": a <= c <= b,
                  "bac": b <= a <= c}
    l, h = RANGES[left]
    lo, hi = RANGES[right]
    return conditions[order] and l <= a <= h and lo <= c <= hi

def envelope(s):
    if s <= 0 or s > 1:
        raise ValueError("envelope defined for 0 < s <= 1")
    which = "d1" if s <= F(1,3) else "d2" if s <= F(1,2) else "d3"
    alpha, beta = ENVELOPES[which]
    return alpha/s + beta - 3*s

def middle(a,c):
    b = 1-a-c
    if not(a>0 and b>0 and c>0):
        raise ValueError("positive gaps a,b,c required")
    return min(a,b)/(a*a*(1-a)**2) + min(b,c)/(c*c*(1-c)**2)

def homogeneous_numerator(order, left, right, vertices):
    """H^9 * 4a²c²(1-a)²(1-c)² * (d(a)+d(c)-M).

    Here H=U+V+W, a=A/H, c=C/H, b=B/H. Every denominator
    multiplied away is strictly positive in the open simplex.
    """
    H = linear((1,1,1))
    A = linear([p[0] for p in vertices])
    C = linear([p[1] for p in vertices])
    B = add(H, scale(A,-1), scale(C,-1))
    L = add(H,scale(A,-1))
    R = add(H,scale(C,-1))
    amin = A if order != "bac" else B
    cmin = C if order == "acb" else B
    al,bl = ENVELOPES[left]
    ar,br = ENVELOPES[right]
    a2,c2,L2,R2,H2,H4 = (power(p,n) for p,n in
         ((A,2),(C,2),(L,2),(R,2),(H,2),(H,4)))
    common = mul(L2,R2)
    p = add(
        scale(mul(A,c2,common,H2),4*al),
        scale(mul(a2,C,common,H2),4*ar),
        scale(mul(a2,c2,common,H),4*(bl+br)),
        scale(mul(power(A,3),c2,common),-12),
        scale(mul(a2,power(C,3),common),-12),
        scale(mul(c2,R2,amin,H4),-4),
        scale(mul(a2,L2,cmin,H4),-4),
    )
    if any(sum(k)!=9 for k in p):
        raise ArithmeticError("certificate numerator is not homogeneous degree 9")
    return p

def certify():
    area = ZERO
    rows = []
    for idx,(order,left,right,raw_vertices) in enumerate(TRIANGLES,1):
        vertices = tuple(tuple(map(F,p)) for p in raw_vertices)
        if not all(in_cell(a,c,order,left,right) for a,c in vertices):
            raise ArithmeticError(f"triangle {idx} is outside its declared cell")
        area += triangle_area(vertices)
        poly = homogeneous_numerator(order,left,right,vertices)
        negative = [(exp,str(coeff)) for exp,coeff in poly.items() if coeff<0]
        if negative:
            raise ArithmeticError(f"negative coefficient on triangle {idx}: {negative[:3]}")
        for weights in ((1,1,1),(1,2,3),(3,2,1)):
            den = sum(weights)
            a = sum(F(weights[k])*vertices[k][0] for k in range(3))/den
            c = sum(F(weights[k])*vertices[k][1] for k in range(3))/den
            rhs = (4*a*a*c*c*(1-a)**2*(1-c)**2
                   *(envelope(a)+envelope(c)-middle(a,c)))
            lhs = evaluate(poly, weights)/den**9
            if lhs != rhs:
                raise ArithmeticError(f"bad numerator algebra on triangle {idx}")
        rows.append({"triangle":idx,"gap_order":order,
                     "endpoint_envelopes":[left,right],
                     "vertices_a_c": [[str(x),str(y)] for x,y in vertices],
                     "nonzero_coefficient_count":len(poly),
                     "negative_coefficient_count":0,
                     "smallest_nonzero_coefficient":str(min(poly.values())),
                     "area":str(triangle_area(vertices)),
                     "coefficients": {f"{i},{j},{k}":str(v)
                                      for (i,j,k),v in sorted(poly.items())}})
    if area != F(1,4):
        raise ArithmeticError(f"triangles do not cover required area 1/4: {area}")
    if len(rows)!=9 or sum(x["nonzero_coefficient_count"] for x in rows)!=447:
        raise ArithmeticError("unexpected P91 certificate counts")
    return {"schema":SCHEMA,"theorem":"3 interior gaps (2 interior nodes): all real positive spacings",
            "three_gap_budget_proved":True,"four_or_more_gaps_proved":False,
            "full_Preissmann_Lemma5_reconstructed":False,
            "full_Preissmann_Lemma6_reproved":False,
            "Riemann_Hypothesis_proved":False,
            "reflection_reduces_to_a_le_c":True,
            "triangles_cover_a_le_c_region_area":"1/4",
            "nonzero_coefficients_total":447,
            "all_coefficients_nonnegative":True,"triangles":rows}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args=p.parse_args()
    result=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("P91: verified",len(result["triangles"]),"triangles,",
          result["nonzero_coefficients_total"],"nonnegative rational coefficients")
    print("Full L6 and RH are NOT proved")

if __name__ == "__main__":
    main()
