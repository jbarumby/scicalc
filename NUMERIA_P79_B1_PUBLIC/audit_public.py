"""NUMERIA P79: conditional public-source replay of P77 B1 with P78 b4 caps.

Preserves P77 source numerics and all earlier certificates; replays P75 in memory.
Do not interpret any finite-window verdict as a global RH proof.
"""
"""Rigorous Arb audit of Chourasiya--Simonic v2 Table 1 B1.

Integral K1: substitute u=x^6; integrate x in [0,2] using acb.integral
and an explicit elementary bound for u>=64. No numerical point samples
are accepted as a certificate.
"""
import argparse
import json
import hashlib
from fractions import Fraction
from pathlib import Path
from flint import acb, arb, ctx
from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q, dyadic
from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table

T0=3*10**12
CUT=64
SOURCE="https://arxiv.org/html/2507.15184v2"

def rat(n,d):
    return arb(n)/arb(d)

def pw(x,e):
    if not x>0:
        raise ArithmeticError("non-positive real power base")
    return (x.log()*e).exp()

def tail_upper(c,c1,c2,L):
    """For u>=64: integrand <=D*exp(-u)*u^3, integrate exactly."""
    if not all(z>0 for z in (c,c1,c2,L)):
        raise ArithmeticError("tail requires positive constants")
    U=arb(CUT)
    D=(1+c/(2*L))**2*(1+c1*c.sqrt()+c2)
    t=(D*(-U).exp()*(U**3+3*U**2+6*U+6)).upper()
    if not t>0 or not t.is_finite():
        raise ArithmeticError("invalid tail")
    return t

def k1_integral(c,c1,c2,L):
    """Rigorous acb.integral enclosure on [0,64] plus strict tail bound."""
    def f(x,analytic):
        v=x**6
        A=acb(1)+(acb(1)+acb(c)*v).log(analytic=analytic)/(2*acb(L))
        B=acb(c1*c.sqrt())*x**3+acb(c2)
        return 6*v*(-v).exp()*A.pow(acb(4)/3,analytic=analytic)*B.pow(acb(2)/3,analytic=analytic)
    I=acb.integral(f,0,2,rel_tol=arb('1e-19'),
                   abs_tol=arb('1e-19'),eval_limit=30000,depth_limit=30)
    if not I.is_finite() or not I.imag.contains(0) or not I.real.lower()>0:
        raise ArithmeticError("quadrature failed or inconclusive")
    tail=tail_upper(c,c1,c2,L)
    return I.real.lower(),(I.real.upper()+tail).upper(),tail

def row_audit(row,cs,b4):
    """Reconstruct every factor in Theorem 1 and (10), (35)--(40)."""
    d,kap,delta=q(row.d),q(row.kappa),q(row.delta)
    s1,s2=rat(row.sigma_start,32),rat(row.sigma_start+1,32)
    L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
    # Choose EXACT real A0=10^printed log10 A0; no unknown optimizer digits.
    A0=(arb(10).log()*q(row.log10_a0)).exp()
    if not arb(10)**10<A0<arb(T0):
        raise ArithmeticError("invalid A0")
    La=A0.log()
    half=(A0/2).log()
    a1=A0/(2*pi*pi)*half**4+q("34.842")*A0*half**3+arb("1.3897e16")
    r=1-arb.const_log2()/La
    a2=r**4/(2*pi*pi)+q("34.842")/La*r**3+arb("9.5079e15")/(A0*La**4)
    if not a1>0 or not a2>0:raise ArithmeticError("invalid a1/a2")
    c=3/(4*kap)
    lower,upper,tail=k1_integral(c,cs["c1"],cs["c2"],L)
    F=pw(2*a2*c.sqrt(),rat(1,3))
    correction=pw(2*a1,rat(1,3))*pw(cs["c2"]+A0*cs["c1"]/T0,rat(2,3))/(pw(arb(T0),rat(1,3))*pw(L,rat(4,3)))
    k1lo=(F*lower+correction).lower()
    k1=(F*upper+correction).upper()
    if not k1lo>0:raise ArithmeticError("K1 lower bound inconclusive")
    two13=pw(arb(2),rat(1,3))
    k2=two13*c.sqrt()*pi.sqrt()/2+two13*(4*kap/3)*pw(A0/arb(T0),arb(3))
    c3=pi*m0*b4/(2*delta)*(2*delta*arb.const_euler()/L).exp()
    c4=b4/(5*delta*delta.exp())+q(".605")*(-2*delta).exp()/L
    k3=c3+c4*(pi.sqrt()/(2*(2*kap).sqrt())+pi*m0/T0)
    l0=(kap*(1+(1/(4*arb(T0))*(1+2*delta/L))**2)).exp()/(1-arb(1)/T0)
    l1=max(arb(1),pw(k1,3*delta/L))*max(arb(1),pw(k3,-3*delta/L))*(3*delta+4*delta*L.log()/L).exp()
    l2=max(arb(1),pw(k2,3*delta/L))*max(arb(1),pw(k3,-3*delta/L))*(3*delta).exp()
    common=max(arb(1),pw(l0,-2*d/L))*max(arb(1),pw(k3,-3*d/L))
    common*=max(pw(l0,2/(2-s1)),pw(l0,2/(2-s2)))
    common*=max(pw(k3,(2*s1-1)/(2-s1)),pw(k3,(2*s2-1)/(2-s2)))
    A1=l1*(d*(3+4*L.log()/L)).exp()*common
    A1*=max(arb(1),pw(k1,3*d/L))
    A1*=max(pw(k1,3*(1-s1)/(2-s1)),pw(k1,3*(1-s2)/(2-s2)))
    A2=l2*(3*d).exp()*common
    A2*=max(arb(1),pw(k2,3*d/L))
    A2*=max(pw(k2,3*(1-s1)/(2-s1)),pw(k2,3*(1-s2)/(2-s2)))
    B1=(A1+A2*pw(L,-6*(1-s2)/(2-s2)))/(2*pi*d)
    upper_B1=B1.upper()
    margin=(q(row.b1)-upper_B1).lower()
    return {
        "sigma":[f"{row.sigma_start}/32",f"{row.sigma_start+1}/32"],
        "printed_B1":row.b1,
        "exact_log10_A0_choice":row.log10_a0,
        "k1_integral_lower_dyadic":dyadic(lower),
        "k1_integral_upper_dyadic":dyadic(upper),
        "k1_tail_upper_dyadic":dyadic(tail),
        "k1_lower_dyadic":dyadic(k1lo),
        "k1_upper_dyadic":dyadic(k1),
        "B1_upper_dyadic":dyadic(upper_B1),
        "margin_lower_dyadic":dyadic(margin),
        "strict":bool(margin>0)
    }


# P79: source-preserving overlay. Imported P77 source above is copied unchanged
# except for explicitly parameterizing the two occurrences of b4 in C3/C4.
P78_CERT = Path(__file__).resolve().parents[1] / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json"
P79_SCHEMA = "NUMERIA.P79.B1.PUBLIC_B4.1"

def exact_from_dyadic(pair):
    """Convert a *directed endpoint* (m, e) into an exact rational m*2**e."""
    m, e = int(pair[0]), int(pair[1])
    return Fraction(m) * (Fraction(2) ** e)

def exact_ball(fr):
    fr = Fraction(fr)
    return arb(fr.numerator) / arb(fr.denominator)

def rounded_up_dyadic(pair, places=6):
    """Return an EXACT decimal string strictly ABOVE a certified upper endpoint."""
    if type(places) is not int or places < 0 or places > 12:
        raise ValueError("unsupported decimal precision")
    x = exact_from_dyadic(pair)
    if x < 0:
        raise ArithmeticError("negative B1 bound")
    scale = 10**places
    scaled = x * scale
    ceil_scaled = (scaled.numerator + scaled.denominator - 1) // scaled.denominator
    # Extra last-place unit gives strict slack even when x lies on the grid.
    n = ceil_scaled + 1
    if places == 0:
        return str(n)
    return f"{n // scale}.{n % scale:0{places}d}"

def load_source():
    raw = P78_CERT.read_bytes()
    cert = json.loads(raw)
    if cert.get("schema") != "NUMERIA.P78.RESTRICTED_B4.1" or len(cert.get("rows",[])) != 16:
        raise ArithmeticError("unrecognized P78 source certificate")
    rows = cert["rows"]
    for idx, (x, row) in enumerate(zip(rows,TABLE),1):
        if x.get("row") != idx or Fraction(x["delta"]) != Fraction(row.delta):
            raise ArithmeticError(f"P78 row or delta mismatch: {idx}")
        b4 = Fraction(x["b4"])
        if not (Fraction("0.54") <= b4 <= Fraction("0.606")):
            raise ArithmeticError("P78 b4 outside audited range")
        if b4 <= Fraction(x["worst_upper"]):
            raise ArithmeticError("P78 rational cap does not dominate source")
        if (idx <= 9) != (b4 == Fraction("0.54")):
            raise ArithmeticError("unexpected P78 b4 switch")
    return rows,hashlib.sha256(raw).hexdigest()

def audit(precision=256):
    if precision not in (192,256,320):
        raise ValueError("precision must be 192, 256 or 320")
    validate_table()
    source_rows, source_sha = load_source()
    with ctx.workprec(precision):
        cs = constants()
        L,pi,m0 = cs["logh"],cs["pi"],cs["m0"]
        cs["c1"] = 6/pi**2 + q("1.048")/L
        cs["c2"] = pi*m0*q(".62")/L + 6*m0/(pi*T0) + pi*m0*q("1.048")/(arb(T0)*L)
        checked = []
        table = []
        for index,row in enumerate(TABLE,1):
            b4_fraction = Fraction(source_rows[index-1]["b4"])
            out = row_audit(row, cs, exact_ball(b4_fraction))
            upper = exact_from_dyadic(out["B1_upper_dyadic"])
            old_ok = bool(out["strict"])
            # Rows 1..9 may preserve the previously certified printed B1.
            # All replaced rows receive outward rational rounding from Arb.
            replacement = row.b1 if index <= 9 else rounded_up_dyadic(out["B1_upper_dyadic"])
            new_margin = (q(replacement)-exact_ball(upper)).lower()
            if not new_margin > 0:
                raise ArithmeticError(f"P79 new B1 not strictly above the Arb endpoint: row {index}")
            if index <= 9 and not old_ok:
                raise ArithmeticError(f"P79 original B1 no longer valid in first nine rows: {index}")
            if index >= 10 and old_ok:
                # Deliberately fail closed if the adversarial change no longer matches the supplied prototype.
                raise ArithmeticError(f"P79 unexpected old B1 survival in row {index}")
            checked.append({
                "row":index,
                "sigma":out["sigma"],
                "b4_source":str(b4_fraction),
                "original_B1":row.b1,
                "old_B1_still_strict":old_ok,
                "B1_upper_dyadic":out["B1_upper_dyadic"],
                "B1_replacement":replacement,
                "replacement_margin_lower_dyadic":dyadic(new_margin),
                "K1_upper_dyadic":out["k1_upper_dyadic"],
                "K1_tail_upper_dyadic":out["k1_tail_upper_dyadic"],
            })
            table.append([replacement,row.b2,row.b3])
        if len(checked)!=16:
            raise ArithmeticError("incomplete P79 table")
        table_json=json.dumps(table,separators=(",",":"),ensure_ascii=True)
        return {
            "schema":P79_SCHEMA,
            "source_P78_certificate_sha256":source_sha,
            "source_P77_method":"NUMERIA_P77_B1.audit_b1.py; two C3/C4 b4 occurrences parametrized",
            "precision_bits":precision,
            "coefficients":checked,
            "table_v2_public_replacement":table,
            "table_sha256":hashlib.sha256(table_json.encode("ascii")).hexdigest(),
            "all_16_strict_new_B1":True,
            "original_printed_B1_fails_rows":list(range(10,17)),
            "external_lemma6_fully_proved":False,
            "global_RH_proved":False,
        }

def replay(start="18",end="79",precision=256,segments=256):
    """No change to P75 sources or their frozen certificates: override in memory."""
    from NUMERIA_P75_CS2 import cs_v2
    from NUMERIA_P75_CS2.certify import run
    coeff = audit(precision)
    original = cs_v2.TABLE_V2
    try:
        cs_v2.TABLE_V2 = tuple(tuple(r) for r in coeff["table_v2_public_replacement"])
        if len(cs_v2.TABLE_V2)!=16:
            raise ArithmeticError("invalid replacement table")
        cover = run(end=end,start=start,count=500,prec=precision,segments=segments,mesh="1/64")
    finally:
        cs_v2.TABLE_V2 = original
    cover["P79_public_B1_table_sha256"]=coeff["table_sha256"]
    cover["P79_source_P78_sha256"]=coeff["source_P78_certificate_sha256"]
    cover["publication_status"]="PRELIMINARY P79: conditional on public-source theorem and all imported analytic gates"
    return cover

def main():
    parser=argparse.ArgumentParser(description="P79 public b4 -> certified outward B1 -> optional P75 finite replay")
    parser.add_argument("--precision",type=int,default=256)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--cover-out",type=Path)
    parser.add_argument("--start",default="18")
    parser.add_argument("--end",default="79")
    parser.add_argument("--segments",type=int,default=256)
    args=parser.parse_args()
    result=audit(args.precision)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"schema":result["schema"],"precision_bits":args.precision,
        "table_sha256":result["table_sha256"],
        "original_printed_B1_fails_rows":result["original_printed_B1_fails_rows"],
        "B1_replacements":[r[0] for r in result["table_v2_public_replacement"]]},indent=2))
    if args.cover_out:
        if args.precision!=256:
            raise ValueError("full replay requires 256-bit precision")
        cover=replay(args.start,args.end,args.precision,args.segments)
        args.cover_out.parent.mkdir(parents=True,exist_ok=True)
        args.cover_out.write_text(json.dumps(cover,indent=2,ensure_ascii=False)+"\n")
        print(json.dumps({k:cover[k] for k in ("status","certified_through_exact",
            "cell_count","first_unresolved_cell","minimum")},indent=2))

if __name__=="__main__":
    main()
