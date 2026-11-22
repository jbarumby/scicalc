"""Restricted, exact-rational Ramaré–Zúñiga source gate for CS v2 Table 1.

NOT a proof of global RH or the unrestricted CS Lemma 6. The bounds from
Ramaré--Zúñiga Alterman arXiv:2603.25961v1 are imported, not re-proven.
"""
import argparse
import json
from fractions import Fraction as F
from math import factorial
from pathlib import Path

DELTAS = tuple(map(F, (".012",".025",".038",".052",".066",".081",".097",
                       ".114",".132",".150",".170",".192",".215",".241",
                       ".270",".303")))
B4_54 = F("0.54")
B4_LIMIT = F("0.606")

def frac(x):
    return f"{x.numerator}/{x.denominator}"

def elementary():
    # For the exponential tail beginning at 1/7!, ratios are <= 1/8.
    e_upper = sum((F(1,factorial(k)) for k in range(7)), F(0))
    e_upper += F(8,7*factorial(7))
    if not e_upper < F(11,4):
        raise ArithmeticError("e upper bound failed")
    for n,threshold in ((28,3*10**12),(36,10**16),(75,10**33)):
        if not F(11,4)**n < threshold:
            raise ArithmeticError("logarithm lower bound failed")
    return frac(e_upper)

def row_result(index,delta):
    if not 0<delta<=F(".303"):
        raise ValueError("scope is 0 < delta <= 0.303")
    ranges=(
        ("3e12 <= X < 1e16",28),
        ("1e16 <= X < 1e33",36),
        ("X >= 1e33",75),
    )
    all_caps=[]
    for label, log_lower in ranges:
        eta_max=2*delta/log_lower
        if log_lower==28:
            if eta_max<=F(1,100):
                s_cap,eps_cap=F(".526"),F(1,100)  # v1, section 6
            elif eta_max<=F(1,50):
                s_cap,eps_cap=F(".5503"),F(1,50)
            elif eta_max<=F(1,25):
                s_cap,eps_cap=F(".593"),F(1,25)
            else:
                raise ArithmeticError("source range exceeded")
        elif log_lower==36:
            s_cap,eps_cap=F(".5104"),F(1,50)
        else:
            s_cap,eps_cap=F(".4794"),F(1,100)
        if eta_max>eps_cap:
            raise ArithmeticError("epsilon not covered by public theorem")
        # zeta(1+eta) <= 1+1/eta, exp(gamma*eta)>=1.
        # Thus eta * exp(-gamma eta) * zeta * S_eta
        # <= S_cap * (1+eta_max), with only rational operations.
        upper=s_cap*(1+eta_max)
        all_caps.append((label, upper, eps_cap))
    worst=max(x[1] for x in all_caps)
    if delta<=F(".132"):
        b4=B4_54
    else:
        scaled=worst*1000
        b4=F(scaled.numerator//scaled.denominator+1,1000)
    if not 0<b4-worst or not b4<=B4_LIMIT:
        raise ArithmeticError("no strict upper bound for b4")
    return {
        "row":index, "delta":frac(delta), "b4":frac(b4),
        "worst_upper":frac(worst), "margin":frac(b4-worst),
    }

def certificate():
    elementary()
    if len(DELTAS)!=16 or list(DELTAS)!=sorted(DELTAS):
        raise ArithmeticError("source table changed")
    rows=[row_result(i,delta) for i,delta in enumerate(DELTAS,1)]
    if [F(x["b4"])==B4_54 for x in rows] != [True]*9+[False]*7:
        raise ArithmeticError("nine-row original bound gate failed")
    return {
        "schema":"NUMERIA.P78.RESTRICTED_B4.1",
        "source":"https://arxiv.org/abs/2603.25961v1",
        "scope":"X >= 3e12, 0 < eta <= 2*delta/log(X)",
        "rows":rows,
        "original_CS_Lemma6_fully_proved":False,
        "B1_recomputed":False,
        "Suzuki_t79_replayed":False,
        "global_RH_proved":False,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    x=parser.add_mutually_exclusive_group(required=True)
    x.add_argument("--out",type=Path)
    x.add_argument("--verify",type=Path)
    args=parser.parse_args()
    data=certificate()
    if args.verify:
        if json.loads(args.verify.read_text(encoding="utf-8"))!=data:
            raise SystemExit("P78 certificate mismatch")
        print("P78 restricted public-source gate: 16 rows verified")
    else:
        args.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
        print("Wrote",args.out)

if __name__=="__main__":
    main()
