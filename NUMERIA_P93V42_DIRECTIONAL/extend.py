"""Extend the fixed N=3309 directional sufficient head; never a global certificate."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from flint import arb, arb_mat, ctx
from NUMERIA_P93V40.assembly import assemble, WIDTH, LOCALIZATION
from NUMERIA_P93V40.edge_schur import (
    edge_constants, complement_budget, _block, _decode, OMEGA)
from NUMERIA_P93V40.edge_direction import directional_cross_majorant
from NUMERIA_P93V40.edge_sinc import endpoint_coefficients
from NUMERIA_P93V40.edge_schur_algebra import _symmetric_enclosure
from NUMERIA_P93V40.prefix import _head_digest
from NUMERIA_P93V40.precondition import (
    _point_factor, verify_triangular_congruence, certify_by_congruence, prove_negative)
from NUMERIA_P93V40.transfer import dyadic

N_INFINITE = 3309
def proof(m):
    n = m.nrows()
    report = certify_by_congruence(m, [n], 0, absolute_bound=True)
    if report["status"] == "inconclusive_congruence":
        detail = certify_by_congruence(m, list(range(1, n+1)), 0)
        report["first_unproved_leading_size"] = next(
            (r["size"] for r in detail["prefixes"] if not r["certified"]), None)
    return report

def radius(m):
    return dyadic(max(x.rad().upper() for x in m.entries()))

def block_extension(matrix, split):
    n = matrix.nrows()
    first, second = range(split), range(split,n)
    a = _block(matrix, first, first)
    b = _block(matrix, first, second)
    d = _block(matrix, second, second)
    row = {"A_sha256": _head_digest(a), "D_sha256": _head_digest(d),
           "A_radius_upper": radius(a), "B_radius_upper": radius(b),
           "D_radius_upper": radius(d)}
    lower, scaling, obstruction = _point_factor(a, arb(0))
    if obstruction is not None:
        row.update(status="inconclusive_A_proposal", A_proof=proof(a))
        return row
    inverse = lower.transpose().inv().mid()
    r = arb_mat(split, split, [(inverse[i,j]*scaling[j]).mid() if i<=j else arb(0)
                               for i in first for j in first])
    verify = verify_triangular_congruence(a,r,[split],0)[0]
    row["A_congruence"] = verify
    if not verify["certified"]:
        row.update(status="inconclusive_A_congruence", A_proof=proof(a))
        return row
    alpha = _decode(verify["row_margin_lower_dyadic"])
    h = _symmetric_enclosure(r.transpose()*a*r)
    f = r.transpose()*b
    solved = None
    try:
        solved = h.solve(f)
        if not all(x.is_finite() for x in solved.entries()):
            solved = None
    except (ArithmeticError, ValueError, ZeroDivisionError):
        pass
    row["interval_solve_finite"] = solved is not None
    if solved is not None:
        x = r*solved
        row["X_radius_upper"] = radius(x)
        s = _symmetric_enclosure(d-f.transpose()*solved)
        row["direct_schur_sha256"] = _head_digest(s)
        row["direct_schur_proof"] = proof(s)
        if row["direct_schur_proof"]["status"] == "certified":
            row["status"] = "certified"
            return row
    # Exact midpoint trial, with a residual correction in the verified basis.
    # S = C - E^T H^-1 E >= C - E^T E/alpha.
    z = solved.mid() if solved is not None else f.mid()
    hz = h*z
    e = f-hz
    c = _symmetric_enclosure(d-f.transpose()*z-z.transpose()*f+z.transpose()*hz)
    s_lower = _symmetric_enclosure(c-(e.transpose()*e)*(1/alpha))
    row["residual_radius_upper"] = radius(e)
    row["lower_schur_sha256"] = _head_digest(s_lower)
    row["lower_schur_proof"] = proof(s_lower)
    if row["lower_schur_proof"]["status"] == "certified":
        row["status"] = "certified"
        return row
    # A negative direction of the sufficient lower bound alone is NOT evidence.
    # Lift proposed y into the original corrected matrix and re-evaluate there.
    for name in ("direct_schur_proof","lower_schur_proof"):
        witness = row.get(name,{}).get("rayleigh_check",{}).get("witness_dyadic")
        if witness:
            y = arb_mat(n-split,1,[_decode(pair) for pair in witness])
            top = (-(r*(z*y))).mid()
            v = arb_mat(n,1,top.entries()+y.entries())
            check = prove_negative(matrix,v,0)
            row.setdefault("lifted_witnesses",[]).append({"origin":name,**check})
            if check["proven_negative"]:
                row["status"] = "negative_direction"
                return row
    row["status"] = "inconclusive"
    return row

def run(modes, split, bits):
    if not 0 < split < modes <= LOCALIZATION:
        raise ValueError("require split < modes <= localization")
    with ctx.workprec(bits):
        sectors = assemble(modes,bits)
        _, epsilon, _ = edge_constants()
        budgets = {q["parity"]:q for q in complement_budget()["results"]
                   if q["modes"] == N_INFINITE}
        endpoints = endpoint_coefficients(modes,width=WIDTH,omega=OMEGA)
        rows = []
        for p, source in enumerate(sectors):
            q = budgets[p]
            if not q["tail_certified_by_bound"]:
                raise ArithmeticError("infinite tail not certified")
            base = arb_mat(source["head"])
            for i in range(modes):
                base[i,i] -= epsilon
            u = arb_mat(modes,1,endpoints[p])
            m = _symmetric_enclosure(
                base - directional_cross_majorant(modes,N_INFINITE,p)*
                _decode(q["cross_majorant_multiplier_upper_dyadic"])
                - (u*u.transpose())*_decode(q["required_head_endpoint_weight_upper_dyadic"]))
            record = block_extension(m,split)
            record.update(parity=p, corrected_matrix_sha256=_head_digest(m),
                          tail_floor_dyadic=q["corrected_tail_reserve_lower_dyadic"],
                          quadrature_entry_error_exact=source["quadrature_entry_error_exact"])
            rows.append(record)
        return {"schema":"NUMERIA.P93V42.DIRECTIONAL_EXTENSION.1",
                "scope":"finite principal compression of fixed directional sufficient head; not a crossing certificate",
                "modes_per_parity":modes,"split_per_parity":split,"infinite_split":N_INFINITE,
                "precision_bits":bits,"sectors":rows,
                "certified":all(q["status"]=="certified" for q in rows)}

def self_test():
    with ctx.workprec(256):
        for d, expected in [(3,"certified"),(0,"negative_direction")]:
            m = arb_mat([[2,1],[1,d]])
            got = block_extension(m,1)["status"]
            if got != expected:
                raise ArithmeticError(f"self-test expected {expected}, got {got}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--modes",type=int,default=400)
    parser.add_argument("--split",type=int,default=200)
    parser.add_argument("--bits",type=int,default=768)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        self_test()
        print("Synthetic positive/negative Schur checks passed")
        return
    if args.out is None:
        parser.error("--out required")
    result=run(args.modes,args.split,args.bits)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"certified":result["certified"],
                      "sectors":[{"parity":x["parity"],"status":x["status"]} for x in result["sectors"]]}))
if __name__=="__main__":
    main()
