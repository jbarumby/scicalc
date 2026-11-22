#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA P78-S3 (prototipo) — Protocolo completo unificado, TODAS as linhas.

Une em um comando a cadeia S1 -> S2 para TODAS as linhas do certificate.json
do P78 (nao apenas 10-16): e a extensao do protocolo as linhas 1-9.

Para cada linha selecionada (--rows, default: todas) e cada precisao da
varredura dyadic (--dyadic, default 256):

  1. K1, K2 rigorosos (acb.integral, uma vez por linha/precisao);
  2. S1: B1_upper em b4 = 0.54 (interno P77) e em b4 publico Ramare-Zuniga;
     margens do b1 impresso atual contra ambos; delta_B1_pct; k3 nos pontos;
  3. MODELO EXATO: e verificado contra e = (2*s2-1)/(2-s2),
     s2 = (sigma_start+1)/32; alpha, beta da linearidade de k3 em b4;
     condicao de validade k3 > 1 (saturacao dos canais max(1, k3^{-x}));
  4. b4*: maior b4 que o b1 impresso atual ainda certifica (bissecao com
     gate certificado), quando o gate passa em 0.54 e falha no publico;
  5. S2: reprint certificado contra b4 publico:
       minimal = menor multiplo de 1/100 com margem estrita > 0;
       robust  = menor multiplo de 1/100 com margem absoluta >= 1e-3;
  6. STATUS da linha:
       ok_keep_printed            — gate passa com b4 publico; sem mudanca;
       broken_only_with_public_b4 — passa em 0.54, falha no publico; reprint;
       broken_at_0_54_and_public  — falha em ambos; reprint (investigar P77);
     (se e2 <= 0, B1 e nao-crescente em b4: registrado, os gates decidem.)

Saidas (em repo.parent):
  p78_s2_full_protocol_results.json — diagnostico completo + SHA256 das entradas;
  p78_s2_table_patch_all.json       — patch da TABLE para TODAS as linhas
                                      (keep_current para as OK; minimal/robust
                                      para as quebradas), pronto para o
                                      p78_s2_recertify.py --proposed.

Checagens cruzadas automaticas (aviso, nao gate) contra artefatos anteriores
em repo.parent: p78_b1_sensitivity_results.json (S1, linhas 10-16),
p78_b1_resolution_results.json e p78_s2_table_patch.json (S2).

NOTAS
  - Prototipo single-run; o certificado OFICIAL e o do pipeline NUMERIA
    (dyadic completo + CI + SHA256). Use --dyadic 192,256,384 para a varredura
    completa antes do CI. Runtime cresce com linhas x precisions (K1/K2 por
    par linha/precisao).
  - Se validate_table() tiver checagens internas sobre b1, atualize-as junto.
  - Linha com b4 publico == 0.54: sem contraste S1 (e_chk/alpha nulos).

Uso:
  python3 p78_s2_full_protocol.py                       # todas as linhas
  python3 p78_s2_full_protocol.py --rows 1-9            # so as novas
  python3 p78_s2_full_protocol.py --dyadic 192,256,384  # varredura completa
"""
import argparse
import hashlib
import json
import math
import platform
import sys
from fractions import Fraction as Fr
from pathlib import Path

SCHEMA = "NUMERIA.P78.S3.FULLPROTOCOL.PROTOTYPE"
B54 = Fr(27, 50)          # b4 interno de referencia do P77 (0.54 exato)
MIN_MARGIN = 1e-3         # folga do reprint robusto (padrao P77)


def find_repo():
    here = Path(__file__).resolve().parent
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    for cand in (Path("/home/z/my-project/scicalc-p78"), Path.cwd()):
        if cand is not None and (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return None


def e_closed(sigma_start):
    s2 = (sigma_start + 1) / 32.0
    return (2.0 * s2 - 1.0) / (2.0 - s2)


def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parse_rows(spec, available):
    """'1-9,12' / '1,5,7' / None -> linhas selecionadas (intersecao com available)."""
    if not spec:
        return list(available)
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            if not (a.strip().isdigit() and b.strip().isdigit()):
                raise SystemExit(f"--rows: trecho invalido {part!r}")
            out |= {r for r in available if int(a) <= r <= int(b)}
        elif part.isdigit():
            out.add(int(part))
        else:
            raise SystemExit(f"--rows: trecho invalido {part!r}")
    sel = sorted(out & set(available))
    if not sel:
        raise SystemExit(f"--rows {spec!r} nao cruza as linhas do certificado "
                         f"{sorted(available)}")
    return sel


def main():
    ap = argparse.ArgumentParser(
        description="Protocolo completo S1+S2 unificado, todas as linhas (prototipo)")
    ap.add_argument("--rows", metavar="SPEC", default=None,
                    help="linhas a processar, ex.: 1-9 ou 1,5,7-12 "
                         "(default: todas as do certificate.json)")
    ap.add_argument("--variant", choices=("minimal", "robust"), default="robust",
                    help="variante padrao do patch_all emitido (default: robust)")
    ap.add_argument("--dyadic", metavar="PRECS", default="256",
                    help="precisoes da varredura dyadic, ex.: 192,256,384 (default: 256)")
    args = ap.parse_args()

    repo = find_repo()
    if repo is None:
        raise SystemExit("repo NUMERIA nao encontrado (NUMERIA_P76_COEFFICIENTS "
                         "ausente); rode dentro do repo")
    sys.path.insert(0, str(repo))
    try:
        import flint as _flint
        from flint import acb, arb, ctx
    except ImportError as exc:
        raise SystemExit(f"python-flint ausente: {exc}")
    from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q
    from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table

    T0 = 3 * 10 ** 12

    # ------------------------------------------------------------------
    # Maquinaria de intervalos (mesmas definicoes do S1/S2)
    # ------------------------------------------------------------------
    def rat(n, d):
        return arb(n) / arb(d)

    def pw(x, e):
        if not x > 0:
            raise ArithmeticError("non-positive real power base")
        return (x.log() * e).exp()

    def row_parts(row, cs):
        d, kap, delta = q(row.d), q(row.kappa), q(row.delta)
        L, pi, m0 = cs["logh"], cs["pi"], cs["m0"]
        A0 = (arb(10).log() * q(row.log10_a0)).exp()
        La = A0.log()
        half = (A0 / 2).log()
        a1 = A0 / (2 * pi * pi) * half ** 4 + q("34.842") * A0 * half ** 3 + arb("1.3897e16")
        r = 1 - arb.const_log2() / La
        a2 = r ** 4 / (2 * pi * pi) + q("34.842") / La * r ** 3 + arb("9.5079e15") / (A0 * La ** 4)
        c = 3 / (4 * kap)
        return d, kap, delta, L, pi, m0, A0, a1, a2, c

    def k1_k2(row, cs, parts):
        d, kap, delta, L, pi, m0, A0, a1, a2, c = parts

        def f(x, analytic):
            v = x ** 6
            A = acb(1) + (acb(1) + acb(c) * v).log(analytic=analytic) / (2 * acb(L))
            B = acb(cs["c1"] * c.sqrt()) * x ** 3 + acb(cs["c2"])
            return 6 * v * (-v).exp() * A.pow(acb(4) / 3, analytic=analytic) \
                * B.pow(acb(2) / 3, analytic=analytic)

        I = acb.integral(f, 0, 2, rel_tol=arb('1e-25'), abs_tol=arb('1e-25'),
                         eval_limit=30000, depth_limit=30)
        if not I.is_finite() or not I.imag.contains(0) or not I.real.lower() > 0:
            raise ArithmeticError("quadrature failed or inconclusive")
        U = arb(64)
        D = (1 + c / (2 * L)) ** 2 * (1 + cs["c1"] * c.sqrt() + cs["c2"])
        tail = (D * (-U).exp() * (U ** 3 + 3 * U ** 2 + 6 * U + 6)).upper()
        lower, upper = I.real.lower(), (I.real.upper() + tail).upper()
        F = pw(2 * a2 * c.sqrt(), rat(1, 3))
        correction = (pw(2 * a1, rat(1, 3)) * pw(cs["c2"] + A0 * cs["c1"] / T0, rat(2, 3))
                      / (pw(arb(T0), rat(1, 3)) * pw(L, rat(4, 3))))
        k1lo = (F * lower + correction).lower()
        k1 = (F * upper + correction).upper()
        if not k1lo > 0:
            raise ArithmeticError("K1 lower bound inconclusive")
        two13 = pw(arb(2), rat(1, 3))
        k2 = two13 * c.sqrt() * pi.sqrt() / 2 + two13 * (4 * kap / 3) * pw(A0 / arb(T0), arb(3))
        return k1, k2

    def b1_upper(row, cs, parts, k1, k2, b4):
        d, kap, delta, L, pi, m0, A0, a1, a2, c = parts
        s1, s2 = rat(row.sigma_start, 32), rat(row.sigma_start + 1, 32)
        c3 = pi * m0 * b4 / (2 * delta) * (2 * delta * arb.const_euler() / L).exp()
        c4 = b4 / (5 * delta * delta.exp()) + q(".605") * (-2 * delta).exp() / L
        k3 = c3 + c4 * (pi.sqrt() / (2 * (2 * kap).sqrt()) + pi * m0 / T0)
        l0 = (kap * (1 + (1 / (4 * arb(T0)) * (1 + 2 * delta / L)) ** 2)).exp() / (1 - arb(1) / T0)
        l1 = max(arb(1), pw(k1, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) \
            * (3 * delta + 4 * delta * L.log() / L).exp()
        l2 = max(arb(1), pw(k2, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) \
            * (3 * delta).exp()
        common = max(arb(1), pw(l0, -2 * d / L)) * max(arb(1), pw(k3, -3 * d / L))
        common *= max(pw(l0, 2 / (2 - s1)), pw(l0, 2 / (2 - s2)))
        common *= max(pw(k3, (2 * s1 - 1) / (2 - s1)), pw(k3, (2 * s2 - 1) / (2 - s2)))
        A1 = l1 * (d * (3 + 4 * L.log() / L)).exp() * common
        A1 *= max(arb(1), pw(k1, 3 * d / L))
        A1 *= max(pw(k1, 3 * (1 - s1) / (2 - s1)), pw(k1, 3 * (1 - s2) / (2 - s2)))
        A2 = l2 * (3 * d).exp() * common
        A2 *= max(arb(1), pw(k2, 3 * d / L))
        A2 *= max(pw(k2, 3 * (1 - s1) / (2 - s1)), pw(k2, 3 * (1 - s2) / (2 - s2)))
        B1 = (A1 + A2 * pw(L, -6 * (1 - s2) / (2 - s2))) / (2 * pi * d)
        return B1.upper(), k3

    def make_constants(prec):
        with ctx.workprec(prec):
            cs = constants()
            cs["c1"] = 6 / cs["pi"] ** 2 + q("1.048") / cs["logh"]
            cs["c2"] = (cs["pi"] * cs["m0"] * q(".62") / cs["logh"]
                        + 6 * cs["m0"] / (cs["pi"] * T0)
                        + cs["pi"] * cs["m0"] * q("1.048") / (arb(T0) * cs["logh"]))
            return cs

    # ------------------------------------------------------------------
    # Plano: linhas do certificado
    # ------------------------------------------------------------------
    cert_path = repo / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json"
    cert = json.loads(cert_path.read_text())
    p78_b4 = {int(r["row"]): Fr(str(r["b4"])) for r in cert["rows"]}
    validate_table()
    row_ids = parse_rows(args.rows, sorted(p78_b4))
    precs = sorted({int(p) for p in args.dyadic.split(",") if p.strip()})
    print(f"[plano] linhas {row_ids} | dyadic precs {precs} | "
          f"variante patch: {args.variant}")

    rows_out = []
    for idx in row_ids:
        row = TABLE[idx - 1]
        bpub = p78_b4[idx]
        old_str = str(TABLE[idx - 1].b1)
        per_prec = {}
        for prec in precs:
            cs = make_constants(prec)
            with ctx.workprec(prec):
                parts = row_parts(row, cs)
                k1, k2 = k1_k2(row, cs, parts)
                b54a = arb(B54.numerator) / arb(B54.denominator)
                bpuba = arb(bpub.numerator) / arb(bpub.denominator)
                up54, k3_54 = b1_upper(row, cs, parts, k1, k2, b54a)
                uppub, k3_pub = b1_upper(row, cs, parts, k1, k2, bpuba)
                printed_old = q(old_str)
                m54 = (printed_old - up54).lower()
                mpub = (printed_old - uppub).lower()
                rec = {
                    "workprec": prec,
                    "k1": float(k1), "k2": float(k2),
                    "k3_at_0_54": float(k3_54), "k3_at_public": float(k3_pub),
                    "B1_upper_0_54": float(up54),
                    "B1_upper_public_lower": float(uppub.lower()),
                    "B1_upper_public_upper": float(uppub.upper()),
                    "margin_old_at_0_54_lower": float(m54),
                    "margin_old_at_public_lower": float(mpub),
                    "gate_old_at_0_54": bool(m54 > 0),
                    "gate_old_at_public": bool(mpub > 0),
                    "delta_B1_upper_pct":
                        100.0 * (float(uppub) - float(up54)) / float(up54),
                }
                if bpub == B54:
                    rec["e_exponent_checked"] = None
                    rec["alpha_b4_linear"] = None
                    rec["beta_b4_linear"] = None
                    rec["note"] = "b4 publico == 0.54: sem contraste S1"
                else:
                    rec["e_exponent_checked"] = math.log(
                        float(uppub) / float(up54)) / math.log(
                        float(k3_pub) / float(k3_54))
                    rec["alpha_b4_linear"] = (float(k3_pub) - float(k3_54)) \
                        / (float(bpub) - float(B54))
                    rec["beta_b4_linear"] = float(k3_54) \
                        - rec["alpha_b4_linear"] * float(B54)
                # reprint certificado contra b4 publico
                n = max(1, int(math.floor(100.0 * float(uppub))) - 2)
                nmin = nrob = None
                mmin = mrob = None
                while n < 1000000:
                    cand = arb(n) / arb(100)
                    m = (cand - uppub).lower()
                    if nmin is None and m > 0:
                        nmin, mmin = n, float(m)
                    if m - arb("0.001") > 0:
                        nrob, mrob = n, float(m)
                        break
                    n += 1
                rec["reprint_minimal"] = f"{nmin / 100:.2f}" if nmin is not None else None
                rec["reprint_robust"] = f"{nrob / 100:.2f}" if nrob is not None else None
                rec["margin_minimal_lower"] = mmin
                rec["margin_robust_lower"] = mrob
                per_prec[str(prec)] = rec

        rec0 = per_prec[str(precs[0])]
        e_c = e_closed(row.sigma_start)
        e_chk = rec0["e_exponent_checked"]
        model_ok = None if e_chk is None else bool(abs(e_chk - e_c) < 1e-9)
        k3_min = min(rec0["k3_at_0_54"], rec0["k3_at_public"])
        sat_ok = bool(k3_min > 1.0)

        # status pela referencia (precs[0]); gates devem ser estaveis em precisao
        g54, gpub = rec0["gate_old_at_0_54"], rec0["gate_old_at_public"]
        if gpub:
            status = "ok_keep_printed"
        elif not g54:
            status = "broken_at_0_54_and_public"
        else:
            status = "broken_only_with_public_b4"

        # b4*: so faz sentido quando 0.54 passa e o publico falha
        b4s = None
        if status == "broken_only_with_public_b4":
            # bissecao com gate certificado na precisao de referencia
            cs0 = make_constants(precs[0])
            with ctx.workprec(precs[0]):
                parts0 = row_parts(row, cs0)
                k1_0, k2_0 = k1_k2(row, cs0, parts0)

                def gate(b4f):
                    up, _k3 = b1_upper(row, cs0, parts0, k1_0, k2_0, b4f)
                    return (q(old_str) - up).lower() > 0

            lo, hi = B54, bpub
            for _ in range(60):
                mid = (lo + hi) / 2
                if gate(arb(mid.numerator) / arb(mid.denominator)):
                    lo = mid
                else:
                    hi = mid
            b4s = {"lower": float(lo), "upper": float(hi)}

        old_f = float(Fr(old_str))
        nmin = rec0["reprint_minimal"]
        rows_out.append({
            "row": idx, "delta": str(row.delta),
            "sigma_start": row.sigma_start,
            "b4_public": str(bpub),
            "b1_old": old_str,
            "status": status,
            "e_exponent_closed": e_c,
            "e_exponent_checked": e_chk,
            "model_consistent": model_ok,
            "k3_min": k3_min, "k3_gt_1_saturation_ok": sat_ok,
            "alpha_b4_linear": rec0["alpha_b4_linear"],
            "beta_b4_linear": rec0["beta_b4_linear"],
            "b4_star_certified": b4s,
            "b4_star_tolerable_float": b4s["upper"] if b4s else None,
            "delta_B1_upper_pct": rec0["delta_B1_upper_pct"],
            "reprint_minimal": nmin, "reprint_robust": rec0["reprint_robust"],
            "margin_minimal_lower": rec0["margin_minimal_lower"],
            "margin_robust_lower": rec0["margin_robust_lower"],
            "growth_vs_old_printed_pct":
                100.0 * (float(nmin) - old_f) / old_f if nmin else None,
            "gate_pass_all_precisions_reprint_minimal": all(
                v["margin_minimal_lower"] is not None and v["margin_minimal_lower"] > 0
                for v in per_prec.values()),
            "dyadic": per_prec,
        })

    # ------------------------------------------------------------------
    # Checagens cruzadas contra artefatos anteriores (aviso, nao gate)
    # ------------------------------------------------------------------
    cross = []
    s1p = repo.parent / "p78_b1_sensitivity_results.json"
    if s1p.exists():
        try:
            s1 = {int(r["row"]): r for r in
                  json.loads(s1p.read_text())["rows"]}
            for r in rows_out:
                if r["row"] in s1 and r["dyadic"].get(str(precs[0])):
                    o = s1[r["row"]]
                    ok_u54 = abs(r["dyadic"][str(precs[0])]["B1_upper_0_54"]
                                 - o["B1_upper_b4_0.54"]) \
                        / o["B1_upper_b4_0.54"] < 1e-8
                    ok_up = abs(r["dyadic"][str(precs[0])]["B1_upper_public_upper"]
                                - o["B1_upper_b4_P78"]) \
                        / o["B1_upper_b4_P78"] < 1e-8
                    cross.append({"vs": "S1", "row": r["row"],
                                  "consistent": bool(ok_u54 and ok_up)})
                    if not (ok_u54 and ok_up):
                        print(f"[aviso] linha {r['row']}: divergencia vs S1")
        except Exception as exc:
            cross.append({"vs": "S1", "error": repr(exc)})
    s2p = repo.parent / "p78_s2_table_patch.json"
    if s2p.exists():
        try:
            s2 = {int(r["row"]): r for r in
                  json.loads(s2p.read_text())["rows"]}
            for r in rows_out:
                if r["row"] in s2 and r["status"] != "ok_keep_printed":
                    okm = (r["reprint_minimal"] == s2[r["row"]]["b1_new_minimal"]
                           and r["reprint_robust"] == s2[r["row"]]["b1_new_robust"])
                    cross.append({"vs": "S2_patch", "row": r["row"],
                                  "consistent": bool(okm)})
                    if not okm:
                        print(f"[aviso] linha {r['row']}: reprint difere do patch S2")
        except Exception as exc:
            cross.append({"vs": "S2_patch", "error": repr(exc)})

    # ------------------------------------------------------------------
    # Patch completo da TABLE (todas as linhas processadas)
    # ------------------------------------------------------------------
    patch_rows = []
    for r in rows_out:
        if r["status"] == "ok_keep_printed":
            patch_rows.append({
                "row": r["row"], "delta": r["delta"],
                "b4_public": r["b4_public"], "b1_old": r["b1_old"],
                "action": "keep_current",
            })
        else:
            patch_rows.append({
                "row": r["row"], "delta": r["delta"],
                "b4_public": r["b4_public"], "b1_old": r["b1_old"],
                "b1_new_minimal": r["reprint_minimal"],
                "b1_new_robust": r["reprint_robust"],
                "expected_upper_public_float":
                    r["dyadic"][str(precs[0])]["B1_upper_public_upper"],
                "expected_margin_minimal_float": r["margin_minimal_lower"],
                "expected_margin_robust_float": r["margin_robust_lower"],
                "e_exponent_closed": r["e_exponent_closed"],
                "growth_vs_old_pct": round(r["growth_vs_old_printed_pct"], 4)
                if r["growth_vs_old_printed_pct"] is not None else None,
                "action": "reprint",
            })
    patch_doc = {
        "schema": "NUMERIA.P78.S2.TABLE.PATCH",
        "based_on": [SCHEMA,
                     "NUMERIA.P78.S1.PROTOTYPE / NUMERIA.P78.S2.REPRINT.PROTOTYPE"],
        "policy": ("default_variant = " + args.variant +
                   "; minimal = margem estrita > 0 (gate P77); robust = margem "
                   "absoluta >= 1e-3. Linhas ok_keep_printed nao mudam."),
        "default_variant": args.variant,
        "verification": ("python3 p78_s2_recertify.py --proposed "
                         "p78_s2_table_patch_all.json [--dyadic 192,256,384]; "
                         "apos editar a TABLE: python3 p78_s2_recertify.py; "
                         "depois certificado oficial dyadic + CI + SHA256 e "
                         "p78_s2_downstream_scan.py."),
        "rows": patch_rows,
    }
    patch_dest = repo.parent / "p78_s2_table_patch_all.json"

    fails_model = [r["row"] for r in rows_out if r["model_consistent"] is False]
    fails_sat = [r["row"] for r in rows_out if not r["k3_gt_1_saturation_ok"]]
    fails_reprint = [r["row"] for r in rows_out
                     if r["status"] != "ok_keep_printed"
                     and not r["gate_pass_all_precisions_reprint_minimal"]]
    summary = {
        "rows_processed": row_ids,
        "status_counts": {
            s: sum(1 for r in rows_out if r["status"] == s)
            for s in ("ok_keep_printed", "broken_only_with_public_b4",
                      "broken_at_0_54_and_public")},
        "rows_needing_reprint": [r["row"] for r in rows_out
                                 if r["status"] != "ok_keep_printed"],
        "rows_ok_keep_printed": [r["row"] for r in rows_out
                                 if r["status"] == "ok_keep_printed"],
        "model_consistent_fail_rows": fails_model,
        "k3_gt_1_fail_rows": fails_sat,
        "reprint_gate_fail_rows": fails_reprint,
        "dyadic_precisions": precs,
        "patch_variant": args.variant,
        "cross_checks": cross,
        "verdict": ("protocolo completo consistente: reprint disponivel para "
                    "todas as linhas quebradas, patch_all emitido; seguir com "
                    "p78_s2_recertify.py --proposed e depois o certificado "
                    "oficial dyadic + CI + SHA256."
                    if not (fails_model or fails_sat or fails_reprint) else
                    "ATENCAO: linhas com problema — model " + str(fails_model) +
                    " / saturacao " + str(fails_sat) +
                    " / reprint " + str(fails_reprint)),
    }

    sha = {}
    for label, p in [("table_py", repo / "NUMERIA_P76_COEFFICIENTS" / "table.py"),
                     ("audit_b23_py", repo / "NUMERIA_P76_COEFFICIENTS" / "audit_b23.py"),
                     ("certificate_json", cert_path)]:
        if p.exists():
            sha[label] = {"path": str(p), "sha256": sha256_file(p)}

    doc = {
        "schema": SCHEMA,
        "rows_spec": args.rows or "todas (certificate.json)",
        "python": platform.python_version(),
        "flint_version": getattr(_flint, "__version__", "unknown"),
        "sha256_inputs": sha,
        "summary": summary,
        "rows": rows_out,
        "notes": [
            "Modelo exato: B1_upper(b4) = D * k3(b4)**e, k3 = alpha*b4 + beta, "
            "e = (2*s2-1)/(2-s2), s2 = (sigma_start+1)/32; valido enquanto k3 > 1.",
            "Status: ok_keep_printed = gate passa com b4 publico (sem mudanca); "
            "broken_* = reprint obrigatorio (minimal/robusto).",
            "b4* reportado apenas para linhas broken_only_with_public_b4 "
            "(bissecao certificada na precisao de referencia).",
            "Prototipo single-run; certificado oficial = pipeline NUMERIA "
            "(dyadic completo + CI + SHA256).",
        ],
    }
    print(json.dumps(doc, indent=1))
    print()
    print(f"{'row':>3} {'delta':>6} {'b4_pub':>9} {'status':>28} {'b4*':>9} "
          f"{'B1_old':>8} {'reprint_min':>11} {'margem_min':>11} {'reprint_rob':>11}")
    for r in rows_out:
        b4s = r["b4_star_tolerable_float"]
        b4s_s = f"{b4s:.6f}" if b4s is not None else "-"
        mmin_s = (f"{r['margin_minimal_lower']:.6f}"
                  if r["margin_minimal_lower"] is not None else "-")
        print(f"{r['row']:>3} {r['delta']:>6} {r['b4_public']:>9} {r['status']:>28} "
              f"{b4s_s:>9} {r['b1_old']:>8} "
              f"{r['reprint_minimal']:>11} {mmin_s:>11} "
              f"{r['reprint_robust']:>11}")
    print("VEREDITO:", summary["verdict"])

    results_dest = repo.parent / "p78_s2_full_protocol_results.json"
    results_dest.write_text(json.dumps(doc, indent=1) + "\n")
    patch_dest.write_text(json.dumps(patch_doc, indent=1) + "\n")
    print(f"[ok] resultados: {results_dest}")
    print(f"[ok] patch completo: {patch_dest}")
    if fails_model or fails_sat or fails_reprint:
        raise SystemExit("S3 GATE FAIL: " + summary["verdict"])


if __name__ == "__main__":
    main()
