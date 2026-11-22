#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA P78-S2 (prototipo) — Diagnostico e RESOLUCAO do gate B1 com b4 publico.

Contexto: o prototipo P78-S1 (p78_b1_sensitivity.py) substituiu o b4 interno
0.54 pelos valores publicos Ramare-Zuniga nas linhas 10-16 e o gate strict_P78
FALHOU em todas as linhas (margens negativas). Este script diagnostica a causa
e entrega a resolucao (etapa S2: re-impressao dos B1 contra b4 publico).

MODELO EXATO (derivado linha a linha de b1_upper do S1/P77):
  B1_upper(b4) = D * k3(b4)**e,     k3(b4) = alpha * b4 + beta   (linear em b4)
  e = (2*s2 - 1)/(2 - s2),  s2 = (sigma_start + 1)/32
Valido enquanto k3 > 1: os canais max(1, k3^{-x}) saturam em 1 e o unico canal
b4-dependente vivo e max(k3^{e1}, k3^{e2}) = k3^{e2}, pois s -> (2s-1)/(2-s) e
crescente e k3 > 1. Nas linhas 10-16, e2 > 0 (s2 >= 26/32 > 1/2). Logo:

  1) B1_upper e ESTRITAMENTE CRESCENTE em b4 neste regime;
  2) os B1 impressos no P77 tinham folga <= 0.045% em b4 e suportam apenas
     b4 <= b4* (~= 0.5400-0.5405); o b4 publico (0.557-0.606) esta muito
     acima => NENHUM re-arranjo dentro desta algebra salva os valores
     impressos atuais (mudar isso exigiria refazer a interpolacao do P77);
  3) a algebra em si CONTINUA VALIDA com b4 publico (k3 >= 4.37 > 1):
     apenas os constantes impressos precisam subir.

RESOLUCAO S2 — re-impressao dos B1 contra b4 publico:
  reprint_minimal = menor multiplo de 1/100 com margem estrita > 0 (gate P77)
  reprint_robust  = menor multiplo de 1/100 com margem absoluta >= 1e-3
                    (mesma ordem das folgas originais do P77: 0.0009-0.0096)
Apos adotar os novos b1 na TABLE: re-certificar no padrao P77 (dyadic + CI +
SHA256) e recalcular as constantes downstream que consomem B1.

MELHORIAS vs S1
  - precisao maior (workprec 256 e tolerancias 1e-25 no modo REPO);
  - b4* (tolerancia maxima de b4 do B1 impresso) calculado por bissecao com
    comparacoes certificadas no modo REPO;
  - selecao CERTIFICADA do reprint (margem via .lower() > 0), nao apenas o gate;
  - checagem do modelo exato contra a forma fechada de e (e sigma_start
    recuperado como inteiro), provando a analise de sensibilidade.

MODOS
  REPO_CERTIFIED            — repo NUMERIA presente (flint + TABLE +
                              certificate.json): recalcula tudo em aritmetica
                              de intervalos e certifica gates/reprints.
  STANDALONE_RECONSTRUCTION — sem repo: reconstrucao a partir do
                              p78_b1_sensitivity_results.json do S1 (fit exato
                              de D, alpha, beta, e; validado contra a forma
                              fechada). Pre-visualizacao em ponto flutuante;
                              valores definitivos exigem REPO_CERTIFIED.

Uso: python3 p78_b1_resolution.py [caminho/p78_b1_sensitivity_results.json]
"""
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
B54 = Fr(27, 50)          # b4 interno de referencia do P77 (0.54 exato)
MIN_MARGIN = 1e-3         # folga alvo do reprint robusto (padrao P77)
NOISE = 1e-9              # guarda contra ruido float no modo standalone
SCHEMA = "NUMERIA.P78.S2.REPRINT.PROTOTYPE"


# ---------------------------------------------------------------------------
# Utilidades comuns
# ---------------------------------------------------------------------------

def find_repo():
    """Sobe na arvore ate achar NUMERIA_P76_COEFFICIENTS; None se ausente."""
    here = HERE
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    for cand in (Path("/home/z/my-project/scicalc-p78"),):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return None


def e_closed(sigma_start):
    """Expoente fechado do canal b4-dependente: e = (2*s2-1)/(2-s2)."""
    s2 = (sigma_start + 1) / 32.0
    return (2.0 * s2 - 1.0) / (2.0 - s2)


def reprint(up, min_margin):
    """Menor n/100 com n/100 - up >= min_margin (min_margin > 0 => estrita)."""
    n = math.ceil(100.0 * (up + min_margin))
    while n / 100.0 - up < min_margin:
        n += 1
    return f"{n / 100:.2f}", n / 100.0 - up


def default_s1():
    cands = [
        Path("/home/z/my-project/upload/p78_b1_sensitivity_results.json"),
        HERE / "p78_b1_sensitivity_results.json",
    ]
    for c in cands:
        if c.exists():
            return c
    raise SystemExit("informe o caminho do p78_b1_sensitivity_results.json")


# ---------------------------------------------------------------------------
# Modo STANDALONE: reconstrucao exata a partir dos resultados do S1
# ---------------------------------------------------------------------------

def fit_from_s1(o):
    k54, kpub = float(o["k3_at_0.54"]), float(o["k3_at_P78"])
    u54, upub = float(o["B1_upper_b4_0.54"]), float(o["B1_upper_b4_P78"])
    bpub = float(Fr(o["b4_P78"]))
    b54 = float(B54)

    # k3 linear em b4 (c3 e c4 lineares): alpha, beta exatos por 2 pontos
    alpha = (kpub - k54) / (bpub - b54)
    beta = k54 - alpha * b54

    # B1 = D * k3^e  =>  e exato pela razao dos dois pontos do S1
    e_fit = math.log(upub / u54) / math.log(kpub / k54)
    s2_rec = (1.0 + 2.0 * e_fit) / (2.0 + e_fit)
    sigma_rec = round(32.0 * s2_rec) - 1
    ec = e_closed(sigma_rec)
    model_ok = (e_fit > 0.0) and (abs(e_fit - ec) < 1e-6) \
        and abs(32.0 * s2_rec - 1.0 - sigma_rec) < 0.01

    D = u54 / (k54 ** e_fit)
    upub_fit = D * (kpub ** e_fit)
    fit_err = abs(upub_fit - upub) / upub

    printed = float(Fr(o["printed_B1"]))
    # b4* : maior b4 que o B1 impresso ainda certifica (B1 = D*k3^e invertido)
    b4_star = ((printed / D) ** (1.0 / e_fit) - beta) / alpha

    rp_min, m_min = reprint(upub, NOISE)
    rp_rob, m_rob = reprint(upub, MIN_MARGIN)

    return {
        "row": int(o["row"]), "delta": o["delta"],
        "sigma_start_recovered": sigma_rec,
        "e_exponent_fit": e_fit, "e_exponent_closed": ec,
        "model_consistent": bool(model_ok),
        "alpha_b4_linear": alpha, "beta_b4_linear": beta,
        "k3_at_0_54": k54, "k3_at_public": kpub,
        "k3_gt_1": bool(min(k54, kpub) > 1.0),
        "b4_public": str(Fr(o["b4_P78"])), "b4_public_dec": bpub,
        "b4_star_tolerable_float": b4_star,
        "b4_star_gap_to_public": bpub - b4_star,
        "printed_B1_old": o["printed_B1"],
        "B1_upper_b4_0_54": u54,
        "B1_upper_b4_public_float": upub,
        "B1_upper_b4_public_refit": upub_fit,
        "upper_refit_rel_err": fit_err,
        "margin_old_at_0_54_float": float(o["margin_0.54"]),
        "margin_old_at_public_float": printed - upub,
        "reprint_minimal": rp_min, "margin_minimal_float": m_min,
        "reprint_robust": rp_rob, "margin_robust_float": m_rob,
        "growth_vs_old_printed_pct": 100.0 * (float(rp_min) - printed) / printed,
    }


def standalone_summary(rows):
    mx = max(r["b4_star_tolerable_float"] for r in rows)
    mn = min(r["b4_public_dec"] for r in rows)
    if mx < mn:
        verdict = (
            f"B1 impresso do P77 NAO e salvavel com b4 publico dentro desta "
            f"algebra: b4* maximo = {mx:.6f} << b4 publico minimo = {mn}. "
            f"A algebra continua valida (k3 > 1 em todas as linhas); adotar o "
            f"reprint S2 (minimal/robusto) e re-certificar no padrao P77.")
    else:
        verdict = "B1 impresso suporta b4 publico em todas as linhas — sem re-impressao."
    return {
        "sigma_start_recovered": {str(r["row"]): r["sigma_start_recovered"] for r in rows},
        "model_consistent_all": all(r["model_consistent"] for r in rows),
        "upper_refit_ok_all": all(r["upper_refit_rel_err"] < 1e-9 for r in rows),
        "k3_gt_1_all": all(r["k3_gt_1"] for r in rows),
        "max_b4_star_tolerable": mx,
        "min_b4_public": mn,
        "old_print_salvageable_with_public_b4": bool(mx >= mn),
        "all_reprints_strict_float": all(r["margin_minimal_float"] > 0 for r in rows),
        "all_reprints_robust_float": all(r["margin_robust_float"] >= MIN_MARGIN for r in rows),
        "reprint_minimal": {str(r["row"]): r["reprint_minimal"] for r in rows},
        "reprint_robust": {str(r["row"]): r["reprint_robust"] for r in rows},
        "max_growth_needed_pct": max(r["growth_vs_old_printed_pct"] for r in rows),
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Modo REPO: recalculo com aritmetica de intervalos (padrao NUMERIA)
# ---------------------------------------------------------------------------

def repo_mode(repo):
    sys.path.insert(0, str(repo))
    from flint import acb, arb, ctx  # noqa: E402
    from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q  # noqa: E402
    from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table  # noqa: E402

    T0 = 3 * 10 ** 12

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

    cert = json.loads((repo / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json").read_text())
    p78_b4 = {r["row"]: Fr(str(r["b4"])) for r in cert["rows"]}
    validate_table()
    rows_out = []
    with ctx.workprec(256):
        cs = constants()
        cs["c1"] = 6 / cs["pi"] ** 2 + q("1.048") / cs["logh"]
        cs["c2"] = (cs["pi"] * cs["m0"] * q(".62") / cs["logh"]
                    + 6 * cs["m0"] / (cs["pi"] * T0)
                    + cs["pi"] * cs["m0"] * q("1.048") / (arb(T0) * cs["logh"]))
        for idx in sorted(p78_b4):
            row = TABLE[idx - 1]
            parts = row_parts(row, cs)
            k1, k2 = k1_k2(row, cs, parts)
            printed = q(row.b1)
            bpub = p78_b4[idx]
            b54a = arb(B54.numerator) / arb(B54.denominator)
            bpuba = arb(bpub.numerator) / arb(bpub.denominator)
            up54, k3_54 = b1_upper(row, cs, parts, k1, k2, b54a)
            up78, k3_78 = b1_upper(row, cs, parts, k1, k2, bpuba)

            def gate(b4a, _row=row, _cs=cs, _parts=parts, _k1=k1, _k2=k2, _p=printed):
                up, _k3 = b1_upper(_row, _cs, _parts, _k1, _k2, b4a)
                return (_p - up).lower() > 0

            e_c = e_closed(row.sigma_start)
            e_chk = math.log(float(up78) / float(up54)) / math.log(float(k3_78) / float(k3_54))
            model_ok = abs(e_chk - e_c) < 1e-9

            g_lo, g_hi = gate(b54a), gate(bpuba)
            b4s, b4s_float = None, None
            if g_hi:
                b4s = {"status": "ok_at_public", "b4_star_geq": float(bpub)}
                b4s_float = float(bpub)
            elif not g_lo:
                b4s = {"status": "broken_even_at_0_54"}
            else:
                lo, hi = B54, bpub
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if gate(arb(mid.numerator) / arb(mid.denominator)):
                        lo = mid
                    else:
                        hi = mid
                b4s = {"status": "bisected", "lower": float(lo), "upper": float(hi)}
                b4s_float = float(hi)

            # reprint certificado a partir do upper com b4 publico
            n = max(1, int(math.floor(100.0 * float(up78))) - 2)
            nmin = nrob = None
            mmin_f = mrob_f = None
            while n < 1000000:
                cand = arb(n) / arb(100)
                m = (cand - up78).lower()
                if nmin is None and m > 0:
                    nmin, mmin_f = n, float(m)
                if m - arb("0.001") > 0:
                    nrob, mrob_f = n, float(m)
                    break
                n += 1

            old = float(Fr(row.b1))
            rows_out.append({
                "row": idx, "delta": str(row.delta),
                "sigma_start": row.sigma_start,
                "e_exponent_closed": e_c, "e_exponent_checked": e_chk,
                "model_consistent": bool(model_ok),
                "k3_at_0_54": float(k3_54), "k3_at_public": float(k3_78),
                "k3_gt_1": bool(min(float(k3_54), float(k3_78)) > 1.0),
                "b4_public": str(bpub),
                "b4_star_tolerable_float": b4s_float,
                "b4_star_certified": b4s,
                "printed_B1_old": row.b1,
                "B1_upper_b4_0_54": float(up54), "B1_upper_b4_public_float": float(up78),
                "margin_old_at_0_54_lower": float((printed - up54).lower()),
                "margin_old_at_public_lower": float((printed - up78).lower()),
                "reprint_minimal": f"{nmin / 100:.2f}" if nmin is not None else None,
                "reprint_robust": f"{nrob / 100:.2f}" if nrob is not None else None,
                "margin_minimal_float": mmin_f,
                "margin_robust_float": mrob_f,
                "growth_vs_old_printed_pct":
                    100.0 * (nmin / 100.0 - old) / old if nmin is not None else None,
            })

    fails = [r["row"] for r in rows_out if not r["margin_old_at_public_lower"] > 0]
    b4s_eff = [r["b4_star_tolerable_float"] for r in rows_out
               if r["b4_star_tolerable_float"] is not None]
    ok_all = all(r["reprint_minimal"] is not None and r["margin_minimal_float"] > 0
                 for r in rows_out)
    summary = {
        "workprec_bits": 256,
        "rows_with_old_print_broken_at_public_b4": fails,
        "model_consistent_all": all(r["model_consistent"] for r in rows_out),
        "k3_gt_1_all": all(r["k3_gt_1"] for r in rows_out),
        "max_b4_star_tolerable": max(b4s_eff) if b4s_eff else None,
        "min_b4_public": min(float(Fr(r["b4_public"])) for r in rows_out),
        "all_reprints_strict_certified": ok_all,
        "all_reprints_robust_certified": all(
            r["reprint_robust"] is not None and r["margin_robust_float"] >= MIN_MARGIN
            for r in rows_out),
        "reprint_minimal": {str(r["row"]): r["reprint_minimal"] for r in rows_out},
        "reprint_robust": {str(r["row"]): r["reprint_robust"] for r in rows_out},
        "verdict": ("reprint S2 certificado disponivel; adotar novos b1 na TABLE, "
                    "re-certificar no padrao P77 (dyadic + CI + SHA256) e "
                    "recalcular as constantes downstream que consomem B1."
                    if ok_all else "HÁ LINHAS SEM REPRINT CERTIFICADO — INVESTIGAR."),
    }
    return rows_out, summary


# ---------------------------------------------------------------------------
# Saida
# ---------------------------------------------------------------------------

def print_table(rows):
    print()
    print(f"{'row':>3} {'delta':>6} {'b4_pub':>7} {'b4*':>9} {'B1_old':>7} "
          f"{'B1_up(pub)':>12} {'reprint_min':>11} {'margem_min':>11} {'reprint_rob':>11}")
    for r in rows:
        b4s = r.get("b4_star_tolerable_float")
        b4s_s = f"{b4s:.6f}" if b4s is not None else "-"
        print(f"{r['row']:>3} {r['delta']:>6} {r['b4_public']:>7} {b4s_s:>9} "
              f"{r['printed_B1_old']:>7} {r.get('B1_upper_b4_public_float', 0):>12.6f} "
              f"{r['reprint_minimal']:>11} {r['margin_minimal_float']:>11.6f} "
              f"{r['reprint_robust']:>11}")


def main():
    repo = find_repo()
    rows_out = summary = None
    mode = None
    if repo is not None:
        try:
            rows_out, summary = repo_mode(repo)
            mode = "REPO_CERTIFIED"
        except Exception as exc:  # fallback seguro
            print(f"[aviso] modo REPO indisponivel: {exc!r}")
            rows_out = summary = None
    if rows_out is None:
        src = Path(sys.argv[1]) if len(sys.argv) > 1 else default_s1()
        s1 = json.loads(Path(src).read_text())
        rows_out = [fit_from_s1(o) for o in s1["rows"]]
        summary = standalone_summary(rows_out)
        mode = "STANDALONE_RECONSTRUCTION"

    doc = {
        "schema": SCHEMA,
        "mode": mode,
        "b4_internal_reference": "27/50 (0.54)",
        "summary": summary,
        "rows": rows_out,
        "notes": [
            "MODELO EXATO: B1_upper(b4) = D * k3(b4)**e com k3 = alpha*b4 + beta; "
            "e = (2*s2-1)/(2-s2), s2 = (sigma_start+1)/32; valido enquanto k3 > 1.",
            "Linhas 10-16: e em [0.526, 1] => B1 estritamente crescente em b4; "
            "canais max(1, k3^{-x}) saturados em 1 (k3 >= 4.37).",
            "Os B1 impressos do P77 suportam b4 <= b4* (~0.5400-0.5405); b4 publico "
            "Ramare-Zuniga = 0.557-0.606 => re-impressao S2 obrigatoria.",
            "Apos adotar os novos b1 na TABLE: re-certificar no padrao P77 (dyadic + "
            "CI + SHA256) e recalcular as constantes downstream que consomem B1.",
            "Linhas 1-9 fora do escopo do S1; se o P78 as afetar, rodar o mesmo "
            "protocolo (o script percorre as linhas do certificate.json no modo REPO).",
            "Modo STANDALONE e pre-visualizacao em ponto flutuante; valores "
            "definitivos exigem o modo REPO_CERTIFIED.",
        ],
    }
    print(json.dumps(doc, indent=1))
    print_table(rows_out)
    print("RESOLUCAO:", summary["verdict"])
    dest = HERE / "p78_b1_resolution_results.json"
    dest.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"[ok] resultados gravados em {dest}")


if __name__ == "__main__":
    main()
