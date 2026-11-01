#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUDITAR_E8K_V18_TESTEMUNHAS.py - Replay independente das testemunhas
publicadas no V17 + controles negativos REAIS (schema
NUMERIA.E8K.V18.WITNESS_INDEPENDENT_REPLAY.1).

Defeitos do V17 que o V18 corrige POR DESENHO:

  1. CONTROLES NEGATIVOS CONSTANTES -> aqui cada mutacao altera o
     RELATORIO PUBLICADO e exige que os DETECTORES rejitem. Mutacao
     nao detectada = selftest reprova. O gate nao pode passar "de
     verdade" sem nunca poder falhar.
  2. CANCELAMENTO DE POLOS AUTO-REFERENTE -> a identidade e testada
     em amostras racionais exatas de X=pi^2 e o ZERO e OBSERVADO
     numericamente com quadratura independente (Simpson), nunca
     A-A=0 simbolico.
  3. CASO RESSONANTE k=5 (Simpson vs Simpson, abs_error 0 por
     construcao) -> identificado, publicado e REBAIXADO a decorativo.
  4. ARTEFATOS V17 SEM MANIFESTO -> o V18 publica os SHA-256 dos tres
     artefatos V17 como ancora (comparavel com --prev-v18 em auditorias
     futuras).

Uso:
    python AUDITAR_E8K_V18_TESTEMUNHAS.py \
        --v17-report RESULTADO_E8K_V17.json \
        --v17-parecer PARECER_TECNICO_V17.md \
        --v17-source AUDITAR_E8K_V17_NORMALIZACAO.py \
        [--e7 ZIP_E7] [--prev-v18 RESULTADO_E8K_V18.json] \
        [--out RESULTADO_E8K_V18.json] [--emit-sums] [--strict]

Codigos de saida: 0 = auditoria limpa; 2 = gate reprovado (em --strict,
pendencia conta como reprovado); 3 = erro de integracao. Stdlib pura.

Limites honestos: o V18 NAO reexecuta o pipeline contra os ZIPs
(E7/Fase4/V16) nem recalcula o LDU da V16; audita o que o V17 PUBLICOU.
Verificacoes float sao sanity; as exatas sao as racionais (suporte,
polos, hashes). Nao prova Guinand-Weil, nao certifica RH.
"""

import argparse
import ast
import copy
import hashlib
import json
import math
import re
import sys
import zipfile
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path

SCHEMA_V18 = "NUMERIA.E8K.V18.WITNESS_INDEPENDENT_REPLAY.1"
SCHEMA_V17 = "NUMERIA.E8K.V17.EXPLICIT_FORMULA_SPECIALIZATION.1"
STATUS_V17 = ("FOURIER_AND_POLE_SPECIALIZATION_AUDITED_FULL_EXPLICIT_"
              "FORMULA_EXTERNAL")
AVISO = ("hashes sem assinacao sao testes de integridade relativa, "
         "nao autenticacao por terceiros")

# ancoras da cadeia V10-V13, revalidadas em chumbo duplo
ANCORA_E7 = ("4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f"
             "387110b02a7")
ANCORA_B5 = ("a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b"
             "7cd53d3df9")

L = Q(9, 20)
DOIS_L = Q(9, 10)

# strings canonicas que o relatorio V17 deve declarar VERBATIM
CANONICO = {
    "sine_transform.w_k": "k*pi/(2L)",
    "sine_transform.F_formula":
        "w_k/sqrt(L)*(exp(-i*r*L)-(-1)^k*exp(+i*r*L))/(w_k^2-r^2), "
        "removable at r=+/-w_k",
    "sine_transform.amplitude_far":
        "k^2*pi^2/L^3 * (cos(L*r)^2 if k odd else sin(L*r)^2)"
        "/(r^2-w_k^2)^2",
    "sine_transform.amplitude_near":
        "4*L*w_k^2*sinc(L*(r-w_k))^2/(r+w_k)^2",
    "sine_transform.amplitude_at_resonance": "L",
    "sine_transform.orthonormality": "int_-L^L f_i(t) f_j(t) dt = delta_ij",
    "pole_moments.identity":
        "M_{k,s}=[2*pi*sqrt(L)*k/(k²*pi²+L²)]*(exp(-s*L/2)-(-1)^k "
        "exp(s*L/2)), s=+/-1",
    "pole_moments.cancel_coefficient":
        "B_{p,k}=-k/p*(p²*pi²+L²)/(k²*pi²+L²), same parity",
    "fourier_convention.forward":
        "F_k(r)=int_-L^L f_k(t) exp(+i*r*t) dt",
    "fourier_convention.inverse":
        "C_ij(u)=1/(2*pi)*int_R H_ij(r) exp(-i*r*u) dr",
    "fourier_convention.spectral_test":
        "H_ij(r)=1/2*(F_i(r)*F_j(-r)+F_j(r)*F_i(-r))",
    "fourier_convention.archimedean":
        "(1/pi)*int_0^infty H_ij(r)*Omega(r) dr (Omega even)",
    "fourier_convention.pole_term": "H_ij(i/2)+H_ij(-i/2)",
    "fourier_convention.prime_term": "-2*(log 2)/sqrt 2*C_ij(log 2)",
    "fourier_convention.regularity":
        "F_i entire and O(r^-2) on each horizontal strip; "
        "H_ij O(r^-4) on each strip",
    "ast_b5_eval_point": "amp2(r, k) * omega(r) / PI",
}


class Gates:
    def __init__(self):
        self.items = {}

    def ok(self, nome, detalhe=""):
        self.items[nome] = {"status": "PASS", "detail": detalhe}

    def fail(self, nome, detalhe=""):
        self.items[nome] = {"status": "FAIL", "detail": detalhe}

    def skip(self, nome, detalhe=""):
        self.items[nome] = {"status": "SKIP", "detail": detalhe}

    def pendencia(self, nome, detalhe, strict):
        if strict:
            self.fail(nome, detalhe + " [strict]")
        else:
            self.items[nome] = {"status": "PENDENCIA", "detail": detalhe}

    def counts(self):
        c = {"PASS": 0, "FAIL": 0, "SKIP": 0, "PENDENCIA": 0}
        for v in self.items.values():
            c[v["status"]] = c.get(v["status"], 0) + 1
        return c


def sha256_arquivo(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for bloco in iter(lambda: f.read(65536), b""):
            h.update(bloco)
    return h.hexdigest()


def para_frac(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, int):
        return Q(v)
    if isinstance(v, str):
        try:
            return Q(v.strip())
        except Exception:
            return None
    return None


def campo(doc, chave):
    """Acessa 'secao.campo' do relatorio V17."""
    secao, sub = chave.split(".", 1)
    obj = doc.get(secao)
    if not isinstance(obj, dict):
        return None
    return obj.get(sub)


# ----------------------------------------------------------------------------
# Reproducao independente: serie atanh com cauda majorante geometrica
# ----------------------------------------------------------------------------

def atanh_interval_indep(q, n):
    """2*atanh(q), q racional em (0,1).

    Soma por ordem DECRESCENTE de j (o V17 soma crescente) + cauda
    geometrica majorante: para j >= n, q^(2j+1)/(2j+1) <=
    q^(2j+1)/(2n+1), logo cauda <= q^(2n+1)/((2n+1)(1-q^2)).
    Raciocinio elementar e independente do V17; resultado exato.
    """
    assert 0 < q < 1
    q2 = q * q
    soma = Q(0)
    for j in range(n - 1, -1, -1):
        soma += 2 * q ** (2 * j + 1) / Q(2 * j + 1)
    cauda = 2 * q ** (2 * n + 1) / (Q(2 * n + 1) * (1 - q2))
    return soma, soma + cauda


# ----------------------------------------------------------------------------
# Detectores: cada um retorna (ok, problemas)
# ----------------------------------------------------------------------------

def detectar_estrutura(doc):
    problemas = []
    gates = doc.get("gates", {})
    if len(gates) != 9 or not all(v is True for v in gates.values()):
        problemas.append("gates: esperados 9, todos True")
    negs = doc.get("negative_controls", {})
    if len(negs) != 6 or not all(v is True for v in negs.values()):
        problemas.append("negative_controls: esperados 6, todos True")
    sn = doc.get("independent_numeric_sanity_only", {})
    casos = (sn.get("fourier_direct_vs_closed", [])
             + sn.get("pole_moments_direct_vs_closed", [])
             + sn.get("orthogonality_direct", []))
    if len(casos) != 25:
        problemas.append("casos de sanity: esperados 25, ha %d" % len(casos))
    if casos and not all(c.get("agrees") is True for c in casos):
        problemas.append("caso de sanity sem agrees=True")
    if sn.get("all_sanity_cases_agree") is not True:
        problemas.append("all_sanity_cases_agree ausente")
    erros = [float(c["abs_error"]) for c in casos if "abs_error" in c]
    if erros:
        if abs(max(erros) - float(sn.get("max_error_float", 0))) > \
                1e-15 * max(erros):
            problemas.append("max_error_float nao e o maximo dos casos")
        pior = max(casos, key=lambda c: float(c.get("abs_error", 0)))
        if not (pior.get("k") == 10 and pior.get("r") == "33"):
            problemas.append("pior caso alegado (k=10, r=33) nao e o pior")
    if doc.get("schema") != SCHEMA_V17:
        problemas.append("schema V17 inesperado")
    if doc.get("status") != STATUS_V17:
        problemas.append("status V17 inesperado")
    if not doc.get("not_certified"):
        problemas.append("not_certified ausente/vazio")
    par = doc.get("parameters", {})
    if par.get("L") != "9/20" or par.get("double_L") != "9/10":
        problemas.append("parametros L/2L divergem do parecer")
    if par.get("support_centered") != ["-9/20", "9/20"]:
        problemas.append("suporte centrado diverge")
    if par.get("N") != [6, 8, 10] or par.get("projection_dimension") != "N-2":
        problemas.append("parametros N/projecao divergem")
    return (not problemas), problemas


def detectar_suporte(doc):
    problemas = []
    w = doc.get("support_witness", {})
    a, b = atanh_interval_indep(Q(1, 3), 70)
    c, d = atanh_interval_indep(Q(1, 2), 100)
    pares = (("ln2_interval[0]", a), ("ln2_interval[1]", b),
             ("ln3_interval[0]", c), ("ln3_interval[1]", d))
    for nome, calc in pares:
        pub = w.get(nome.split("[")[0])
        idx = int(nome.split("[")[1][:-1]) if "[" in nome else None
        pub_valor = pub[idx] if isinstance(pub, list) and idx is not None \
            else None
        if pub_valor is None or str(calc) != str(pub_valor):
            problemas.append("testemunha %s nao reproduzida bit a bit"
                             % nome)
    dois_l = para_frac(w.get("support_double_L", ""))
    if dois_l != DOIS_L:
        problemas.append("support_double_L != 9/10")
    if not (b < dois_l < c):
        problemas.append("desigualdade estrita log2 < 2L < log3 falhou")
    if w.get("strict_inequality") is not True:
        problemas.append("strict_inequality nao declarada True")
    mid2 = float((a + b) / 2)
    mid3 = float((c + d) / 2)
    if abs(mid2 - math.log(2)) > 1e-12 or abs(mid3 - math.log(3)) > 1e-12:
        problemas.append("intervalos nao cercam ln2/ln3 (sanity float)")
    return (not problemas), problemas


def detectar_polos(doc):
    """Identidade exata em amostras racionais de X = pi^2 + colunas."""
    problemas = []
    amostras = [Q(2), Q(7, 3), Q(5), Q(13, 2), Q(10), Q(47, 5),
                Q(3, 7), Q(11, 13)]
    for k in range(3, 11):
        p = 1 if k % 2 else 2
        if k % 2 != p % 2:
            problemas.append("pivo de paridade errado (k=%d)" % k)
        for x in amostras:
            den_k = k * k * x + L * L
            den_p = p * p * x + L * L
            beta = -Q(k, p) * den_p / den_k
            soma = Q(k) / den_k + beta * Q(p) / den_p
            if soma != 0:
                problemas.append("identidade C_k+beta*C_p != 0 (k=%d, X=%s)"
                                 % (k, x))
    linhas = doc.get("pole_moments", {}).get("two_pole_cancellations", [])
    if len(linhas) != 8:
        problemas.append("esperadas 8 colunas de cancelamento, ha %d"
                         % len(linhas))
    for r in linhas:
        k = r.get("mode")
        p_esp = 1 if k % 2 else 2
        if (r.get("pivot") != p_esp
                or r.get("parity") != ("odd" if k % 2 else "even")
                or r.get("numerator_coeff_X") != "0"
                or r.get("numerator_constant") != "0"
                or r.get("plus_and_minus_moments_cancel_exactly")
                is not True):
            problemas.append("coluna de polo k=%s inconsistente" % k)
    return (not problemas), problemas


def quadratura_simpson(f, a, b, n=4096):
    if n % 2:
        raise ValueError("n deve ser par")
    dx = (b - a) / n
    return dx / 3 * (f(a) + f(b)
                     + 4 * sum(f(a + j * dx) for j in range(1, n, 2))
                     + 2 * sum(f(a + j * dx) for j in range(2, n, 2)))


def momento(k, s):
    """Integral de f_k(t)*exp(s*t/2) em [-L, L]: fechada e por Simpson.

    Fechada: derivada por integracao elementar (parecer, secao 5);
    Simpson: caminho independente da formula fechada.
    """
    lf = 9.0 / 20.0
    w = k * math.pi / (2 * lf)
    fechado = (w / math.sqrt(lf)) / (w * w + 0.25) * (
        math.exp(-s * lf / 2) - ((-1) ** k) * math.exp(s * lf / 2))
    base = (lambda t:
            math.sin(k * math.pi * (t + lf) / (2 * lf))
            / math.sqrt(lf) * math.exp(s * t / 2))
    quad = quadratura_simpson(base, -lf, lf)
    return fechado, quad


def detectar_cancelamento_numerico():
    """O ZERO do cancelamento, OBSERVADO com dados independentes.

    Nada de A-A simbolico: avalia M_{k,s} (fechada e por quadratura)
    e exige M_k + beta*M_p ~ 0 relativo a escala, para s=+-1.
    """
    problemas = []
    lf = 9.0 / 20.0
    pior_fech = 0.0
    pior_quad = 0.0
    for k in range(3, 11):
        p = 1 if k % 2 else 2
        beta = -(k / p) * ((p * p * math.pi ** 2 + lf * lf)
                           / (k * k * math.pi ** 2 + lf * lf))
        for s in (-1, 1):
            fech_k, quad_k = momento(k, s)
            fech_p, quad_p = momento(p, s)
            if abs(fech_k - quad_k) > 1e-9 * max(abs(fech_k), 1e-30):
                problemas.append("momento fechado != quadratura (k=%d,s=%+d)"
                                 % (k, s))
            escala = max(abs(fech_k), abs(quad_k), 1e-30)
            res_fech = abs(fech_k + beta * fech_p) / escala
            res_quad = abs(quad_k + beta * quad_p) / escala
            pior_fech = max(pior_fech, res_fech)
            pior_quad = max(pior_quad, res_quad)
            if res_fech > 1e-9 or res_quad > 1e-7:
                problemas.append("cancelamento k=%d s=%+d nao observado "
                                 "(fech=%.2e, quad=%.2e)"
                                 % (k, s, res_fech, res_quad))
    info = {"pior_residuo_fechado": "%.3e" % pior_fech,
            "pior_residuo_quadratura": "%.3e" % pior_quad}
    return (not problemas), problemas, info


def detectar_amplitude():
    """Tres formas da amplitude + valor removivel na ressonancia."""
    problemas = []
    lf = 0.45
    pior = 0.0
    for k in (1, 2, 3, 4, 7, 10):
        for r in (0.1, 0.5, 1.3, 5.2, 20.1, 33.7):
            w = k * math.pi / (2 * lf)
            trig = math.cos(lf * r) if k % 2 else math.sin(lf * r)
            far = (k * k * math.pi ** 2 / lf ** 3 * trig ** 2
                   / (r * r - w * w) ** 2)
            x = lf * (r - w)
            near = (4 * lf * w * w * (math.sin(x) / x) ** 2
                    / (r + w) ** 2)
            direta = ((w * w / lf)
                      * (2 - 2 * ((-1) ** k) * math.cos(2 * r * lf))
                      / (w * w - r * r) ** 2)
            rel = max(abs(far - near), abs(far - direta),
                      abs(near - direta)) / max(abs(direta), 1e-300)
            pior = max(pior, rel)
    if pior > 1e-9:
        problemas.append("formas far/near/direta divergem (rel=%.2e)" % pior)
    w5 = 5 * math.pi / (2 * lf)
    r5 = w5 * (1 + 1e-9)
    x5 = lf * (r5 - w5)
    near5 = 4 * lf * w5 * w5 * (math.sin(x5) / x5) ** 2 / (r5 + w5) ** 2
    if abs(near5 - lf) > 1e-6:
        problemas.append("amplitude removivel na ressonancia != L")
    return (not problemas), problemas, {"max_desvio_relativo": "%.3e" % pior}


def detectar_canonicos(doc):
    """Toda formula publica do V17 deve casar verbatim com o canonico."""
    problemas = []
    for chave, esperado in CANONICO.items():
        if chave == "ast_b5_eval_point":
            if doc.get("ast_b5_eval_point") != esperado:
                problemas.append("ast_b5_eval_point diverge do canonico")
            continue
        valor = campo(doc, chave)
        if valor != esperado:
            problemas.append("formula %s diverge do canonico" % chave)
    return (not problemas), problemas


def detectar_condicionais(doc):
    problemas = []
    conds = doc.get("v16_inherited_finite_conditionals", [])
    ns = sorted(c.get("N") for c in conds)
    if ns != [6, 8, 10]:
        problemas.append("N herdados != [6, 8, 10]")
    for c in conds:
        n = c.get("N")
        for exato_c, aprox_c in (("delta_exact", "delta_approx"),
                                 ("min_shifted_pivot_exact",
                                  "min_shifted_pivot_approx")):
            exato = para_frac(c.get(exato_c, ""))
            try:
                aprox = float(c.get(aprox_c))
            except (TypeError, ValueError):
                aprox = None
            if exato is None or aprox is None:
                problemas.append("N=%s: %s/%s ilegivel" % (n, exato_c,
                                                           aprox_c))
                continue
            if abs(float(exato) - aprox) > 1e-12:
                problemas.append("N=%s: %s inconsistente com %s"
                                 % (n, exato_c, aprox_c))
        if c.get("positive_definite_if_analytic_enclosures_hold") is not True:
            problemas.append("N=%s: condicional de positividade ausente" % n)
        if c.get("projected_dimension") != (n or 0) - 2:
            problemas.append("N=%s: projecao != N-2" % n)
    return (not problemas), problemas


def detectar_hashes_fonte(doc, fonte_v17):
    """Constantes *_SHA do codigo V17 vs relatorio vs ancoras da cadeia."""
    problemas = []
    padrao = re.compile(r"^([A-Z0-9_]+)_SHA\s*=\s*['\"]([0-9a-f]{64})['\"]",
                        re.M)
    consts = dict(padrao.findall(fonte_v17))
    mapa = {"E7": "e7_zip_sha256", "PHASE4": "phase4_zip_sha256",
            "B5": "b5_source_sha256", "PHASE3": "phase3_source_sha256",
            "V16_REPORT": "v16_report_sha256",
            "V16_SOURCE": "v16_source_sha256"}
    for k, campo_hash in mapa.items():
        if k not in consts:
            problemas.append("constante %s_SHA ausente do codigo V17" % k)
            continue
        if doc.get("source_hashes", {}).get(campo_hash) != consts[k]:
            problemas.append("hash %s diverge entre codigo e relatorio"
                             % campo_hash)
    if consts.get("E7") != ANCORA_E7:
        problemas.append("ancora E7 diverge da cadeia V10-V13")
    if consts.get("B5") != ANCORA_B5:
        problemas.append("ancora B5 diverge da cadeia V10-V13")
    return (not problemas), problemas


def detectar_parecer(texto, doc):
    problemas = []
    if "RH permanece não demonstrada" not in texto:
        problemas.append("parecer sem a clausula de nao-demonstracao")
    if STATUS_V17 not in texto:
        problemas.append("status do relatorio ausente do parecer")
    if ("2,86e−13" not in texto) and ("2,86e-13" not in texto):
        problemas.append("parecer nao cita o erro maximo 2.86e-13")
    for frase in ("25 verificações", "8/8 cancelamentos", "6/6 mutações",
                  "Guinand", "Chirre", "Bombieri"):
        if frase not in texto:
            problemas.append("parecer nao cita: %r" % frase)
    return (not problemas), problemas


def casos_decorativos(doc):
    """Casos em que 'direto' e 'fechado' sao o MESMO caminho de Simpson
    (r == w_k: F_closed retorna a propria quadratura). abs_error 0 por
    construcao — evidencia decorativa, nao verificacao."""
    lf = 0.45
    saida = []
    for caso in doc.get("independent_numeric_sanity_only", {}).get(
            "fourier_direct_vs_closed", []):
        try:
            k = int(caso["k"])
            r = float(caso["r"])
        except (KeyError, TypeError, ValueError):
            continue
        w = k * math.pi / (2 * lf)
        if abs(r - w) < 1e-9 and caso.get("abs_error") == "0.000000e+00":
            saida.append({"k": k, "r": caso.get("r"),
                          "motivo": "F_closed usa Simpson quando "
                                    "|r-w_k|<1e-10: compara consigo mesmo"})
    return saida


# ----------------------------------------------------------------------------
# Selftest: mutacoes REAIS sobre o relatorio publicado
# ----------------------------------------------------------------------------

def detectar_doc(doc, fonte_v17):
    """Todos os detectores dependentes do relatorio, agrupados."""
    problemas = []
    for nome, det in (("estrutura", detectar_estrutura),
                      ("suporte", detectar_suporte),
                      ("polos", detectar_polos),
                      ("canonicos", detectar_canonicos),
                      ("condicionais", detectar_condicionais),
                      ("hashes_fonte", lambda d:
                       detectar_hashes_fonte(d, fonte_v17))):
        ok, probs = det(doc)
        problemas.extend("%s: %s" % (nome, p) for p in probs)
    return problemas


def troca_digito(s):
    """Incrementa o PRIMEIRO digito do numerador (muda o valor ~10%)."""
    num, den = s.split("/")
    primeiro = num[0]
    novo = "0" if primeiro == "9" else str(int(primeiro) + 1)
    return novo + num[1:] + "/" + den


def selftest_mutacoes(doc, fonte_v17):
    """Cada mutacao altera o relatorio e DEVE ser detectada.
    Mutacao nao detectada = defeito no detector = selftest reprova."""
    casos = []

    def mut(nome, aplicar):
        m = copy.deepcopy(doc)
        aplicar(m)
        detectado = len(detectar_doc(m, fonte_v17)) > 0
        casos.append({"mutacao": nome, "detectada": detectado})

    mut("testemunha ln2 adulterada (primeiro digito)",
        lambda m: m["support_witness"]["ln2_interval"].__setitem__(
            0, troca_digito(m["support_witness"]["ln2_interval"][0])))
    mut("desigualdade de suporte quebrada (ln2 recebe fração de ln3)",
        lambda m: m["support_witness"]["ln2_interval"].__setitem__(
            1, m["support_witness"]["ln3_interval"][0]))
    mut("strict_inequality removida",
        lambda m: m["support_witness"].__setitem__("strict_inequality",
                                                   False))
    mut("beta com sinal invertido no parecer publicado",
        lambda m: m["pole_moments"].__setitem__(
            "cancel_coefficient",
            m["pole_moments"]["cancel_coefficient"].replace("-k/p",
                                                            "+k/p")))
    mut("fator primo -2 vira -1 (o controle falso do V17, agora real)",
        lambda m: m["fourier_convention"].__setitem__(
            "prime_term",
            m["fourier_convention"]["prime_term"].replace("-2*(log 2)",
                                                          "-1*(log 2)")))
    mut("denominador da regularização trocado ((r+w_k) -> (r-w_k))",
        lambda m: m["sine_transform"].__setitem__(
            "amplitude_near",
            m["sine_transform"]["amplitude_near"].replace("r+w_k", "r-w_k")))
    mut("paridade do pivô violada (k=3 com pivô 2)",
        lambda m: m["pole_moments"]["two_pole_cancellations"][0].__setitem__(
            "pivot", 2))
    mut("numerador de cancelamento adulterado ('0' -> '1')",
        lambda m: m["pole_moments"]["two_pole_cancellations"][0].__setitem__(
            "numerator_coeff_X", "1"))
    mut("delta exato da V16 adulterado (primeiro digito)",
        lambda m: m["v16_inherited_finite_conditionals"][0].__setitem__(
            "delta_exact",
            troca_digito(
                m["v16_inherited_finite_conditionals"][0]["delta_exact"])))
    mut("max_error_float subestimado",
        lambda m: m["independent_numeric_sanity_only"].__setitem__(
            "max_error_float", 1e-20))
    mut("status inflado para RIEMANN_HYPOTHESIS_CERTIFIED",
        lambda m: m.__setitem__("status", "RIEMANN_HYPOTHESIS_CERTIFIED"))
    mut("not_certified esvaziado",
        lambda m: m.__setitem__("not_certified", []))
    mut("hash de fonte E7 adulterado no relatorio",
        lambda m: m["source_hashes"].__setitem__(
            "e7_zip_sha256",
            ("5" if m["source_hashes"]["e7_zip_sha256"][0] != "5" else "6")
            + m["source_hashes"]["e7_zip_sha256"][1:]))
    return casos


# ----------------------------------------------------------------------------
# Opcional: AST da fonte B5, se o ZIP E7 for fornecido
# ----------------------------------------------------------------------------

def ast_eval_point(fonte_bytes):
    try:
        arvore = ast.parse(fonte_bytes.decode("utf-8"))
    except SyntaxError as exc:
        return False, "ilegivel: %r" % exc
    funcs = {n.name: n for n in arvore.body
             if isinstance(n, ast.FunctionDef)}
    if "eval_point" not in funcs:
        return False, "funcao eval_point ausente"
    rets = [n for n in funcs["eval_point"].body
            if isinstance(n, ast.Return)]
    if len(rets) != 1:
        return False, "esperado um unico return em eval_point"
    return ast.unparse(rets[0].value) == CANONICO["ast_b5_eval_point"], \
        ast.unparse(rets[0].value)


# ----------------------------------------------------------------------------
# Principal
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="V18: replay independente das testemunhas V17 "
                    "+ mutacoes reais")
    ap.add_argument("--v17-report", type=Path, required=True)
    ap.add_argument("--v17-parecer", type=Path, required=True)
    ap.add_argument("--v17-source", type=Path, required=True)
    ap.add_argument("--e7", type=Path, default=None,
                    help="ZIP E7 (opcional; habilita o gate AST do B5)")
    ap.add_argument("--prev-v18", type=Path, default=None,
                    help="relatorio V18 anterior (cotejo de artefatos)")
    ap.add_argument("--out", type=Path,
                    default=Path("RESULTADO_E8K_V18.json"))
    ap.add_argument("--emit-sums", action="store_true")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    inicio = datetime.now(timezone.utc).isoformat()
    gates = Gates()
    try:
        doc = json.loads(args.v17_report.read_text(encoding="utf-8"))
        parecer = args.v17_parecer.read_text(encoding="utf-8")
        fonte_v17 = args.v17_source.read_text(encoding="utf-8")
    except Exception as exc:
        print("ERRO DE INTEGRACAO: artefatos V17 ilegiveis: %r" % (exc,),
              file=sys.stderr)
        return 3

    # ---- gates dependentes do relatorio ----
    for nome, det, detalhe_ok in (
            ("estrutura_relatorio_v17", detectar_estrutura,
             "9 gates, 6 negativos, 25 sanity, atribuicao do erro maximo"),
            ("testemunhas_suporte_reproduzidas", detectar_suporte,
             "ln2/ln3 bit a bit; log2 < 9/10 < log3 exato"),
            ("polos_identidade_e_colunas", detectar_polos,
             "C_k+beta*C_p=0 em 8 amostras de X por k; 8 colunas"),
            ("formulas_canonicas_declaradas", detectar_canonicos,
             "%d strings verbatim" % (len(CANONICO) - 1)),
            ("condicionais_v16_consistentes", detectar_condicionais,
             "N=6,8,10: racionais vs floats; condicional preservada"),
            ("hashes_codigo_vs_relatorio", lambda d:
             detectar_hashes_fonte(d, fonte_v17),
             "6 constantes *_SHA conferem; ancoras E7/B5 da cadeia")):
        ok, probs = det(doc)
        if ok:
            gates.ok(nome, detalhe_ok)
        else:
            gates.fail(nome, "; ".join(probs[:6]))

    # ---- gates independentes do relatorio ----
    ok, probs, info_canc = detectar_cancelamento_numerico()
    if ok:
        gates.ok("cancelamento_zero_observado",
                 "M_k+beta*M_p~0 (fechada %.0e, quadratura %.0e)"
                 % (float(info_canc["pior_residuo_fechado"]),
                    float(info_canc["pior_residuo_quadratura"])))
    else:
        gates.fail("cancelamento_zero_observado", "; ".join(probs[:6]))

    ok, probs, info_amp = detectar_amplitude()
    if ok:
        gates.ok("amplitude_formas_concordam",
                 "far/near/direta ate %s; ressonancia -> L"
                 % info_amp["max_desvio_relativo"])
    else:
        gates.fail("amplitude_formas_concordam", "; ".join(probs[:6]))

    ok, probs = detectar_parecer(parecer, doc)
    if ok:
        gates.ok("parecer_consistente",
                 "classificacao, erro maximo e clausulas citadas conferem")
    else:
        gates.fail("parecer_consistente", "; ".join(probs[:6]))

    # ---- caso ressonante decorativo: achado honesto, nao reprovacao ----
    decorativos = casos_decorativos(doc)
    if decorativos:
        gates.pendencia("ressonancia_declarada_decorativa",
                        "%d caso(s) k=%s com abs_error 0 por construcao "
                        "(Simpson vs Simpson); rebaixado no relatorio V18"
                        % (len(decorativos),
                           ",".join(str(c["k"]) for c in decorativos)),
                        args.strict)
    else:
        gates.ok("ressonancia_declarada_decorativa",
                 "nenhum caso direto==fechado por construcao")

    # ---- selftest: mutacoes reais ----
    mutacoes = selftest_mutacoes(doc, fonte_v17)
    n_det = sum(1 for c in mutacoes if c["detectada"])
    if n_det == len(mutacoes) and mutacoes:
        gates.ok("selftest_mutacoes_reais",
                 "%d/%d mutacoes detectadas pelos detectores"
                 % (n_det, len(mutacoes)))
    else:
        gates.fail("selftest_mutacoes_reais",
                   "nao detectadas: %s" % ", ".join(
                       c["mutacao"] for c in mutacoes
                       if not c["detectada"]))

    # ---- manifesto: ancora publica dos artefatos V17 ----
    artefatos = {
        "v17_report_sha256": sha256_arquivo(args.v17_report),
        "v17_parecer_sha256": sha256_arquivo(args.v17_parecer),
        "v17_source_sha256": sha256_arquivo(args.v17_source),
    }
    if args.prev_v18:
        try:
            prev = json.loads(args.prev_v18.read_text(encoding="utf-8"))
            prev_hashes = prev.get("v17_artifact_sha256", {})
            diverg = [k for k, v in artefatos.items()
                      if prev_hashes.get(k) not in (None, v)]
            if diverg:
                gates.pendencia("manifesto_v17",
                               "artefatos mudaram desde a auditoria "
                               "anterior: %s" % ", ".join(diverg),
                               args.strict)
            else:
                gates.ok("manifesto_v17",
                         "3 artefatos inalterados desde a auditoria "
                         "anterior")
        except Exception as exc:
            gates.pendencia("manifesto_v17",
                            "relatorio anterior ilegivel: %r" % exc,
                            args.strict)
    else:
        gates.ok("manifesto_v17",
                 "3 artefatos V17 ancorados por hash (primeira publicacao)")

    # ---- AST opcional da fonte B5 ----
    if args.e7:
        try:
            with zipfile.ZipFile(args.e7) as zf:
                b5 = zf.read("fontes/AUDITOR_PONTOS_B5.py")
            if sha256_arquivo_bytes(b5) != ANCORA_B5:
                gates.fail("b5_ast_opcional", "SHA da fonte B5 divergente")
            else:
                ok, como = ast_eval_point(b5)
                if ok:
                    gates.ok("b5_ast_opcional", "eval_point: %s" % como)
                else:
                    gates.fail("b5_ast_opcional", "AST: %s" % como)
        except Exception as exc:
            gates.pendencia("b5_ast_opcional",
                            "ZIP E7 ilegivel: %r" % exc, args.strict)
    else:
        gates.skip("b5_ast_opcional",
                   "sem --e7; o gate AST do B5 fica por conta do V17")

    # ---- monta relatorio ----
    reprovado = any(v["status"] == "FAIL"
                    for v in gates.items.values())
    saida = {
        "schema": SCHEMA_V18,
        "utc": inicio,
        "rounding": "outward (racionais exatos; floats sao sanity)",
        "warning": AVISO,
        "status": ("AUDIT_FAILED" if reprovado else
                   "V17_WITNESSES_REPLAYED_MUTATIONS_DETECTED"),
        "gates": gates.items,
        "v17_artifact_sha256": artefatos,
        "mutations": mutacoes,
        "decorative_cases": decorativos,
        "cancelamento_numerico": info_canc,
        "amplitude": info_amp,
        "not_certified": [
            "os ZIPs E7/Fase4/V16 nao foram reexecutados; o V18 audita o "
            "relatorio publicado, nao a fonte",
            "verificacoes float sao sanity; as exatas sao as racionais "
            "(suporte, polos, hashes)",
            "LDU da V16 nao foi recalculado (sem as matrizes G/P/B)",
            "Guinand-Weil permanece teorema externo assumido",
            "mutacoes provam que os DETECTORES respondem a adulteracao "
            "do relatorio; nao substituem a execucao do V17 contra os "
            "ZIPs",
            "Hipotese de Riemann nao certificada",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(saida, indent=2, ensure_ascii=False, sort_keys=True)
        + "\n", encoding="utf-8")

    if args.emit_sums:
        linhas = ["%s  %s" % (artefatos["v17_report_sha256"],
                              args.v17_report.name),
                  "%s  %s" % (artefatos["v17_parecer_sha256"],
                              args.v17_parecer.name),
                  "%s  %s" % (artefatos["v17_source_sha256"],
                              args.v17_source.name)]
        Path("SHA256SUMS_V18.txt").write_text("\n".join(linhas) + "\n",
                                              encoding="utf-8")

    print("Gates: %s" % json.dumps(gates.counts()))
    for nome, v in gates.items.items():
        print("  [%s] %s: %s" % (v["status"], nome, v["detail"][:90]))
    print("Relatorio: %s" % args.out)
    return 2 if reprovado else 0


def sha256_arquivo_bytes(b):
    return hashlib.sha256(b).hexdigest()


if __name__ == "__main__":
    sys.exit(main())