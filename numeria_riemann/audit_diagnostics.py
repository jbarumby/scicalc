"""Diagnóstico de resultados NUMERIA-RIEMANN (sem assumir campos opcionais).

Uso:
    python -m numeria_riemann.audit_diagnostics RESULTADO.json [OUTRO.json ...] [--out RELATORIO.md | --out-dir PASTA]

Lê relatórios schema 1 (LDL intervalar) e schemas 2-3 (LDL + teste de margem tipo Rump).
Termina com código 1 se algum relatório for internamente inconsistente
(ex.: status "certificado" com pivô de limite inferior não positivo).
"""
import argparse
import json
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

CERT = "certified_positive_definite_finite_matrix"


def _numbers(values):
    """Converte para Decimal; None para valores inválidos ou não finitos (NaN/Inf)."""
    result = []
    for value in values or []:
        try:
            d = Decimal(str(value))
            result.append(d if d.is_finite() else None)
        except (InvalidOperation, ValueError, TypeError):
            result.append(None)
    return result


def _fmt(v):
    return "inválido" if v is None else str(v)


def check_consistency(report):
    """Lista de inconsistências internas do relatório (vazia se coerente)."""
    problems = []
    status = str(report.get("status", ""))
    dim = int(report.get("projected_dimension") or 0)
    lowers = _numbers(report.get("interval_ldl_pivot_lower_bounds"))
    certs = report.get("certificates")
    if status == CERT:
        if certs is not None:
            if not any((c or {}).get("passed") for c in certs.values()):
                problems.append("status certificado, mas nenhum certificado em `certificates` passou.")
        else:  # schema 1: só o LDL intervalar certifica
            if len(lowers) < dim or any(v is None or v <= 0 for v in lowers):
                problems.append("status certificado, mas os pivôs não são todos estritamente positivos.")
    elif status.startswith("inconclusive") and certs is not None:
        if any((c or {}).get("passed") for c in certs.values()):
            problems.append("status inconclusivo, mas há certificado que passou.")
    return problems


def summarize(report):
    status = report.get("status", "não informado")
    parameters = report.get("parameters") or {}
    lowers = _numbers(report.get("interval_ldl_pivot_lower_bounds"))
    uppers = _numbers(report.get("interval_ldl_pivot_upper_bounds"))
    pivots = report.get("interval_ldl_pivots") or []
    dim = int(report.get("projected_dimension") or 0)
    certs = report.get("certificates") or {}
    lines = [
        "# Diagnóstico de auditoria intervalar", "",
        f"- Status declarado: `{status}`",
        f"- N: `{parameters.get('N', 'não informado')}`",
        f"- Integração: `{parameters.get('archimedean', 'rangesum (schema 1)')}`",
        f"- Células de quadratura: `{parameters.get('quadrature_cells') or 'não se aplica'}`",
        f"- Dimensão projetada: `{report.get('projected_dimension', 'não informada')}`",
        f"- Pivôs intervalares registrados: `{len(pivots)}` de `{dim}`",
    ]
    if certs:
        lines += ["", "## Certificados", ""]
        for name, info in certs.items():
            info = info or {}
            extra = f", delta={info['delta_upper']}" if info.get("delta_upper") else ""
            fail = f", falhou no pivô {info['failed_at_pivot']}" if info.get("failed_at_pivot") else ""
            lines.append(f"- `{name}`: {'passou' if info.get('passed') else 'não passou'}{extra}{fail}")
    lines += ["", "## Pivôs (limite inferior / superior / largura)", ""]
    if lowers:
        for i, lo in enumerate(lowers, 1):
            up = uppers[i - 1] if i - 1 < len(uppers) else None
            width = (up - lo) if (up is not None and lo is not None) else None
            lines.append(f"- Pivô {i}: {_fmt(lo)} / {_fmt(up) if uppers else 'ausente'} / "
                         f"{_fmt(width) if width is not None else 'n/d'}")
    else:
        lines.append("- Não fornecidos.")
    if not uppers:
        lines.append("- Campo `interval_ldl_pivot_upper_bounds` ausente; larguras dos pivôs indisponíveis.")
    lines += ["", "## Interpretação", ""]
    if str(status).startswith("inconclusive"):
        lines.append("A execução é **inconclusiva**. Nenhuma certificação de positividade decorre deste relatório.")
    elif status == CERT:
        lines.append("O relatório declara certificação para a **matriz finita restrita** dos parâmetros acima; "
                     "este diagnóstico não revalida a aritmética de intervalos.")
    else:
        lines.append("Este relatório reproduz o status declarado pelo arquivo; não valida independentemente a certificação.")
    if any(v is not None and v <= 0 for v in lowers):
        lines.append("Há limite inferior não positivo; o LDL intervalar não certificou todos os pivôs como estritamente positivos.")
    if len(lowers) < dim:
        lines.append("A lista de limites inferiores é menor que a dimensão projetada; o LDL parou antes do fim.")
    diag = (report.get("diagnostic_only_midpoint_eigenvalues") or {}).get("values")
    if diag:
        lines.append(f"Diagnóstico NÃO rigoroso (autovalores do ponto médio): {', '.join(diag)}.")
    problems = check_consistency(report)
    if problems:
        lines += ["", "## ⚠ Inconsistências", ""] + [f"- {p}" for p in problems]
    lines += ["", "Escopo: matrizes finitas; isto não constitui demonstração da hipótese de Riemann.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("reports", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, help="arquivo de saída (apenas com um relatório)")
    parser.add_argument("--out-dir", type=Path, help="pasta de saída: um .md por relatório")
    args = parser.parse_args(argv)
    if args.out and len(args.reports) > 1:
        parser.error("--out aceita um único relatório; use --out-dir")
    exit_code = 0
    for path in args.reports:
        report = json.loads(path.read_text(encoding="utf-8"))
        content = summarize(report)
        if check_consistency(report):
            exit_code = 1
        target = args.out or (args.out_dir / f"diagnostics-{path.stem}.md" if args.out_dir else None)
        if target:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
            print(f"Relatório salvo: {target}")
        else:
            print(content)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
