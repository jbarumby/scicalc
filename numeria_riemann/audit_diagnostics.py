"""Diagnóstico de resultados NUMERIA-RIEMANN (sem assumir campos opcionais).

Uso: python -m numeria_riemann.audit_diagnostics RESULTADO.json --out RELATORIO.md
"""
import argparse
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path


def _numbers(values):
    result = []
    for value in values or []:
        try:
            result.append(Decimal(str(value)))
        except (InvalidOperation, ValueError, TypeError):
            result.append(None)
    return result


def summarize(report):
    status = report.get('status', 'não informado')
    parameters = report.get('parameters') or {}
    lowers = _numbers(report.get('interval_ldl_pivot_lower_bounds'))
    uppers = _numbers(report.get('interval_ldl_pivot_upper_bounds'))
    pivots = report.get('interval_ldl_pivots') or []
    lines = [
        '# Diagnóstico de auditoria intervalar', '',
        f"- Status declarado: `{status}`",
        f"- N: `{parameters.get('N', 'não informado')}`",
        f"- Células de quadratura: `{parameters.get('quadrature_cells', 'não informado')}`",
        f"- Dimensão projetada: `{report.get('projected_dimension', 'não informada')}`",
        f"- Pivôs intervalares registrados: `{len(pivots)}`",
        '', '## Limites inferiores dos pivôs', ''
    ]
    if lowers:
        for i, v in enumerate(lowers, 1):
            lines.append(f'- Pivô {i}: {v if v is not None else "inválido"}')
    else:
        lines.append('- Não fornecidos.')
    lines += ['', '## Limites superiores dos pivôs', '']
    if uppers:
        for i, v in enumerate(uppers, 1):
            lines.append(f'- Pivô {i}: {v if v is not None else "inválido"}')
    else:
        lines.append('- Campo `interval_ldl_pivot_upper_bounds` ausente; não é possível informar larguras dos pivôs.')
    lines += ['', '## Interpretação', '']
    if status.startswith('inconclusive'):
        lines.append('A execução é **inconclusiva**. Nenhuma certificação de positividade decorre deste relatório.')
    else:
        lines.append('Este relatório reproduz o status declarado pelo arquivo; não valida independentemente a certificação.')
    if any(v is not None and v <= 0 for v in lowers):
        lines.append('Há limite inferior não positivo; o teste LDL intervalar não certificou todos os pivôs registrados como estritamente positivos.')
    if len(lowers) < int(report.get('projected_dimension') or 0):
        lines.append('A lista de limites inferiores é menor que a dimensão projetada; pode representar interrupção antecipada da fatoração.')
    lines += ['', 'Escopo: matrizes finitas; isto não constitui demonstração da hipótese de Riemann.', '']
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    content = summarize(json.loads(args.report.read_text(encoding='utf-8')))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(content, encoding='utf-8')
        print(f'Relatório salvo: {args.out}')
    else:
        print(content)


if __name__ == '__main__':
    main()
