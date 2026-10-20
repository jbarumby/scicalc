"""Offline diagnostics for JSON reports produced by weil_interval_audit.py.

Does not require python-flint.  Parses certified directed decimal endpoint strings.
"""
import argparse
from decimal import Decimal
import json
from pathlib import Path


def summarize(report):
    lower = [Decimal(s) for s in report['interval_ldl_pivot_lower_bounds']]
    upper = [Decimal(s) for s in report['interval_ldl_pivot_upper_bounds']]
    assert len(lower) == len(upper)
    lines = [f"# N={report['parameters']['N']} — interval LDL diagnostics", "",
             f"- Status: `{report['status']}`",
             f"- Dimension: {report['projected_dimension']}",
             f"- Cells: {report['parameters']['quadrature_cells']}",
             f"- Computed pivots: {len(lower)}", "",
             "| Pivot | Lower bound | Upper bound | Width | Relative width |",
             "|---:|---:|---:|---:|---:|"]
    for i, (lo, hi) in enumerate(zip(lower, upper), 1):
        mid = (lo + hi) / 2
        width = hi - lo
        rel = width / abs(mid) if mid else 'undefined'
        lines.append(f"| {i} | {lo:.9f} | {hi:.9f} | {width:.9f} | {rel if isinstance(rel, str) else format(rel, '.6g')} |")
    if len(lower) < report['projected_dimension']:
        lines += ["", f"Stopped before pivot {len(lower)+1}; this is not a negative-eigenvalue certificate."]
    if len(lower) == report['projected_dimension'] and all(x > 0 for x in lower):
        lines += ["", "All reported LDL pivot lower bounds are positive (conditional on correctness of the underlying interval model)."]
    else:
        lines += ["", "The reported LDL bounds do not certify positive definiteness."]
    return '\n'.join(lines) + '\n'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('reports', nargs='+', type=Path)
    p.add_argument('--out', type=Path)
    a = p.parse_args()
    content = '\n'.join(summarize(json.loads(p.read_text(encoding='utf8'))) for p in a.reports)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(content, encoding='utf8')
    print(content)


if __name__ == '__main__':
    main()
