#!/usr/bin/env python3
"""Conservative composite Gauss-Legendre/Bernstein bound for NUMERIA P227.

The analytic bound |f_nm(z)| < 100000 on every rho=3 panel ellipse is proved
in README_ES8P229.md. All final error bounds below are exact Fractions.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext

RHO = 3
GL_ORDER = 60
M_MAJORANT = F(100000)


def total_error_bound(T: int) -> F:
    """Composite GL bound over 2*T panels of physical width 1/2.

    Trefethen's analytic-ellipse bound is
      |I-Q_n| <= 64 M / (15 (rho^2-1) rho^(2n)).
    The map from a panel to [-1,1] contributes 1/4. With rho=3,n=60,
    replacing 3^120 by 3^120-1 only enlarges the bound.
    """
    if T not in (150, 200):
        raise ValueError('only the audited T=150 or T=200 window is supported')
    panels = 2*T
    panel_scale = F(1,4)
    per_reference_interval = F(64)*M_MAJORANT / (F(15)*(RHO*RHO-1)*(RHO**(2*GL_ORDER)-1))
    return panels*panel_scale*per_reference_interval


def decimal_bound(x: F, digits=12) -> str:
    with localcontext() as ctx:
        ctx.prec=digits+8
        y=Decimal(x.numerator)/Decimal(x.denominator)
        return f'{y:.{digits}E}'


def result(T: int, p193_E: F | None = None):
    E=total_error_bound(T)
    out={'T':T,'panel_count':2*T,'panel_width':'1/2','gauss_order':GL_ORDER,
         'rho':RHO,'uniform_integrand_majorant_exact':str(M_MAJORANT),
         'E_GL_exact':str(E),'E_GL_decimal':decimal_bound(E),
         'theorem':'Trefethen 2008, Theorem 4.5, ellipse-analytic Gauss bound',
         'analytic_bound':'|f_nm(z)| < 100000 on every mapped rho=3 ellipse; proof in README_ES8P229.md'}
    if p193_E is not None:
        out['P193_E_GL_exact']=str(p193_E)
        out['P193_E_GL_decimal']=decimal_bound(p193_E)
        out['new_bound_below_P193']=E<p193_E
        out['factor_tighter_than_P193_decimal']=decimal_bound(p193_E/E,8)
    return out

if __name__=='__main__':
    import json
    from pathlib import Path
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--p193',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    sectors=[]
    for name,T in [('even',200),('odd',150)]:
        E193=None
        if a.p193:
            import json as _json
            c=_json.loads((a.p193/f'P193_CERT_{name}.json').read_text())
            E193=F(c['E_GL_total']['exact'])
        sectors.append(result(T,E193))
    a.out.write_text(json.dumps({'stage':'NUMERIA ES8-P229 independent GL60 error bounds','sectors':sectors},indent=2)+'\n')
    print(json.dumps({'stage':'NUMERIA ES8-P229 independent GL60 error bounds','sectors':sectors},indent=2))
