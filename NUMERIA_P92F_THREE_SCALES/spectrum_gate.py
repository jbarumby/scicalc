"""P92-F scale-spectrum Schur gate with exact rational-square ratios.

Universal theorem for arbitrary real nearest gaps is in SCALE_SPECTRUM_GATE.md.
Exact Fraction here handles rational-square spectra; zero float.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

ZETA2_UP=F(5,3)
PI_LOWER=F(25,8)
SQRT6_5_LOWER=F(109,100)
TARGET_LOWER=PI_LOWER**2/F(3)*SQRT6_5_LOWER
SCHEMA='NUMERIA.P92F.SCALE_SPECTRUM.RATIONAL_SCHUR.1'

def spectrum_rows(positive_square_roots):
    ss=tuple(sorted(set(map(F,positive_square_roots))))
    if not ss or any(s<=0 for s in ss):
        raise ValueError('positive roots required')
    rows=[]
    for a in ss:
        small=sum(((1-b/a)**2*(1+(b/a)**2) for b in ss if b<a),F(0))
        big=sum((ZETA2_UP*(1-a/b)**2 for b in ss if b>a),F(0))
        rows.append({'root_scale':str(a),'gap_scale':str(a*a),
                     'lower_scale_contribution':str(small),
                     'higher_scale_contribution_zeta_upper':str(big),
                     'rational_upper':str(small+big)})
    return rows

def certify():
    choices=(
        ('single',(1,)),
        ('three_squares',(1,2,3)),
        ('four_squares',(1,2,3,4)),
        ('five_squares',(1,2,3,4,5)),
        ('five_scaled_squares',(F(3,2),F(3),F(9,2),F(6),F(15,2))),
        ('five_close',(10,11,12,13,14)),
        ('six_close',(10,11,12,13,14,15)),
        ('four_extreme',(1,10,100,1000)),
    )
    out=[]
    for name,root_scales in choices:
        rows=spectrum_rows(root_scales)
        ceiling=max(F(row['rational_upper']) for row in rows)
        out.append({'name':name,'scale_count':len(rows),
                    'max_rational_schur_ceiling':str(ceiling),
                    'P88_target_proved_by_rational_lower':ceiling<TARGET_LOWER,
                    'rows':rows})
    five=next(row for row in out if row['name']=='five_squares')
    five_scaled=next(row for row in out if row['name']=='five_scaled_squares')
    if F(five['max_rational_schur_ceiling'])!=F(6829,2160):
        raise ArithmeticError('unexpected five-square exact spectrum bound')
    if five['max_rational_schur_ceiling']!=five_scaled['max_rational_schur_ceiling']:
        raise ArithmeticError('scaling invariance lost')
    if not (F(6829,2160)<TARGET_LOWER):
        raise ArithmeticError('five-scale residual target not certified')
    extreme=next(row for row in out if row['name']=='four_extreme')
    if extreme['P88_target_proved_by_rational_lower']:
        raise ArithmeticError('extreme ratio unsupported example certified')
    return {'schema':SCHEMA,
            'general_spectrum_schur_bound_proved_analytically':True,
            'five_scales_c_times_1_4_9_16_25_unlimited_points_proved':True,
            'all_arbitrary_scale_spectra_global_Delta_proved':False,
            'full_Preissmann_Lemma6_proved':False,
            'Riemann_Hypothesis_proved':False,
            'pi_lower':str(PI_LOWER),'sqrt6_over5_lower':str(SQRT6_5_LOWER),
            'P88_target_lower':str(TARGET_LOWER),
            'five_scale_rational_norm_upper':'6829/2160',
            'cases':out}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    d=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('P92F SPECTRUM PASS: five-scale rational Schur ceiling 6829/2160')
    print('Open: arbitrary spectra without a residual certificate; full L6 and RH')

if __name__=='__main__':
    main()
