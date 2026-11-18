#!/usr/bin/env python3
"""P238: exact-input Arb audit of the finite margin and Schur bridge threshold.

Certifies only the finite-block margin derived from P237 and emits an explicit
sufficient inequality for a hypothetical infinite tail block. It does NOT
certify tail coercivity, coupling bounds, or the Weil/Zhu identification.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'p237'))
from flint import arb,arb_mat,ctx
from equilibrated_ldl import N,ball,load_matrix
from p237_center_perturbation import interval_ldl

COMMON_MARGIN=arb(1)/(arb(10)**18)

def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()

def certify_sector(data:Path,sector:str,bits:int):
    input_path=data/'P237_FINITE_LDL_2048.json.gz'
    center_path=data/'P237_CENTER_WEYL_2048.json.gz'
    raw=input_path.read_bytes();src=json.loads(gzip.decompress(raw))
    existing=json.loads(gzip.decompress(center_path.read_bytes()))
    if src.get('schema')!='NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1' or src.get('full_quadrature_coverage_verified') is not True or src.get('all_142506_interval_sums_shadow_audited') is not True or src.get('global_GL60_error_added_once_per_matrix_entry') is not True:
        raise ValueError('P237 full assembly flags missing')
    if existing.get('input_artifact_sha256')!=sha(input_path) or existing.get('both_377_mode_interval_matrices_certified_spd') is not True:
        raise ValueError('P237 center-Weyl certificate not bound to input')
    audit=data/'P237_SHADOW_FULL_AUDIT.json'
    if sha(audit)!=src.get('independent_shadow_exact_audit_sha256'):
        raise ValueError('shadow audit hash mismatch')
    full_c=data/f'P237_FULL_C_{sector}.json.gz'
    if sha(full_c)!=src['sectors'][sector]['source_full_C_sha256']:
        raise ValueError(f'{sector} full-C hash mismatch')
    R=load_matrix(src['sectors'][sector]['assembled_full_R_upper_triangle'])
    C=[[R[i][j].mid() for j in range(N)] for i in range(N)]
    ok,count,pivots,failed=interval_ldl(C)
    if not ok or count!=N:raise ArithmeticError(f'{sector} center LDL not complete: {count}, {failed}')
    inv=arb_mat(C).inv();nu=arb(0)
    for i in range(N):
        for j in range(N):
            x=inv[i,j]
            if not x.is_finite():raise ArithmeticError('nonfinite inverse enclosure')
            u=x.abs_upper();nu+=u*u
    delta=arb(0)
    for i in range(N):
        row=arb(0)
        for j in range(N):row+=R[i][j].rad().upper()
        if row.upper()>delta.upper():delta=row
    delta=delta.upper();nu=nu.upper()
    margin=1/nu.sqrt()-delta
    if not (margin.lower()>COMMON_MARGIN):
        raise ArithmeticError(f'{sector} margin does not exceed 1e-18: {margin!r}')
    ratio=delta*delta*nu
    if not (ratio.upper()<1):raise ArithmeticError('Weyl transfer criterion failed')
    return {'dimension':N,'center_ldl_positive_pivots':count,
      'center_spd':True,'inverse_frobenius_squared_upper_ball':ball(nu),
      'matrix_box_infinity_perturbation_upper_ball':ball(delta),
      'transferred_lambda_lower_ball':ball(margin),
      'common_margin_lower_exact':'1/1000000000000000000',
      'common_margin_strictly_exceeded':True,
      'squared_weyl_ratio_upper_ball':ball(ratio),
      'full_interval_box_spd':True,
      'source_full_C_sha256':sha(full_c)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--bits',type=int,default=2048);a=p.parse_args()
    if a.bits<1152:raise ValueError('precision below source enclosure precision')
    ctx.prec=a.bits
    data=a.data
    results={s:certify_sector(data,s,a.bits) for s in ('even','odd')}
    input_path=data/'P237_FINITE_LDL_2048.json.gz';cert_path=data/'P237_CENTER_WEYL_2048.json.gz';audit=data/'P237_SHADOW_FULL_AUDIT.json'
    result={'schema':'NUMERIA.ES8P238.P237_SCHUR_MARGIN_AUDIT.1','precision_bits':a.bits,
      'P237_matrix_artifact_sha256':sha(input_path),'P237_center_Weyl_certificate_sha256':sha(cert_path),
      'P237_shadow_audit_sha256':sha(audit),'sectors':results,
      'common_finite_block_margin_lower_exact':'1/1000000000000000000',
      'conditional_global_block_sufficient_condition':'For an identified symmetric block operator [[A,B],[B*,D]], if A >= (m0-epsilon_bridge) I, D >= eta I with eta>0, ||B|| <= beta, and epsilon_bridge + beta^2/eta < m0, then the full block is strictly positive.',
      'conditional_normalized_coupling_target_when_epsilon_zero':'beta/sqrt(eta) < 1/1000000000',
      'finite_block_margin_rederived':True,'tail_coercivity_proved':False,
      'coupling_bound_proved':False,'finite_matrix_identified_with_global_Weil_compression':False,
      'infinite_Weil_positivity_proved':False,'RH_proved':False}
    if not all(results[s]['common_margin_strictly_exceeded'] and results[s]['full_interval_box_spd'] for s in ('even','odd')):raise ArithmeticError('finite margin audit incomplete')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n')
    print('P238 finite lower margin 1e-18 verified in both sectors; global bridge remains conditional.')

if __name__=='__main__':main()
