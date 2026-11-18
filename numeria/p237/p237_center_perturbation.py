#!/usr/bin/env python3
"""Certify P237 full interval matrices via midpoint SPD + Weyl perturbation.

The midpoint matrix is certified by interval LDL.  A rigorous upper bound on
its inverse Frobenius norm and a row-sum bound on the interval radii then prove
that every matrix in the assembled Arb box is positive definite.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
from flint import arb,arb_mat,ctx
from equilibrated_ldl import N,UPPER,dyadic,ball,load_matrix

def interval_ldl(A):
    L=[[arb(0) for _ in range(N)] for _ in range(N)]
    D=[];piv=[]
    for k in range(N):
        p=A[k][k]
        for j in range(k):p-=L[k][j]*L[k][j]*D[j]
        if not p.is_finite() or p.lower()<=0:
            return False,k,piv,ball(p) if p.is_finite() else None
        D.append(p);piv.append(ball(p))
        for i in range(k+1,N):
            t=A[i][k]
            for j in range(k):t-=L[i][j]*L[k][j]*D[j]
            L[i][k]=t/p
        if (k+1)%25==0:print('CENTER LDL positive pivots',k+1,'/',N,flush=True)
    return True,N,piv,None

def certify(R):
    C=[[R[i][j].mid() for j in range(N)] for i in range(N)]
    ok,count,pivots,failed=interval_ldl(C)
    if not ok:return {'center_spd':False,'center_positive_pivots':count,'center_failed_pivot':failed}
    # Arb inverse enclosure is rigorous; sum of squares of entrywise absolute
    # upper bounds encloses ||C^{-1}||_F^2.
    inv=arb_mat(C).inv()
    sum_sq=arb(0)
    for i in range(N):
        for j in range(N):
            x=inv[i,j]
            if not x.is_finite():raise ArithmeticError(f'nonfinite inverse at {(i,j)}')
            u=x.abs_upper()
            sum_sq+=u*u
    # Symmetry makes ||E||_2 <= ||E||_infinity.  The row sum of Arb radii
    # encloses ||E||_infinity for every matrix in the interval box.
    delta=arb(0)
    for i in range(N):
        row=arb(0)
        for j in range(N):row+=R[i][j].rad().upper()
        if row.upper()>delta.upper():delta=row
    delta=delta.upper()
    criterion=delta*delta*sum_sq
    transfer=criterion<1
    return {'center_spd':True,'center_positive_pivots':N,'center_pivots':pivots,
      'inverse_frobenius_squared_upper':repr(sum_sq),
      'interval_box_infinity_norm_perturbation_upper':repr(delta),
      'squared_weyl_ratio_upper':repr(criterion),
      'perturbation_strictly_below_center_lambda_lower':bool(transfer),
      'all_matrices_in_interval_box_spd':bool(transfer)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--bits',type=int,default=2048);a=p.parse_args()
    if a.bits<1152:raise ValueError('precision below prior assembled artifact precision')
    ctx.prec=a.bits
    raw=a.input.read_bytes();src=json.loads(gzip.decompress(raw))
    if src.get('schema')!='NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1' or src.get('full_quadrature_coverage_verified') is not True or src.get('all_142506_interval_sums_shadow_audited') is not True or src.get('global_GL60_error_added_once_per_matrix_entry') is not True:raise ValueError('input artifact not full/shadow-gated')
    sectors={}
    for sec in ('even','odd'):
        old=src['sectors'][sec];R=load_matrix(old['assembled_full_R_upper_triangle']);res=certify(R)
        res['source_full_C_sha256']=old['source_full_C_sha256'];sectors[sec]=res
        print('CENTER+WEYL RESULT',sec,res.get('center_spd'),res.get('all_matrices_in_interval_box_spd'),flush=True)
    result={'schema':'NUMERIA.ES8P237.CENTER_PLUS_WEYL.1','precision_bits':a.bits,
      'input_artifact_sha256':hashlib.sha256(raw).hexdigest(),
      'input_shadow_audit_sha256':src['independent_shadow_exact_audit_sha256'],
      'same_full_matrices_GL60_beta_and_pole_enclosures':True,
      'symmetry_basis':'upper triangle mirrored exactly; bilinear form and diagonal terms are symmetric',
      'method':'Arb interval LDL of exact dyadic midpoint, Arb inverse-Frobenius bound, and symmetric row-sum Weyl perturbation bound',
      'both_377_mode_interval_matrices_certified_spd':all(sectors[s].get('all_matrices_in_interval_box_spd',False) for s in ('even','odd')),
      'infinite_Weil_positivity_proved':False,'RH_proved':False,'sectors':sectors}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),compresslevel=9,mtime=0))
    print('CERTIFICATE',result['both_377_mode_interval_matrices_certified_spd'],'sha256',hashlib.sha256(a.out.read_bytes()).hexdigest(),flush=True)
if __name__=='__main__':main()
