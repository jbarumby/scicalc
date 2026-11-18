#!/usr/bin/env python3
"""Independent exact-rational audit of two directed-MPFR LDL/inverse traces.

Inputs are the *previously assembled* P237 symmetric interval matrices.  No
assumption that their entries correctly represent a Weil quadratic form is made.
The MPFR executable computes centers and inverse enclosures without Arb.
This exact Python verifier audits the independently emitted dyadic results and
sums all 142,506 original radius entries into symmetric row-norm bounds.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
from fractions import Fraction as F

N=377
SHA_ASSEMBLY='2695768a3cca11237dd63cb40c0b6c34427b390329e5ff2ec089e31111c7898a'
SHA_CENTER='272caafee875c2d5301782ed56bcb4d895df98ee0551a2e878d99577f44156d3'
SHA_SHADOW='e51a28fe0180c1a6215be23b0f268d570bed3c748210dcca0a0d9ff3f59bea26'
SHA_C={'even':'b5419ac9af938915e646fa34e6990efb78c5f64aca62519720569320c8cd8c80',
       'odd':'0232b9ca923654ff89330efa0d1a8945a37085e9bbd4ca7ce218fdc7879f195c'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dyad(m,e):
    m,e=int(m),int(e)
    return F(m<<e) if e>=0 else F(m,1<<-e)

def ball_bounds(b):
    if not (isinstance(b,list) and len(b)==2):raise ValueError('invalid ball')
    mid=dyad(*b[0]);rad=dyad(*b[1])
    if rad<0:raise ValueError('negative radius')
    return mid-rad,mid+rad

def load_mpfr_trace(path):
    rows=[];nu=None
    for line in Path(path).read_text().splitlines():
        if not line:continue
        if line.startswith('# inverse_frobenius_squared_upper '):
            if nu is not None:raise ValueError('duplicate nu')
            a=line.split();nu=dyad(a[-2],a[-1])
        elif line.startswith('#'):continue
        else:
            a=line.split()
            if len(a)!=5 or int(a[0])!=len(rows):raise ValueError('bad MPFR trace entry')
            lo,hi=dyad(a[1],a[2]),dyad(a[3],a[4])
            if not (0<lo<=hi):raise ValueError('nonpositive independent LDL pivot')
            rows.append([lo,hi])
    if len(rows)!=N or nu is None or nu<=0:raise ValueError('trace incomplete')
    return rows,nu

def ceil_fraction_decimal_scientific(q, digits=12):
    # Human-readable approximation only; no use in proof inequalities.
    if q==0:return '0'
    import decimal
    with decimal.localcontext() as ctx:
        ctx.prec=digits+7
        return f'{decimal.Decimal(q.numerator)/decimal.Decimal(q.denominator):.{digits-1}e}'

def audit(files,trace_dir):
    for filename,expected in [('P237_FINITE_LDL_2048.json.gz',SHA_ASSEMBLY),
      ('P237_CENTER_WEYL_2048.json.gz',SHA_CENTER),
      ('P237_SHADOW_FULL_AUDIT.json',SHA_SHADOW),
      *[(f'P237_FULL_C_{sec}.json.gz',digest) for sec,digest in SHA_C.items()]]:
        if sha(files/filename)!=expected:raise ValueError(f'SHA-256 mismatch {filename}')
    matrix=json.loads(gzip.decompress((files/'P237_FINITE_LDL_2048.json.gz').read_bytes()))
    center=json.loads(gzip.decompress((files/'P237_CENTER_WEYL_2048.json.gz').read_bytes()))
    if matrix['schema']!='NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1':raise ValueError('matrix schema mismatch')
    if center['input_artifact_sha256']!=SHA_ASSEMBLY:raise ValueError('center not bound')
    if matrix['independent_shadow_exact_audit_sha256']!=SHA_SHADOW:raise ValueError('shadow not bound')
    if not (matrix['full_quadrature_coverage_verified'] and matrix['all_142506_interval_sums_shadow_audited'] and matrix['global_GL60_error_added_once_per_matrix_entry']):
        raise ValueError('not full matrix')
    report={'schema':'NUMERIA.ES8P242.INDEPENDENT_MPFR_1024_LDL_WEYL.1',
        'source_full_matrix_sha256':SHA_ASSEMBLY,'source_center_sha256':SHA_CENTER,
        'input_pivot_algorithm':'independent C MPFR directed interval LDL on exact dyadic centers (Schur rank-one updates)',
        'input_inverse_algorithm':'independent C MPFR directed triangular solves for C inverse',
        'exact_verification':'Python stdlib integer Fraction; no floats for inequalities',
        'sectors':{},'total_independent_positive_pivots':0,'RH_proved':False,
        'full_Weil_Zhu_bridge_verified':False,'source_interval_matrix_assembly_independently_reproduced':False}
    for sector in ['even','odd']:
        mpfr_pivots,nu=load_mpfr_trace(trace_dir/f'P242_PIVOTS_{sector}_1024.txt')
        original_pivots=center['sectors'][sector]['center_pivots']
        if len(original_pivots)!=N:raise ValueError('original not 377 pivots')
        mismatches=[]
        for k,(mine,orig) in enumerate(zip(mpfr_pivots,original_pivots)):
            ol,oh=ball_bounds(orig)
            if max(mine[0],ol)>min(mine[1],oh):mismatches.append(k)
        if mismatches:raise ValueError(f'{sector} independent pivot bounds do not intersect original at {mismatches[:5]}')
        arr=matrix['sectors'][sector]['assembled_full_R_upper_triangle']
        if len(arr)!=N*(N+1)//2:raise ValueError('matrix triangle size mismatch')
        sums=[F(0) for _ in range(N)]
        maximum_rad=F(0)
        idx=0
        for i in range(N):
            for j in range(i,N):
                ii,jj,b=arr[idx];idx+=1
                if (i,j)!=(ii,jj):raise ValueError(f'{sector} matrix order mismatch')
                mid=dyad(*b[0]);rad=dyad(*b[1]);
                if rad<0:raise ValueError('bad radius')
                sums[i]+=rad
                if i!=j:sums[j]+=rad
                if maximum_rad<rad:maximum_rad=rad
        delta=max(sums)
        # Exact Weyl tests, without calculating a square root:
        #    delta^2 * nu < 1   => every matrix inside box SPD
        #    (delta + m)^2 * nu < 1 => every matrix >= m I
        criterion=delta*delta*nu
        m=F(1,10**18)
        criterion_margin=(delta+m)*(delta+m)*nu
        if not criterion<1:raise ArithmeticError(f'{sector} full box Weyl failed')
        if not criterion_margin<1:raise ArithmeticError(f'{sector} margin 1e-18 failed')
        if not delta<F(4,10**47):raise ArithmeticError(f'{sector} radius too large')
        if nu>F(10**36):raise ArithmeticError(f'{sector} inverse norm unexpectedly large')
        report['total_independent_positive_pivots']+=N
        report['sectors'][sector]={
            'dimension':N,
            'independent_mpfr_pivot_count':N,
            'independent_mpfr_pivot_enclosures_all_strictly_positive':True,
            'all_original_Arb_pivot_balls_intersect_independent_mpfr_balls':True,
            'original_interval_matrix_entries_rechecked':N*(N+1)//2,
            'symmetric_row_radii_summed_as_exact_fractions':True,
            'row_radii_max_exact_comparison':'delta < 4/10^47',
            'row_radii_max_approx':ceil_fraction_decimal_scientific(delta,16),
            'inverse_frobenius_squared_upper_approx':ceil_fraction_decimal_scientific(nu,16),
            'full_box_weyl_criterion_upper_approx':ceil_fraction_decimal_scientific(criterion,16),
            'criterion_delta_squared_nu_strictly_below_1':True,
            'criterion_delta_plus_1e-18_squared_nu_strictly_below_1':True,
            'proved_lower_eigenvalue_for_every_matrix_in_box':'1/10^18',
            'independent_trace_sha256':sha(trace_dir/f'P242_PIVOTS_{sector}_1024.txt')}
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--files',type=Path,default=Path('.'))
    p.add_argument('--traces',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'P242_AUDIT_CERT.json')
    a=p.parse_args()
    report=audit(a.files,a.traces)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    for sec,v in report['sectors'].items():print(sec,v['independent_mpfr_pivot_count'],'delta',v['row_radii_max_approx'],'nu',v['inverse_frobenius_squared_upper_approx'])
    print('P242 CERTIFIED complete centers, inverse bounds and boxes with margin 1e-18')

if __name__=='__main__':main()
