#!/usr/bin/env python3
"""P246: independent rational audit of uploaded P237/P245D JSONs.

This verifier does not execute the P245D GL60 integrand or authenticate missing
P237 chunk blobs. It validates the archived inclusions, exact rational source
balls and independently reconstructs all 754 Bernstein/Bessel row-budget values.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path

PINS = {
 'test_p237_analytic_certificates.py':'3d5e8ca9c8e1905b1a344fc2c92f69a7a1e9afd5e885436166b1221e770c727d',
 'test_p245d_additional_chunk_replay.py':'87de36698024ac435936e470666f2eb5ea7008e3b22ade3f1564c834c6de6260',
 'p245d_additional_chunk_replay.py':'bdbc99b9e267bc54d488a987ba5b5a9a386b2d3f984f04d45c55cc2ee50abe85',
 'P245D_SELECTED_975_1000.json':'719af5bf1937fa0bb584c1db436f1501a27196e5ecbeb8b9229e0ddb40140e53',
 'P245D_SELECTED_500_525.json':'597ea1038dc99f0a71074f33f1be27922f2a5e15011f5e5e447f656f604e1cd2',
 'P245D_SELECTED_25_50.json':'40ba13381884c5e11db58aa95014053a2affaf11a4b1eb76122dba02ffbe037b',
 'P237_SYMBOL_E3_RATIONAL_GATE.json':'3e9a33fb6e762f87f905cb85b3b9db7c1a895968446fd6ab1b77ba8153b79ae1',
 'P237_GL60_ANALYTIC_CERT.json':'d2daefbc570d6b43fc64206bad71dcd82b94519e52900f3b2233b977628c033d',
 'P237_ANALYTIC_GATES.json':'e27a32503d1d01196dc2173f156c4a7107a2040ebe44dcb9645133672925f5ec',
}
PAIRS=((0,0),(0,1),(0,2),(0,3),(1,1),(1,2),(2,2),(3,3),
       (0,16),(16,16),(0,376),(1,376),(376,376))
REF_PINS={
(25,'even'):'3903e507a9b7c492cdab22b8a290481d43593d69d9ca0f8fde9fb80cb6d73da8',
(25,'odd'):'ff0f5256cd92893ee2b27c2753b9c97e8a17c43e55c92c8c76c908beb899eeda',
(500,'even'):'b448bcaad18577b10a056970c0e277873f316b9b8ff35dffb6f09e4547505d71',
(500,'odd'):'0380e8986f374d00fcf97f6b4d851cffbd596db89a518ae71b7eb983bf484607',
(975,'even'):'7a7d4f3290c0ff4a14a728af37af322bc723cb73bd9c92f5f012b12d7c49ed2a',
(975,'odd'):'d26488fdf2af0b26fa3bd7c0e7b767824c202432265b13e697979143cdc344be'}
ROOT_SHA='bf3a05df56b98e7a947eed516e25e8342a473d3c9aa01f3d194a582580742218'

def need(test,msg):
    if not test:raise ArithmeticError(msg)

def dyad(m,e):
    m,e=int(m),int(e)
    return F(m << e) if e>=0 else F(m,1<<(-e))

def approx(x):
    with localcontext() as ctx:
        ctx.prec=18
        return f'{Decimal(x.numerator)/Decimal(x.denominator):.12e}'

def verify(data:Path,check_hashes=True):
    if check_hashes:
        for filename,digest in PINS.items():
            path=data/filename
            need(sha256(path.read_bytes()).hexdigest()==digest,'source hash mismatch: '+filename)
    gates=json.loads((data/'P237_ANALYTIC_GATES.json').read_text())
    symbol=json.loads((data/'P237_SYMBOL_E3_RATIONAL_GATE.json').read_text())
    gl60=json.loads((data/'P237_GL60_ANALYTIC_CERT.json').read_text())
    L=F(161,200)
    rem=F(2*1007+1,4)+F(3+F(1,3),8)
    imm=F(3-F(1,3),8)
    need(rem==F(3025,6) and imm==F(1,3),'ellipse geometry independent')
    need(F(gates['t_real_abs_upper'])==rem,'P237 analytic geom mismatch')
    need(F(gates['t_imag_abs_upper'])==imm,'P237 imag geom mismatch')
    need(gates['panels']==1008 and gates['Tsharp']==504,'P237 size mismatch')
    need(F(gates['bessel_Lt_modulus_squared_upper'])==L*L*(rem*rem+imm*imm),'Bessel magnitude mismatch')
    need(L*L*(rem*rem+imm*imm)<406*406,'Bessel magnitude not <406')
    need(F(symbol['symbol_minus_beta_bound_strict_less'])==F(274,7),'symbol bound mismatch')
    need(symbol['beta_rational_lo_gt_zero'] and symbol['beta_rational_hi_lt_two'],'beta gate unavailable')
    need(F(symbol['cosh_2over3_bound'])==F(9,7),'cosh gate mismatch')
    need(gl60['Gauss_legendre_nodes']==60 and gl60['degree_of_exactness']==119,'GL60 mismatch')
    need(F(gl60['symbol_complex_bound_strict_less'])==F(274,7),'symbol/gl60 mismatch')
    need(gl60['per_panel_error_coefficient_strict_less']=='50 / 3^120','old error factor mismatch')

    coeff_old=F(1008*50,3**120)
    coeff_new=F(1008,3**120)*F(274,7)
    need(coeff_new/coeff_old==F(137,175),'ratio mismatch')
    amplitudes={0:[],1:[]}; odd_df=1; xpow=1
    for n in range(754):
        nr=2*L*(2*n+1)
        ceilroot=isqrt(nr.numerator//nr.denominator)
        if F(ceilroot*ceilroot)<nr:ceilroot+=1
        amplitudes[n&1].append(ceilroot*min(F(2),F(2*xpow,odd_df)))
        xpow*=406
        odd_df*=2*n+3
    rows={}
    for p,sector in ((0,'even'),(1,'odd')):
        amp=amplitudes[p]
        original=gl60['sectors'][sector]
        need(len(amp)==377,'mode count')
        need(str(max(amp))==original['max_basis_bound_exact'],'maximum amplitude mismatch')
        need(str(sum(amp))==original['sum_basis_bounds_exact'],'sum amplitude mismatch')
        row_old=coeff_old*max(amp)*sum(amp)
        row_new=coeff_new*max(amp)*sum(amp)
        need(row_old==F(original['row_sum_GL60_error_exact']),'source GL60 row budget mismatch')
        need(row_new<F(4,10**47),'strengthened row budget over threshold')
        rows[sector]={'modes':377,'original_row_budget_approx':approx(row_old),
          'strengthened_row_budget_approx':approx(row_new),
          'strengthened_below_4e-47':True,
          'ratio_new_old_exact':'137/175'}

    incl=0;coarse=0;per_block=[]
    for start in (25,500,975):
        cert=json.loads((data/f'P245D_SELECTED_{start}_{start+25}.json').read_text())
        need(cert['parameters']['panel_range']==[start,start+25],'wrong panel range')
        need(cert['root_brackets_sha256']==ROOT_SHA,'root provenance mismatch')
        need(cert['parameters']['GL60_nodes_per_panel']==60,'bad node count')
        need(cert['parameters']['fixed_point_bits']==256,'bad precision metadata')
        need(cert['panel_pair_evaluations']==39000,'invalid work count')
        need(cert['verified_chunk_sum_enclosures']==26,'invalid claim count')
        need(len(cert['comparisons'])==2,'sector count incorrect')
        for comparison in cert['comparisons']:
            sector=comparison['sector'];need(sector in ('even','odd'),'invalid sector')
            need(comparison['source']['sha256']==REF_PINS[start,sector],'chunk SHA claim mismatch')
            need(cert['reference_chunks'][sector]['sha256']==REF_PINS[start,sector],'chunk metadata mismatch')
            recs=comparison['selected_entry_comparisons']
            need(len(recs)==13,'wrong comparison length')
            need([tuple(r['pair']) for r in recs]==list(PAIRS),'wrong pairs/order')
            for r in recs:
                (m,e),(rad,k)=r['source_Arb_ball_dyadic']
                midpoint=dyad(m,e);rad=dyad(rad,k)
                need(rad>=0,'negative source radius')
                source_lower,source_upper=midpoint-rad,midpoint+rad
                need(F(r['source_interval_lower_fraction'])==source_lower and F(r['source_interval_upper_fraction'])==source_upper,
                     'archived source ball/fraction mismatch')
                interval=r['replay_interval']
                need(interval['scale_power_of_two']==256,'bad scale')
                lower,upper=F(int(interval['lower']),2**256),F(int(interval['upper']),2**256)
                need(lower<=source_lower<=source_upper<=upper,'stored interval does not enclose source')
                if 376 in r['pair']:coarse+=1
                incl+=1
        per_block.append({'panel_range':[start,start+25], 'archived_inclusions':26,
                          'original_chunk_sha256_pinned_but_bytes_unavailable':True})
    return {'schema':'NUMERIA.P246.INDEPENDENT_UPLOADED_CERT_AUDIT.1',
       'input_artifacts_hash_pinned':check_hashes,'input_artifacts_checked':len(PINS),
       'archived_p245d_inclusions_verified_rationally':incl,
       'archived_high_degree_coarse_cases':coarse,
       'p245d_gl60_contributions_reexecuted_from_integrand':False,
       'p237_analytic_amplitudes_independently_reconstructed':754,
       'improved_symbol_bound':'274/7','old_to_new_error_ratio':'137/175',
       'GL60_row_errors':rows,'blocks':per_block,
       'p237_chunk_reference_blobs_available':False,
       'full_142506_entries_reconstructed':False,
       'Weil_Zhu_identity_proved':False,'RH_proved':False}


def main():
    pa=argparse.ArgumentParser();pa.add_argument('--data',type=Path,required=True)
    pa.add_argument('--out',type=Path,required=True)
    a=pa.parse_args()
    report=verify(a.data)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('PASS P246 exact stored inclusions:',report['archived_p245d_inclusions_verified_rationally'])
    print('PASS P246 independently generated basis bounds:',report['p237_analytic_amplitudes_independently_reconstructed'])
    print('New GL60 ratio:',report['old_to_new_error_ratio'])

if __name__=='__main__':main()
