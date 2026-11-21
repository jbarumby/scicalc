"""NUMERIA P259 Gate B: independently check fraction arithmetic in a certified
negative witness and transport it below the prime-7 threshold.

The analytic assumptions (P258/P239 Fourier multiplier and pole estimates)
are listed explicitly, not asserted as numerically re-proved here.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import file_digest
import json
from pathlib import Path

EXPECTED_GATE_B_SHA256 = "09b5a27700fa90f9c268882dbff4459fef5e6d72bd0b528f8e6e9a52d628b69a"
L_TEST = Q('0.972955074527656652552676371721590')
L_BELOW = Q('0.9729550745276566525526763717215898648185')
OLD_DILATION_LIPSCHITZ = 21250
ADDED_STEP_LIPSCHITZ = 426
NEW_DILATION_LIPSCHITZ = OLD_DILATION_LIPSCHITZ + ADDED_STEP_LIPSCHITZ


def need(value, message):
    if not value:
        raise ArithmeticError(message)


def log7_bounds(terms=180):
    x=Q(3,4)
    lower=2*sum((x**(2*n+1)/Q(2*n+1) for n in range(terms)),Q(0))
    upper=lower+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
    return lower,upper


def verify_derivative_budget():
    # For g supported [-1,1] with ||g||_2=1, Fourier dilation gives:
    # ||dF_L/dL||_{|t|<=T} <= 1/(2L)+T <=1+T for L>=0.95.
    # d/dL ||P_<=T F_L||^2 <= 2*(T+1), so 1/4 highpass
    # contributes at most (T+1)/2 <= 426 at T=850.
    p258_part=2*(504+1)*21+36
    banded_part=Q(1,4)*2*(850+1)
    need(p258_part <= OLD_DILATION_LIPSCHITZ, "P258 book-keeping exceeds claimed budget")
    need(banded_part <= ADDED_STEP_LIPSCHITZ, "new step budget exceeded")
    return {'P258_existing_band_and_poles_upper_exact':str(p258_part),
            'P259_new_850_highpass_derivative_upper_exact':str(banded_part),
            'used_P258_lipschitz':OLD_DILATION_LIPSCHITZ,
            'used_new_highpass_lipschitz':ADDED_STEP_LIPSCHITZ,
            'total_lipschitz':NEW_DILATION_LIPSCHITZ,
            'analytic_premises':[
                'P239/P258 |Psi4(t)-1/4|<21 on |t|<=504',
                'P258 uniform derivative budget of pole part <=36',
                'unitary Fourier dilation, Plancherel, support [-1,1]',
                'P259 highpass step = (1/4)*P_{|t|>850}'
            ]}


def replay_result(result):
    records=result['checked_prefixes']
    need(result['status']=='CERTIFIED_NEGATIVE_DIRECTION',"Gate B not certified negative")
    need(result['source_P258_sha256']==
         '1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c',
         'input P258 SHA mismatch')
    need(Q(result['proven_floor']['floors'][0])==Q(1,4) and
         Q(result['proven_floor']['floors'][1])==Q(1,2), 'incorrect step floors')
    need([a['dimension'] for a in records]==[29,36,48,64],"unexpected scan")
    need(all(a['certified_positive'] for a in records[:3]),'prior prefix not certified')
    bad=records[-1]
    need(bad['first_inconclusive_pivot_index']==56 and bad['positive_pivot_count']==56,
         'negative witness pivot changed')
    witness=bad['negative_witness']
    need(witness is not None and witness['coordinate_count']==57,
         'expected exact 57-coordinate witness')
    r=witness['exact_rational_replay']
    need(r['witness_strict_sign']=='negative' and r['dimension']==57,
         'fraction replay did not certify negative sign')
    vec=[Q(value) for i,value in r['nonzero_witness_exact']]
    need(len(vec)==57 and all(int(pair[0])==i for i,pair in enumerate(r['nonzero_witness_exact'])),
         'witness entries are missing or misordered')
    norm2=sum((v*v for v in vec),Q(0))
    need(norm2==Q(r['norm_squared_exact']) and norm2>0,'norm mismatch')
    err=sum((Q(v) for v in r['operator_norm_errors_exact'].values()),Q(0))
    need(err==Q(r['total_operator_norm_error_exact']) and err>=0,'GL60 error mismatch')
    components=r['matrix_components']
    need(set(components)=={'P259_banded'},'unexpected replay components')
    lo,hi=map(Q,components['P259_banded']['quadratic_interval_exact'])
    ul,uh=map(Q,r['unnormalized_interval_exact'])
    need(ul==lo-err*norm2 and uh==hi+err*norm2,'operator error not applied correctly')
    rl,rh=map(Q,r['rayleigh_interval_exact'])
    need(rl==ul/norm2 and rh==uh/norm2 and rl<=rh<0,'normalized Rayleigh arithmetic invalid')
    return rl,rh,{'dimension':57,'norm_squared_exact':str(norm2),'error_norm_upper_exact':str(err),
                  'sign_certified_with_fraction_arithmetic':True}


def audit(path):
    with Path(path).open('rb') as source:
        digest=file_digest(source,'sha256').hexdigest()
    need(digest==EXPECTED_GATE_B_SHA256,'Gate B certificate SHA256 mismatch')
    doc=json.loads(Path(path).read_text())
    oldlo,oldhi,arith=replay_result(doc)
    low,high=log7_bounds()
    need(L_BELOW<low/2<high/2<L_TEST, 'below-threshold support inequality not established')
    budget=verify_derivative_budget()
    loss=Q(NEW_DILATION_LIPSCHITZ)*(L_TEST-L_BELOW)
    transferred_upper=oldhi+loss
    need(transferred_upper<0,'negative witness cannot be certified below threshold')
    return {'schema':'NUMERIA.ES8P259.GATE_B.TRANSPORT.1',
            'source_Gate_B_sha256':digest,
            'exact_original_rayleigh_lower':str(oldlo),
            'exact_original_rayleigh_upper':str(oldhi),
            'original_upper_display':f'{float(oldhi):.17e}',
            'L_test_exact':str(L_TEST),'L_below_exact':str(L_BELOW),
            'verified_L_below_lt_log7_half_lt_L_test':True,
            'dilation_derivative_budget':budget,
            'loss_exact':str(loss),'loss_display':f'{float(loss):.17e}',
            'transported_upper_exact':str(transferred_upper),
            'transported_upper_display':f'{float(transferred_upper):.17e}',
            'rational_internal_replay':arith,
            'same_fixed_piecewise_four_prime_auxiliary_not_positive_for_any_L_ge_log7_half':True,
            'scope':'conditional on the P239/P258 analytic input bounds listed above',
            'full_Weil_form_positivity_or_RH_proved':False}


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--gate-b',type=Path,required=True)
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args()
    j=audit(a.gate_b)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(j,ensure_ascii=False,indent=2)+'\n')
    print('P259 Gate B verified negative at below log7/2',j['transported_upper_display'])

if __name__=='__main__':main()
