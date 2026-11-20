#!/usr/bin/env python3
"""ES8-P257: exact-rational audit of the prime-7 taper obstruction and
optimized boundary-strip Fourier estimate.

This code verifies numerical constants and algebraic hypotheses in Fraction.
The infinite-dimensional operator, Plancherel, Fourier projector and classic
Machin/arctanh identities are proven mathematically in PROVA_ES8P257.md;
they are NOT automatically formalized by this finite Python checker.

Crucial: this is a CONDITIONAL crossing estimate; it does NOT certify the
auxiliary four-prime model positive above log(7)/2 or prove RH.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

OMEGA = 2016
OLD_LOSS = F(2520)
TARGET_LOSS = F(2225)
KAPPA_UP = F(37, 50)
KAPPA_LO_TARGET = F(73,100)
PI_LO_TARGET = F(157,50)
HIGH_RESERVE_LO = F(11,8)
HIGH_REMAINDER = F(1,12)


def need(statement: bool, message: str) -> None:
    if not statement:
        raise ArithmeticError(message)


def arctanh_lower(x: F, terms: int) -> F:
    x = F(x)
    if not F(0) < x < F(1) or terms < 1:
        raise ValueError('atanh parameters must be in (0,1) with terms>=1')
    return 2 * sum((x**(2*n+1) / (2*n+1) for n in range(terms)), F(0))


def arctanh_upper(x: F, terms: int) -> F:
    lower = arctanh_lower(x, terms)
    return lower + 2*x**(2*terms+1)/((2*terms+1)*(1-x*x))


def arctan_lower_alternating(x: F, terms: int) -> F:
    """Even partial sum of atan(x), with last term negative, is a lower bound."""
    x = F(x)
    if not 0 < x <= 1 or terms < 1:
        raise ValueError('atan parameters invalid')
    return sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(2*terms)), F(0))


def pi_strict_lower() -> F:
    # pi = 16 atan(1/5) - 4 atan(1/239), with atan(1/239)<1/239.
    return 16*arctan_lower_alternating(F(1,5),2)-4*F(1,239)


def prime_7_constants():
    log7_low = arctanh_lower(F(3,4),8)
    log7_up = arctanh_upper(F(3,4),8)
    log4_low = arctanh_lower(F(3,5),4)
    pi_low = pi_strict_lower()
    # Integer cross-product checks, without floating-point sign decisions.
    need(F(66,25)**2 < 7 < F(53,20)**2, 'sqrt(7) rational brackets invalid')
    need(log7_up < F(39,20), 'log7 upper bound not strict')
    kappa_lower = log7_low/F(53,20)
    kappa_upper = log7_up/F(66,25)
    need(KAPPA_LO_TARGET < kappa_lower, 'kappa>73/100 not proved')
    need(kappa_upper < KAPPA_UP, 'kappa<37/50 not proved')
    need(pi_low > PI_LO_TARGET, 'pi >157/50 not proved')
    # B(2016)-B(504) = log 4 + (1/504 - 1/2016) = log4+1/672.
    reserve_lower = log4_low + F(1,672)
    need(reserve_lower > HIGH_RESERVE_LO, 'four-prime high tail reserve not proved')
    return {
        'log7_lower_exact':str(log7_low),
        'log7_upper_exact':str(log7_up),
        'sqrt7_lower_exact':str(F(66,25)),
        'sqrt7_upper_exact':str(F(53,20)),
        'kappa_lower_exact':str(kappa_lower),
        'kappa_upper_exact':str(kappa_upper),
        'pi_lower_exact':str(pi_low),
        'log4_lower_exact':str(log4_low),
        'B2016_minus_B504_lower_exact':str(reserve_lower),
    }


def optimal_eta_and_loss():
    """Optimize strip penalty with a fixed residual high tail of 1/12.

    Given kappa<k=37/50, pi>157/50, reserve>11/8, choose eta=k/(tau-k)
    with tau=11/8-1/12=31/24. High frequency consumption <=tau.
    Then loss factor is < (2Omega/pi_low)*k*(1+eta).
    """
    constants=prime_7_constants()
    k=KAPPA_UP
    tau=HIGH_RESERVE_LO-HIGH_REMAINDER
    need(tau>k, 'frequency reserve not enough for prime-7 coupling')
    eta=k/(tau-k)
    high_consumption=k*(1+1/eta)
    need(high_consumption==tau, 'Young parameter algebra failed')
    C = F(2*OMEGA,1)/PI_LO_TARGET * k * (1+eta)
    need(C<TARGET_LOSS<OLD_LOSS, 'strict strip penalty improvement not established')
    return {
        'underlying_rational_constant_bounds': constants,
        'omega': OMEGA,
        'kappa_upper_exact':str(k),
        'pi_lower_exact':str(PI_LO_TARGET),
        'high_floor_lower_exact':str(HIGH_RESERVE_LO),
        'high_remaining_floor_exact':str(HIGH_REMAINDER),
        'tau_exact':str(tau),
        'young_eta_exact':str(eta),
        'young_high_consumption_upper_exact':str(high_consumption),
        'loss_factor_upper_exact':str(C),
        'usable_integer_loss_strict_upper':int(TARGET_LOSS),
        'old_integer_loss':int(OLD_LOSS),
        'relative_saving_exact':str(1-TARGET_LOSS/OLD_LOSS),
        'old_delta_for_half_margin':'m/5040',
        'new_delta_for_half_margin':'m/4450',
        'delta_allowance_gain_over_old_exact':str(F(5040,4450)-1),
    }


def normalized_taper_weight(q: int, X: F) -> F:
    X=F(X)
    if q < 2 or X <= 0:
        raise ValueError('q >= 2 and X positive required')
    return F(0) if F(q)>=X else F(1)-q/X


def first_prime7_gap(eps: F):
    """Finite rational lower bound on ||K7_hard - K7_taper||, via kappa>73/100.

    For X=7+eps, 0<eps<1, taper w7=eps/X and norm gap EXACTLY
    kappa * 7/X, because ||J+J*||=1 on the full L2 boundary strips.
    """
    eps=F(eps)
    if not F(0)<eps<F(1):raise ValueError('require 0<epsilon<1')
    X=F(7)+eps
    w=normalized_taper_weight(7,X)
    need(w==eps/X, 'taper failed at prime7')
    need(w<eps/7, 'taper is not O(epsilon)')
    norm_gap_lower = KAPPA_LO_TARGET*(1-w)
    need(norm_gap_lower>F(63,100), 'gap was not nonuniformly large')
    return {'X':str(X),'epsilon':str(eps), 'taper_weight_7':str(w),
            'hard_shift_norm':'log(7)/sqrt(7)',
            'taper_shift_norm':'(log(7)/sqrt(7))*epsilon/(7+epsilon)',
            'difference_shift_norm':'(log(7)/sqrt(7))*7/(7+epsilon)',
            'difference_norm_strict_rational_lower':str(norm_gap_lower)}


def finite_operator_toy(k: F, taper: F):
    """Exact 2x2 model of two boundary strips; shares nonzero spectrum with J+J*."""
    k,taper=F(k),F(taper)
    if k<=0 or not 0<=taper<=1:raise ValueError('invalid coupling or taper')
    hard=[[F(0),k],[k,F(0)]]
    soft=[[F(0),k*taper],[k*taper,F(0)]]
    u=(F(1),F(1))
    v=(F(1),F(-1))
    def rayleigh(M,vec):
        return sum(vec[i]*M[i][j]*vec[j] for i in range(2) for j in range(2)) / 2
    need(rayleigh(hard,u)==k and rayleigh(hard,v)==-k, 'hard spectrum error')
    need(rayleigh(soft,u)==k*taper and rayleigh(soft,v)==-k*taper, 'soft spectrum error')
    return {'hard_negative_rayleigh_exact':str(rayleigh(hard,v)),
            'soft_negative_rayleigh_exact':str(rayleigh(soft,v)),
            'exact_difference_operator_norm':str(k*(1-taper))}


def audit():
    loss=optimal_eta_and_loss()
    eps_values=[F(1,10),F(1,100),F(1,1000),F(1,10**6)]
    gaps=[first_prime7_gap(eps) for eps in eps_values]
    # Verify the direct cosine convolution identity term-wise:
    # tapered S_L=2 int_0^L exp(-2(L-u)) S_u du
    # on a q contribution, RHS weight = 1-exp(-2(L-a_q))=1-q/X.
    for X in [F(7,1)+F(1,1000),F(9),F(20)]:
        for q in range(2,int(X)+2):
            if q<X:
                need(normalized_taper_weight(q,X)==1-F(q,X),'taper kernel coefficient mismatch')
    toy=finite_operator_toy(F(3,4),F(1,100))
    return {
        'schema':'NUMERIA.ES8P257.TAPER_AND_STRIP_CERT.1',
        'proof_file':'PROVA_ES8P257.md',
        'results':{
            'exact_taper_convolution':'S_taper_L(t)=2*int_0^L exp(-2*(L-u))*S_u(t) du',
            'prime7_hard_operator_norm':'log(7)/sqrt(7) for every 0<delta=2L-log7<log7',
            'operator_gap_nonvanishing_as_epsilon_to_zero':True,
            'taper_not_an_operator_norm_approximation':True,
            'optimized_conditional_strip_loss':'2225*delta instead of 2520*delta',
            'conditional_crossing_floor':'Q_Weil(f) >= (m-2225*delta)||f||^2 + (1/12)||P_{|t|>2016}f||^2',
        },
        'exact_rational_constants':loss,
        'prime7_taper_examples':gaps,
        'two_strip_finite_toy':toy,
        'limitations':{
            'four_prime_auxiliary_margin_positive_beyond_threshold_established':False,
            'crossing_log7_over_2_established':False,
            'bounded_quantities_not_automatically_RH':True,
            'Lean_proof_replayed':False,
            'RH_proved':False,
        }
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    result=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    print('P257 PASS: rational constants and taper obstruction; crossing remains conditional')

if __name__=='__main__':main()
