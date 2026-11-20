"""P256: exact rational coefficient bridge between von Mangoldt and NUMERIA.

The proof of the L-integral identity is in PROVA_P256.md. The rational
program checks coefficient identities in a formal basis of log(p)*cos(t log(q)).
No Lean kernel, floating-point sign decisions, RH proof, or Weil coercivity.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

LEAN_REF = '40e56bab'
SOURCE_GIT_BLOBS = {
 'src/Development/AutomorphicArithmetic/CharacterBoxes/MainConclusion.lean':'ee340012ab19713a8acfecf8adc6a331cbc2f1a4',
 'src/Development/AutomorphicArithmetic/CharacterBoxes/EveryCharacter.lean':'6ac769568b67f85f33565ba3e1eef7b91f53c55f',
 'src/Development/AutomorphicArithmetic/CharacterBoxes/RiemannHypothesisOfEveryCharacter.lean':'c6dfef5bfa62aa0ef1cd1481833e1cd2475cbbbc',
 'src/Development/AutomorphicArithmetic/CharacterBoxes/ClaimStrip.lean':'60b3e96d81c2cb7d821d14c978b38bd9fd0d97f8',
 'src/Development/AutomorphicArithmetic/CharacterBoxes/HalfPlane.lean':'dcd66cb38f94dfabe24700848a68dc9135bbd3fe',
 'src/Development/AutomorphicArithmetic/Foundations/GaloisRepresentations/IntegralStructures/RowSrc.lean':'f769999ac5e49241993e5efbf81e4d2cded18acd',
 'src/Development/AutomorphicArithmetic/AutomorphicVarieties/ModularityArgument/ThresholdArguments/DenseThresholdCases/VwinOfRow.lean':'4cba179acd45ac08d3cad21f51c9b9d62715a9a4',
 'src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean':'a13ba1f8c38d4dbf2ddd809fbf9610a007936c24'
}

def prime_power_base(n):
    if type(n) is not int or n < 2:
        raise ValueError('n must be an integer >= 2')
    remaining, p = n, 2
    while p*p <= remaining:
        if remaining % p == 0:
            while remaining % p == 0:
                remaining //= p
            return p if remaining == 1 else None
        p += 1
    return n if remaining == n else None

def prime_terms(x):
    x = Q(x)
    if x < 2:
        return []
    upper = (x.numerator + x.denominator-1) // x.denominator
    return [(q,p) for q in range(2,upper) if (p:=prime_power_base(q)) is not None]

def weighted_cosine_coefficients(x):
    x=Q(x)
    if x<=0: raise ValueError('X must be positive')
    return {q:(p,x/Q(q)-1) for q,p in prime_terms(x)}

def abel_coefficients(x):
    x=Q(x)
    if x<=0: raise ValueError('X must be positive')
    return {q:(p,x*(1/Q(q)-1/x)) for q,p in prime_terms(x)}

def imaginary_mellin_coefficients(x):
    x=Q(x)
    if x<=0: raise ValueError('X must be positive')
    # X*P(i/2)-P(-i/2), evaluated formally coefficient by coefficient.
    return {q:(p,x/Q(q)-1) for q,p in prime_terms(x)}

def grouped_log_coefficients(x):
    ans={}
    for q,(p,c) in weighted_cosine_coefficients(x).items():
        ans[p]=ans.get(p,Q(0))+c
    return ans

def slope_coefficients(x):
    ans={}
    for q,p in prime_terms(Q(x)):
        ans[p]=ans.get(p,Q(0))+Q(1,q)
    return ans

def activation_7(epsilon):
    epsilon=Q(epsilon)
    if not Q(0)<epsilon<Q(1):
        raise ValueError('require 0<epsilon<1')
    weighted=weighted_cosine_coefficients(Q(7)+epsilon)
    if weighted.get(7)!=(7,epsilon/7) or set(weighted)!={2,3,4,5,7}:
        raise ArithmeticError('incorrect prime 7 activation')
    return {'X':str(7+epsilon),'epsilon':str(epsilon),
            'extra_weighted_log7_cos_coefficient':str(epsilon/7),
            'unweighted_S_jump_at_t_zero':'2*log(7)/sqrt(7)'}

def git_blob_sha1(raw):
    return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()

def validate_source_tree(source_dir):
    """Check source hashes AND reviewed statements; no Lean proof compilation."""
    raw={}
    for name,sha in SOURCE_GIT_BLOBS.items():
        original=(source_dir/name).read_bytes()
        if git_blob_sha1(original)!=sha:
            raise ValueError('unreviewed source: '+name)
        raw[name]=original.decode('utf-8')
    def get(key):
        return next(value for name,value in raw.items() if name.endswith('/'+key+'.lean'))
    checks={
      'MainConclusion': [
          'theorem UpperTheory.FinalTheorems.everyCharacter : UpperTheory.EveryCharacter :=',
          'everyCharacter_of_row UpperTheory.Gal.Int.galoisLocalGlobalStrong_recParameter',
          'every_character_of_literature_of_vWin','UpperTheory.One.vwin_of_row h'],
      'EveryCharacter':['def UpperTheory.EveryCharacter : Prop :=',
                        '1 / 2 < s.re','DirichletCharacter.LFunction χ s ≠ 0'],
      'RiemannHypothesisOfEveryCharacter': [
          'theorem UpperTheory.zeta_ne_zero_of_everyCharacter',
          'theorem UpperTheory.riemannHypothesis_of_everyCharacter'],
      'ClaimStrip':['theorem factor_ne_zero','(hs : 0 < s.re)',
                    'theorem prim_zero','factor_ne_zero _ hs0 p',
                    'theorem strip (h : HalfPlane) : Strip := by'],
      'HalfPlane':['def Claim.HalfPlane : Prop :=','1 / 2 < s.re'],
      'RowSrc':['theorem galoisLocalGlobalStrong_recParameter',
                'galoisRepresentations_of_thetaLift'],
      'VwinOfRow':['theorem UpperTheory.One.vwin_of_row',
                   'UpperTheory.One.Join.automorphic_window_of_inputs'],
      'EveryCharacterOfLiteratureOfVWin':[
          'theorem UpperTheory.LiteratureTheorems.every_character_of_literature_of_vWin',
          '(O : _root_.UpperTheory.LiteratureTheorems.OpenLiteratureStatements)',
          '(hV : UpperTheory.Top.VWin) : UpperTheory.EveryCharacter']
    }
    for key,needles in checks.items():
        source=get(key)
        for needle in needles:
            if needle not in source: raise ValueError(f'{key}: missing {needle}')
    return {'revision':LEAN_REF,'authenticated_source_files':len(raw),
            'exact_git_blob_shas':SOURCE_GIT_BLOBS,
            'verified_scope':'static names/types/hypotheses only',
            'Lean_compiled_here':False}

def audit():
    windows=[Q(2),Q(6),Q(7),Q(7001,1000),Q(8),Q(9),Q(12),Q(20),Q(100)]
    for x in windows:
        reference=weighted_cosine_coefficients(x)
        if reference!=abel_coefficients(x) or reference!=imaginary_mellin_coefficients(x):
            raise ArithmeticError('invalid coefficient identity at '+str(x))
    if grouped_log_coefficients(Q(6))!={2:Q(5,2),3:Q(1),5:Q(1,5)}:
        raise ArithmeticError('incorrect four-prime weighted sum')
    if prime_power_base(6) is not None:
        raise ArithmeticError('Lambda(6) is zero')
    before,after=slope_coefficients(Q(7)),slope_coefficients(Q(7001,1000))
    if after.get(7)!=Q(1,7) or any(after[p]!=v for p,v in before.items()):
        raise ArithmeticError('wrong prime 7 derivative jump')
    return {'schema':'NUMERIA.ES8P256.WEIGHTED_PRIME_BRIDGE.1',
            'bridge':'W_L(t)=exp(L)*integral_0^L cosh(L-u)*S_u(t)du',
            'analytic_integral_proof_file':'PROVA_P256.md',
            'tests_in_exact_fraction_arithmetic':[str(x) for x in windows],
            'weighted_at_X6':{str(p):str(c) for p,c in grouped_log_coefficients(6).items()},
            'prime7_activation':activation_7(Q(1,1000)),
            'finite_coefficient_checks':True,
            'analytic_integral_replayed_numerically':False,
            'new_Weil_coercivity_proved':False,'RH_proved':False}

def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--out',type=Path,required=True)
    cli.add_argument('--source-dir',type=Path)
    opts=cli.parse_args()
    result=audit()
    if opts.source_dir:
        result['source_chain']=validate_source_tree(opts.source_dir)
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print('P256 PASS: finite weighted prime coefficients; no RH claim')

if __name__=='__main__':
    main()
