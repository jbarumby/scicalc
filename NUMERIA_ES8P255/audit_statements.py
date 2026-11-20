#!/usr/bin/env python3
"""P255: source-pinned static check of the five Alpöge Lean headlines.

SECURITY/LOGIC: This does NOT invoke the Lean kernel. Hashes pin reviewed
source bytes; theorem declarations and definitions alone are NOT proofs.
The MIT challenge includes sorryAx by design and must not be treated as
an independent proof replay. This checker never makes such a claim.
"""
from __future__ import annotations
import argparse
from hashlib import sha256, sha1
import json
from pathlib import Path
import re

PINNED = {
    'Main.lean': '017645ca448f2fd1f6fd5c4c2cb2c6e96e472725be13f6e180ca6e1574360939',
    'Statements.lean': '79ef7ca776e42abe1f46caee13ea09bd19f1e2d0b1620d8d5d200258f33b47bf',
    'MainTheorems.lean': '1e734200b9b6a6d2a968743bd0c433d614bd6606cf538ac6d4f0ff70044dde53',
}
# Git blob SHA returned for the pinned Mathlib revision
MATHLIB_GIT_BLOB_SHA1 = '3c92195fce20df7fda41cfd823b3266f83cf3070'
MATHLIB_PINNED_REF = '0df444a360eaa60ab8c11dca51a86af692955474'
EXPECTED = {
    'riemannHypothesis': ('RiemannHypothesis', 'rhOfAll'),
    'riemannZeta_eq_zero_re': ('RiemannZetaZerosOnLine', 'rhOfAll'),
    'LFunction_ne_zero_of_half_lt_re': ('DirichletLNonvanishingRightOfHalf', 'allChi'),
    'LFunction_eq_zero_re': ('DirichletLZerosOnLine', 'claimAll.2.1'),
    'LFunction_eq_zero_re_of_isPrimitive': ('PrimitiveDirichletLZerosOnLine', 'critAll'),
}
KNOWN_DEFINITIONS = set(v[0] for v in EXPECTED.values()) - {'RiemannHypothesis'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_source(name: str, data: bytes):
    digest = sha256(data).hexdigest()
    require(digest == PINNED[name], f'{name}: SHA-256 mismatch; unreviewed source')
    return data.decode('utf-8')


def check_mathlib(data: bytes):
    header = b'blob '+str(len(data)).encode()+b'\0'
    git_blob_sha = sha1(header+data).hexdigest()
    require(git_blob_sha == MATHLIB_GIT_BLOB_SHA1, 'Mathlib git blob mismatch')
    needle = (
        'def RiemannHypothesis : Prop :=\n'
        '  ∀ (s : ℂ) (_ : riemannZeta s = 0) '
        '(_ : ¬∃ n : ℕ, s = -2 * (n + 1)) '
        '(_ : s ≠ 1), s.re = 1 / 2'
    )
    require(needle in data.decode('utf-8'), 'Mathlib RH definition mismatch')
    return {'mathlib_revision':MATHLIB_PINNED_REF, 'git_blob_sha1':git_blob_sha, 'rh_definition_matched':True}


def check_headlines(main_text, statements_text, kernel_text):
    # Exact hash must be verified BEFORE these text/regex checks: extraction
    # by regex alone cannot account for Lean macro elaboration or shadowing.
    found = re.findall(r'^theorem\s+(\w+)\s*:\s*(\w+)\s*:=\s*([^\n]+)$', main_text, re.M)
    require(len(found)==5, 'expected exactly five headline theorems')
    by_name = {n:(t,rhs.strip()) for n,t,rhs in found}
    require(len(by_name)==5 and by_name==EXPECTED, 'headline types or dependencies differ')
    definitions = set(re.findall(r'^def\s+(\w+)\s*:\s*Prop\s*:=', statements_text, re.M))
    require(definitions==KNOWN_DEFINITIONS, 'Dirichlet/RIemann statement declarations differ')
    require('theorem allChi : UpperTheory.EveryCharacter := UpperTheory.FinalTheorems.everyCharacter' in kernel_text.splitlines(),
            'source of the all-character nonvanishing result changed')
    require('theorem rhOfAll : RiemannHypothesis := Formalism.VW.rh_of_every allChi' in kernel_text.splitlines(),
            'RH reduction source changed')
    require('import Mathlib' in statements_text, 'Mathlib import missing')
    return [{'name':n,'type':t,'depends_immediately_on':rhs} for n,(t,rhs) in by_name.items()]


def audit(folder: Path, mathlib:Path|None=None):
    contents = {name: (folder/name).read_bytes() for name in PINNED}
    texts = {name: check_source(name, raw) for name,raw in contents.items()}
    headline_list = check_headlines(texts['Main.lean'], texts['Statements.lean'], texts['MainTheorems.lean'])
    return {
        'schema':'NUMERIA.ES8P255.ALPOGE_STATEMENT_AUDIT.1',
        'source':'Alpöge/Claude pinned original; MIT challenge audited in NUMERIA PR58',
        'verified_input_sha256':dict(PINNED),
        'five_headline_types_verified_from_exact_pinned_sources':headline_list,
        'official_mathlib_rh_definition':
            'forall complex s, zeta(s)=0, s not negative even and s!=1 => Re(s)=1/2',
        'mathlib_revision':MATHLIB_PINNED_REF,
        'mathlib_file_locally_verified':check_mathlib(mathlib.read_bytes()) if mathlib else None,
        'declared_grh_scope':'Dirichlet L-functions; NOT general automorphic/Dedekind L-functions',
        'solution_full_export_replayed_in_this_audit':False,
        'Lean_kernel_invoked_in_this_audit':False,
        'theorem_statement_verification_is_not_a_proof':True,
        'MIT_full_export_independent_kernel_report':'nanoda accepted per MIT status 2026-11-20 19:45 UTC',
        'RH_proved_independently_by_NUMERIA':False,
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--input-dir',type=Path,default=Path(__file__).with_name('fixtures'))
    cli.add_argument('--mathlib',type=Path,default=None,help='optional pinned Mathlib RiemannZeta.lean')
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args()
    result=audit(a.input_dir,a.mathlib)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('P255 verified pinned source/type declarations (NOT Lean proof replay)',len(result['five_headline_types_verified_from_exact_pinned_sources']))

if __name__=='__main__':main()
