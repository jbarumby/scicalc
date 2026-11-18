#!/usr/bin/env python3
"""Fail-closed provenance link from independent P237 shadow audit to finite LDL.

Checks bytes SHA-256 of BOTH full C files against the separately generated
shadow audit. Validates the complete source range manifest hashes, but does
not itself establish Arb numerical correctness or Weil/Zhu identities.
"""
import hashlib,json
from pathlib import Path


def check_shadow(audit_file,full_folder):
    d=json.loads(Path(audit_file).read_text())
    if d.get('schema')!='NUMERIA.ES8P237.SHADOW_EXACT_AUDIT.1' or d.get('both_sectors_complete') is not True or d.get('independent_source_bound_enclosures_confirmed') is not True or d.get('does_not_verify_analytic_bridge') is not True or d.get('RH_proved') is not False:
        raise ArithmeticError('shadow audit full gate failed')
    arr=d.get('sectors')
    if not isinstance(arr,list) or len(arr)!=2 or [v.get('sector') for v in arr]!=['even','odd']:
        raise ValueError('shadow audit sectors do not match')
    digest={}
    for s in arr:
        sec=s['sector']
        file=f'P237_FULL_C_{sec}.json.gz'
        if s.get('merged_file')!=file or s.get('panels_recomputed')!=1008 or s.get('triangle_entries_exactly_rechecked')!=71253 or s.get('enclosure_contains_independent_unrounded_exact_sum') is not True or s.get('addition_used_only_python_integers') is not True or s.get('all_panel_ranges_contiguous') is not True or s.get('has_full_quadrature') is not True:
            raise ArithmeticError(f'shadow audit incomplete for {sec}')
        files=s.get('source_hashes_verified')
        if type(files)!=list or len(files)!=41 or files[0].get('sha256') not in ({'even':'6cfc2d3080389751b08d62dbdcdabdbd16a8be64ea5b23993dac0f27f61383f5', 'odd':'76a7c147bf805fca4b67f65cbfd01464cf95134b30e990d3b5f24226a4de1f3d'}[sec],):
            raise ValueError('missing P236 locked baseline in shadow audit')
        for f in files:
            if len(f.get('sha256',''))!=64 or f.get('entries')!=71253:raise ValueError('bad source hash/count')
        path=Path(full_folder)/file
        dg=hashlib.sha256(path.read_bytes()).hexdigest()
        if dg!=s.get('merged_sha256'):raise ArithmeticError(f'full matrix SHA256 changed after audit: {sec}')
        digest[sec]=dg
    return digest
