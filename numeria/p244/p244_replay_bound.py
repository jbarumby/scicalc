#!/usr/bin/env python3
"""P244 replay wrapper with explicit exporter and certificate pins.

The original fail-closed replay is preserved unchanged under package/. This
wrapper additionally authenticates export_centers.py and P242_AUDIT_CERT.json,
then labels hash/byte matching separately from whether MPFR was rerun now.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent / 'package'
sys.path.insert(0, str(PACKAGE))
from p244_replay_fail_closed import run as base_run  # noqa: E402

EXPECTED_EXPORTER_SHA256 = '58be39fcfd6b424cdec9af9b20b5a3b9d78e6fa3f0a40d4127916f8c97b92fac'
EXPECTED_P242_CERT_SHA256 = '9f89c17339c9daca7d1e2e020dbfdbecc83fc47a79301b6ebde544dad0099447'


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_static_pins(cert: Path) -> dict[str, str]:
    exporter_hash = sha256(PACKAGE / 'export_centers.py')
    cert_hash = sha256(cert)
    if exporter_hash != EXPECTED_EXPORTER_SHA256:
        raise ArithmeticError('P244 refuses modified center exporter')
    if cert_hash != EXPECTED_P242_CERT_SHA256:
        raise ArithmeticError('P244 refuses unpinned P242 certificate')
    return {
        'verified_exporter_sha256': exporter_hash,
        'verified_p242_certificate_sha256': cert_hash,
    }


def run_bound(data: Path, traces: Path, cert: Path, work: Path, recompute: bool) -> dict:
    pins = verify_static_pins(cert)
    result = base_run(data, traces, cert, work, recompute)
    result.update(pins)
    for sector in ('even', 'odd'):
        sector_result = result['MPFR_traces'][sector]
        # The base wrapper's name implies a fresh execution even in check-only
        # mode; this field states only what its digest and byte comparisons show.
        sector_result['trace_digest_and_bytes_match'] = sector_result.pop('identical_reexecution')
    result['trace_check_mode'] = 'regenerated_this_run' if recompute else 'preexisting_trace_check_only'
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--traces', type=Path, required=True)
    parser.add_argument('--cert', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--check-existing', action='store_true', help='check existing pinned traces without compiling/running MPFR')
    args = parser.parse_args()
    result = run_bound(args.data, args.traces, args.cert, args.work, not args.check_existing)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print('P244 bound replay: PASS, 754 pivots; mode=' + result['trace_check_mode'])


if __name__ == '__main__':
    main()
