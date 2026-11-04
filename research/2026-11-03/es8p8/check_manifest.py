#!/usr/bin/env python3
"""Verify versioned SHA-256 content against the immutable local manifest."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
p=json.loads((ROOT/'MANIFESTO_SHA256.json').read_text(encoding='utf-8'))
assert p['schema']=='NUMERIA.ES8P8.SHA256.1'
for name,hashsum in p['files'].items():
    path=ROOT/name
    assert path.is_file(),f'missing {name}'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==hashsum, f'SHA256 mismatch {name}'
print('PASS: SHA256 of',len(p['files']),'ES8-P8 artifacts')
