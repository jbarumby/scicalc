#!/usr/bin/env python3
"""Verify SHA-256 of files listed in a versioned ES8-P6 manifest."""
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
m=json.loads((root/'MANIFESTO_SHA256.json').read_text(encoding='utf-8'))
for path,sha in m['sha256'].items():
    f=root/path
    assert f.is_file(),f'Missing {f}'
    assert hashlib.sha256(f.read_bytes()).hexdigest()==sha,f'Hash mismatch {path}'
assert len(m['sha256'])>=7
print('ES8-P6 manifest PASS:',len(m['sha256']),'files')
