#!/usr/bin/env python3
"""Byte-level SHA-256 integrity of ES8-P9 research files."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
v=json.loads((root/'MANIFESTO_SHA256.json').read_text(encoding='utf-8'))
assert v['schema']=='NUMERIA.ES8P9.SHA256.1'
for name,digest in v['files'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
print('ES8-P9 SHA-256 verified:',len(v['files']))
