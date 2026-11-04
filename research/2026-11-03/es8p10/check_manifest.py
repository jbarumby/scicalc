#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
items=json.loads((root/'MANIFESTO_SHA256.json').read_text())['files']
for name,sha in items.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha, name
print('ES8-P10:',len(items),'source hashes match')
