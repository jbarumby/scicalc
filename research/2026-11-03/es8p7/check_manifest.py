from pathlib import Path
import hashlib, json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFESTO_SHA256.json').read_text(encoding='utf-8'))
assert manifest['schema']=='NUMERIA.ES8P7.MANIFEST.1'
for path,digest in manifest['files'].items():
    assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest,path
print('SHA256 MANIFEST: PASS',len(manifest['files']))
