#!/usr/bin/env python3
"""Verifica o ZIP E7 exato e extrai somente a fonte B5 após conferir hashes."""
import argparse, hashlib, json
from pathlib import Path
from zipfile import ZipFile

EXPECTED_ZIP = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
EXPECTED_B5 = 'a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
MEMBER = 'fontes/AUDITOR_PONTOS_B5.py'

def h(b): return hashlib.sha256(b).hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--zip', type=Path, required=True)
    p.add_argument('--extract-to', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    raw = a.zip.read_bytes()
    if h(raw) != EXPECTED_ZIP:
        raise SystemExit('ERRO: SHA-256 do ZIP E7 diverge do valor aprovado')
    with ZipFile(a.zip) as z:
        src = z.read(MEMBER)
    if h(src) != EXPECTED_B5:
        raise SystemExit('ERRO: SHA-256 da fonte B5 diverge do valor aprovado')
    a.extract_to.parent.mkdir(parents=True, exist_ok=True)
    a.extract_to.write_bytes(src)
    a.report.parent.mkdir(parents=True,exist_ok=True)
    a.report.write_text(json.dumps({
        'schema':'NUMERIA.E8K.CI.SOURCE.1',
        'source_zip_sha256':h(raw), 'b5_source_sha256':h(src),
        'source_member': MEMBER, 'verified': True}, indent=2)+'\n')
    print('ZIP e fonte B5: SHA-256 verificados')

if __name__ == '__main__': main()
