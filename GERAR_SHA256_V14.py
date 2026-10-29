#!/usr/bin/env python3
"""Manifesta arquivos V14; exclui o proprio SHA256SUMS e o relatorio derivado de conferencia.
Arquivo manifesto cobre todos os inputs, fontes, outputs reprodutiveis e metadados.
"""
import argparse,hashlib
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));a=p.parse_args();r=a.root.resolve()
 out=[]
 for f in sorted(r.rglob('*')):
  if not f.is_file() or f.is_symlink():continue
  rel=f.relative_to(r).as_posix()
  if rel in ('SHA256SUMS.txt','CONFERENCIA_E8K_V14.json','CONFERENCIA_E8K_V14_REPETIDA.json') or '/__pycache__/' in '/'+rel or rel.endswith('.pyc'):continue
  out.append(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {rel}')
 (r/'SHA256SUMS.txt').write_text('\n'.join(out)+'\n',encoding='utf8')
 print('sha entries',len(out))
if __name__=='__main__':main()
