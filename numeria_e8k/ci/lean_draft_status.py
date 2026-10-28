#!/usr/bin/env python3
"""Inventário do rascunho Lean; não compila e não valida os teoremas."""
import argparse, json, re
from pathlib import Path
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a=p.parse_args()
    lines=a.source.read_text(encoding='utf-8').splitlines()
    outstanding=[]
    in_block=False
    for i,line in enumerate(lines,1):
        code=line
        if '/-' in code:
            in_block=True
        if not in_block and re.search(r'\bsorry\b',code.split('--',1)[0]):
            outstanding.append(i)
        if '-/' in code:
            in_block=False
    result={'schema':'NUMERIA.E8K.LEAN.DRAFT.1',
            'compilation':'NOT_RUN','formal_proof_complete':False,
            'remaining_sorry_lines':outstanding,
            'statement':'A presença de sorry impede classificar os lemas como provas fechadas.'}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Lean: rascunho não compilado; sorry nas linhas',outstanding)
if __name__=='__main__':main()
