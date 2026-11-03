"""Rank the user-supplied Scopus/LeapSpace bibliographic entries.

Classification is scope-based, not evidence that the paper proves RH or Weil
positivity. No fetching of abstracts or online content is required.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

DIRECT = {
    '10.1016/S0377-0427(00)00336-8': ('ALTA - auditoria de zeta', 'Algoritmos independentes para avaliar zeta e validar oráculos numéricos; não reconstrói E7.'),
    '10.1016/j.jmaa.2017.12.045': ('ALTA - linha ξ derivada', 'Zeros na reta crítica de combinações de derivadas de ξ; estudo de proporções, não prova RH.'),
}
INDIRECT = {
    '10.1016/j.jmaa.2020.124541': ('MÉDIA - Dirichlet', 'Relações entre séries de potências e séries de Dirichlet; investigar aplicabilidade, sem fórmula de Weil demonstrada.'),
    '10.1016/j.jmaa.2019.04.019': ('MÉDIA - Mellin', 'Sequências de Appell, Mellin e Fourier; matemática auxiliar, não certifica matrizes E7.'),
    '10.1007/s10474-023-01363-1': ('MÉDIA - Möbius', 'Somas parciais de Möbius generalizado, potencial rota Nyman/Mellin; não resolve projeção V17.'),
    '10.2140/pjm.1977.68.277': ('MÉDIA - Möbius', 'Definições de Möbius de ordem k; útil como contexto da rota alternativa.'),
    '10.1016/j.jmaa.2005.01.020': ('BAIXA - funções especiais', 'Polinômios Apostol-Bernoulli/Euler; nenhuma identidade direta com a zeta de Riemann no escopo resumido.'),
}

def triage(paths):
    raw = []
    for path in paths:
        with open(path,encoding='utf-8-sig',newline='') as f:
            rows=list(csv.DictReader(f))
        for row in rows:
            raw.append({'filename':Path(path).name,**row})
    unique={}
    for row in raw:
        doi=row.get('DOI','').strip().lower()
        key=doi if doi else 'title:'+row.get('Title','').strip().casefold()
        if key not in unique:
            unique[key]={**row, 'copies':1,'original_references':[f"{row['filename']}#{row['Reference']}"]}
        else:
            unique[key]['copies']+=1
            unique[key]['original_references'].append(f"{row['filename']}#{row['Reference']}")
    sorted_records=[]
    for key,row in unique.items():
        doi=row.get('DOI','').strip()
        bucket,why=(DIRECT.get(doi) or INDIRECT.get(doi) or
                    ('BAIXA - indireta', 'Amostra bibliográfica: título e resumo não estabelecem vínculo direto com o bloqueio ES8-P4.'))
        sorted_records.append({
            'reference':row['Reference'], 'doi':doi,
            'title':row['Title'],'year':(row.get('Publication Date') or '')[:4],
            'level':bucket,'rationale':why,'occurrences':row['copies'],
            'original_references':'; '.join(row['original_references']),
            'source_type':row.get('Source Type',''),
        })
    sorted_records.sort(key=lambda r:({'ALTA':0,'MÉDIA':1,'BAIXA':2}.get(r['level'].split()[0],3),r['title'].casefold()))
    return {'input_records':len(raw),'unique_records':len(sorted_records),
            'duplicate_excess':len(raw)-len(sorted_records),
            'records':sorted_records,'status':'TITLES_ABSTRACTS_ONLY_NO_MATH_PROOF_VERIFIED'}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--base',type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args()
    data=triage([a.base/'References_LeapSpace (14).csv',a.base/'References_LeapSpace (15).csv'])
    dest=a.base/'REFERENCIAS_TRIADAS_ES8P4.csv'
    with dest.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data['records'][0]));w.writeheader();w.writerows(data['records'])
    (a.base/'RESUMO_REFERENCIAS_ES8P4.json').write_text(json.dumps({k:v for k,v in data.items() if k!='records'},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f"References {data['input_records']} rows, {data['unique_records']} distinct, {data['duplicate_excess']} duplicate copies")
if __name__=='__main__':main()
