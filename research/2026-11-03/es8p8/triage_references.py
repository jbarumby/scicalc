#!/usr/bin/env python3
"""Scope limits of the two user-supplied bibliographic records."""
import csv
from pathlib import Path
PATH=Path(__file__).resolve().parent
USES={
'10.1016/j.amc.2006.08.061':(
    'Apoio para representacoes da digama e log-gama por contornos de Hankel',
    'Nao implementado como uma terceira avaliacao intervalar de todos os nos',
    'https://www.sciencedirect.com/science/article/abs/pii/S0096300306010836'),
'10.1007/pl00000117':(
    'Contexto teorico de continuacao analitica para a zeta de Hurwitz-Lerch',
    'Nao fornece por si so a cota nodal da forma de Weil',
    'https://keio.elsevierpure.com/en/publications/on-the-hurwitz-lerch-zeta-function/')
}
def run():
    with (PATH/'REFERENCIAS_ORIGINAIS_17.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    out=[]
    for row in rows:
        doi=row['DOI'].strip().lower()
        if doi not in USES:raise ValueError('Unexpected DOI not triaged: '+doi)
        use,lim,link=USES[doi]
        out.append({'Title':row['Title'],'DOI':doi,'Researchers':row['Researchers'],
                    'Aplicacao':use,'Limitacao':lim,'Fonte_verificada':link})
    with (PATH/'REFERENCIAS_TRIADAS_ES8P8.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    assert len(out)==2
    print('Bibliography: two records examined and triaged')
if __name__=='__main__':run()
