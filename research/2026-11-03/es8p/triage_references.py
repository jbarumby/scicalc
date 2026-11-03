"""Classify user-supplied reference CSV as metadata, not as verified full papers."""
import csv
import argparse
from pathlib import Path

LABELS={
'378':('BAIXA','Contexto de caos clássico; não fixa coeficientes da fórmula explícita','Sem operador da zeta nem prova de Weil estabelecida'),
'379':('EXPLORATÓRIA','Analogia Majorana/espalhamento que pode informar ES3','Ligação com zeros condicional e física especulativa'),
'380':('MÉDIA','Potenciais com zeros como polos redundantes; controle adversarial espectral','Codificar zeros em polos não implica autoadjunção nem RH'),
'381':('MÉDIA','Fase da zeta e analogias de traço; comparação com ES4/ES5','Não certifica entradas finitas de Weil'),
'382':('MÉDIA','Modelo de espalhamento hiperbólico; ressonâncias e atrasos','Analogia espectral; não substitui demonstração de equivalência'),
'383':('ALTA (ES ESPECTRAL)','Sierra: H=xp com interação; função de Jost e consequências analíticas','Resonâncias seguem posições médias; espectro exato de zeta não provado'),
'384':('DUPLICATA','Mesmo DOI e título do registro 380','Não contar como artigo novo'),
'385':('MÉDIA-ALTA (ADÉLICO)','S-matrizes p-ádicas e produto adélico; possível cotejo local-global','Não fornece diretamente a matriz de Weil E7/V24'),
'386':('ALTA (ES ESPECTRAL)','Sierra: reconstrução quântica e função de Jost para sigma>1','Caso sigma=1/2 não é construção global certificada'),
'387':('BAIXA','Modelo de vácuo quântico com analogias de zeta','Distante da certificação de Fourier, polos e primos'),
'388':('MÉDIA-ALTA (ES ESPECTRAL)','LeClair: fluxo espectral baseado no produto de Euler; teste de hipóteses','Cuidado com extensão sigma→1/2 e auto-adjunção não demonstrada'),
'389':('MÉDIA-ALTA (DOMÍNIOS)','Bender/Brody/Müller: operador com condição de contorno ligada a zeros','O próprio artigo deixa autoadjunção como etapa heurística')}

def triage(input_path,output_path):
    with open(input_path,encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    counts={};out=[]
    for row in rows:
        id=row['Reference'];doi=row['DOI'].strip().lower()
        priority,use,caution=LABELS.get(id,('A CLASSIFICAR','Sem triagem','Verificar fonte'))
        if doi:
            counts[doi]=counts.get(doi,0)+1
        title=row['Title'].strip()
        if id=='385' and not title:title='[título ausente no CSV] Scattering on p-adic and adelic symmetric spaces (identificação provisória pelo resumo)'
        out.append({'id_original':id,'titulo_informado_ou_inferido':title,'doi':doi,
          'ano':row['Publication Date'][:4], 'triagem':priority,'utilidade_numera':use,
          'limite_interpretacao':caution,'url_fonte':'https://doi.org/'+doi if doi else '',
          'duplicado_do':'380' if id=='384' else '', 'nivel_verificacao':'resumo CSV; títulos e informações centrais de artigos prioritários confrontados com fontes públicas; não revisão do texto integral'})
    dest=Path(output_path);dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=out[0].keys(),lineterminator='\n');w.writeheader();w.writerows(out)
    return {'total_registros':len(rows),'distintos_por_doi':len(counts),'duplicatas':len(rows)-len(counts),'mais_proximos_de_ES8P':[385], 'artigos_que_constituem_prova_de_RH':0}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(triage(a.source,a.out))
