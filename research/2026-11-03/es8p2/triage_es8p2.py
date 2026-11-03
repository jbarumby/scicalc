"""Deterministic deduplication of uploaded LeapSpace metadata.

The decisions are research triage, not verification of full paper claims.
"""
from __future__ import annotations
import argparse, csv, json, hashlib
from collections import defaultdict
from pathlib import Path

# Curation informed by the supplied abstracts and spot-checked publisher pages.
RANK={
 '10.1007/s000290050042':('alta','fórmula de traço, adeles e interpretação de zeros','Verificar sinais e normalização pelo artigo de Connes; nao substitui auditoria V17/V24.'),
 '10.5565/publmat6321906':('alta','transformada de Fourier e funções de suporte limitado','Ferramentas de extremalidade; teoremas de cotas em parte assumem RH.'),
 '10.1016/j.aim.2019.106926':('alta','programação semidefinida e seleção de testes','Métodos SDP são úteis no NUMERIA; resultados de correlação não certificam Weil e algumas consequências assumem RH.'),
 '10.1016/0022-1236(79)90049-1':('alta','operadores integrais em teoria dos números','Investigar operadores e kernels originais; metadados insuficientes para extrair lema aplicável.'),
 '10.1090/s0273-0979-99-00766-1':('média','simetrias de espectros zeta e geometria','Panorama espectral; não fixa matriz finita do projeto.'),
 '10.1088/1751-8121/ad4c2d':('média','candidato Hamiltoniano Hilbert–Pólya','Auditar domínio, autoadjunção e dependência da condição de contorno; não pressupor prova de RH.'),
 '10.1103/physrevlett.118.130201':('média','Hamiltoniano e condições de contorno','Rota espectral separada da certificação da fórmula de Weil.'),
 '10.1016/j.nuclphysb.2007.03.049':('média','modelo xp e funções de Jost','Explorar somente como controle para modelos espectrais.'),
 '10.1103/revmodphys.83.307':('média','panorama física da hipótese de Riemann','Revisão contextual, sem certificados novos para a matriz.'),
 '10.1016/j.jnt.2023.09.003':('média','cotas explícitas de zeta na reta crítica','Potencial controle de caudas ou diagnósticos; ainda sem vinculação a E7/V24.'),
 '10.1016/j.jnt.2011.10.002':('média','mollificadores e proporção de zeros na reta','Mostra apenas fração de zeros na linha, não a totalidade.'),
 '10.1016/j.jnt.2016.01.026':('média','zeros na reta crítica','Ferramentas analíticas, mas não identidade das matrizes.'),
 '10.1016/j.na.2010.07.037':('média','inversão de Fourier e equação funcional','Conferir convenções exatas de transformada antes de aplicar.'),
 '10.1016/j.jnt.2016.09.022':('média','mollificadores e zeta perturbada','Controle adicional de distribuições de zeros, não Weil local.'),
 '10.1016/j.jnt.2021.02.007':('média','argumento da função zeta','Potencial uso em limites espectrais; não prova positividade.'),
 '10.1016/j.jnt.2023.07.004':('média','estimativas de zeta na faixa crítica','Sem conexão direta com kernel finito V17/V24.'),
 '10.1007/s00601-026-02032-z':('contextual','ressonâncias de função de Jost multicanal','Modelo físico de espalhamento, sem relação matemática demonstrada com zeros zeta.'),
 '10.1016/j.aop.2020.168274':('baixa','dinâmica caótica','Contexto físico distante do gate de Weil.'),
}

def classify(record):
    doi=(record.get('DOI') or '').strip().lower()
    title=(record.get('Title') or '').lower()
    if doi in RANK:return RANK[doi]
    if any(w in title for w in ('scattering','jost','resonance','coulomb','regge','gamow','multi-channel','coupled channel')):
        return ('contextual','métodos de espalhamento e funções de Jost','Não identifica matrizes E7/V24 nem a zeta sem uma nova construção.')
    if any(w in title for w in ('riemann','zeta','fourier','gamma function','special functions','dirichlet')):
        return ('contextual','funções especiais / zeta / Fourier','Ler proposições específicas antes de alegar aplicação ao gate da ES8-P2.')
    return ('baixa','tema distante do gate ES8-P2','Sem contribuição direta comprovada pelos metadados enviados.')

def triage(root:Path):
    groups=defaultdict(list)
    for n in (10,11,12):
        path=root/f'References_LeapSpace ({n}).csv'
        with path.open(encoding='utf-8-sig',newline='') as fh:
            for row in csv.DictReader(fh):
                row['_source_file']=path.name
                key=(row.get('DOI') or '').lower().strip()
                if not key:key='title:'+row.get('Title','').lower().strip()
                groups[key].append(row)
    order={'alta':0,'média':1,'contextual':2,'baixa':3}
    results=[]
    for doi,recs in groups.items():
        ref=recs[0];rank,theme,comment=classify(ref)
        authors=ref.get('Researchers','')
        correction=''
        if doi=='10.1007/s000290050042':
            authors='Alain Connes'
            correction='Autoria corrigida: CSV informava Cormes; verificado no site do autor.'
        results.append({'Prioridade':rank,'Título':ref.get('Title',''),
            'Autores (metadados/correções)':authors,'DOI':ref.get('DOI',''),
            'Ano':ref.get('Publication Date','')[:4],
            'Referências originais':'; '.join(r['Reference'] for r in recs),
            'Arquivos originais':'; '.join(r['_source_file'] for r in recs),
            'Tema relevante':theme,'Avaliação para ES8-P2 / RH':comment,
            'Correção de metadados':correction,
            'Nível de verificação':'CSV resumo; conferência externa seletiva, não leitura integral',
            'Fonte':('https://doi.org/'+ref.get('DOI','')) if ref.get('DOI') else ''})
    results.sort(key=lambda r:(order[r['Prioridade']],r['Título'].casefold()))
    output=root/'REFERENCIAS_TRIADAS_ES8P2.csv'
    with output.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    summary={'raw':sum(map(len,groups.values())),'unique':len(results),
        'duplicates':sum(len(v)-1 for v in groups.values()),
        'priorities':{tier:sum(r['Prioridade']==tier for r in results) for tier in order},
        'notes':['Metadata error Connes/Cormes corrected for display; original preserved.',
                 'Several SDP/pair-correlation results invoke RH as an assumption; do not use circularly.',
                 'No reviewed abstract proves RH or V17/V24 matrix identity.']}
    (root/'RESUMO_REFERENCIAS_ES8P2.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'))
    args=p.parse_args();print(json.dumps(triage(args.root),ensure_ascii=False))
