#!/usr/bin/env python3
"""Consolidate TWO supplied Scopus CSVs; do not infer full text availability."""
import argparse
import csv
import hashlib
from pathlib import Path

HIGH = {
 '314':('espectral','alta','Critérios espectrais/torais equivalentes; não produz RH automaticamente.'),
 '315':('espectral','media','Modelo de Dirac com ajuste espectral; hipótese física condicional.'),
 '316':('espectral','alta','Fórmula de traço, geometria não comutativa, lacuna de identificação.'),
 '320':('espectral','media','Estatística dos zeros para testes secundários.'),
 '321':('espectral','media','Assintótica espectral; controle de compatibilidade.'),
 '328':('nao_linear','media','Aproximações por ambos lados para operadores não lineares; precisa hipóteses concretas.'),
 '329':('dirac_numerico','baixa','Métodos espectrais para dinâmica Hamiltoniana não linear; não constrói operador zeta.'),
 '330':('nao_linear','media','Iterações de Newton de alta ordem para descoberta, sem certificado de existência.'),
 '341':('zeta_computacional','alta','Algoritmos rápidos para avaliação da zeta em grandes alturas, não prova de zeros.'),
 '342':('zeta_computacional','alta','Métodos e escolhas de precisão para avaliações da zeta.'),
 '343':('zeta_computacional','media','Múltiplas avaliações, possível aceleração exploratória.'),
 '344':('krawczyk','altissima','Localização certificada de zeros de funções analíticas; exemplos com zeta.'),
 '345':('localizacao_heuristica','media','Localização multiescala sob ruído; precede certificação intervalar.'),
 '349':('nao_linear','media','Galerkin/colocação para equações integrais não lineares; exige correspondência analítica.'),
}
MEDIUM_SPECTRAL={'313','317','318','319','322','323','324','325','326','327'}
LOW_DIRAC={'332','333','334','335','336','337','338','339','340','353','354','355','356'}
LOW_OTHER={'331','350','351','352'}
PERIPHERAL={'346','347','348'}


def classify(r):
    id_=r['Reference']
    if id_ in HIGH:
        return HIGH[id_]
    if id_ in MEDIUM_SPECTRAL:
        return ('espectral','media','Relaciona zeta a espectros/Dirac, sem resolver o obstáculo de Weil por si só.')
    if id_ in LOW_DIRAC:
        return ('dirac_numerico','baixa','EDPs de Dirac temporais: úteis só se surgir operador concreto e condições de conservação.')
    if id_ in LOW_OTHER:
        return ('numerico_indireto','baixa','Método geral; a aplicabilidade exige equações do NUMERIA explicitadas.')
    if id_ in PERIPHERAL:
        return ('periferico','baixa','Não aborda diretamente zeta, Krawczyk ou matriz de Weil.')
    return ('nao_classificado','baixa','Exige leitura e conferência bibliográfica adicional.')


def build(paths, output):
    files=[]; bydoi={}
    for path in paths:
        p=Path(path)
        filehash=hashlib.sha256(p.read_bytes()).hexdigest()
        files.append({'filename':p.name,'sha256':filehash})
        with p.open(newline='',encoding='utf-8-sig') as f:
            for r in csv.DictReader(f):
                key=(r['DOI'] or r['Title']).casefold().strip()
                if key in bydoi:
                    bydoi[key]['arquivos_origem'] += ';'+p.name
                    continue
                cat,pri,note=classify(r)
                bydoi[key]={
                    'numero_referencia':r['Reference'],
                    'titulo':r['Title'],'autoria_csv':r['Researchers'],
                    'ano':r['Publication Date'][:4],
                    'doi':r['DOI'], 'url':('https://doi.org/'+r['DOI']) if r['DOI'] else '',
                    'categoria':cat,'prioridade':pri,'aproveitamento_e_cautela':note,
                    'arquivos_origem':p.name,
                    'tipo_conteudo_registrado':r['Reference Type'],
                    'nota_correcao_autoria':'DOI 10.1007/s000290050042: Alain Connes (e não A. Cormes no CSV).' if r['Reference']=='316' else '',
                }
    rows=sorted(bydoi.values(),key=lambda r:int(r['numero_referencia']))
    with Path(output).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return files,rows


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--files',nargs='+',required=True)
    a.add_argument('--out',required=True)
    args=a.parse_args()
    fi,rows=build(args.files,args.out)
    print('bibliographic_entries',len(rows),'source_files',len(fi),'priority_high+',sum(z['prioridade'] in ('alta','altissima') for z in rows))

if __name__=='__main__':
    main()
