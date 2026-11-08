#!/usr/bin/env python3
"""ES8-P94: quadratura GL VALIDADA para o termo arquimediano, com Arb.

Inovações em relação à P89:
  (a) nós/pesos GL de Legendre isolados por troca de sinal RACIONAL,
      usando a rotina fornecida pelo usuário (cópia em originais/);
  (b) erro GL de ordem 12 por Cauchy, com M=1024 e raio=1/4;
  (c) painéis racionais comprimento <=1/5 (13 por segmento de 5/2).

Momento B(r) mantido como integração Arb validada de b_h(x)x^k exp(i*c*x),
com série de Taylor grau 14 nos segmentos. P89 usa mesmo tipo de momentização,
portanto NÃO é uma implementação totalmente independente da integração interna.

Este programa, e somente sua execução em python-flint com os 240 segmentos,
pode produzir um novo enclosure para o termo arquimediano [0,600].
"""
from __future__ import annotations
import argparse,sys,os,json,platform,math
from fractions import Fraction as Q
from pathlib import Path
from p94_bounds import (R,SEGMENTS,SEG_LEN,H,L,K,M,panel_lengths_for_segment,
                        segment_gl_error,moment_taylor_error,full_gl_error,proven_preconditions)

NAMES=(1,2,3,5,7,11,13,17)
ORDER=12

def parse_pairs(values:str)->list[tuple[int,int]]:
    if values.strip().lower()=='all':return [(i,j) for i in range(8) for j in range(i,8)]
    parsed=set()
    for spec in values.split(';'):
        parts=spec.split(',')
        if len(parts)!=2:raise ValueError('pares separados por ; e números separados por vírgula')
        a,b=map(int,parts)
        if a not in NAMES or b not in NAMES:raise ValueError('centro não pertence aos oito centros')
        i,j=sorted((NAMES.index(a),NAMES.index(b)))
        parsed.add((i,j))
    if not parsed:raise ValueError('nenhum par')
    return sorted(parsed)


def run(shard:int,shards:int=8,bits:int=256,pilot_segments:int=0,pairspec:str='1,2',out:Path|None=None):
    if not 0<=shard<shards<=240:raise ValueError('shard inválido')
    if bits<192:raise ValueError('mínimo 192 bits, preferível 256')
    if pilot_segments<0:raise ValueError('pilot_segments não pode ser negativo')
    proof=proven_preconditions()
    from flint import arb,acb,ctx
    import flint
    sys.path.insert(0,str(Path(__file__).parent/'originais'))
    from gauss_legendre_intervalar_flint import gauss_legendre_rule_arb
    ctx.prec=bits
    pairs=parse_pairs(pairspec)
    h=arb(2)/5
    pi=arb.pi()
    logpi=pi.log()
    from flint import fmpq
    def qarb(x:Q)->arb:
        return arb(fmpq(x.numerator,x.denominator))
    taylor_q=moment_taylor_error()
    taylor_err=arb(0,fmpq(taylor_q.numerator,taylor_q.denominator))
    segment_q=segment_gl_error(ORDER)
    segment_err=arb(0,fmpq(segment_q.numerator,segment_q.denominator))
    nodes_and_weights=gauss_legendre_rule_arb(ORDER)
    xs,ws=nodes_and_weights
    if not len(xs)==len(ws)==ORDER:raise ArithmeticError('ordem GL inválida')
    if not bool(sum(ws,arb(0)).contains(arb(2))):
        raise ArithmeticError('pesos não encerram soma exata 2')
    # Verificação redundante de exatidão polinomial GL até grau 2n−1.
    # São fatos algébricos; o teste apenas protege bugs na avaliação dos nós/pesos.
    for k in range(2*ORDER):
        expected=arb(0) if k%2 else arb(2)/(k+1)
        moment=sum((w*x**k for x,w in zip(xs,ws)),arb(0))
        if not bool(moment.contains(expected)):
            raise ArithmeticError(f'erro na regra GL, monômio {k}')
    centers=[arb(v).log() for v in NAMES]
    def bump(z,analytic):
        re=z.real
        if bool(re.lower()>-h) and bool(re.upper()<h):
            return (1-acb(1)/(1-(z/h)**2)).exp()
        if bool(re.upper()<-h) or bool(re.lower()>h):
            return acb(0)
        if analytic:return acb(arb('nan'))
        if not bool(z.imag==0):return acb(arb('nan'))
        return acb(arb('0.5 +/- 0.5'))
    def binary(x):
        return [str(y) for v in (x.mid(),x.rad()) for y in v.man_exp()]
    total={p:arb(0) for p in pairs}
    processed=[]
    for s in range(shard,SEGMENTS,shards):
        if pilot_segments and len(processed)>=pilot_segments:break
        c=arb(5)*s/2+arb(5)/4
        moments=[]
        for k in range(K+1):
            def g(z,analytic,k=k):
                return bump(z,analytic)*z**k*(acb(0,1)*c*z).exp()
            v=acb.integral(g,-h,h,abs_tol=arb('1e-13'),rel_tol=arb('1e-24'),
                           deg_limit=180,eval_limit=300000,depth_limit=90)
            if not v.is_finite():
                raise ArithmeticError(f'momento não finito, segmento={s},k={k}')
            if not bool(v.real.rad().upper()<arb('2e-12') and v.imag.rad().upper()<arb('2e-12')):
                raise ArithmeticError(f'momento impreciso, segmento={s},k={k}: {v}')
            moments.append(v)
        # 12 painéis de comprimento 1/5 e o último de comprimento 1/10
        pa=Q(5*s,2)
        for length in panel_lengths_for_segment():
            pb=pa+length
            mid=qarb((pa+pb)/2)
            half=qarb(length/2)
            for x,w in zip(xs,ws):
                r=mid+half*x
                u=r-c
                F=acb(0)
                term=acb(1)
                for k in range(K+1):
                    F+=term*moments[k]
                    term*=acb(0,u)/(k+1)
                B=F.real+taylor_err
                argument=acb(0,r/2)
                m=((acb(arb(1)/4)+argument).digamma()
                  +(acb(arb(1)/4)-argument).digamma()).real/2-logpi
                base=w*half/pi*m*B*B
                if not base.is_finite():raise ArithmeticError(f'integrando não finito, segmento {s}')
                for i,j in pairs:
                    total[(i,j)]+=base*((centers[j]-centers[i])*r).cos()
            pa=pb
        if pa!=Q(5*(s+1),2):raise ArithmeticError('buraco no suporte de Fourier')
        for p in pairs:total[p]+=segment_err
        processed.append(s)
        print('SEGMENT',s,'OF',SEGMENTS,'PAIR_COUNT',len(pairs),flush=True)
    complete=(pilot_segments==0)
    recs=[{'i':i,'j':j,'n_i':NAMES[i],'n_j':NAMES[j],
           'arch_finite':binary(total[(i,j)])} for i,j in pairs]
    report={'status':'RIGOROUS_FINITE_GL_ARCH_0_TO_600' if complete else 'PILOT_PARTIAL_SEGMENTS_ONLY',
        'algorithm':'ES8-P94-rational-GL12-Cauchy-Taylor14',
        'bits':bits,'python_flint':flint.__version__,'python':sys.version,
        'shard':shard,'shards':shards,'segments':processed,
        'centers_integers':NAMES,'h':'2/5','integral_range':'0..600',
        'gauss_legendre_order':ORDER,'max_panel_length':'1/5',
        'proof':proof,
        'records':recs,
        'warnings':['apenas termo arquimediano FINITO: cauda r>600 deve ser adicionada',
                    'momento interno compartilha método Arb acb.integral com P89',
                    'positividade de matriz finita não prova RH']}
    if out:
        if out.exists():raise FileExistsError('não sobrescrever evidência: '+str(out))
        out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--shard',type=int,default=0)
    p.add_argument('--shards',type=int,default=8)
    p.add_argument('--prec',type=int,default=256)
    p.add_argument('--pairs',default='1,2',help='centros 1,2 ou all ou 1,2;1,1')
    p.add_argument('--pilot-segments',type=int,default=0)
    p.add_argument('--out',type=Path,default=Path('ARCH_P94_SHARD.json'))
    a=p.parse_args()
    try:
        r=run(a.shard,a.shards,a.prec,a.pilot_segments,a.pairs,a.out)
        print(r['status'],'segments',len(r['segments']))
    except Exception as e:
        print('INCONCLUSIVO_ARB_P94',type(e).__name__,str(e),file=sys.stderr)
        raise SystemExit(2)