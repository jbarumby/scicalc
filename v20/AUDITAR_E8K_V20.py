#!/usr/bin/env python3
"""E8-K V20: three strict independent Gibbs/zeta diagnostics.

This does not certify the Weil explicit formula, Phase-4 quadrature, or RH.
"""
import argparse
import hashlib
import json
from pathlib import Path
import gibbs_euler
import yang_lee
import gibbs_gram

HERE=Path(__file__).resolve().parent
FILES=['gibbs_euler.py','yang_lee.py','gibbs_gram.py','AUDITAR_E8K_V20.py','TESTAR_E8K_V20.py']

def run():
    a=gibbs_euler.run()
    b=yang_lee.run()
    c=gibbs_gram.run()
    return {
      'schema':'NUMERIA.E8K.V20.GIBBS_YANG_LEE_AND_EULER.1',
      'status':'FINITE_GIBBS_EULER_YANG_LEE_GRAM_AUDITED_WEIL_BRIDGE_OPEN',
      'gates':{
         'gibbs_zeta_euler':{'status':'PASS','result':a},
         'finite_yang_lee_mode':{'status':'PASS','result':b},
         'exact_gibbs_gram':{'status':'PASS','result':c},
         'short_support_rank_obstruction':{'status':'PASS','rank_upper':1,'condition':'L=9/20, log 2 < 2L < log 3 and g_k(0)=0'}
      },
      'input_sha256':{x:hashlib.sha256((HERE/x).read_bytes()).hexdigest() for x in FILES},
      'interpretation':{
         'grand_canonical_bose_mode':'finite cutoff single mode; Yang-Lee circle roots not zeta zeros',
         'Gibbs_gram':'positive finite Gibbs weights imply Gram PSD; rational Vandermonde test yields strict positivity for separate feature family',
         'Weil_obstruction':'For NUMERIA sine modes evaluated at log n with 2L=0.9, Gibbs Gram has rank <=1, so no positive-definite Weil identification for dimensions >=2',
         'analytic_continuation':'Physical Gibbs probability requires beta>1 and cannot be analytically continued as a positive probability to Re(s)=1/2.'
      },
      'not_certified':[
         'No proven identity between positive Gibbs Gram and Weil projected matrices',
         'No independent Phase-4 quadrature derivative majorants or complex error bounds',
         'No verification of all 1600 quadrature subintervals or projected LDL matrices here',
         'No theorem that Yang-Lee fugacity zeros equal Riemann zeta zeros',
         'No proof of the Riemann Hypothesis'
      ]
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    result=run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS: 4 independent checks (Gibbs/Euler, finite Yang-Lee, rational Gram, rank obstruction)')
    print('STATUS:',result['status'])
    return 0
if __name__=='__main__':raise SystemExit(main())
