"""NUMERIA ES7: documentary normalization gate. NOT a mathematical proof of matrix equivalence."""
from pathlib import Path
import argparse, hashlib, json

def inspect(v24:Path,es6:Path,output:Path):
    v=v24.read_text(encoding='utf-8')
    e=es6.read_text(encoding='utf-8')
    checks={
      'v24_uses_digamma_arch_kernel': 'psi(1/4+ir/2)' in v,
      'v24_arch_matrix_weight_is_Fi_times_conjugate_Fj': 'F_i(r)\\overline{F_j(r)}' in v,
      'v24_arch_integral_has_1_over_pi': 'pi^{-1}' in v,
      'es6_gaussian_weight': 'h_a(t)=exp(-a*t^2)' in e,
      'es6_integral_kernel': 'Re psi(1/4+it/2)' in e,
    }
    # Expressions are present in the README and ES6 implementation; the
    # consistency check is lexical and does not establish analytic equality.
    if not all(checks.values()):
        raise ValueError('Documentation normalization anchor missing; re-audit formulas')
    result={'schema':'NUMERIA.ES7.V24_NORMALIZATION_GATE.1',
      'source_sha256':{str(v24):hashlib.sha256(v24.read_bytes()).hexdigest(),str(es6):hashlib.sha256(es6.read_bytes()).hexdigest()},
      'lexical_anchors':checks,
      'agreed_arch_kernel':'Omega(r)=Re digamma(1/4+i*r/2)-log(pi)',
      'agreed_even_normalization':'1/pi integral_[0,infinity] Omega(r)*weight(r) dr',
      'weight_ES7':'h_a(r)=exp(-a*r*r), scalar Gaussian test',
      'weight_V24':'F_i(r)*conjugate(F_j(r)), matrix built from compactly supported test modes',
      'matrix_equivalence':'NOT_ESTABLISHED',
      'basis_projection':'NOT_AUDITED',
      'pole_and_prime_sign_transport':'NOT_AUDITED',
      'claim':'Same documented archimedean kernel and even-normalization; no V24 entry or LDL certificate transferred'}
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--v24',type=Path,required=True)
    p.add_argument('--es6',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('COTEJO_ES7_V24.json'))
    args=p.parse_args()
    print('Gate:',inspect(args.v24,args.es6,args.out)['matrix_equivalence'])
