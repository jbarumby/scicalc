#!/usr/bin/env python3
"""Independent rational fixed-point interval reconstruction B^T(G+P)B.
No Arb and no floating point enters any certified calculation.
Only accepts validated Phase 3 alternative bounds as numerical input.
"""
from __future__ import annotations
import math,json,sys
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
import hashlib

S=10**90

def floordiv(x,y):return x//y

def ceildiv(x,y):return -((-x)//y)

class I:
  __slots__=('lo','hi')
  def __init__(self, a=0,b=None,raw=False):
    if raw:self.lo,self.hi=int(a),int(b);return
    q=F(a);r=F(b) if b is not None else q
    self.lo=floordiv(q.numerator*S,q.denominator)
    self.hi=ceildiv(r.numerator*S,r.denominator)
  @classmethod
  def raw(cls,lo,hi):return cls(lo,hi,raw=True)
  def __add__(a,b):
    if not isinstance(b,I):b=I(b)
    return I.raw(a.lo+b.lo,a.hi+b.hi)
  __radd__=__add__
  def __neg__(a):return I.raw(-a.hi,-a.lo)
  def __sub__(a,b):return a + (-b if isinstance(b,I) else -I(b))
  def __rsub__(a,b):return I(b)-a
  def __mul__(a,b):
    if not isinstance(b,I):b=I(b)
    pp=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi]
    return I.raw(floordiv(min(pp),S),ceildiv(max(pp),S))
  __rmul__=__mul__
  def __truediv__(a,b):
    if not isinstance(b,I):b=I(b)
    if b.lo<=0<=b.hi:raise ZeroDivisionError(f'{b}')
    inv=I.raw(floordiv(S*S,b.hi),ceildiv(S*S,b.lo))
    return a*inv
  def __rtruediv__(a,b):return I(b)/a
  def __pow__(a,n):
    if n<0:return I(1)/a**(-n)
    v=I(1);b=a
    while n:
      if n%2:v=v*b
      b=b*b;n//=2
    return v
  def min(self):return F(self.lo,S)
  def max(self):return F(self.hi,S)
  def mid(self):return F(self.lo+self.hi,2*S)
  def width(self):return F(self.hi-self.lo,S)
  def __repr__(a):return f'[{float(a.min()):.14g}, {float(a.max()):.14g}]'

def atan_frac(q,K=110):
    x=F(1,q);s=F(0)
    for n in range(K):s+=((-1)**n)*x**(2*n+1)/F(2*n+1)
    t=x**(2*K+1)/F(2*K+1)
    return I(s,s+t) if K%2==0 else I(s-t,s)

def ln_via_atanh(q,K=120):
    z=F(1,q)
    s=2*sum((z**(2*n+1)/F(2*n+1) for n in range(K)),F(0))
    t=2*z**(2*K+1)/(F(2*K+1)*(1-z*z))
    return I(s,s+t)

PI=16*atan_frac(5)-4*atan_frac(239)
LN2=ln_via_atanh(3)
SQ2=I.raw(math.isqrt(2*S*S),math.isqrt(2*S*S)+1)
L=I(F(9,20))

# The argument is reduced using an integer multiple of 2*pi; the interval test
# red subset [-4,4] is exact and the Taylor error uses exact |sin^(n)|<=1.
E_SIN=I(F(4**120,math.factorial(120)))

def reduce(x):
  nearest=round(float(x.mid()/(2*PI.mid()))) # selects q only; bounds checked below
  y=x-2*nearest*PI
  if not (-4*S <= y.lo <= y.hi <=4*S):raise ArithmeticError('argument not in [-4,4]')
  return y

def sin_i(x):
  y=reduce(x); z=y*y;term=y; s=term
  for k in range(1,60):
    term= -term*z/((2*k)*(2*k+1))
    s+=term
  return s+I(-E_SIN.max(),E_SIN.max())

def cos_i(x):
  y=reduce(x);z=y*y;term=I(1);s=term
  for k in range(1,60):
    term=-term*z/((2*k-1)*(2*k))
    s+=term
  return s+I(-E_SIN.max(),E_SIN.max())

def overlap(a,i,j):
  b=2*L-a
  def directed(k,m):
    wk=k*PI/(2*L);wm=m*PI/(2*L)
    first_phase=wk*a
    if k==m: first=b*cos_i(first_phase)
    else: first=(sin_i((wk-wm)*b+first_phase)-sin_i(first_phase))/(wk-wm)
    second=(sin_i((wk+wm)*b+first_phase)-sin_i(first_phase))/(wk+wm)
    return (first-second)/(2*L)
  return (directed(i,j)+directed(j,i))/2

def exact_ldl_pivots(M):
  n=len(M)
  Ld=[[F(0) for _ in range(n)] for _ in range(n)]
  diag=[]
  for i in range(n):
    Ld[i][i]=F(1)
    pivot=M[i][i]-sum((Ld[i][k]**2*diag[k] for k in range(i)), F(0))
    diag.append(pivot)
    if pivot<=0:return diag,i
    for j in range(i+1,n):
      Ld[j][i]=(M[j][i]-sum((Ld[j][k]*Ld[i][k]*diag[k] for k in range(i)), F(0)))/pivot
  return diag,None


def check_reports(root):
  from quadratura_independente_fase3 import majorante_analitico_racional,produto_maximo_conservador
  LOG480=9*LN2-ln_via_atanh(31)  # 480=512*(15/16); ln(15/16)=-2atanh(1/31)
  h=F(1,10);rho=F(3,10);m=16;Ncells=1600
  prod=produto_maximo_conservador(m)
  fact=h*prod*(h/rho)**m
  erro=sum((majorante_analitico_racional(F(9,20),z*h+h/2,h,rho)*fact for z in range(Ncells)),F(0))
  pairs_expected={(i,j) for i in range(1,11) for j in range(i,11) if (i-j)%2==0}
  datas=[]
  for bits in (192,256):
    f=root/f'FASE3_G_N10_{bits}bits_30entradas.json'
    data=json.loads(f.read_text(encoding='utf-8'))
    mpath=root/f'MANIFESTO_FASE3_N10_{bits}bits.json'
    man=json.loads(mpath.read_text(encoding='utf-8'))
    assert man['result_file']==f.name
    assert man['bits']==bits and man['N']==10 and man['entradas_distintas']==30
    assert man['source_sha256']==hashlib.sha256((root/'quadratura_independente_fase3.py').read_bytes()).hexdigest()
    assert data['status']=='PASSOU' and len(data['entries'])==30
    assert {(e['i'],e['j']) for e in data['entries']}==pairs_expected
    for e in data['entries']:
      assert (e['L'],e['R'],e['precision_bits'],e['cells'],e['nodes_per_cell'],e['cell_width'],e['cauchy_radius'])==('9/20','160',bits,1600,16,'1/10','3/10')
      assert not e['uses_acb_integral'] and e['reference_consistent'] and e['reference_intersects_independent']
      assert F(e['interpolation_error_bound_exact'])==erro
      assert F(e['cauchy_error_factor'])==fact
      assert F(e['interpolation_product_bound'])==prod
      assert F(e['integral_low']) <= F(e['integral_high'])
      tail=I(4*e['i']*e['j'])*PI/(L*L*L)*(8+LOG480+I(F(1,3)))/(3*160**3)
      if F(e['tail_upper'])<tail.max():raise ArithmeticError(f'Incorrect tail bound {bits},{e["i"]},{e["j"]}')
    datas.append(data)
  intersections=0
  m0={(e['i'],e['j']):e for e in datas[0]['entries']}
  m1={(e['i'],e['j']):e for e in datas[1]['entries']}
  for k,a in m0.items():
    b=m1[k]
    alo,ahi=F(a['integral_low'])-F(a['tail_upper']),F(a['integral_high'])+F(a['tail_upper'])
    blo,bhi=F(b['integral_low'])-F(b['tail_upper']),F(b['integral_high'])+F(b['tail_upper'])
    if not (alo<=bhi and blo<=ahi):raise ArithmeticError(f'disjoint 192-256 {k}')
    intersections+=1
  return {'intervalos_integrais_conferidos':60,'cotas_de_cauda_validadas':60,
          'cotas_de_interpolacao_validadas':60,
          'pares_192x256_com_intersecao':intersections,
          'erro_total_interpolacao_exato':str(erro),
          'erro_total_interpolacao_decimal':float(erro),
          'largura_integral_192_max':float(max(F(e['integral_high'])-F(e['integral_low']) for e in datas[0]['entries']))}


def matrix(N=10,bits=192):
  root=Path(__file__).parent
  data=json.loads((root/f'FASE3_G_N10_{bits}bits_30entradas.json').read_text())
  pairs={(e['i'],e['j']) for e in data['entries']}
  expected={(i,j) for i in range(1,11) for j in range(i,11) if (i-j)%2==0}
  assert data['status']=='PASSOU' and pairs==expected and len(data['entries'])==30
  G=[[I(0) for _ in range(N)] for __ in range(N)]
  for e in data['entries']:
    i,j=e['i']-1,e['j']-1
    if max(i,j)>=N: continue
    assert e['precision_bits']==bits and (e['L'],e['R'])==('9/20','160')
    q=I(F(e['integral_low'])-F(e['tail_upper']),F(e['integral_high'])+F(e['tail_upper']))
    G[i][j]=G[j][i]=q
  P=[[I(0) for _ in range(N)] for __ in range(N)]
  factor=-2*LN2/SQ2
  for i in range(1,N+1):
    for j in range(i,N+1):
      # The symmetrized overlap of unlike-parity basis elements is 0.
      if (i-j)%2:continue
      v=factor*overlap(LN2,i,j)
      P[i-1][j-1]=P[j-1][i-1]=v
  B=[[I(0) for _ in range(N-2)] for __ in range(N)]
  p2=PI*PI
  for c,k in enumerate(range(3,N+1)):
    p=1 if k%2 else 2
    B[p-1][c]= -(I(F(k,p)) * (p*p*p2+L*L)/(k*k*p2+L*L))
    B[k-1][c]=I(1)
  W=[[G[i][j]+P[i][j] for j in range(N)] for i in range(N)]
  D=N-2
  WB=[[sum((W[i][k]*B[k][j] for k in range(N)),I(0)) for j in range(D)]for i in range(N)]
  A=[[sum((B[k][i]*WB[k][j] for k in range(N)),I(0)) for j in range(D)]for i in range(D)]
  for i in range(D):
    for j in range(i+1,D):
      z=I(min(A[i][j].min(),A[j][i].min()),max(A[i][j].max(),A[j][i].max()))
      A[i][j]=A[j][i]=z
  M=[[(A[i][j].min()+A[i][j].max())/2 for j in range(D)]for i in range(D)]
  delta=max(sum((a.width()/2 for a in row),F(0))for row in A)
  shift=[[M[i][j]-(delta if i==j else 0) for j in range(D)]for i in range(D)]
  pivots,failed=exact_ldl_pivots(shift)
  pmin=min(pivots)
  pivot_lower=F((pmin*10**6).__floor__(),10**6)
  # ceil(delta*1e9), done directly with rational arithmetic (no binary floats).
  delta_upper=F(-((-delta*10**9).__floor__()),10**9)
  assert delta<=delta_upper and pivot_lower<=pmin
  return dict(N=N,bits=bits,verified=failed is None,
      delta=float(delta),delta_upper_rational=str(delta_upper),
      min_pivot=float(pmin),min_pivot_lower_rational=str(pivot_lower),
      failed_at=None if failed is None else failed+1,diag_first=float(M[0][0]),
      max_matrix_width=float(max(a.width() for row in A for a in row)),
      pi_width=float(PI.width()),ln2_width=float(LN2.width()),
      Gmaxwidth=float(max(g.width()for row in G for g in row)))

if __name__=='__main__':
 root=Path(__file__).parent
 print('audit checks:',json.dumps(check_reports(root),ensure_ascii=False),flush=True)
 results=[]
 for bits in (192,256):
  for N in (6,8,10):
   ans=matrix(N,bits)
   results.append(ans)
   print(json.dumps(ans,ensure_ascii=False),flush=True)
 summary={'metodo':'quadratura polinomial + cota racional + projecao racional exata',
          'observacao':'Validacao condicional as bolas nodais Arb e a prova da cota analitica; nao demonstra RH',
          'source_sha256':hashlib.sha256((root/'quadratura_independente_fase3.py').read_bytes()).hexdigest(),
          'integrals':check_reports(root),'projected_certificates':results,
          'status':'PASSOU' if all(r['verified'] for r in results) else 'INCONCLUSIVO'}
 (root/'RESULTADO_AUDITORIA_PROJECAO_RACIONAL.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 assert summary['status']=='PASSOU'

