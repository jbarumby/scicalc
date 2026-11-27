"""Independent integer replay of adaptive weighted row bounds (no optimizer)."""
from fractions import Fraction as F
from math import lcm
import csv,json
from collections import Counter

def need(ok,msg):
 if not ok: raise ArithmeticError(msg)
with open('deficit_2000000.csv') as f:
 rows=[tuple(map(int,r)) for r in list(csv.reader(f))[1:]]
parts=[];prev=-1
for l,r,d in rows:
 need(r==l+1 and l>=prev and d>0,'cell order/shape');prev=r
 if l<4500*16:continue
 if parts and parts[-1][1]==l:parts[-1][1]=r;parts[-1][2]+=d;parts[-1][3]=max(parts[-1][3],d)
 else:parts.append([l,r,d,d])
x=json.load(open('hybrid39.json'));ks=x['ks'];w=x['weights'];need(len(parts)==len(w)==len(ks),'component lengths');need(all(k==(14 if m>=16000000 else 2) for k,(_,_,m,_) in zip(ks,parts)),'factor rule');need(min(w)>0,'weights')
scale=2**20;den=lcm(*[k-1 for k in ks]);cap=(1387*scale+999)//1000
v=[m*wj*k*(den//(k-1)) for (_,_,m,_),wj,k in zip(parts,w,ks)]
mx=F(0)
for i,(l,r,_,_) in enumerate(parts):
 tot=0
 for (L,R,_,_),vv in zip(parts,v):
  dist=max(0,L-r,l-R);kernel=cap if not dist else min(cap,(16*scale+dist-1)//dist)
  kernel+=(16*scale+l+L-1)//(l+L);tot+=kernel*vv
 mx=max(mx,F(tot*106,den*scale*16000000*333*w[i]))
d=max(F(d*k,10**6*(k-1)) for (_,_,_,d),k in zip(parts,ks));need(mx==F(x['c']) and d==F(x['d']),'bounds replay')
lam=F(x['lambda']);s=lam*d/((1-lam)*mx);f1=mx/lam*s/(1+s);need(s==F(x['s']) and f1==F(x['f1']) and f1<2,'harmonic budget')
# New Binet upper: |remainder|<=1/(3t), Re log z<=log(t/2)+1/(8t^2).
need(F(9944,1000)+F(1,3*4500)+F(1,8*4500**2)-2<8,'Binet constant')
print(json.dumps({'components':len(parts),'factors':dict(Counter(ks)),'c_exact':str(mx),'d_exact':str(d),'f1_exact':str(f1),'weighted_row_checks':len(parts),'kernel_terms':len(parts)**2,'Binet_constant_verified':True,'scope':'replay of inherited deficit table, not independent regeneration'},indent=2))
