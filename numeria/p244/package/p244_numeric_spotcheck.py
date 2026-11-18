#!/usr/bin/env python3
"""Independent float64 GL60 spot-check. Diagnostic only, NEVER a certificate."""
import gzip,json,math
from pathlib import Path
import numpy as np
from scipy.special import digamma,spherical_jn
from numpy.polynomial.legendre import leggauss

L=161/200; T=504
pvals=np.arange(1008,dtype=np.float64)
x,w=leggauss(60)
t=((2*pvals+1)/4)[:,None]+x[None,:]/4
weights=np.broadcast_to(w[None,:]/4,t.shape)
A=2*np.log(2)/np.sqrt(2)+2*np.log(3)/np.sqrt(3)+np.log(2)+2*np.log(5)/np.sqrt(5)
beta=np.log(T/(2*np.pi))-1/T-A
psi=np.real(digamma(0.25+0.5j*t))-np.log(np.pi)
for a,log_n in [(2*np.log(2)/np.sqrt(2),np.log(2)),(2*np.log(3)/np.sqrt(3),np.log(3)),(np.log(2),2*np.log(2)),(2*np.log(5)/np.sqrt(5),np.log(5))]:
    psi-=a*np.cos(t*log_n)
common=(weights/np.pi)*(psi-beta)
out={'diagnostic_only_float64':True,'beta_float64':beta,'sectors':{}}
for parity,fn in [('even','P237_FULL_C_even.json.gz'),('odd','P237_FULL_C_odd.json.gz')]:
    data=json.loads(gzip.decompress((Path('/mnt/data')/fn).read_bytes()))
    rows=data['C_quadrature_outer_upper_triangle']
    modes=[parity=='odd',2+(parity=='odd'),4+(parity=='odd'),6+(parity=='odd')]
    U=[]
    for i,n in enumerate(modes):
        U.append((-1)**i*np.sqrt(2*L*(2*n+1))*spherical_jn(n,L*t))
    checks=[]
    for i,j in [(0,0),(0,1),(0,2),(1,1),(1,2),(3,3)]:
        # index in row-major upper triangle 377x377
        k=i*377-i*(i-1)//2+(j-i)
        ii,jj,b=rows[k]
        assert [ii,jj]==[i,j]
        mid=int(b[0][0])*(2.0**int(b[0][1])) if abs(int(b[0][1]))<1000 else math.ldexp(float(b[0][0]),int(b[0][1]))
        num=np.sum(U[i]*U[j]*common)
        checks.append({'indices':[i,j],'degrees':[modes[i],modes[j]],'float64_gl60':float(num),'p237_C_midpoint_approx':float(mid),'abs_diff':float(abs(num-mid))})
    out['sectors'][parity]=checks
Path('/mnt/data/NUMERIA_ES8P244_AUDITORIA/P244_NUMERICAL_SPOTCHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print('P244 independent scipy float64 GL60 diagnostic:')
for sec,checks in out['sectors'].items():
    print(sec,[(v['indices'],f"{v['abs_diff']:.3e}") for v in checks])
