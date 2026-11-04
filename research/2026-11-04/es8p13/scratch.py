from fractions import Fraction as F
from decimal import Decimal
from AUDITAR_FASE3_PROJECAO import I,PI,LN2,SQ2,L,overlap
from es8p6_independent_audit import sharp_tail_upper
from es8p10_lemmas import new_intervals

def matrix(N,records):
 G=[[I(0) for _ in range(N)]for _ in range(N)];P=[[I(0)for _ in range(N)]for _ in range(N)];
 for i in range(1,N+1):
  for j in range(i,N+1):
   if (i-j)%2:continue
   rec=records[f'{i},{j}'];lo,hi=map(lambda s:F(Decimal(s)),rec['independent_finite_with_cauchy']);t=sharp_tail_upper(i,j)
   G[i-1][j-1]=G[j-1][i-1]=I(lo-t,hi+t);P[i-1][j-1]=P[j-1][i-1]=(-2*LN2/SQ2)*overlap(LN2,i,j)
 B=[[I(0) for _ in range(N-2)]for _ in range(N)];pi2=PI*PI
 for col,k in enumerate(range(3,N+1)):
  p=1 if k%2 else 2
  B[p-1][col]=-(I(F(k,p))*(p*p*pi2+L*L)/(k*k*pi2+L*L));B[k-1][col]=I(1)
 W=[[G[i][j]+P[i][j]for j in range(N)]for i in range(N)];d=N-2
 WB=[[sum((W[i][k]*B[k][j] for k in range(N)),I(0))for j in range(d)]for i in range(N)]
 A=[[sum((B[k][i]*WB[k][j]for k in range(N)),I(0))for j in range(d)]for i in range(d)]
 for i in range(d):
  for j in range(i+1,d):
   q=I(min(A[i][j].min(),A[j][i].min()),max(A[i][j].max(),A[j][i].max()));A[i][j]=A[j][i]=q
 return A

if __name__=='__main__':
 data=new_intervals()
 for n in(6,8,10):
  A=matrix(n,data);d=n-2
  print('N',n)
  for i in range(d):
   mid=A[i][i].mid();rest=sum(abs(A[i][j].mid())for j in range(d)if j!=i)
   print(' '.join('%+.5f'%float(A[i][j].mid()) for j in range(d)), 'diagonal dominance minus .05',float(mid-rest-F(1,20)))
