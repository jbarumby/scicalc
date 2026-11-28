"""Exact algebraic identity and finite assembly for adaptive hybrid. No full numerical head claimed."""
import sympy as sp
import json

def assemble(A,B,H,C,D,lam,theta=sp.Integer(1)):
 """A: full head A0; B: leading compression sum k_i K_i; H: compression of QHQ on trailing block.
 All matrices expressed in orthonormal parity basis. Returns congruence without inverse.
 Requires 0<=H<=I, 0<lam<1, 0<=theta<=1, C,D>0. """
 n=A.rows;m=B.rows;r=n-m
 if A.shape!=(n,n) or B.shape!=(m,m) or H.shape!=(r,r):raise ValueError('dimension')
 if not(0<lam<1 and 0<=theta<=1 and C>0 and D>0):raise ValueError('parameters')
 Q=sp.diag(sp.zeros(m),sp.eye(r));BP=sp.diag(B,sp.zeros(r));s=lam*D/((1-lam)*C)
 J=sp.diag(sp.eye(m),(sp.eye(r)+s*H)/(1+s));R=J*(A-BP-((1-theta)*C+theta*C/lam)*Q)*J+theta*C/(lam*(1+s))*J*Q
 return R,J

def test():
 F=sp.Rational;A=sp.Matrix([[2,F(1,5),F(1,6)],[F(1,5),3,F(1,7)],[F(1,6),F(1,7),4]]);B=sp.Matrix([[F(1,10)]]);H=sp.Matrix([[F(1,3),F(1,7)],[F(1,7),F(2,3)]]);C=F(19,10);D=F(21,5);lam=F(91,100);s=lam*D/((1-lam)*C);Q=sp.diag(sp.zeros(1),sp.eye(2));Harm=C/lam*(sp.eye(2)-(sp.eye(2)+s*H).inv());BP=sp.diag(B,sp.zeros(2));rows=[]
 for t in [F(0),F(1,2),F(1)]:
  R,J=assemble(A,B,H,C,D,lam,t);M=A-BP-(1-t)*C*Q-t*sp.diag(sp.zeros(1),Harm)
  if R-J*M*J!=sp.zeros(3):raise ArithmeticError('congruence identity')
  rows.append({'theta':str(t),'exact_identity':True})
 print(json.dumps({'test':'exact rational 3x3 algebraic example, not a NUMERIA head','results':rows},indent=2))
if __name__=='__main__':test()
