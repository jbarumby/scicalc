# NUMERIA ES8-P102 — analytical audit of P92 2x2 (research only)

**Reference:** Bombieri–Lagarias, *Complements to Li's Criterion for the Riemann Hypothesis*, section 3, "Explicit Formula", manuscript p.6:
https://websites.umich.edu/~lagarias/doc/bombieri.pdf

The calculations below examine only the 2-centre test subspace,
**not** universal Weil positivity or RH.

## Test family and admissibility

h=2/5; b(u)=exp(-1/(1-(u/h)^2)) for |u|<h, zero outside;
R(s)=integral_R b(y+s)b(y) dy is even and compactly supported
in [-4/5,4/5]; delta=log(j/i), i,j in {1,2}.
Define H_delta(t)=R(t+delta), F_delta(t)=[R(delta+t)+R(delta-t)]/2,
f_delta(x)=x^(-1/2) H_delta(log x) for x>0.
The f_delta are smooth with support compact inside (0,infty), hence
satisfy the regularity/decay hypotheses of the published explicit formula.

## From Bombieri–Lagarias to the implementation

The published arithmetic formula is

 W(f) = integral_0^infty (f(x)+f_tilde(x)) dx
      - sum_{n>=2} Lambda(n) [f(n)+f_tilde(n)]
      - (log(pi)+EulerGamma) f(1)
      - integral_1^infty [f(x)+f_tilde(x)-2 f(1)/x^2] x/(x^2-1) dx,

where f_tilde(x)=x^(-1) f(1/x).

Set x=e^t. Then f_delta(e^t)=e^(-t/2) H_delta(t),
f_tilde_delta(e^t)=e^(-t/2) H_delta(-t), f(1)=R(delta).
The final integral becomes

 I_BL(delta)=integral_0^infty
 [e^(t/2)*F_delta(t) - e^(-t)*R(delta)]/sinh(t) dt.

The code uses instead

 I_P92(delta)=integral_0^infty
 [e^(t/2)*F_delta(t) - R(delta)]/sinh(t) dt.

The identity

 integral_0^infty (1-e^(-t))/sinh(t) dt
  = integral_0^infty 2/(e^t+1) dt = log(4)

gives I_BL=I_P92+log(4)*R(delta), proving EXACTLY

 -(log(pi)+gamma)R(delta)-I_BL
  = -(log(4*pi)+gamma)R(delta)-I_P92.

Thus the apparently different archimedean constants agree.

## Poles, primes, support

A=integral_R b(y)e^(y/2)dy. Because b is even,

  P_ij=A^2(e^(delta/2)+e^(-delta/2)).

The prime-power term is

 S_ij=sum_{n>=2} Lambda(n)/sqrt(n)
  [R(log(n)+delta)+R(log(n)-delta)].

For delta in {0,log2}, only n in {2,3,4} can contribute:
for n>=5, log(n/2)>=log(5/2)>6/7>4/5, where the
log lower bound follows from log(u)=2atanh((u-1)/(u+1))
>2(u-1)/(u+1) for u>1.
Lambda(4)=log2. This matches the hard-coded list.

## Origin and tail

K(t)=(exp(t/2)-1)/sinh(t), with y=exp(-t/2), has
K(t)=2y/[(1+y)(1+y^2)] and K(0)=1/2.

||b'||_1=2/e, ||b'||_infty<=8/(e^2*h),
||R''||_infty<=16/(e^3*h)<2/h.
The symmetric Taylor difference therefore satisfies
|F_delta(t)-R(delta)|<=t^2/h, so for 0<epsilon<=1,

 | integral_0^epsilon exp(t/2)*(F_delta(t)-R(delta))/sinh(t) dt |
 <= epsilon^2/h.

For P101 epsilon=1/2048: upper bound 5/8388608.
For the user's separate P92 computation epsilon=1/65536:
upper bound 5/8589934592.

For T=2h+|delta|, F_delta(t)=0 for t>=T, and

  integral_T^infty -R(delta)/sinh(t) dt
    = R(delta)*log(tanh(T/2)),

matching the code's compact-support tail.

## Status and outstanding obligation

**Derived/checked:** the published explicit identity, constant,
poles, prime-power support, regularization, origin error and tail.

**Still requires independent review:** the interval extension of
every nested callback for all real balls passed by acb.integral,
including balls crossing support boundaries; scope of analytic=True
and outward rounding of constructed interval balls. The smoke tests
in this PR are NOT a mathematical proof of all nested enclosures.

The exact rational determinant audit was successful for ES8-P101
at 128 and 256 bits, conditional on these numerical premises:
https://github.com/jbarumby/scicalc/actions/runs/42449573849

DO NOT merge to main without independent review. No RH claim.
