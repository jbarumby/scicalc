# ES8-P96: certified small-t corrective remainder (P92 normalization)

Let b_h(x)=exp(-1/(1-(x/h)^2)) for |x|<h and 0 outside, h=2/5.
Let R(s)=integral b_h(y+s)b_h(y) dy and
F_delta(t)=[R(delta+t)+R(delta-t)]/2.

**Claim:** for 0<epsilon<1,

| integral_0^epsilon e^(t/2)(F_delta(t)-R(delta))/sinh(t) dt |
  <= epsilon^2 / h.

**Proof.** With u=x/h and z=1/(1-u^2)>=1,

|b'_h(x)| = (2|u|/h) z^2 exp(-z) <= 8/(e^2*h),

since z^2 exp(-z) attains its maximum at z=2. The bump is even
and strictly decreasing on (0,h), so ||b'_h||_1=2/e.
Integration by parts gives R''(s)=-integral b'_h(y+s)b'_h(y) dy.
Therefore ||R''||_inf <= 16/(e^3*h) < 2/h (since e>2).

Taylor's symmetric second-difference bound yields
|F_delta(t)-R(delta)| <= t^2 ||R''||_inf/2 <= t^2/h.
For 0<t<epsilon<1, sinh(t)>=t and exp(t/2)<2.
The corrective integrand is bounded in magnitude by 2t/h,
hence the integral is bounded by epsilon^2/h. QED.

With epsilon=1/65536 and h=2/5 the rational bound is
exactly 5/8589934592 ~= 5.820766091346741e-10.
This is 12 times sharper than the previous 12*epsilon^2/h.

Scope: only the small-t correction. No conclusion on the complete
nested integral, the sign of W12, or the Riemann hypothesis.
