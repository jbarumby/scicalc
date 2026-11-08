from flint import acb,arb,ctx
import flint,json
ctx.prec=256
h=arb(2)/5
def b(z,analytic):
    re=z.real
    if bool(re.lower()>-h) and bool(re.upper()<h):
        return (1-acb(1)/(1-(z/h)**2)).exp()
    if bool(re.upper() < -h) or bool(re.lower()>h):
        return acb(0)
    if analytic:return acb(arb("nan"))
    if not bool(z.imag==0):return acb(arb("nan"))
    return acb(arb("0.5 +/- 0.5"))
for c in (arb(5)/4,arb(401)*arb(5)/4):
    for k in (0,1,14):
        def f(z,analytic):
            return b(z,analytic)*z**k*(acb(0,1)*c*z).exp()
        out=acb.integral(f,-h,h,abs_tol=arb("1e-13"),rel_tol=arb("1e-24"),
                         deg_limit=180,eval_limit=300000,depth_limit=90)
        if not out.is_finite() or not bool(out.real.rad().upper()<arb("2e-12")):
            raise ArithmeticError(f"pilot failed center={c} k={k}: {out}")
        print("PILOT_OK",str(c),k,str(out),flush=True)
print("FLINT_FOURIER_PILOT_PASS",flint.__version__,ctx.prec)
