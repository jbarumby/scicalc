#!/usr/bin/env python3
"""P104: independently exercise Arb interval-parametric nested integration.

These tests DO NOT prove the full P92 formula. They test exact primitives,
real parameter balls and interval-valued integration endpoints at 192 bits.
"""
import json, platform, unittest
from fractions import Fraction as F
from importlib import metadata
from pathlib import Path
import flint
from flint import arb, acb, ctx
ctx.prec=192

def exact(q):
    f=q.fmpq()
    return F(int(f.p), int(f.q))

def bounds(z):
    re=z.real if isinstance(z,acb) else z
    return exact(re.lower()), exact(re.upper())

def contains(z,low,high):
    l,h=bounds(z)
    return l <= low and h >= high

S=acb(arb("0.2 +/- 0.01"))
A=acb(arb("0 +/- 0.01"))
B=acb(arb("1 +/- 0.01"))

class ProofObligations(unittest.TestCase):
    def test_parameter_ball_polynomial_integral(self):
        # integral_0^1 (y+s)^2 dy = 1/3+s+s^2; increasing in 0.19 <= s <= 0.21
        out=acb.integral(lambda y,analytic: (y+S)**2,acb(0),acb(1),
                         abs_tol=arb("1e-18"),rel_tol=arb("1e-18"),
                         eval_limit=8000,depth_limit=36)
        minimum=F(1,3)+F(19,100)+F(19,100)**2
        maximum=F(1,3)+F(21,100)+F(21,100)**2
        self.assertTrue(contains(out,minimum,maximum),str(out))

    def test_interval_endpoints_constant(self):
        # any a in [-0.01,0.01], b in [0.99,1.01]
        out=acb.integral(lambda y,analytic: acb(1),A,B,
                         abs_tol=arb("1e-18"),rel_tol=arb("1e-18"),
                         eval_limit=8000,depth_limit=36)
        self.assertTrue(contains(out,F(98,100),F(102,100)),str(out))

    def test_interval_endpoints_and_parameter(self):
        # integrate (y+s)dy = (b-a)*s + (b^2-a^2)/2
        # monotonically increasing in s and b, decreasing in a
        out=acb.integral(lambda y,analytic: y+S,A,B,
                         abs_tol=arb("1e-18"),rel_tol=arb("1e-18"),
                         eval_limit=8000,depth_limit=36)
        def primitive(a,b,s):
            return (b*b-a*a)/2+s*(b-a)
        minimum=primitive(F(1,100), F(99,100), F(19,100))
        maximum=primitive(F(-1,100), F(101,100), F(21,100))
        self.assertTrue(contains(out,minimum,maximum),str(out))

    def test_analytic_denial_keeps_valid_direct_enclosure(self):
        # No analytic extension is supplied for absolute value; direct enclosures only.
        def g(y,analytic):
            if analytic:
                return acb(arb.nan(),arb.nan())
            return acb(abs(y.real))
        out=acb.integral(g, acb(-1), acb(1),
                         abs_tol=arb("1e-7"),rel_tol=arb("1e-7"),
                         eval_limit=15000,depth_limit=44)
        self.assertTrue(contains(out,F(1),F(1)),str(out))

if __name__=="__main__":
    path=Path("P104_ENVIRONMENT.json")
    version_keys={k:str(getattr(flint,k)) for k in dir(flint)
                  if "version" in k.lower() and len(k)<60}
    payload={"python_flint_distribution": metadata.version("python-flint"),
             "python":platform.python_version(),
             "platform":platform.platform(),
             "flint_version_attributes":version_keys,
             "note":"FLINT linked C library version may differ from the source revision cited in the audit."}
    path.write_text(json.dumps(payload,indent=2)+"\n")
    print("ENVIRONMENT",json.dumps(payload),flush=True)
    unittest.main(verbosity=2)
