#!/usr/bin/env python3
"""Exact V24 perturbation tests: domain, sensitivity and output provenance."""
import hashlib,json,subprocess,sys,tempfile
from fractions import Fraction as F
from pathlib import Path
import AUDITAR_V24_CONTINUIDADE as v
ROOT=Path(__file__).resolve().parents[1]
EXPECTED_V23={
 "AUDITAR_V23_DERIVADAS.py":"76761b42c3f95142665d92718dba5a7c8d0b49ccfe2299f7e070a23ebba0cb75",
 "AUDITAR_V23_SUAVIZACAO.py":"4f6aea9747b6c87061e343322d54b6e4ed0a23b3278be5f7fb99a517acc08383"
}
for name,expected in EXPECTED_V23.items():
 actual=hashlib.sha256((ROOT/"v23"/name).read_bytes()).hexdigest()
 assert actual==expected, (name,actual)
assert v.LOG2_HI<2*v.L<v.log_enclosure(3)[0]
assert v.STEP_D1==8 and v.STEP_D2==128
assert v.ARCH_TAIL_INTEGRAL==16 and v.ARCH_SMALL_WEIGHT==12

for k,delta in ((1,F(1,1000)),(2,F(1,10000)),(10,F(1,1000000))):
 p=v.parameters(k,delta)
 assert p["l1"]>0 and p["Dtaper"]>p["Dbase"]>0
 a=v.components(p,p)
 assert a["arch_0_1"]>0 and a["arch_1_infty"]>0
 assert a["total"]==a["arch_total"]+a["poles_s0_s1"]+a["prime_two"]
 assert a["arch_total"]==a["arch_0_1"]+a["arch_1_infty"]
 b=v.parameters(k,delta/10)
 assert p["l1"]==100*b["l1"]
 assert a["total"]==100*v.components(b,b)["total"]
# Asymmetric cross-term still gives symmetric entry bound.
a=v.parameters(2,F(1,10000));b=v.parameters(9,F(1,10000))
assert v.components(a,b)==v.components(b,a)
# Failure conditions and proof against inadvertent broadening of domain.
for k,d in ((0,F(1,100)),(-1,F(1,100)),(1,F(0)),(1,-F(1,50)),(1,v.L),(1,v.L+F(1,100))):
 try:v.parameters(k,d)
 except ValueError:pass
 else:raise AssertionError("invalid taper/mode was accepted")
# Omitting derivative transition terms invalidates the 1/r² taper envelope.
p=v.parameters(2,F(1,100000))
assert p["Dtaper"] > 2*v.L*(F(2)*v.PI_HI/(2*v.L))**2/v.SQLO
# No deterministic JSON may silently change between runs.
with tempfile.TemporaryDirectory() as d:
 paths=[Path(d)/"first.json",Path(d)/"second.json"]
 for p in paths:
  run=subprocess.run([sys.executable,str(Path(__file__).with_name("AUDITAR_V24_CONTINUIDADE.py")),"--out",str(p)],
      capture_output=True,text=True,timeout=60)
  assert run.returncode==0,run.stderr+run.stdout
 assert paths[0].read_bytes()==paths[1].read_bytes()
 saved=ROOT/"v24/RESULTADO_V24_CONTINUIDADE.json"
 if saved.is_file():assert saved.read_bytes()==paths[0].read_bytes(),"committed witness changed"
 j=json.loads(paths[0].read_text())
 assert j["schema"]=="NUMERIA.E8K.V24.GLOBAL_TAPER_PERTURBATION.1"
 assert len(j["test_deltas"])==4 and j["modes"]==10
 assert j["status"].endswith("CONDITIONAL")
 assert "Riemann Hypothesis" in " ".join(j["not_certified"])
print("PASS: V24 global r-range, exact prime/pole errors, symmetric 10-mode perturbations, provenance and two deterministic replays")
