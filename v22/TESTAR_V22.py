#!/usr/bin/env python3
"""V22 exact, adversarial and deterministic checks; uses stdlib only."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
from fractions import Fraction as F
import AUDITAR_V22_QUADRATURA as quadrature
import AUDITAR_V22_FROBENIUS as frob
import AUDITAR_V22_SUPORTE as sup
import AUDITAR_V22_SCREW as screw
ROOT=Path(__file__).resolve().parents[1]
E7=ROOT/"NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip"
assert hashlib.sha256(E7.read_bytes()).hexdigest()==quadrature.E7_SHA
q=quadrature.audit(E7)
assert q["b5_records"]==15754
assert q["components"]["G11"]["historical_underbounds"]==0
assert q["components"]["G22"]["historical_underbounds"]==1
assert q["components"]["G22"]["correction"]==str(F(1,10**45))
assert sup.audit()["sharp_sine_C_infinity"]==False
assert frob.audit()["points"]=={"F5":8,"F25":32}
assert screw.audit()["principal_minors"]==31
with tempfile.TemporaryDirectory() as d:
 temp=Path(d)
 bad=temp/"bad-e7.zip";original=bytearray(E7.read_bytes());original[123]^=1;bad.write_bytes(original)
 try:quadrature.audit(bad)
 except AssertionError:pass
 else:raise AssertionError("Tampered E7 accepted")
 for target in ("QUADRATURA","FROBENIUS","SUPORTE","SCREW"):
  s=ROOT/"v22"/("AUDITAR_V22_"+target+".py")
  p1,p2=temp/(target+"1.json"),temp/(target+"2.json")
  for out in (p1,p2):
   args=[sys.executable,str(s)]
   if target=="QUADRATURA":args+=["--e7",str(E7)]
   args+=["--out",str(out)]
   r=subprocess.run(args,text=True,capture_output=True,timeout=60)
   assert r.returncode==0,r.stdout+r.stderr
  assert p1.read_bytes()==p2.read_bytes(),target+" non-deterministic"
  saved=ROOT/"v22"/("RESULTADO_V22_"+target+".json")
  if saved.is_file():assert saved.read_bytes()==p1.read_bytes(),target+" saved result differs"
print("PASS: V22 11 B3 bands, 15754 B5 records, G22 one-tick fix, F5/F25, Burnol domain, 31 toy screw minors, tampered E7 and 4 deterministic JSON pairs")
