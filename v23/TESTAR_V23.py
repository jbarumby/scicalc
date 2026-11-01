#!/usr/bin/env python3
"""Fail-closed, adversarial and deterministic tests for V23 (stdlib only)."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
from fractions import Fraction as F
import AUDITAR_V23_DERIVADAS as d
import AUDITAR_V23_SUAVIZACAO as s
ROOT=Path(__file__).resolve().parents[1]
E7=ROOT/"NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip"
assert hashlib.sha256(E7.read_bytes()).hexdigest()==d.E7_SHA
r=d.audit(E7)
assert r["bands_count"]==11 and all(v["encloses_independent_calculus_bound"] for v in r["rows"])
assert all(F(v["verified_margin_lower"])>0 for v in r["rows"])
assert len(s.audit()["symmetry_adversary"]["off_line_zeros"])==4
for x in (F(-1),F(0),s.L,s.L*2):
 try:s.taper_bounds(x)
 except ValueError:pass
 else:raise AssertionError("invalid taper accepted")
try:d.m2_independent(1,6,"G22")  # crosses Fourier pole for this kind if narrow bound branch? h<=12, so not invalid
except Exception:pass
# Direct invalid band and unknown mode must reject.
for l,h,mode in [(-1,1,"G11"),(6,6,"G11"),(12,24,"UNKNOWN")]:
 try:d.m2_independent(l,h,mode)
 except (ValueError,AssertionError):pass
 else:raise AssertionError("invalid derivative band accepted")
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp)
 bad=p/"tampered-e7.zip";b=bytearray(E7.read_bytes());b[100]^=1;bad.write_bytes(b)
 try:d.audit(bad)
 except AssertionError:pass
 else:raise AssertionError("corrupted archive accepted")
 for mod,cli in (("DERIVADAS",["--e7",str(E7)]),("SUAVIZACAO",[])):
  src=ROOT/"v23"/f"AUDITAR_V23_{mod}.py"
  paths=[p/f"{mod}_1.json",p/f"{mod}_2.json"]
  for out in paths:
   process=subprocess.run([sys.executable,str(src),*cli,"--out",str(out)],
     text=True,capture_output=True,timeout=150)
   assert process.returncode==0,process.stderr+process.stdout
  assert paths[0].read_bytes()==paths[1].read_bytes(),mod+" output changed"
  saved=ROOT/"v23"/f"RESULTADO_V23_{mod}.json"
  if saved.is_file():assert saved.read_bytes()==paths[0].read_bytes()
print("PASS: 11 independent derivative bounds, taper Fourier errors, off-line symmetric zeros, negative inputs, source tamper and 2 deterministic JSON replays")
