#!/usr/bin/env python3
"""Compare the complete regenerated rational certificate byte-for-byte."""
from pathlib import Path
p=Path(__file__).resolve().parent
assert (p/'RESULTADOS_ES8P10_REPLAY.json').read_bytes()==(p/'RESULTADOS_ES8P10.json').read_bytes()
print('ES8-P10 REPLAY: byte-for-byte identical exact rational output')
