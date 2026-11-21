"""Verbatim decimal inputs from Chourasiya--Simonic, arXiv:2507.15184v2, Table 1.

All decimal strings represent *exact decimal rationals* in the input.
Do not silently treat an input printed decimal as an outward-rounded bound.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class TableRow:
    sigma_start: int  # numerator / 32
    b1: str
    b2: str
    b3: str
    d: str
    kappa: str
    delta: str
    log10_a0: str


_RAW = (
    ("5.069", "9.461", "167.8", "0.28", "0.34", "0.012", "11.25"),
    ("6.046", "9.168", "162.6", "0.28", "0.32", "0.025", "11.26"),
    ("7.077", "8.875", "157.4", "0.28", "0.30", "0.038", "11.27"),
    ("8.185", "8.582", "152.2", "0.28", "0.29", "0.052", "11.27"),
    ("9.387", "8.290", "147.0", "0.28", "0.27", "0.066", "11.28"),
    ("10.70", "7.997", "141.9", "0.28", "0.25", "0.081", "11.29"),
    ("12.11", "7.704", "136.7", "0.28", "0.23", "0.097", "11.30"),
    ("13.66", "7.411", "131.5", "0.28", "0.21", "0.114", "11.31"),
    ("15.35", "7.118", "126.3", "0.28", "0.20", "0.132", "11.32"),
    ("17.21", "6.826", "121.1", "0.28", "0.18", "0.150", "11.33"),
    ("19.31", "6.533", "115.9", "0.28", "0.16", "0.170", "11.35"),
    ("21.77", "6.240", "110.7", "0.28", "0.14", "0.192", "11.36"),
    ("24.83", "5.745", "101.9", "0.29", "0.13", "0.215", "11.37"),
    ("29.01", "5.463", "96.87", "0.29", "0.11", "0.241", "11.40"),
    ("35.35", "5.180", "91.86", "0.29", "0.09", "0.270", "11.42"),
    ("46.06", "4.897", "86.84", "0.29", "0.07", "0.303", "11.46"),
)
TABLE = tuple(TableRow(16 + i, *data) for i, data in enumerate(_RAW))


def validate_table():
    """Fail closed if transcription, interval count or endpoint conventions change."""
    from fractions import Fraction
    if len(TABLE) != 16 or [r.sigma_start for r in TABLE] != list(range(16, 32)):
        raise ArithmeticError("Table 1 must contain the 16 consecutive intervals n/32")
    if not all(Fraction(r.d) == Fraction(28,100) for r in TABLE[:12]):
        raise ArithmeticError("expected d=0.28 in the first twelve rows")
    if not all(Fraction(r.d) == Fraction(29,100) for r in TABLE[12:]):
        raise ArithmeticError("expected d=0.29 in the last four rows")
    if not all(Fraction(r.b1) > 0 and Fraction(r.b2) > 0 and Fraction(r.b3) > 0
               and Fraction(r.kappa) > 0 and Fraction(r.delta) > 0
               for r in TABLE):
        raise ArithmeticError("published coefficients and parameters must be positive")
    return True


validate_table()
