"""Verify the inherited prime-deficit enclosure, including all omitted cells."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from struct import Struct

from flint import acb, arb, ctx


START = 4500
TAIL_START = 10**6
SCALE = 16 * 2**14
SUBCELLS = SCALE // 16
BETA = 2
FREQUENCIES = ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7),
               (8, 2), (9, 3), (11, 11), (13, 13))
CSV_PATH = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39" / "deficit_2000000.csv"
CSV_SHA256 = "2c9111dfbe2c77d2d3775c9d4187e658a6c784b3e416c417b3c59dd907dfb8cc"
TUPLE = Struct(">QQI")


def dyadic(value):
    mantissa, exponent = value.man_exp()
    return [str(mantissa), int(exponent)]


class NineTermSymbol:
    def __init__(self):
        self.twice_pi = 2 * arb.pi()
        self.log_pi = arb.pi().log()
        self.frequencies = tuple((arb(q).log(), 2 * arb(p).log() / arb(q).sqrt())
                                 for q, p in FREQUENCIES)
        self.mass = sum((c for _, c in self.frequencies), arb(0))
        self.curvature = sum((c * w**2 for w, c in self.frequencies), arb(0)).upper()
        self.slopes = tuple((w, c, c*w) for w, c in self.frequencies)
        self.largest_first = tuple(sorted(self.frequencies, key=lambda z: z[1], reverse=True))

    def arch_lower(self, left):
        """Expanded Binet floor, increasing for t >= 4500."""
        return ((left / self.twice_pi).log() - 1 / (2*left**2) - 2 / (15*left**3)).lower()

    def cell_margin(self, left, right, cap, choose_direct):
        """Lower bound on Psi_9(t) - beta + cap for the entire closed cell."""
        mid, rad = (left + right) / 2, (right - left) / 2
        base = self.arch_lower(left) - BETA + cap
        if choose_direct:
            # Ignoring each not-yet-evaluated term uses its exact maximum 1.
            gain = base - self.mass
            interval = mid + arb(0, rad)
            one = arb(1)
            for frequency, mass in self.largest_first:
                gain += mass * (one - min(one, (frequency * interval).cos().upper()))
                if gain > 0:
                    return gain.lower()
        value, slope = arb(0), arb(0)
        for frequency, amplitude, derivative_amplitude in self.slopes:
            sine, cosine = (frequency * mid).sin_cos()
            value += amplitude * cosine
            slope += derivative_amplitude * sine
        return (base - value - abs(slope) * rad - self.curvature * rad**2 / 2).lower()

    def point_margin(self, t, cap):
        arch = acb(arb(1) / 4, t / 2).digamma().real - self.log_pi
        trig = sum((c * (w * t).cos() for w, c in self.frequencies), arb(0))
        return arch - trig - BETA + cap


def read_rows(path):
    if sha256(path.read_bytes()).hexdigest() != CSV_SHA256:
        raise ValueError("inherited deficit CSV SHA-256 mismatch")
    previous = 0
    with path.open(newline="", encoding="ascii") as stream:
        rows = csv.reader(stream)
        if next(rows) != ["left16", "right16", "deficit_micro"]:
            raise ValueError("unexpected deficit CSV columns")
        for lft, rgt, cap in rows:
            lft, rgt, cap = map(int, (lft, rgt, cap))
            if not (0 <= previous <= lft and rgt == lft + 1 and 0 < cap < 2**32):
                raise ValueError("invalid deficit CSV support")
            previous = rgt
            yield lft, rgt, cap


def source_cells(source, start, stop, omitted_width=32):
    """Cover [start,stop] once, applying zero deficit on omitted cells."""
    cursor, limit = start * 16, stop * 16
    for left, right, micro in read_rows(source):
        if right <= cursor:
            continue
        if left >= limit:
            break
        while cursor < left:
            end = min(left, cursor + omitted_width)
            yield cursor * SUBCELLS, end * SUBCELLS, 0
            cursor = end
        yield left * SUBCELLS, right * SUBCELLS, micro
        cursor = right
    while cursor < limit:
        end = min(limit, cursor + omitted_width)
        yield cursor * SUBCELLS, end * SUBCELLS, 0
        cursor = end


def certify(source=CSV_PATH, *, start=START, stop=TAIL_START, precision=128):
    if type(start) is not int or type(stop) is not int or not START <= start < stop <= TAIL_START:
        raise ValueError("expected integer START <= start < stop <= TAIL_START")
    if type(precision) is not int or precision < 80:
        raise ValueError("use integer Arb precision of at least 80 bits")
    digest = sha256()
    table_cells = zero_chunks = table_leaves = zero_leaves = subdivisions = 0
    min_table = min_zero = None
    with ctx.workprec(precision):
        sym = NineTermSymbol()
        step = arb(1) / SCALE
        final_floor = (sym.arch_lower(arb(TAIL_START)) - BETA - sym.mass).lower()
        if not final_floor > 0:
            raise ArithmeticError("the high-frequency Binet floor did not close")

        def check(lft, rgt, micro, cap):
            nonlocal table_leaves, zero_leaves, subdivisions, min_table, min_zero
            margin = sym.cell_margin(lft*step, rgt*step, cap, micro == 0 and lft >= 50000*SCALE)
            if margin > 0:
                digest.update(TUPLE.pack(lft, rgt, micro))
                if micro:
                    table_leaves += 1
                    if min_table is None or margin < min_table[0]:
                        min_table = margin, lft, rgt, micro
                else:
                    zero_leaves += 1
                    if min_zero is None or margin < min_zero[0]:
                        min_zero = margin, lft, rgt, micro
                return None
            mid = (lft + rgt) // 2
            if mid == lft:
                point = sym.point_margin((lft + rgt)*step/2, cap)
                return {
                    "cell": [[lft, SCALE], [rgt, SCALE]],
                    "cap_micro": micro,
                    "cell_margin_lower_dyadic": dyadic(margin),
                    "point_margin_upper_dyadic": dyadic(point.upper()),
                    "point_proves_invalid_cap": bool(point < 0),
                }
            subdivisions += 1
            return check(lft, mid, micro, cap) or check(mid, rgt, micro, cap)

        unresolved = None
        for lft, rgt, micro in source_cells(source, start, stop):
            if micro:
                table_cells += 1
            else:
                zero_chunks += 1
            unresolved = check(lft, rgt, micro, arb(micro) / 10**6)
            if unresolved is not None:
                break

        def summary(record):
            if record is None:
                return None
            bound, left, right, micro = record
            return {"margin_lower_dyadic": dyadic(bound), "cell": [[left, SCALE], [right, SCALE]],
                    "cap_micro": micro}

        return {
            "schema": "NUMERIA.P93V40.NINE_TERM_DEFICIT_REPLAY.1",
            "source_sha256": CSV_SHA256,
            "precision_bits": precision,
            "start": start,
            "stop": stop,
            "status": "certified" if unresolved is None else "inconclusive",
            "scope_has_no_gap_above_4500": unresolved is None and start == START and stop == TAIL_START,
            "ordered_leaf_partition_sha256": digest.hexdigest(),
            "listed_cells_attempted": table_cells,
            "omitted_chunks_attempted": zero_chunks,
            "listed_accepted_leaves": table_leaves,
            "omitted_accepted_leaves": zero_leaves,
            "subdivisions": subdivisions,
            "minimum_listed": summary(min_table),
            "minimum_omitted": summary(min_zero),
            "first_unresolved": unresolved,
            "total_amplitude_upper_dyadic": dyadic(sym.mass.upper()),
            "binet_tail_floor_margin_lower_dyadic": dyadic(final_floor),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    parser.add_argument("--start", type=int, default=START)
    parser.add_argument("--stop", type=int, default=TAIL_START)
    parser.add_argument("--precision", type=int, default=128)
    args = parser.parse_args()
    if args.verify is not None:
        expected = json.loads(args.verify.read_text())
        result = certify(start=expected["start"], stop=expected["stop"],
                         precision=expected["precision_bits"])
        if json.dumps(expected, sort_keys=True) != json.dumps(result, sort_keys=True):
            raise ValueError("spectral replay does not match the supplied certificate")
        print(f"Verified {args.verify}")
    else:
        result = certify(start=args.start, stop=args.stop, precision=args.precision)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({"status": result["status"], "first_unresolved": result["first_unresolved"]}))
    if result["status"] != "certified":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
