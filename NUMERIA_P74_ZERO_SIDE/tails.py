"""Bound the potentially off-line part of Suzuki's zero expansion."""

from flint import arb

from NUMERIA_ES8P259.suzuki import _dyadic


HEIGHT = 3 * 10**12

# Kadiri–Lumley–Ng, Table 1. The earlier rows are redundant at this
# height after bounding D(sigma) by half the total tail mass.
KLN_ROWS = (
    (80, "6.918", "3.997"),
    (85, "8.975", "3.588"),
    (86, "9.441", "3.514"),
    (87, "9.926", "3.430"),
    (88, "10.431", "3.346"),
    (89, "10.955", "3.262"),
    (90, "11.499", "3.186"),
    (91, "12.063", "3.102"),
    (92, "12.646", "3.017"),
    (93, "13.250", "2.941"),
    (94, "13.872", "2.856"),
    (95, "14.513", "2.772"),
    (96, "15.173", "2.694"),
    (97, "15.850", "2.609"),
    (98, "16.544", "2.532"),
    (99, "17.253", "2.446"),
)


def count_tail_upper(height):
    """Upper dyadic for the reciprocal-square mass above height >= e."""
    height = arb(height)
    if not height >= arb(1).exp():
        raise ValueError("Bellotti–Wong requires height >= e")
    ell = height.log()
    two_pi = 2 * arb.pi()
    a, b, c = arb("0.10076"), arb("0.24460"), arb("8.08344")
    error = a * ell + b * ell.log() + c
    slope = a + b / ell
    # E(u) <= E(height) + slope*log(u/height), from log concavity in log(u).
    return (((ell - two_pi.log() + 1) / (two_pi * height))
            + (2 * error + slope / 2) / height**2).upper()


class TailBounds:
    """Enclose upper bounds for the loss from off-line zero quartets."""

    def __init__(self, profile="graded"):
        if profile not in ("coarse", "graded"):
            raise ValueError("tail_profile must be coarse or graded")
        self.profile = profile
        height = arb(HEIGHT)
        log_height = height.log()
        self.all_mass_upper = count_tail_upper(height)

        # D(sigma) <= 2*integral_{HEIGHT}^infty N(sigma,u)*u^(-3) du.
        # Substitute u=exp(v) in the KLN power term. The logarithmic term
        # has a closed form from integrating v^2 exp(-2v) by parts.
        log_mass = (log_height**2 + log_height + arb(1) / 2) / height**2
        self.thresholds = tuple(row[0] for row in KLN_ROWS)
        mass_upper = (self.all_mass_upper / 2).upper()
        masses = []
        for hundredths, a, b in KLN_ROWS:
            sigma = arb(hundredths) / 100
            shape = 6 - 2 * sigma
            decay = 2 - arb(8) / 3 * (1 - sigma)
            first = (decay * log_height).gamma_upper(shape) / decay**shape
            bound = (2 * arb(a) * first + arb(b) * log_mass).upper()
            mass_upper = min(mass_upper, bound)
            masses.append(mass_upper)
        self.band_mass_uppers = tuple(masses)

        edges = tuple(arb(sigma - 50) / 100 for sigma in self.thresholds) + (arb(1) / 2,)
        self.base_theta = edges[0]
        self.increments = tuple(((after + before) / 2, (after - before) / 2)
                                for before, after in zip(edges, edges[1:]))

    def coarse(self, t):
        """Use the critical strip and the all-zero reciprocal-square tail."""
        return 2 * self.all_mass_upper * (t / 2).cosh()

    def graded(self, t):
        """Use the cumulative KLN bounds and positive increments of cosh."""
        # Charge every zero the baseline. A right-hand zero above a
        # threshold adds the corresponding increment for its whole quartet.
        result = 2 * self.all_mass_upper * (self.base_theta * t).cosh()
        for mass, (mean, half_difference) in zip(self.band_mass_uppers, self.increments, strict=True):
            # cosh(after*t) - cosh(before*t), in a form stable near t=0.
            result += 8 * mass * (mean * t).sinh() * (half_difference * t).sinh()
        return result

    def bound(self, t):
        """Use a monotone pointwise majorant on the selected profile."""
        coarse = self.coarse(t)
        if self.profile == "coarse":
            return coarse.upper()
        return min(coarse.upper(), self.graded(t).upper())

    def metadata(self):
        """Return the dyadic upper endpoints and their analytic provenance."""
        sources = {
            "finite_height": "https://doi.org/10.1112/blms.12460",
            "zero_count": "https://doi.org/10.1090/mcom/4133",
            "zero_density_kln": "https://doi.org/10.1016/j.jmaa.2018.04.071",
        }
        return {
            "verified_critical_line_height": HEIGHT,
            "all_positive_ordinate_mass_upper_dyadic": _dyadic(self.all_mass_upper),
            "density_bounds": [
                {
                    "sigma": f"{hundredths}/100",
                    "mass_upper_dyadic": _dyadic(mass),
                }
                for hundredths, mass in zip(self.thresholds, self.band_mass_uppers, strict=True)
            ],
            "sources": sources,
        }
