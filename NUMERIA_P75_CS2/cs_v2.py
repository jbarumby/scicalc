"""Directed Arb majorants from Chourasiya--Simonic arXiv:2507.15184v2.

This implements Corollary 1, Table 1 of the *v2* preprint, never v1.
The source is a PREPRINT. Certification of interval arithmetic is separate
from a mathematical audit of the source theorem and its hypotheses.
"""
from fractions import Fraction
from flint import arb
from NUMERIA_P74_ZERO_SIDE.tails import HEIGHT, count_tail_upper
from NUMERIA_ES8P259.suzuki import _dyadic

# The 16 consecutive CLOSED strips [k/32,(k+1)/32], k=16,...,31.
# Literal printed decimal constants, as in arXiv:2507.15184v2 table 1.
TABLE_V2 = (
 ('5.069','9.461','167.8'), ('6.046','9.168','162.6'),
 ('7.077','8.875','157.4'), ('8.185','8.582','152.2'),
 ('9.387','8.290','147.0'), ('10.70','7.997','141.9'),
 ('12.11','7.704','136.7'), ('13.66','7.411','131.5'),
 ('15.35','7.118','126.3'), ('17.21','6.826','121.1'),
 ('19.31','6.533','115.9'), ('21.77','6.240','110.7'),
 ('24.83','5.745','101.9'), ('29.01','5.463','96.87'),
 ('35.35','5.180','91.86'), ('46.06','4.897','86.84'),
)


def as_ball(q: Fraction) -> arb:
    return arb(q.numerator) / q.denominator


def row_indices(sigma: Fraction):
    """Both closed-interval rows are valid at a shared endpoint."""
    if not Fraction(1, 2) <= sigma <= 1:
        raise ValueError('sigma must be in [1/2,1]')
    return tuple(i for i in range(16)
                 if Fraction(16+i,32) <= sigma <= Fraction(17+i,32))


def density_integral_upper(sigma: Fraction) -> arb:
    """Upper bound for D_sigma(H) = sum_{gamma>H, beta>=sigma} gamma^-2.

    For sigma>1/2, the verified finite RH height gives N(sigma,H)=0.
    Integration by parts yields 2*int_H^inf N(sigma,u)*u^-3 du.
    Integrate the explicit Corollary 1 v2 majorant; use outward Arb.
    """
    sigma = Fraction(sigma)
    if sigma <= Fraction(1, 2) or sigma > 1:
        raise ValueError('density integral is only used at 1/2<sigma<=1')
    if sigma == 1: return arb(0)
    s = as_ball(sigma)
    logh = arb(HEIGHT).log()
    alpha = 3 * (1-s) / (2-s)
    power = (7-5*s) / (2-s)
    lam = 2-alpha
    if not lam > 0:
        raise ArithmeticError('lambda must be positive')
    shape = power + 1
    x = lam * logh
    # x.gamma_upper(shape) is the rigorous upper incomplete gamma integral
    #  Gamma(shape,x) = integral_x^infinity exp(-v) v^(shape-1) dv.
    gamma_piece = (x.gamma_upper(shape) / lam**shape).upper()
    quadratic = ((logh**2+logh+arb(1)/2) / arb(HEIGHT)**2).upper()
    linear = ((logh+arb(1)/2) / arb(HEIGHT)**2).upper()
    estimates = []
    for idx in row_indices(sigma):
        b1,b2,b3 = [arb(y) for y in TABLE_V2[idx]]
        estimates.append((2*b1*gamma_piece+b2*quadratic+b3*linear).upper())
    if not estimates:
        raise ArithmeticError('no valid Table 1 row')
    return min(estimates)


class ChourasiyaV2Tail:
    """A monotone upper bound on the potentially off-line zero contribution.

    The mass W counts positive ordinates, with multiplicity. At sigma>1/2,
    the upper half strip has mass at most W/2 by reflection beta<->1-beta.
    The layer-cake loss is at most:
       U(t) = 2W + 4 sum_j D(s_j+) (cosh(theta_{j+1}t)-cosh(theta_j t)).
    The first band uses D(1/2+) <= W/2, not N(1/2,H).
    """
    def __init__(self, segments=128):
        if type(segments) is not int or segments<16 or segments>512 or segments & (segments-1):
            raise ValueError('segments must be a power of two between 16 and 512')
        self.segments=segments
        self.total_mass=count_tail_upper(HEIGHT).upper()
        sigma = [Fraction(1,2)+Fraction(i,2*segments) for i in range(segments+1)]
        cap=(self.total_mass/2).upper()
        self.mass_caps=[]
        for s in sigma[:-1]:
            if s>Fraction(1,2):
                cap=min(cap, density_integral_upper(s))
            self.mass_caps.append(cap)
        self.theta=[as_ball(s-Fraction(1,2)) for s in sigma]

    def bound(self,t):
        if not t >= 0:
            raise ValueError('t must be nonnegative')
        result=2*self.total_mass
        # cosh(b)-cosh(a) = 2*sinh((a+b)/2)*sinh((b-a)/2),
        # each factor >=0 for nonnegative t. No cancellation issue.
        for j,cap in enumerate(self.mass_caps):
            mean=(self.theta[j+1]+self.theta[j])*t/2
            half=(self.theta[j+1]-self.theta[j])*t/2
            increment=(2*mean.sinh()*half.sinh()).upper()
            result += 4*cap*increment
        # The all-strip bound is always a valid alternative.
        coarse=(2*self.total_mass*(t/2).cosh()).upper()
        return min(result.upper(),coarse)

    def metadata(self):
        return {
          'source':'https://arxiv.org/html/2507.15184v2',
          'source_type':'preprint_v2_not_peer_review_gate',
          'source_formula':'Corollary 1 (3), Table 1; exponent (7-5*sigma)/(2-sigma)',
          'finite_rh_height':HEIGHT,
          'all_positive_ordinate_mass_upper_dyadic':_dyadic(self.total_mass),
          'segments':self.segments,
          'row_coefficients_v2':[list(r) for r in TABLE_V2],
          'mass_caps_dyadic':[_dyadic(v) for v in self.mass_caps],
        }

# P84: explicit S2 public-source B1 table, opt-in only.  The original
# TABLE_V2 above is retained literally for reproducibility of P75/P77.
from NUMERIA_P76_COEFFICIENTS.table import TABLE_S2_PUBLIC_ROBUST as _TABLE_S2
TABLE_V2_S2_PUBLIC_ROBUST = tuple((r.b1, r.b2, r.b3) for r in _TABLE_S2)
if len(TABLE_V2_S2_PUBLIC_ROBUST) != 16:
    raise ArithmeticError("S2 B1 table requires all sixteen strips")

