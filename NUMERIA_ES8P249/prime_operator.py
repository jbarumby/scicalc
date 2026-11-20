"""Matrix-free finite prime operators on a uniform grid of orthonormal cells."""

from fractions import Fraction
import math
import operator

import numpy as np


P236_HALF_WIDTH = Fraction(161, 200)
P236_PRIME_POWERS = ((2, 2), (3, 3), (4, 2), (5, 5))  # (prime power, prime)


def _cell_width(cells, half_width):
    cells = operator.index(cells)
    half_width = float(half_width)
    if cells < 1 or not math.isfinite(half_width) or half_width <= 0:
        raise ValueError('cells and half_width must be positive')
    width = 2 * half_width / cells
    if width <= 0 or not math.isfinite(width):
        raise ValueError('cell width must be representable in binary64')
    return width


def _prime_parameters(powers):
    parameters = []
    for power, prime in powers:
        power, prime = operator.index(power), operator.index(prime)
        if power < 2 or prime < 2:
            raise ValueError('prime and prime power must be at least two')
        parameters.append((math.log(power), 2 * math.log(prime) / math.sqrt(power)))
    return parameters


def _frequencies(values):
    frequencies = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(frequencies)):
        raise ValueError('frequencies must be finite')
    return frequencies


def _evaluate_symbol(frequencies, parameters):
    result = np.zeros_like(frequencies)
    for shift, weight in parameters:
        result += weight * np.cos(frequencies * shift)
    return result


def prime_symbol(frequencies, powers=P236_PRIME_POWERS):
    """Evaluate S(t) = sum_q a_q cos(t log q) with the P236 sign before subtraction."""
    return _evaluate_symbol(_frequencies(frequencies), _prime_parameters(powers))


def dyadic_prime_components(frequencies, max_level, powers=P236_PRIME_POWERS):
    """Return the coarse approximation, dyadic details, and final residual.

    The level-j approximation is first-order Taylor on cells of width 2**-j.
    Its analytic error is at most 2**(-2*j) * sum_q(a_q * log(q)**2) / 8.
    The arrays are binary64 evaluations of the real decomposition.
    """
    max_level = operator.index(max_level)
    if not 0 <= max_level <= 52:
        raise ValueError('max_level must be between 0 and 52 for binary64')
    frequencies = _frequencies(frequencies)
    parameters = _prime_parameters(powers)
    previous = np.zeros_like(frequencies)
    components = []
    for level in range(max_level + 1):
        width = math.ldexp(1.0, -level)
        centers = (np.floor(frequencies / width) + 0.5) * width
        offset = frequencies - centers
        current = np.zeros_like(frequencies)
        for shift, weight in parameters:
            phase = centers * shift
            current += weight * (np.cos(phase) - offset * shift * np.sin(phase))
        components.append(current - previous)
        previous = current
    residual = _evaluate_symbol(frequencies, parameters) - previous
    return (*components, residual)


class ToeplitzFFT:
    """Real symmetric Toeplitz matvec with a zero-padded circulant embedding."""

    def __init__(self, lags):
        lags = np.array(lags, dtype=np.float64, copy=True)
        if lags.ndim != 1 or not lags.size or not np.all(np.isfinite(lags)):
            raise ValueError('lags must be a nonempty finite vector')
        lags.flags.writeable = False
        self.lags = lags
        self._size = lags.size
        self._fft_size = 1 << (2 * self._size - 2).bit_length()
        embedding = np.zeros(self._fft_size)
        embedding[:self._size] = lags
        if self._size > 1:
            embedding[-(self._size - 1):] = lags[:0:-1]
        self._spectrum = np.fft.fft(embedding)

    def matvec(self, vector):
        """Apply the Toeplitz operator to a real or complex coefficient vector."""
        vector = np.asarray(vector)
        if vector.shape != (self._size,) or vector.dtype.kind not in 'biufc':
            raise ValueError('vector must be numeric and match the number of lags')
        if not np.all(np.isfinite(vector)):
            raise ValueError('vector must be finite')
        transformed = np.fft.fft(vector, n=self._fft_size)
        product = np.fft.ifft(self._spectrum * transformed)[:self._size]
        return product if np.iscomplexobj(vector) else product.real


class PrimeShiftOperator(ToeplitzFFT):
    """Compression of sum_q a_q(T_{log(q)} + T_{-log(q)})/2 to the cell basis.

    Cell functions are extended by zero. This is the full-frequency prime
    term S(t) before the subtraction in P236, with no spectral cutoff.
    """

    def __init__(self, cells, half_width=P236_HALF_WIDTH, powers=P236_PRIME_POWERS):
        width = _cell_width(cells, half_width)
        lags = np.zeros(cells)
        for shift, weight in _prime_parameters(powers):
            if shift >= cells * width:
                continue
            distance = shift / width
            lag = math.floor(distance)
            fraction = distance - lag
            if lag < cells:
                # The two shift directions coincide on lag zero.
                lags[lag] += weight * (1 - fraction) / (1 if lag == 0 else 2)
            if lag + 1 < cells:
                lags[lag + 1] += weight * fraction / 2
        super().__init__(lags)


def fourier_cell_lags(cells, half_width, frequencies, weights):
    """Lags of a finite real Fourier quadrature on the same cell basis.

    Weights include quadrature normalization and any signed multiplier.
    Integrating each cell's Fourier phase supplies the squared sinc factor.
    """
    width = _cell_width(cells, half_width)
    frequencies = _frequencies(frequencies)
    weights = np.asarray(weights, dtype=np.float64)
    if frequencies.ndim != 1 or weights.shape != frequencies.shape:
        raise ValueError('frequencies and weights must be equally sized vectors')
    if not np.all(np.isfinite(weights)):
        raise ValueError('weights must be finite')
    distances = width * np.arange(cells)
    lags = np.zeros(cells)
    for t, weight in zip(frequencies, weights):
        sinc = np.sinc(t * width / (2 * np.pi))
        lags += (width * weight * sinc**2) * np.cos(t * distances)
    return lags
