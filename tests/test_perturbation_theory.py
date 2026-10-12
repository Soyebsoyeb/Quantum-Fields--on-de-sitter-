import numpy as np
import pytest

from desitter_qft.in_in import (conformal_contact_correlator, contact_correlator, contact_integral_mp,
                                contact_integral_real_axis_abel, massless_minimal_log_coefficient,
                                schwinger_keldysh_branches)
from desitter_qft.perturbation_theory import (anharmonic_ground_state_cumulant, qm_first_order_slope,
                                              qm_inin_prediction, secular_ratio)

KS = [1.0, 1.3, 0.7, 1.1]


def test_inin_formula_against_exact_diagonalisation_sign_and_normalisation():
    """Independent check: 0+1-dim anharmonic oscillator."""
    for w in (1.0, 2.0):
        assert abs(qm_first_order_slope(w) / qm_inin_prediction(1.0, w) - 1) < 1e-4
    # negative control: wrong sign would be caught
    assert qm_first_order_slope(1.0) * qm_inin_prediction(1.0, 1.0) > 0


def test_qm_cumulant_second_order_scaling():
    g = np.array([0.01, 0.02])
    k = [anharmonic_ground_state_cumulant(x, 1.0) for x in g]
    dev = [abs(ki - qm_inin_prediction(x, 1.0)) for ki, x in zip(k, g)]
    assert 3.2 < dev[1] / dev[0] < 4.8  # O(g^2) remainder


@pytest.mark.parametrize("eta", [-3.0, -0.5, -0.02])
def test_conformal_contact_numeric_equals_analytic(eta):
    assert abs(contact_correlator(eta, KS, 0.5) / conformal_contact_correlator(eta, KS, 1.0) - 1) < 1e-12


def test_conformal_invariance_trispectrum_scales_as_eta4():
    a = contact_correlator(-1.0, KS, 0.5)
    b = contact_correlator(-0.1, KS, 0.5)
    assert abs(b / a - 1e-4) < 1e-12


def test_contour_rotation_agrees_with_real_axis_abel_limit():
    I = complex(contact_integral_mp(-0.5, KS, 1.5))
    eps = np.array([0.02, 0.01, 0.005])
    vals = np.array([contact_integral_real_axis_abel(-0.5, KS, 1.5, e) for e in eps])
    extrap = 2 * vals[2] - vals[1]  # linear Richardson in epsilon
    assert abs(extrap - I) < 5e-4 * abs(I)


def test_scipy_and_mpmath_backends_agree_for_generic_nu():
    a = contact_correlator(-0.05, KS, 1.0, backend="scipy")
    b = contact_correlator(-0.05, KS, 1.0, backend="mp")
    assert abs(a / b - 1) < 1e-8


def test_massless_minimal_late_time_logarithm_has_derived_coefficient():
    e1, e2 = -1e-4, -1e-5
    slope = (contact_correlator(e2, KS, 1.5) - contact_correlator(e1, KS, 1.5)) / (np.log(-e2) - np.log(-e1))
    assert abs(slope / massless_minimal_log_coefficient(KS, 1.0) - 1) < 1e-4


def test_secular_growth_pattern():
    N = np.array([2.0, 4.0, 6.0]); eta = -np.exp(-N)
    conf = secular_ratio(eta, [1.0] * 4, 0.5)
    assert np.allclose(conf, conf[0], rtol=1e-9)             # conformal invariance: no growth
    m = secular_ratio(eta, [1.0] * 4, 1.5)
    assert (m[1] - m[0]) < 0 and (m[2] - m[1]) < 0           # monotone growth in magnitude
    assert abs((m[2] - m[1]) - (m[1] - m[0])) < 0.15 * abs(m[2] - m[1])   # ~linear in N
    s = secular_ratio(eta, [1.0] * 4, 1.0)
    assert abs(s[2] - s[1]) < 0.25 * abs(s[1] - s[0])        # light massive: saturates


def test_schwinger_keldysh_branches_sum_to_commutator_form():
    plus, minus = schwinger_keldysh_branches(-0.5, KS, 0.5)
    assert abs(plus - np.conj(minus)) < 1e-14 or abs(plus.imag + minus.imag) < 1e-12
    expected = -1.0 / (8 * sum(KS) * np.prod(KS))
    assert abs((plus + minus).real - expected) < 1e-12 and abs((plus + minus).imag) < 1e-12
