import sympy as sp

from desitter_qft import backgrounds as bg

H = sp.Symbol("H", positive=True)


def test_ricci_scalar_in_three_charts_is_positive_12H2():
    for metric in (bg.flat_chart_metric_cosmic, bg.flat_chart_metric_conformal, bg.static_patch_metric):
        g, X = metric()
        assert sp.simplify(bg.ricci_scalar(g, X) - 12 * H**2) == 0


def test_ricci_tensor_is_proportional_to_metric():
    g, X = bg.flat_chart_metric_conformal()
    Ric = bg.ricci_tensor(g, X)
    assert sp.simplify(Ric - 3 * H**2 * g) == sp.zeros(4, 4)  # R_mn = Lambda g_mn, Lambda = 3H^2


def test_generic_frw_formula_R_equals_6_app_over_a3():
    g, X, a = bg.frw_conformal_metric()
    R = bg.ricci_scalar(g, X)
    assert sp.simplify(R - 6 * sp.diff(a, X[0], 2) / a**3) == 0


def test_sign_convention_negative_control():
    """With the opposite Riemann sign convention R would be -12 H^2: the test discriminates."""
    g, X = bg.flat_chart_metric_cosmic()
    assert sp.simplify(bg.ricci_scalar(g, X) + 12 * H**2) != 0
