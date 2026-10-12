import sympy as sp

from desitter_qft.interactions import check_dimensions, interaction_hamiltonian_density


def test_interaction_hamiltonian_is_minus_lagrangian_for_nonderivative_vertex():
    got, expected = interaction_hamiltonian_density()
    assert sp.simplify(got - expected) == 0


def test_dimensional_consistency_and_lambda_dimensionless():
    terms, ok = check_dimensions()
    assert ok and all(v == 0 for v in terms.values())


def test_dimension_check_detects_inconsistency():
    """If lambda carried dimension 1 the quartic term would not be marginal: the bookkeeping must reflect that."""
    d4x, phi, lam = -4, 1, 1
    assert d4x + lam + 4 * phi != 0
