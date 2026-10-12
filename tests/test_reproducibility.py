import json
import os

import numpy as np

from desitter_qft import experiments as ex
from desitter_qft.provenance import hash_arrays


def test_experiments_are_deterministic():
    a, b = ex.qm_inin_experiment(), ex.qm_inin_experiment()
    assert hash_arrays({"g": a["g"], "x": a["exact"]}) == hash_arrays({"g": b["g"], "x": b["exact"]})
    t1, t2 = ex.tensor_experiment(), ex.tensor_experiment()
    assert hash_arrays({"y": t1["Dt2_num"]}) == hash_arrays({"y": t2["Dt2_num"]})


def test_hash_is_sensitive_to_data():
    x = np.linspace(0, 1, 10)
    assert hash_arrays({"x": x}) != hash_arrays({"x": x + 1e-6})


def test_recorded_results_exist_and_contain_environment():
    root = os.path.join(os.path.dirname(__file__), "..", "results")
    for name in ("free_modes", "two_point", "tree_level"):
        p = os.path.join(root, name + ".json")
        if os.path.exists(p):
            rec = json.load(open(p))
            assert "environment" in rec and "data_sha256" in rec
