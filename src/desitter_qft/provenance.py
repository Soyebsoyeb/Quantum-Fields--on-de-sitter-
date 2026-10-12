"""Run metadata and deterministic hashing for reproducibility records."""
from __future__ import annotations

import hashlib
import json
import platform
import sys

import numpy as np


def environment_info() -> dict:
    import matplotlib
    import mpmath
    import scipy
    import sympy

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "sympy": sympy.__version__,
        "mpmath": mpmath.__version__,
        "matplotlib": matplotlib.__version__,
        "float_precision": "IEEE-754 binary64 (mpmath used where stated)",
    }


def hash_arrays(data: dict, decimals: int = 10) -> str:
    """SHA-256 of rounded array content, insensitive to last-bit noise (rounded to `decimals`)."""
    h = hashlib.sha256()
    for key in sorted(data):
        val = data[key]
        h.update(key.encode())
        if isinstance(val, dict):
            h.update(hash_arrays(val, decimals).encode())
            continue
        arr = np.asarray(val)
        if np.iscomplexobj(arr):
            arr = np.stack([arr.real, arr.imag])
        if arr.dtype.kind in "fc" or arr.dtype.kind == "f":
            arr = np.round(arr.astype(float), decimals) + 0.0
        h.update(np.ascontiguousarray(arr).tobytes() if arr.dtype.kind in "fiub" else str(arr).encode())
    return h.hexdigest()


def write_record(path, record: dict) -> None:
    record = {**record, "environment": environment_info()}
    with open(path, "w") as f:
        json.dump(record, f, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
