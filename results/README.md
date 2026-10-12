# results/

Machine-readable records written by `desitter-qft` (do not edit by hand):

| File | Written by | Content |
|---|---|---|
| `free_modes.json` | `scripts/reproduce_free_modes.py` | parameters, solver settings, Wronskian reports, tolerance sweep, data hash, environment |
| `two_point.json` | `scripts/reproduce_two_point.py` | mode-sum vs closed-form errors, propagator diagnostics |
| `tree_level.json` | `scripts/reproduce_tree_correlator.py` | tree-level correlator, in-in vs exact-diagonalisation, tensor spectrum |
| `validation.json` | `desitter-qft validate` | all programmatic checks with tolerances |

`data_sha256` is a SHA-256 of the output arrays rounded to 10 decimals (insensitive to last-bit noise). `reports/technical_report.md` section 17 is generated from these files by `scripts/build_report.py`.
