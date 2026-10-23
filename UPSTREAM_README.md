# Yildirim Invariant Geometry - Reproducibility and Submission Package

Version: 1.0 (12 August 2026)

This package contains a corrected, invariant-centered reformulation of the legacy recursive geometry together with exact and synthetic falsification experiments.

## Core conclusion

The stated recursion

`Conv(P x {-1,0,+1})`

is exactly `P x [-1,1]` for convex P. Starting from a point, the construction produces the standard hypercube family. The package therefore does not retain the earlier octagonal/new-polytope interpretation. It studies instead which structural, arithmetic-lattice, polarized-toric, and discriminant information survives declared transformation quotients and how discriminative those invariants are in controlled synthetic experiments.

## Contents

- `manuscript/` - main article (Markdown, DOCX, PDF, LaTeX when generated)
- `supplement/` - supplementary appendix
- `code/run_all_experiments.py` - full deterministic experiment driver
- `code/requirements.txt` - minimal dependencies
- `results/results_summary.json` - main machine-readable output
- `results/tables/` - detailed CSV/JSON experiment results
- `results/figures/` - figures used in the manuscript
- `docs/CLAIM_BOUNDARIES.md` - statements that are supported vs. explicitly not claimed
- `docs/REVISION_MAP.md` - mapping from the legacy formulation to the corrected invariant formulation

## Reproduce

From the package root:

```bash
python code/run_all_experiments.py --atlas-n 2000
```

Master seed: `20260812`.

The final reference run used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Exact integer/Fraction arithmetic is used for discrete claims; SciPy is used only for an independent numerical convex-hull check.

## Interpretation

`R=1` for an invariant in the random atlas means collision-free/sufficient behavior for the explicitly defined finite operational target. It is not a universal completeness theorem. The manuscript states this boundary explicitly.

## Submission note

The documents are prepared in a high-rigor journal/preprint style, but no quartile placement or acceptance is guaranteed. Journal-specific formatting, word limits, disclosure language and data/code repository links should be adapted after a target journal is selected.
