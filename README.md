# Selected code and numerical data: finite-noise lattice-surgery decoding

Niharika Verma — Indian Institute of Technology, Patna.

A deliberately small research companion containing only the code and final numerical data needed to inspect the central lattice-surgery experiments. This is **not** the full working directory. The manuscript, figures, experiment plans, audit records, drafts, research reports and historical snapshots are excluded.

## Contents

- `simulate.py`: sample a selected circuit and compare ordinary or correlated matching with the specified local decisions.
- `predict.py`: compute the all-local-orders / at-most-two-background-fault likelihood approximation.
- `physical_channels.py`: the physical-channel propagation helper used by that calculation.
- `analyze.py`: recompute paired effects and confidence intervals from aggregate counts.
- `circuits/`: 37 serialized Stim circuits for the original noise sweep, matched intervention, mixed-bias comparison, and temporal-boundary extension.
- `data/`: eight selected numerical datasets: condition parameters, paired counts, likelihood forecasts, joint likelihood results, intervention channels, leading coefficients, correlated comparisons, and the fixed-gate reset/readout results.

## Usage

Use Python 3.12 and the pinned dependencies:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python analyze.py
python simulate.py mixed_00 --shots 10000
python simulate.py mixed_00 --shots 10000 --correlated
python predict.py mixed_00
```

`--original-shots` explicitly requests the original full sampling budget (10 or 20 million shots per condition). The default 10,000 shots is a smoke run, not a statistically adequate confirmation of ppm-scale differences. The saved condition seeds and 50,000-shot chunking describe the original paired simulations. Reproducing seeded results depends on retaining the stated software versions. No complete full-budget replay is claimed for this curated export.

## Reading the data

Probabilities and effects are stored as fractions, not percentages or ppm. Multiply an effect by 1e6 to obtain ppm. A positive paired effect means the modified policy makes fewer parity errors. Counts are `[ordinary_errors, modified_errors, repairs, introduced_errors]`.

Each row retains its simultaneous confidence-family size. Intervals use conservative Clopper–Pearson marginal bounds. Separate experiments do not share one global 95% statement. Correlated matching reuses the mixed-bias records and is not additional independent evidence. Exact identical policies must not be interpreted as a measured nonzero improvement.

The likelihood approximation retains mutually exclusive outcomes at each physical location and does not renormalize discarded background mass. Its omitted-mass bound is too broad for a general sign guarantee. The selected local decisions can sacrifice correction of individual physical faults. This code is not a new general near-ML decoder or a full logical-instrument benchmark.

The joint-likelihood and reset/readout tables preserve final all-order estimates and intervals. Their original large estimator checkpoints are intentionally excluded; this small package does not by itself reproduce every historical likelihood-estimator interval or every result in the manuscript.

## Provenance and visibility

The paired circuits were exported from the study and their SHA-256 hashes checked against the saved experiment records. The compact scripts are a curated extraction of the sampling, inference and interval calculations, not byte-identical copies of the original experiment drivers. They contain no plans, audit reports or figure-generation code.

Circuit construction used FTPrimitiveBench revision `b75fbd1c514bc4fb7fcca75618a191e6b58fa7ed` (https://github.com/ShuwenKan/FTPrimitiveBench). No third-party source tree is bundled. Stim and PyMatching are installed as dependencies.

This repository remains **private** as requested. Its selected contents are being prepared for a possible public release; no public deposit, public preregistration, submission, or DOI is claimed. No new reuse license has been assigned.
