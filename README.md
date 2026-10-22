# Competing Fault Histories and Finite-Noise Trade-offs in Lattice-Surgery Decoding


This repository contains code and numerical data for a study of competing fault histories in lattice-surgery decoding. The question is simple: when different faults produce the same detector record, can correcting a particular single fault increase the average logical error?

The files here cover the main sampling experiments and likelihood calculations. The manuscript, figures, experiment plans, internal reports, and earlier working files are kept separately.

## What is included

- `simulate.py` runs a selected circuit and compares ordinary or correlated matching with the local decoding changes.
- `predict.py` calculates logical-label probabilities using all local fault orders and up to two background faults.
- `physical_channels.py` provides the physical-fault propagation routines used in that calculation.
- `analyze.py` recomputes paired effects and confidence intervals from the saved counts.
- `circuits/` contains 52 Stim circuits for the original and extended noise sweeps, matched interventions, mixed-bias tests, and preparation–recovery extension.
- `data/` contains ten datasets covering the experimental conditions, counts, predictions, intervention channels, expansion coefficients, and reset/readout comparisons.

## Running the code

Use Python 3.12 with the dependencies in `requirements.txt`:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

Start by analyzing the saved results:

```sh
python analyze.py
```

To sample a circuit or evaluate its likelihood prediction:

```sh
python simulate.py mixed_00 --shots 10000
python simulate.py mixed_00 --shots 10000 --correlated
python predict.py mixed_00
python simulate.py dense_00 --shots 10000
```

A 10,000-shot run is useful for checking the installation. It cannot resolve the small error-rate differences reported in the study.

The `--original-shots` option uses the original budget of 10 or 20 million shots per condition. Reproducing a seeded run requires the recorded seed, the original 50,000-shot batching, and the specified software versions. The exported package has not been replayed in full at those budgets.

## Reading the data

Probabilities are stored as fractions. Multiply by \(10^6\) to convert an effect to parts per million.

A positive effect means that the modified decoder makes fewer parity errors. The count arrays use the order

```text
[ordinary_errors, modified_errors, repairs, introduced_errors]
```

The intervals are built from Clopper–Pearson bounds, with simultaneous coverage calculated separately for each experimental family. Each row records the family size used. There is no single 95% confidence statement covering the whole study.

`data/dense_results.json` contains all 27 conditions used in the final dense-sweep figures: 17 noise scales, nine targeted attenuation fractions, and one matched suppression. Each has 20 million shots. These rows reuse the original counts where applicable; they are not additional independent experiments. Their displayed intervals use a retrospective family of 270 marginals (four decoder-event counts and six syndrome/label counts per condition). The earlier intervals in `paired_results.json` are retained with their original family sizes.

`data/attenuation_predictions.json` gives the nine deterministic predictions overlaid on the attenuation sweep. These predictions were calculated after the experiment, without fitting to the observations. They are separate from the prospective mixed-noise and circuit-extension tests. `predict.py` can recompute them from the corresponding condition identifiers.

Correlated matching is evaluated on the same mixed-bias records as ordinary matching. Those comparisons are not independent replications. Some policies also make exactly the same decisions; their zero difference is an identity, not an estimated improvement.

## What the package can reproduce

The scripts support circuit sampling, the selected likelihood approximation, and analysis of the saved aggregate counts. `analyze.py` also checks the circuit hashes and reconstructs all 270 marginal intervals used in the dense figures, along with the target-versus-control contrasts. `simulate.py` reports the complete-syndrome label counts as well as the paired decoder outcomes. They do not reproduce every calculation in the paper.

In particular, the joint-likelihood and reset/readout datasets contain the final estimates and intervals, but their large estimator checkpoints are not included. The scripts were extracted from the research code and are not byte-for-byte copies of the original experiment drivers.

The likelihood approximation preserves mutually exclusive outcomes at each fault location and leaves discarded background probability unnormalized. Its omitted-mass bound does not generally certify the predicted sign.

The decoding changes studied here can improve average parity error while failing on selected individual faults. The results concern that trade-off in the tested circuits. They do not establish a general near-maximum-likelihood decoder or describe every error in the logical operation.

## Circuit provenance

The circuits were constructed using [FTPrimitiveBench](https://github.com/ShuwenKan/FTPrimitiveBench), revision

```text
b75fbd1c514bc4fb7fcca75618a191e6b58fa7ed
```

Their SHA-256 hashes were checked against the saved experiment records. No third-party source tree is included; Stim and PyMatching are installed as dependencies.
