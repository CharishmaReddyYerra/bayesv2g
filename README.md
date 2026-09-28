# BayesV2G

BayesV2G is a reproducible Bayesian variant-to-gene prioritization engine for
LDL cholesterol. It ranks candidate genes within a fine-mapped GWAS signal and
reports posterior uncertainty and a held-out calibration assessment rather than
an unscaled score.

This repository uses real, versioned data throughout:

- Open Targets Locus-to-Gene 26.06 train/test data: 115,852 gene-locus rows
  across 7,231 loci, pinned to a Hugging Face commit.
- Open Targets Platform 26.09: fine-mapped credible sets and functional evidence
  at ten real LDL signals from study `GCST90019512`.
- GWAS Catalog `GCST90018961`: 5,772 real summary-statistic variants in the
  rs12740374/SORT1 region, downloaded from its indexed GRCh38 file.
- eQTL Catalogue: 8,631 GTEx liver eQTL rows and 3,892 Sun et al. plasma pQTL
  rows for regional, allele-matched colocalization at the SORT1 locus.
- GTEx v10 liver and Ensembl: lead-variant eQTL, gene, VEP, and LD annotations
  cached with checksums for the independent SORT1 audit.

The ready-to-read output is [the HTML report](results/ldl_demo/report.html).
Machine-readable rankings are in
[`prioritized_genes.csv`](results/ldl_demo/prioritized_genes.csv), and the exact
run provenance is in [`run_manifest.json`](results/ldl_demo/run_manifest.json).
The interactive [Evidence Workbench](https://bayesv2g-evidence-workbench.sibunnu.chatgpt.site)
provides a polished locus explorer, validation view, and independent SORT1 audit
over those same checked-in outputs. Its source and local setup instructions are
in [`dashboard/`](dashboard/README.md).

## Main result

On 1,048 held-out labelled loci, BayesV2G has higher observed ranking metrics
than an equal-probability baseline and the distance-only comparator:

| Method | Top-1 accuracy | Recall@3 | Mean reciprocal rank |
|---|---:|---:|---:|
| BayesV2G | 0.731 | 0.860 | 0.816 |
| Distance-only | 0.658 | 0.803 | 0.744 |
| Uniform | 0.114 | 0.327 | 0.294 |

The paired top-1 difference versus distance-only is 0.073, with a 95% locus-
bootstrap interval of 0.046–0.102 (10,000 fixed-seed resamples).

Several real loci illustrate the intended behavior: ABCG8, HMGCR, and APOE rank
above closer candidates when the integrated evidence supports them. The SORT1
region also exposes genuine ambiguity: rs12740374 is a strong liver eQTL for
SORT1, CELSR2, and PSRC1, and the model ranks CELSR2 first under the selected
features. The project reports that result rather than forcing a familiar answer.

## Quick start

Python 3.11 or newer is required.

```bash
cd bayesv2g
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,data,download]'
```

For the exact recorded analysis/test versions on Python 3.12, install
`requirements.lock` and then the package without dependency resolution:

```bash
python -m pip install -r requirements.lock
python -m pip install -e . --no-deps
```

Download the pinned Open Targets artifacts with the Hugging Face CLI:

```bash
hf download opentargets/locus_to_gene_26.06 \
  --revision 6abe9cc490222986925e1acca696d48c9b2437a4 \
  --include 'train.parquet' 'test.parquet' 'config.json' 'README.md' \
  --local-dir data/raw/open_targets_l2g_26.06
```

The checked-in snapshots are sufficient for an offline, deterministic run:

```bash
bayesv2g run
pytest
```

To refresh the ten loci from the live Open Targets API:

```bash
bayesv2g run --refresh
```

To rebuild the independent SORT1 raw-data snapshot before running the model:

```bash
python scripts/fetch_sort1_case_study.py
bayesv2g run
```

Some institutional TLS proxies require the fetch script's explicit `--insecure`
option. The resulting manifest records whether certificate verification was
disabled. Never use that option silently in an automated production workflow.

## Statistical model

For gene `g` at locus `l`, the conditional model is

```text
eta_lg = beta' x_lg
P(g | l, surfaced candidates) = exp(eta_lg) / sum_h exp(eta_lh)
beta_k ~ Normal(0, 0.75^2)
```

The six pre-specified features are distance-to-TSS, eQTL colocalization, pQTL
colocalization, sQTL colocalization, VEP consequence, and enhancer-to-gene
evidence. A MAP estimate and analytic Hessian define a Laplace posterior; draws
from that posterior produce 90% probability intervals.

The categorical likelihood makes candidates compete within a locus and avoids
the incoherent independent-Bernoulli result in which several genes can each have
high probability. To acknowledge genes omitted by candidate ascertainment and
incomplete credible-set coverage, reported absolute probabilities are

```text
P_absolute(g) = (1 - rho) * credible_set_coverage * P_conditional(g)
P_unresolved  = 1 - (1 - rho) * credible_set_coverage
```

with `rho = 0.20` in the demonstration. This is a transparent sensitivity
parameter, not an empirically identified biological rate.

The raw SORT1 audit separately applies single-signal Wakefield approximate Bayes
factors to beta/SE pairs and runs the classic H0–H4 colocalization calculation on
2,400+ allele-matched SNPs per molecular trait. It is deliberately labelled
approximate because it ignores LD and cannot replace multi-signal fine-mapping or
SuSiE-coloc. The primary ranking uses Open Targets SuSiE-inf credible sets and
precomputed feature summaries instead.

## Outputs

`bayesv2g run` writes:

- candidate evidence, credible-set variants, and locus metadata;
- posterior gene rankings and coefficient intervals;
- held-out metrics, calibration, and sensitivity analyses;
- the fitted model and a checksum-rich run manifest;
- a `SHA256SUMS` file covering every generated artifact;
- a portable HTML report with linked publication-ready PNG figures;
- the complete and top-25 SORT1 ABF fine-mapping tables;
- regional SORT1 eQTL/pQTL colocalization summaries and variant-level audit data.

See [Methods](docs/METHODS.md), [Data dictionary](docs/DATA_DICTIONARY.md),
[Data provenance](data/README.md), and [Limitations](docs/LIMITATIONS.md) for
the details needed to reproduce and defend the analysis.

## Repository layout

```text
config/                 analysis configuration
data/raw/               pinned Open Targets parquet inputs
data/snapshots/         immutable API and regional-source snapshots
docs/                   methods, limitations, and output definitions
dashboard/              interactive evidence workbench website
results/ldl_demo/       completed analysis and HTML report
scripts/                auditable data-acquisition scripts
src/bayesv2g/           package source
tests/                   statistical and data-transformation tests
requirements.lock       exact reference analysis/test environment
```

## Scope

This is a complete research prototype, not a clinical decision system. The
learned ten-locus path consumes Open Targets' precomputed fine-mapping and
functional features; the independent SORT1 path reconstructs regional ABF
fine-mapping and colocalization directly from raw public summary statistics.
Probabilities are conditional on the model, labels, data releases, candidate
set, and unresolved-mass assumption. Open Targets' external L2G score is used
only to surface candidates and as a comparison field; it is not a model feature.

MIT licensed. Source datasets retain their original terms; see
[THIRD_PARTY_DATA_NOTICES.md](THIRD_PARTY_DATA_NOTICES.md) and
[data/README.md](data/README.md).
