# Methods

## Scientific target

For each independently fine-mapped LDL signal, BayesV2G estimates a probability
distribution over the candidate genes surfaced for that signal. The estimand is
therefore conditional:

```text
P(gene mediates locus | candidate set, evidence, labels, model, data releases)
```

It is not the unconditional probability that a gene is biologically causal.

## Data layers

The completed analysis has two intentionally separate paths.

The primary ten-locus analysis uses Open Targets Platform release 26.09 for
`GCST90019512`, an LDL biomarker GWAS. Its GraphQL response supplies SuSiE-inf
credible-set variants, posterior inclusion probabilities, locus metadata, and
gene-level evidence. The response is frozen as JSON before transformation.

The model is fitted and evaluated using Open Targets'
`locus_to_gene_26.06` dataset at commit
`6abe9cc490222986925e1acca696d48c9b2437a4`. The official train set contains
99,048 rows at 6,183 loci and the test set contains 16,804 rows at 1,048 loci.
Any configured demonstration loci found in the train partition are removed
before fitting; this removed 49 rows at three loci in the recorded run.

The independent SORT1 audit starts from GWAS Catalog `GCST90018961`, the Sakaue
et al. cross-population LDL meta-analysis. Its remote Tabix index is used to
retrieve only GRCh38 chromosome 1 positions 108,774,968–109,774,968. The local
slice contains 5,772 variants. Regional eQTL Catalogue release 7 summary data are
cached for reprocessed GTEx v8 liver `QTD000266` and Sun et al. plasma pQTL
`QTD000584`, alongside
GTEx v10 and Ensembl annotations. This audit is not used to fit the model.

## Evidence matrix

The primary model uses six pre-specified Open Targets features:

| Feature | Interpretation |
|---|---|
| `distanceTssMean` | PIP-weighted distance-to-TSS evidence |
| `eQtlColocH4Maximum` | maximum cis-eQTL shared-causal-variant support |
| `pQtlColocH4Maximum` | maximum cis-pQTL shared-causal-variant support |
| `sQtlColocH4Maximum` | maximum cis-sQTL shared-causal-variant support |
| `vepMean` | PIP-weighted VEP consequence evidence |
| `e2gMean` | PIP-weighted enhancer-to-gene evidence |

These fields are already bounded evidence features produced by the Open Targets
pipeline. The pinned matrices encode absent evidence as zero; the model rejects
missing or non-finite numeric inputs rather than silently imputing them. The
upstream zero convention does not imply that missing and biologically null
evidence are scientifically identical.

The API's external L2G score is retained for comparison and candidate provenance
but is excluded from the design matrix. The API exposes candidates at an
external score threshold of 0.05, so candidate ascertainment still depends on
the upstream system even though its aggregate score is not fitted as a feature.

## Conditional Bayesian model

Let `x_lg` be the six-feature vector for gene `g` at locus `l`. A locus has one
unit of target mass. When a training locus has multiple curated positives, that
mass is divided equally among them. The model is

```text
eta_lg = x_lg' beta
p_lg = exp(eta_lg) / sum_h exp(eta_lh)
beta_k ~ Normal(0, sigma^2), sigma = 0.75
```

The likelihood is the grouped categorical cross-entropy. There is no intercept,
because a common intercept cancels within a locus softmax and is unidentifiable.
All loci receive equal total likelihood weight regardless of candidate count.

L-BFGS-B finds the maximum a posteriori coefficient vector. The analytic Hessian
of the negative log posterior is inverted to form a multivariate-normal Laplace
approximation. The reported 5th and 95th percentiles come from 2,000 fixed-seed
coefficient draws in the configured analysis (the checked run manifest records
the exact override used for that run).

## Candidate-conditional and absolute probability

Conditional softmax probabilities sum to one over surfaced candidates. The
report additionally reserves mass for an omitted or unresolved gene:

```text
resolvable_l = (1 - rho) * C_l
P_absolute(g) = resolvable_l * p_lg
P_unresolved(l) = 1 - resolvable_l
```

where `C_l` is the summed PIP of the 95% credible-set records returned for the
locus, clipped to `[0, 1]`, and `rho = 0.20`. The value of `rho` is a declared
modeling choice. Both the conditional and absolute probabilities are exported.

## Held-out evaluation

Evaluation is performed on all 1,048 loci in the official test partition. The
metrics are top-1 accuracy, recall at three, mean reciprocal rank, multiclass
log loss, Brier score, and entropy. Tied scores use exact expectation over a
random ordering of the tied block; this prevents the uniform baseline from
benefiting from input row order.

Pre-specified comparisons are uniform allocation, distance only, eQTL-coloc
only, and VEP only. Scalar baseline scores are transformed into within-locus
probabilities with a fixed temperature of 12. No temperature is tuned on the
test set.

The observed BayesV2G-minus-distance top-1 difference is evaluated with a paired
locus bootstrap. The recorded analysis uses 10,000 fixed-seed resamples, keeping
all candidate rows from a sampled locus together.

## Sensitivity analysis

The complete model is refitted with prior standard deviations 0.5 and 1.5. It
is also refitted six times, omitting one evidence feature each time. These
comparisons use MAP probabilities to avoid posterior-draw Monte Carlo noise.
Unresolved masses of 0.10 and 0.40 are evaluated without refitting. For every
scenario and demonstration locus, the export records the top gene, whether it
changed, and the largest absolute probability change across candidates.

## Independent SORT1 fine-mapping audit

For each of the 5,772 GWAS variants with finite beta and standard error,
BayesV2G computes the Wakefield approximate Bayes factor for a quantitative
trait with prior variance `W = 0.04`:

```text
V_i = SE_i^2
r_i = W / (V_i + W)
log(ABF_i) = 0.5 * [log(1 - r_i) + r_i * (beta_i / SE_i)^2]
PIP_i = ABF_i / sum_j ABF_j
```

The smallest descending-PIP set reaching 95% is reported. This is a transparent
single-causal-variant approximation and ignores LD; it is not presented as a
replacement for SuSiE or conditional multi-signal fine-mapping.

For colocalization, biallelic SNVs are matched by GRCh38 position plus the
unordered allele pair. Duplicate QTL rows caused by multiple rsID mappings are
collapsed only when beta and standard error agree. For each molecular trait,
the classic five configurations are computed in log space:

```text
H0: neither trait is associated
H1: GWAS only
H2: QTL only
H3: both traits, different causal variants
H4: both traits, one shared causal variant
```

The per-variant priors are `p1 = 1e-4`, `p2 = 1e-4`, and `p12 = 1e-5`. The
primary true-effect prior variance is `W = 0.04` for both traits. The analysis
also exports PP(H4) under `W = 0.01` and `W = 0.16`; these are scale-sensitivity
checks, not evidence that one prior is biologically correct. Outputs include
PP(H0)–PP(H4), the top shared variant, and its conditional-H4 PIP. Effect signs
are not used by the ABF and therefore do not establish a directional mediation
relationship.

## Reproducibility

The random seed is fixed. All source snapshots and parquet inputs are identified
by SHA-256 in `run_manifest.json`, which also records package versions, platform,
model fit diagnostics, exclusions, and source release identifiers. Network calls
occur only when a refresh is explicitly requested. The run manifest also stores
per-file and aggregate SHA-256 hashes for the exact analysis source/configuration,
and `requirements.lock` captures the recorded Python 3.12 reference environment.
