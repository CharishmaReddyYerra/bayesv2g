# Output data dictionary

All tables are written to `results/ldl_demo/`. Identifiers and genomic positions
use GRCh38 unless stated otherwise.

## `prioritized_genes.csv`

The primary deliverable: one row per surfaced gene-locus pair.

| Column or family | Definition |
|---|---|
| `studyLocusId` | Stable Open Targets study-locus identifier |
| `locus_label` | Human-readable label from the analysis configuration |
| `lead_variant`, `chromosome`, `position` | Fine-mapped signal lead variant and coordinate |
| `geneId`, `gene_symbol`, `gene_name`, `gene_biotype` | Ensembl target identity and description |
| `gene_tss`, `distance_to_tss_bp` | Strand-aware TSS and absolute lead-to-TSS distance |
| `l2g_score_external` | Open Targets aggregate L2G score, retained only for comparison |
| `candidate_source` | Rule by which the candidate entered the table |
| `*ColocClpp*`, `*ColocH4*` | Molecular-QTL colocalization evidence from Open Targets |
| `distance*` | Distance evidence features from Open Targets |
| `vep*` | PIP-aggregated predicted-consequence evidence |
| `e2g*` | PIP-aggregated enhancer-to-gene evidence |
| `conditional_posterior_*` | Mean and 5th/95th percentiles conditional on surfaced genes |
| `pip_coverage_95` | Sum of PIPs flagged as the 95% credible set |
| `posterior_*` | Conditional posterior scaled by resolvable locus mass |
| `unresolved_probability` | Probability reserved outside surfaced candidates |
| `surfaced_candidate_count` | Number of externally prefiltered candidates at the locus |
| `candidate_set_is_singleton` | True when no surfaced gene competes with this candidate |
| `probability_scope` | Explicit statement of the candidate-conditional probability scope |
| `rank` | Descending absolute posterior within locus |
| `evidence_coverage` | Fraction of five non-distance evidence classes with value > 0 |
| `locus_entropy` | Entropy of the conditional mean probability vector |

## `candidate_evidence.csv`

The transformed Open Targets candidate feature matrix before model predictions.
It contains the identity, distance, candidate-provenance, and evidence columns
defined above, plus `locus_order` for report ordering.

## `loci.csv`

| Column | Definition |
|---|---|
| `study_id`, `trait` | Source GWAS accession and reported trait |
| `lead_variant_id` | Open Targets chromosome-position-allele identifier |
| `p_value`, `beta`, `standard_error`, `sample_size` | Source association metadata |
| `finemapping_method` | Fine-mapping method reported by Open Targets |
| `credible_set_confidence` | Open Targets quality/confidence label |
| `credible_set_size_95`, `credible_set_size_99` | Variant counts in flagged sets |
| `pip_coverage_95`, `pip_coverage_99` | PIP sum across flagged set records |
| `candidate_count` | Number of L2G candidates surfaced by the API |
| `publication_title`, `pubmed_id` | Study provenance |

## `credible_set_variants.csv`

One row per Open Targets fine-mapped variant. It records alleles, rsID, PIP,
95%/99% membership, log Bayes factor, and most severe predicted consequence.

## `posterior_coefficients.csv`

| Column | Definition |
|---|---|
| `feature` | Model feature name |
| `posterior_mean` | MAP coefficient, also the Laplace posterior mean |
| `posterior_sd` | Square root of the Laplace covariance diagonal |
| `posterior_lower`, `posterior_upper` | Normal-approximation 5th/95th percentiles |
| `odds_ratio` | `exp(posterior_mean)` for a full one-unit feature contrast |

Because features are correlated and differently distributed, these coefficients
must not be interpreted as causal effects.

## `validation_metrics.csv`

One row for BayesV2G and each baseline. `n_loci` is the number of held-out choice
sets. Ranking columns are `top1_accuracy`, `recall_at_3`, and
`mean_reciprocal_rank`; probability-quality columns are
`multiclass_log_loss`, `brier_score`, and `mean_entropy`.

## `paired_top1_comparison.csv`

The BayesV2G-minus-distance top-1 difference and its paired 95% locus-bootstrap
interval, plus counts of improved, worsened, and tied loci and the number of
bootstrap resamples.

## `calibration.csv`

Top-prediction calibration bins with `mean_confidence`, `observed_accuracy`, and
the number of loci `n`. Empty bins are not emitted.

## `sensitivity.csv`

One row per analysis scenario and demonstration locus. `scenario_type`
distinguishes MAP model refits from probability-assumption changes. `top_gene_changed`
compares each refit with the main model; `max_probability_change` is the largest
candidate-level absolute shift at that locus.

## `sort1_finemapping_all.csv.gz` and `sort1_finemapping_top25.csv`

The raw GWAS columns are followed by `log_approximate_bayes_factor`,
`approximate_pip`, `in_95_credible_set`, and `pip_rank`. Extremely small source
p-values are preserved as text; beta and standard error drive the ABF.

## `sort1_gene_evidence.csv`

GTEx v10 liver association fields for SORT1, CELSR2, and PSRC1 joined to their
BayesV2G rank, posterior, external comparison score, and selected model evidence.

## `sort1_colocalization.csv`

One row per regional molecular trait. It records dataset, evidence type, gene,
number of allele-matched SNVs, PP(H0)–PP(H4), PP(H4) sensitivity at effect-prior
variances 0.01 and 0.16, the top shared variant, its PIP conditional on H4, and
the 95% shared-variant credible-set size.

## `sort1_coloc_shared_variants.csv.gz`

Variant-level inputs and posteriors for every regional coloc analysis: aligned
GWAS/QTL beta and SE, single-signal trait PIPs, shared-variant PIP conditional on
H4, unconditional shared posterior, and shared 95% credible-set membership.

## `model.json`

A portable fitted artifact containing model type/version, feature order,
hyperparameters, MAP coefficients, Hessian, covariance, and optimizer metadata.

## `run_manifest.json`

The reproducibility record: UTC timestamp, resolved configuration, software and
platform versions, data releases, SHA-256 input hashes, training exclusions,
fit diagnostics, and raw SORT1 audit statistics.

## `SHA256SUMS`

Standard SHA-256 lines for every generated table, model, figure, report, and run
manifest, using paths relative to the output directory. Verify from that
directory with `sha256sum -c SHA256SUMS`.
