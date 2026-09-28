# Limitations and responsible interpretation

BayesV2G is a transparent research prototype. Its posterior probabilities are
useful for prioritization, but they are conditional on assumptions that should
remain visible in any presentation or manuscript.

## Candidate ascertainment

The ten-locus candidate set is supplied by Open Targets at an external L2G score
threshold of 0.05. The aggregate score is not a model feature, but it influences
which genes are eligible. Genes absent from that set cannot win. The explicit
unresolved category acknowledges this problem but cannot identify which omitted
gene is responsible. The train/test matrices contain about 16 candidates per
locus on average, whereas the demonstration exposes 2.6 (range 1–5). Conditional
softmax probabilities are not invariant to that shift. APOB and LDLR each have
one surfaced candidate, so their conditional probability of 1.0 is structural,
not evidence-derived certainty.

## Reference-label circularity

The public gold-standard labels reflect curated biological knowledge and may
partly use coding, proximity, QTL, or related evidence. A model trained and
tested with those evidence classes can therefore look better than it would
against independent perturbation experiments. The official held-out partition
prevents row leakage; it does not eliminate conceptual label circularity.

## Evidence dependence and missingness

Distance, VEP, enhancer-to-gene, and molecular-QTL summaries are correlated.
Coefficient magnitudes are predictive weights, not causal effect sizes or an
independent measure of biological importance. Maximum QTL values collapse many
tissues and studies. A zero can mean no signal, no relevant tissue, low power,
or unavailable data, yet the current matrix cannot distinguish those states.
The Sun et al. pQTL uses an aptamer assay; protein-altering variants can change
aptamer binding rather than circulating protein abundance.

## Fine-mapping and ancestry

The primary credible sets are Open Targets SuSiE-inf products based on external
LD. Their accuracy depends on LD reference, ancestry match, imputation, and model
assumptions. The independent `GCST90018961` audit combines European- and East
Asian-ancestry cohorts; a single-population LD matrix would not fully represent
that meta-analysis. The audit's Wakefield calculation assumes one causal variant
and ignores LD, so its PIP should be treated as a sanity check only.
The regional coloc audit shares the same one-causal-variant assumption. Very high
PP(H4) at this locus can reflect a common sharp lead signal and should not be
read as proof that all three eQTL genes mediate LDL biology.
The shared-effect prior is not learned from these data; exported values at three
prior variances show numerical sensitivity but do not solve trait-scale choice.

## One primary target distribution

The categorical likelihood distributes one unit of probability over a locus.
This is coherent for ranking, but biology can involve multiple mediating genes,
pleiotropic mechanisms, or context-specific effects. Multiple accepted training
genes share the target mass, which is a pragmatic labeling rule rather than a
mechanistic claim.

## Posterior approximation

The Laplace approximation assumes that the coefficient posterior is adequately
represented by a Gaussian near its MAP. It is fast and reproducible but does not
capture skew, heavy tails, or multimodality. A future hierarchical implementation
could use MCMC or variational inference and model tissue-level effects directly.

## Unresolved mass

The outside-candidate parameter `rho = 0.20` is selected, not learned. Absolute
posterior values will change with this choice and with credible-set coverage.
Rank order within a locus does not change because the same scale factor is
applied to every surfaced candidate.

## Validation scope

The current benchmark is cross-trait Open Targets L2G data, not prospective wet
lab validation. The ten LDL examples are illustrative and not an independent
accuracy estimate. The ANGPTL3 locus, where DOCK7 ranks above ANGPTL3, is a useful
face-validity warning: the model may emphasize learned feature patterns over a
well-known gene. Disagreements should trigger evidence review, not be hidden.

## Clinical scope

The output must not be used for diagnosis, treatment, or patient-level decisions.
It does not establish that perturbing a gene is safe, efficacious, or directionally
appropriate. Therapeutic translation needs independent functional, phenotypic,
and safety evidence.
