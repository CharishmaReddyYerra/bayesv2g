# Project summary

## Research question

Given a fine-mapped LDL-cholesterol GWAS signal, which nearby gene is most likely
to mediate the association after integrating proximity, molecular-QTL
colocalization, predicted consequence, and enhancer-to-gene evidence?

## Contribution

BayesV2G turns six evidence streams into a coherent within-locus posterior
distribution. It adds three elements that are often absent from demonstration
projects: an explicit unresolved category, held-out grouped evaluation against
honest baselines, and a complete provenance trail down to source-file hashes.

The implementation is intentionally lightweight. Conditional Bayesian softmax
regression, analytic derivatives, and a Laplace posterior make the full 99,000-row
fit practical on a laptop without hiding the statistical assumptions behind a
large probabilistic-programming stack.

## Empirical result

On 1,048 official held-out loci, top-1 accuracy is 0.731 versus 0.658 for the
distance-only baseline. The paired difference is 0.073 (95% locus-bootstrap
interval 0.046–0.102). Recall@3 is 0.860 versus 0.803. These numbers show useful
incremental ranking information, not proof of causal-gene discovery, because the
curated labels may share evidence with the predictors.

Ten real LDL signals provide interpretable case studies. ABCG8, HMGCR, LPA,
NPC1L1, LDLR, APOE, APOB, and PCSK9-related candidates illustrate different
evidence combinations. The SORT1 locus demonstrates why the problem is hard:
the lead variant has strong liver eQTL effects for multiple genes. The ANGPTL3
locus supplies an equally important negative example because the learned model
ranks DOCK7 first despite strong prior biological knowledge about ANGPTL3.

## Defensible presentation framing

The strongest accurate claim is:

> BayesV2G is a reproducible probabilistic prioritization model that integrates
> real fine-mapped and functional-genomics evidence, improves held-out ranking
> over proximity alone, and exposes rather than suppresses uncertainty and model
> failures.

Avoid claiming that the posterior is a direct probability of biological truth,
that held-out labels are independent experimental validation, or that the raw
SORT1 ABF calculation is full LD-aware fine-mapping.

## Logical next extensions

1. Replace maximum QTL summaries with tissue-specific latent effects and partial
   pooling across liver, adipose, artery, heart, and blood.
2. Add explicit missingness indicators and measurement-error models for each
   evidence source.
3. Learn the unresolved mass from orthogonal perturbation data instead of fixing
   it globally.
4. Validate prospectively against CRISPR perturbation, rare-variant burden, and
   independent drug-target evidence withheld from label creation.
5. Use ancestry-matched, multi-signal fine-mapping and propagate PIP uncertainty
   into the gene model rather than consuming fixed aggregate features.
