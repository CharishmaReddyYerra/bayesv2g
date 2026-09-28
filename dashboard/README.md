# BayesV2G Evidence Workbench

Interactive dashboard for the Bayesian variant-to-gene prioritization project. It renders the checked-in LDL demonstration outputs from `../results/ldl_demo` and does not recompute model results in the browser.

## Included views

- ten-locus explorer with absolute and within-candidate posterior modes
- ranked gene hypotheses with 90% intervals and explicit unresolved mass
- gene-level evidence inspector and credible-set variant track
- held-out model validation, calibration, and paired bootstrap comparison
- independent SORT1 fine-mapping and colocalization audit
- source ledger, probability construction, and interpretive caveats

## Real-data snapshot

The displayed outputs use Open Targets L2G, GWAS Catalog, eQTL Catalogue/GTEx, and Ensembl-derived evidence. The analysis run manifest and checksums live in `../results/ldl_demo`.

## Local development

```sh
npm ci
npm run dev
```

The local dashboard runs at `http://localhost:5173/` by default.

## Refresh dashboard data

After regenerating the BayesV2G result tables:

```sh
node scripts/generate-dashboard-data.mjs
```

Then verify the site with:

```sh
npm run lint
npm run build
```

This dashboard is for research prioritization only and is not clinical guidance.
