# Data provenance and retrieval

Only compact, redistributable snapshots and the pinned Open Targets L2G parquet
files are kept locally. Top-level inputs have SHA-256 digests in
`results/ldl_demo/run_manifest.json`; the hashed SORT1 source manifest contains
and the pipeline verifies the digest of every regional source file.

## Open Targets L2G 26.06

Repository: `opentargets/locus_to_gene_26.06` on Hugging Face
Pinned revision: `6abe9cc490222986925e1acca696d48c9b2437a4`

```bash
hf download opentargets/locus_to_gene_26.06 \
  --revision 6abe9cc490222986925e1acca696d48c9b2437a4 \
  --include 'train.parquet' 'test.parquet' 'config.json' 'README.md' \
  --local-dir data/raw/open_targets_l2g_26.06
```

The analysis expects the two parquet files at that exact path. The source
repository README and configuration remain beside them for provenance. These
data contain Open Targets' curated L2G labels and aggregate evidence features.

## Open Targets Platform LDL snapshot

`snapshots/open_targets_26.09_ldl_loci.json` is the raw GraphQL response for ten
configured `GCST90019512` study-locus IDs, captured from Platform release 26.09.
It includes every returned credible-set variant and candidate-gene feature, not
only the top-ranked records. The source is the UK Biobank study “Genetics of 35
blood and urine biomarkers,” with a reported sample size of 355,197. Refresh
explicitly with:

```bash
bayesv2g fetch --output data/snapshots/open_targets_26.09_ldl_loci.json
```

The live API can evolve even at a named release. A refreshed file will therefore
have a new hash and may produce different results.

## Independent SORT1 audit

`snapshots/sort1_case/` combines official public sources:

- GWAS Catalog `GCST90018961`, Sakaue et al. 2021 LDL meta-analysis;
- eQTL Catalogue release 7, uniformly reprocessed GTEx v8 liver eQTL
  `QTD000266` (208 samples);
- eQTL Catalogue Sun et al. plasma pQTL `QTD000584` (3,301 participants);
- GTEx v10 liver dynamic eQTL API records;
- Ensembl gene overlap, VEP, and 1000 Genomes phase 3 EUR LD endpoints.

The retrieval script reads the remote GWAS Tabix index and downloads only
GRCh38 `1:108774968-109774968`, yielding 5,772 rows rather than copying the
roughly 431 MB whole-genome file. Run:

```bash
python scripts/fetch_sort1_case_study.py
```

The source study combines 343,621 European-ancestry and 72,866 East
Asian-ancestry participants. The local audit does not claim that EUR-only LD is
an ancestry-matched reference for the meta-analysis.

The two QTL snapshots retain only SORT1, CELSR2, and PSRC1 rows in the same 1 Mb
window: 8,631 liver eQTL records and 3,892 plasma pQTL records. The pQTL dataset
contains a SORT1 aptamer but no CELSR2 or PSRC1 aptamer in this region.

The manifest records source URLs, retrieval time, build, row count, TLS status,
file sizes, and SHA-256 hashes. Donor-level GTEx expression/genotype arrays are
deliberately discarded; only public summary associations are stored.

## Builds, identifiers, and terms

- Genomic coordinates in this project are GRCh38 and 1-based.
- Candidate genes use Ensembl identifiers; version suffixes are not used in the
  primary L2G tables.
- The GWAS Catalog harmonized summary statistics are released under CC0.
- eQTL Catalogue summary statistics are available under CC BY 4.0.
- GTEx summary API use remains subject to GTEx portal terms.
- Ensembl REST data remain subject to Ensembl terms and attribution guidance.
- Open Targets and Hugging Face-hosted artifacts retain their upstream licenses
  and attribution requirements.

The MIT license in this repository covers the project code, not third-party
datasets. See [`THIRD_PARTY_DATA_NOTICES.md`](../THIRD_PARTY_DATA_NOTICES.md)
for file-level sources, licenses, acknowledgments, and transformations.
