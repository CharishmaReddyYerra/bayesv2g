# Third-party data notices

The repository's MIT license applies to BayesV2G code and original
documentation only. The files below retain their upstream licenses, terms,
and attribution requirements. Source snapshots are included to make the
published analysis reproducible; they are not relicensed by this project.

## Open Targets

- `data/snapshots/open_targets_26.09_ldl_loci.json` is an Open Targets
  Platform 26.09 GraphQL response captured on 2026-09-26. Platform data are
  made available under [CC0 1.0](https://platform-docs.opentargets.org/licence),
  subject to the [Open Targets terms of use](https://platform-docs.opentargets.org/licence/terms-of-use).
  Citation guidance is available from [Open Targets](https://platform-docs.opentargets.org/citation).
- `data/raw/open_targets_l2g_26.06/` is from the Open Targets
  [`locus_to_gene_26.06`](https://huggingface.co/opentargets/locus_to_gene_26.06)
  Hugging Face model repository at revision
  `6abe9cc490222986925e1acca696d48c9b2437a4`, whose repository card identifies
  the license as MIT. The bundled upstream README is retained.

Please cite Buniello et al., *Nucleic Acids Research* (2025),
[doi:10.1093/nar/gkae1128](https://doi.org/10.1093/nar/gkae1128), and Ghoussaini
et al., *Nature Genetics* 53, 1527–1533 (2021),
[doi:10.1038/s41588-021-00945-5](https://doi.org/10.1038/s41588-021-00945-5),
as appropriate.

## NHGRI-EBI GWAS Catalog

`data/snapshots/sort1_case/gwas_gcst90018961_*` contains a regional,
harmonised subset and metadata for study `GCST90018961`, downloaded on
2026-09-26 and restricted to GRCh38 `1:108774968-109774968`. The study's
[official API record](https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90018961)
identifies its terms as CC0. See the [GWAS Catalog citation and reuse guidance](https://www.ebi.ac.uk/gwas/docs/faq/).

Please cite Cerezo et al., *Nucleic Acids Research* 53(D1), D998–D1005
(2025), [doi:10.1093/nar/gkae1070](https://doi.org/10.1093/nar/gkae1070), and
Sakaue et al., *Nature Genetics* 53, 1415–1424 (2021),
[doi:10.1038/s41588-021-00931-x](https://doi.org/10.1038/s41588-021-00931-x).

## eQTL Catalogue

`eqtl_catalogue_qtd000266_liver_sort1_genes.tsv.gz` and
`eqtl_catalogue_qtd000584_plasma_sort1_genes.tsv.gz` are from eQTL Catalogue
release 7 datasets `QTD000266` and `QTD000584`. They are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) according to the
[eQTL Catalogue license](https://www.ebi.ac.uk/eqtl/License/). The upstream
files were modified here by retaining only GRCh38
`1:108774968-109774968` and rows for SORT1, CELSR2, and PSRC1. Study metadata
are available from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Studies/).

Please cite Kerimov et al., *Nature Genetics* 53, 1290–1299 (2021),
[doi:10.1038/s41588-021-00924-w](https://doi.org/10.1038/s41588-021-00924-w),
and Kerimov et al., *PLOS Genetics* 19, e1010932 (2023),
[doi:10.1371/journal.pgen.1010932](https://doi.org/10.1371/journal.pgen.1010932).
For the underlying studies, cite the GTEx Consortium, *Science* 369,
1318–1330 (2020), [doi:10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776),
and Sun et al., *Nature* 558, 73–79 (2018),
[doi:10.1038/s41586-018-0175-2](https://doi.org/10.1038/s41586-018-0175-2).

## GTEx Portal

`gtex_variant_lookup_rs12740374.json` and
`gtex_v10_liver_rs12740374_eqtl.json` are frozen GTEx v10 API responses
captured on 2026-09-26. Donor-level expression and genotype arrays were not
retained. Redistribution is governed by the [GTEx Portal license and data-use notice](https://gtexportal.org/home/license).

> The data used for the analyses described in this repository were obtained
> from the GTEx Portal on 09/26/2026.

These frozen responses may not reflect the most current or accurate data
available from the GTEx Portal. GTEx does not endorse this project.

## Ensembl and 1000 Genomes

`ensembl_protein_coding_genes.json`, `ensembl_vep_rs12740374.json`, and
`ensembl_1000g_eur_ld_rs12740374.json` were retrieved through the Ensembl REST
API on 2026-09-26 using Ensembl release 116-era services. They remain subject
to the [Ensembl disclaimer and data terms](https://www.ensembl.org/info/about/legal/disclaimer.html)
and [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use/).
Third-party `colocated_variants` records returned by VEP were deliberately
omitted because they can include redistribution-restricted HGMD or COSMIC
content; the project uses only Ensembl consequence annotations. See Ensembl's
[citation guidance](https://jun2026.archive.ensembl.org/info/about/publications.html).

Please cite Dyer et al., “Ensembl 2025,” *Nucleic Acids Research* 53(D1),
D948–D957, [doi:10.1093/nar/gkae1071](https://doi.org/10.1093/nar/gkae1071).
The LD response derives from 1000 Genomes phase 3; see the
[IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer/) and
cite Auton et al., *Nature* 526, 68–74 (2015),
[doi:10.1038/nature15393](https://doi.org/10.1038/nature15393).
