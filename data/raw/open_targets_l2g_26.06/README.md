---
language: en
library_name: sklearn
license: mit
tags:
- sklearn
- tabular-classification
- genomics
- gwas
- gene-prioritization
---
# Locus-to-Gene (L2G) Model

The locus-to-gene (L2G) model prioritises likely causal genes at each GWAS locus based on genetic and functional genomics features.

## Model Description

This is a **Gradient Boosting Classifier** (XGBoost) trained to predict causal genes at GWAS loci.

Limited to protein-coding genes with available feature data.

**Key Features:**
- **Distance**: proximity from credible set variants to gene
- **Molecular QTL Colocalization**: evidence from expression/protein QTL studies
- **Variant Pathogenicity**: VEP (Variant Effect Predictor) scores

## Usage

```python
from gentropy.method.l2g.model import LocusToGeneModel
from gentropy.common.session import Session

# Load model from Hugging Face Hub
session = Session()
model = LocusToGeneModel.load_from_hub(
    session=session,
    hf_model_id="opentargets/locus_to_gene"
)

# Make predictions on your L2G feature matrix
predictions = model.predict(your_feature_matrix, session)
```

## Training

- **Architecture**: XGBoost Gradient Boosting Classifier
- **Training Data**: Curated positive/negative gene-locus pairs from Open Targets
- **Evaluation Metric**: Area under precision-recall curve (AUCPR)

## Citation

If you use this model, please cite:

```bibtex
@article{ghoussaini2021open,
title={Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics},
author={Ghoussaini, Maya and Mountjoy, Edward and Carmona, Maria and others},
journal={Nature Genetics},
volume={53},
pages={1527--1533},
year={2021},
doi={10.1038/s41588-021-00945-5}
}
```

## More Information

- **Repository**: [opentargets/gentropy](https://github.com/opentargets/gentropy)
- **Documentation**: [L2G Method Docs](https://opentargets.github.io/gentropy/python_api/methods/l2g/_l2g/)
- **Developer**: Open Targets
