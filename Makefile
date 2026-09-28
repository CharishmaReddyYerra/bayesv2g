.PHONY: install data fetch-sort1 run test verify clean-results

PYTHON ?= python3
HF_REVISION := 6abe9cc490222986925e1acca696d48c9b2437a4

install:
	$(PYTHON) -m pip install -e '.[dev,data,download]'

data:
	hf download opentargets/locus_to_gene_26.06 --repo-type dataset --revision $(HF_REVISION) --include 'train.parquet' 'test.parquet' 'config.json' 'README.md' --local-dir data/raw/open_targets_l2g_26.06

fetch-sort1:
	$(PYTHON) scripts/fetch_sort1_case_study.py

run:
	$(PYTHON) -m bayesv2g run

test:
	$(PYTHON) -m pytest

verify: test run

clean-results:
	@echo "Results are provenance-bearing artifacts; remove results/ldl_demo manually if replacement is intended."
