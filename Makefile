.PHONY: install test smoke features evaluate clean lint

install:
	pip install -e ".[dev,full]"

test:
	python -m pytest tests/ -q

smoke:
	python scripts/00_smoke_test.py

features:
	python scripts/01_extract_features.py --samples config/samples.tsv \
		--outdir results/features --fasta data/ref/hg38.fa

evaluate:
	python scripts/02_train_evaluate.py --features results/features

lint:
	ruff check src/ scripts/ tests/

clean:
	rm -rf results/smoke results/cache data/simulated .pytest_cache
	find . -name __pycache__ -type d -exec rm -rf {} +
