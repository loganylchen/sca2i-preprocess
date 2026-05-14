.PHONY: help dryrun lint test clean clean-all envs dag report

SNAKEMAKE ?= snakemake
PROFILE   ?= profiles/default
JOBS      ?= 4

help:
	@echo "sca2i-preprocess — Snakemake workflow"
	@echo ""
	@echo "Targets:"
	@echo "  make dryrun     - snakemake -n (dry-run, prints jobs)"
	@echo "  make lint       - snakemake --lint"
	@echo "  make test       - snakemake --use-conda -j\$$JOBS (default 4) on test config"
	@echo "  make envs       - pre-create all conda envs without running rules"
	@echo "  make dag        - render DAG to dag.svg"
	@echo "  make report     - generate snakemake report.html"
	@echo "  make clean      - remove results/ logs/ benchmarks/"
	@echo "  make clean-all  - clean + remove .snakemake/ resources/"
	@echo ""
	@echo "Variables: SNAKEMAKE=$(SNAKEMAKE) PROFILE=$(PROFILE) JOBS=$(JOBS)"

dryrun:
	$(SNAKEMAKE) -n --profile $(PROFILE)

lint:
	$(SNAKEMAKE) --lint

envs:
	$(SNAKEMAKE) --use-conda --conda-create-envs-only -j$(JOBS)

test:
	$(SNAKEMAKE) --use-conda -j$(JOBS) --profile $(PROFILE)

dag:
	$(SNAKEMAKE) --dag | dot -Tsvg > dag.svg
	@echo "Wrote dag.svg"

report:
	$(SNAKEMAKE) --report report.html

clean:
	rm -rf results/ logs/ benchmarks/

clean-all: clean
	rm -rf .snakemake/ resources/
