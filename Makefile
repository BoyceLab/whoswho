.PHONY: fetch build test release
fetch:
	python -m gene_spine.fetch
build:
	python -m gene_spine.build
test:
	python -m pytest -q tests
release: fetch build
	@echo "See outputs/release.json"
