.PHONY: fetch build lists test release
fetch:
	python -m gene_spine.fetch
build:
	python -m gene_spine.build
# The derived lists read the spine the build wrote, so they run after it, never instead of it.
lists:
	python -m gene_spine.epilepsy_list
test:
	python -m pytest -q tests
release: fetch build lists
	@echo "See outputs/release.json"
