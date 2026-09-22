.DEFAULT_GOAL := help

PYTHON ?= python3
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))

.PHONY: help test quality check release-check tag release major minor patch

help:
	@echo "PyFiberModes development commands"
	@echo ""
	@echo "  make test                  Run the test suite"
	@echo "  make quality               Run static quality checks"
	@echo "  make check                 Run tests, quality, and release checks"
	@echo "  make release-check         Check version metadata consistency"
	@echo "  make tag VERSION=vX.Y.Z    Create a release commit and annotated tag"
	@echo "  make release patch         Create and push the next patch release"
	@echo "  make release minor         Create and push the next minor release"
	@echo "  make release major         Create and push the next major release"

ifneq ($(filter tag,$(MAKECMDGOALS)),)
ifneq ($(strip $(TAG_VERSION)),)
.PHONY: $(TAG_VERSION)
$(TAG_VERSION):
	@:
endif
endif

test:
	MPLBACKEND=Agg $(PYTHON) -m pytest

quality:
	$(PYTHON) -m ruff check PyFiberModes tests tools

release-check:
	$(PYTHON) tools/check_release.py $(if $(VERSION),--version $(VERSION),)

check: quality test release-check

tag:
	$(PYTHON) tools/release_tag.py "$(TAG_VERSION)"

release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	@set -eu; release_tag="$$($(PYTHON) tools/next_release_version.py $(RELEASE_KIND))"; \
	$(PYTHON) tools/release_tag.py "$$release_tag"; \
	git push origin HEAD "refs/tags/$$release_tag"

major minor patch:
	@:
