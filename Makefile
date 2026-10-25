.DEFAULT_GOAL := help

PYTHON ?= python3
SPHINXOPTS ?=
DOCS_BUILD_DIR ?= docs/build
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))

.PHONY: help docs docs-install test quality check release-check tag release major minor patch

help:
	@echo "PyFiberModes development commands"
	@echo ""
	@echo "  make docs                  Build documentation and execute gallery examples"
	@echo "  make docs-install          Install the documentation dependencies"
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

docs-install:
	$(PYTHON) -m pip install -e ".[documentation]"

docs:
	MPLBACKEND=Agg $(PYTHON) -m sphinx -M html docs/source "$(DOCS_BUILD_DIR)" $(SPHINXOPTS)

test:
	MPLBACKEND=Agg $(PYTHON) -m pytest

quality:
	$(PYTHON) -m ruff check PyFiberModes tests tools
	$(PYTHON) -m mypy

release-check:
	$(PYTHON) -m martinpdes_maintenance.cli check $(if $(VERSION),--version $(VERSION),)

check: quality test release-check

tag:
	$(PYTHON) -m martinpdes_maintenance.cli tag "$(TAG_VERSION)"

release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	@set -eu; release_tag="$$($(PYTHON) -m martinpdes_maintenance.cli next $(RELEASE_KIND))"; \
	$(PYTHON) -m martinpdes_maintenance.cli tag "$$release_tag"; \
	git push origin HEAD "refs/tags/$$release_tag"

major minor patch:
	@:
