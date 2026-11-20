# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
#
# SCUBAS developer entry points.  `make` (or `make help`) lists them.
# Override the interpreter with  make PY=/path/to/python  (default: the conda env named "scubas").
# Ruff is a dev tool: RUFF defaults to "$(PY) -m ruff"; point it at any ruff binary with  make RUFF=/path/to/ruff lint.

ENV      ?= scubas
CONDA    ?= conda
PY       ?= $(shell $(CONDA) run -n $(ENV) which python 2>/dev/null || echo python)
RUFF     ?= $(PY) -m ruff
export RUFF
SRC      = scubas tests scripts examples

ifeq ($(NO_COLOR),)
  C0 := \033[0m
  CB := \033[1m
  CG := \033[38;5;78m
  CC := \033[38;5;45m
  CD := \033[2m
  CO := \033[38;5;214m
endif

.DEFAULT_GOAL := help
.PHONY: env build check logo codecov-on codecov-off codecov-status publish-on publish-off publish-status publish-tag rtd-on rtd-off rtd-status help install dev test cov cov-xml cov-floors lint lint-ratchet fmt format fmt-check spdx release-check changelog-check compat-baseline pre-commit-install pre-commit qa docs docs-serve clean

help: ## Show this list
	@printf "$(CB)SCUBAS$(C0)  $(CD)Submarine Cables Upset by Auroral Streams$(C0)\n\n"
	@printf "$(CB)Usage:$(C0) make $(CC)<target>$(C0)\n\n"
	@awk 'BEGIN {FS = ":.*##"} /^[a-zA-Z0-9_-]+:.*##/ {printf "  \033[38;5;45m%-16s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@printf "\n$(CD)python: $(PY)$(C0)\n"

env: ## Create the conda environment "scubas" (environment.yml; then `make dev`)
	@printf "$(CG)creating conda env '$(ENV)'$(C0)\n"
	$(CONDA) env create -f environment.yml -y

install: ## Install the package (runtime dependencies only)
	$(PY) -m pip install .

dev: ## Editable install with the dev extras (pytest, coverage, ruff, pre-commit, build, twine)
	$(PY) -m pip install -e ".[dev]"

test: ## Run the tests (the compatibility suite in tests/compat is part of them)
	@printf "$(CG)pytest$(C0)\n"
	$(PY) -m pytest

cov: ## Tests with a coverage report
	$(PY) -m pytest --cov=scubas --cov-report=term-missing

cov-xml: ## Tests with coverage.xml, the overall gate, then the per-module coverage floors
	$(PY) -m pytest --cov=scubas --cov-report=xml --cov-fail-under=90 -q
	$(PY) scripts/check_coverage_floors.py

cov-floors: ## Check the per-module coverage floors against an existing coverage.xml
	$(PY) scripts/check_coverage_floors.py

lint: ## Lint with ruff (the legacy debt in ruff.toml is allowed, anything new is not)
	@printf "$(CC)ruff check$(C0)\n"
	$(RUFF) check $(SRC)

lint-ratchet: ## The legacy lint debt is exactly what ruff.toml lists (it may only shrink)
	$(PY) scripts/ruff_ratchet.py --check

fmt: ## Auto-fix what ruff can fix in files that are not legacy debt (review the diff!)
	$(RUFF) check --fix $(SRC)

format: ## Format the code with ruff format (width 88, legacy exclusions in ruff.toml)
	$(RUFF) format $(SRC)

fmt-check: ## Check the formatting without changing files
	@printf "$(CC)ruff format --check$(C0)\n"
	$(RUFF) format --check $(SRC)

spdx: ## Every Python file has the SPDX header, except the listed legacy files
	$(PY) scripts/check_spdx.py

release-check: ## The version in the code, CITATION.cff (when present) and CHANGELOG.md agree (add TAG=vX.Y.Z to compare a tag)
	$(PY) scripts/check_release.py $(if $(TAG),--tag $(TAG),)

changelog-check: ## CHANGELOG.md must change with scubas/ (staged files)
	$(PY) scripts/check_changelog.py

compat-baseline: ## Regenerate tests/compat/baseline/ (ONLY on unchanged legacy code, or for an approved change of a number or signature)
	$(PY) scripts/gen_compat_baseline.py

pre-commit-install: ## Install the git hooks (pre-commit and pre-push)
	$(PY) -m pre_commit install --hook-type pre-commit --hook-type pre-push

pre-commit: ## Run every pre-commit hook on all tracked and staged files
	$(PY) -m pre_commit run --all-files

build: ## Build the sdist and the wheel into dist/
	$(PY) -m build

check: build ## Check the built distributions (twine, and that they hold what they should)
	$(PY) -m twine check dist/*
	$(PY) scripts/check_dist.py

docs: ## Build the documentation (Sphinx, dirhtml, warnings are errors)
	$(PY) -m sphinx -W --keep-going -b dirhtml docs docs/_build/dirhtml

REPO     ?= shibaji7/SCUBAS
VARAPI   = repos/$(REPO)/actions/variables

# setvar NAME VALUE: create or update a GitHub Actions repository variable (needs gh and admin rights; these targets change the repository's settings, so run them on purpose)
setvar = gh api -X PATCH $(VARAPI)/$(1) -f value=$(2) >/dev/null 2>&1 || gh api -X POST $(VARAPI) -f name=$(1) -f value=$(2) >/dev/null
# getvar NAME: the variable's value, or "not set"
getvar = $$(gh api $(VARAPI)/$(1) --jq .value 2>/dev/null || echo "not set")

codecov-off: ## Turn the Codecov upload in CI off (variable CODECOV_UPLOAD=false)
	@$(call setvar,CODECOV_UPLOAD,false)
	@printf "$(CO)CODECOV_UPLOAD=false$(C0) on $(REPO)\n"

codecov-on: ## Turn the Codecov upload in CI on again (variable CODECOV_UPLOAD=true; unset also means on)
	@$(call setvar,CODECOV_UPLOAD,true)
	@printf "$(CG)CODECOV_UPLOAD=true$(C0) on $(REPO)\n"

codecov-status: ## Show whether the Codecov upload is on (unset or anything but false: on)
	@echo "CODECOV_UPLOAD=$(call getvar,CODECOV_UPLOAD)  (unset or anything but false: the upload runs)"

publish-on: ## Let the release workflow publish tags to TestPyPI, PyPI and GitHub Releases (variable PUBLISH_ENABLED=true)
	@$(call setvar,PUBLISH_ENABLED,true)
	@printf "$(CG)PUBLISH_ENABLED=true$(C0) on $(REPO): a pushed tag vX.Y.Z now publishes\n"

publish-off: ## Stop the release workflow from publishing (tags are still checked and built, nothing is uploaded)
	@$(call setvar,PUBLISH_ENABLED,false)
	@printf "$(CO)PUBLISH_ENABLED=false$(C0) on $(REPO)\n"

publish-status: ## Show whether publishing is on (the variable PUBLISH_ENABLED)
	@echo "PUBLISH_ENABLED=$(call getvar,PUBLISH_ENABLED)  (unset or false: nothing is published)"

publish-tag: ## Publish one existing tag that was held back: make publish-tag TAG=v0.2.0 (run publish-on first)
	@test -n "$(TAG)" || { echo "usage: make publish-tag TAG=vX.Y.Z"; exit 1; }
	@git rev-parse -q --verify refs/tags/$(TAG) >/dev/null || { echo "no local tag $(TAG)"; exit 1; }
	gh workflow run release.yml --ref $(TAG) -R $(REPO)

logo: ## Regenerate the logo SVGs in docs/_static/ (needs fonttools and matplotlib)
	$(PY) scripts/make_logo.py

rtd-off: ## Stop Read the Docs from publishing: creates .rtd-skip (builds are cancelled after the next push)
	@touch .rtd-skip
	@printf "$(CC)Read the Docs builds off$(C0): .rtd-skip created, commit and push it\n"

rtd-on: ## Let Read the Docs publish again: removes .rtd-skip (takes effect with the next push)
	@rm -f .rtd-skip
	@printf "$(CG)Read the Docs builds on$(C0): .rtd-skip removed, commit and push\n"

rtd-status: ## Show whether Read the Docs builds are switched off (.rtd-skip present)
	@if [ -f .rtd-skip ]; then echo "Read the Docs builds are OFF (.rtd-skip exists)"; else echo "Read the Docs builds are ON"; fi

docs-serve: docs ## Build the documentation and serve it on http://localhost:8000
	$(PY) -m http.server 8000 --directory docs/_build/dirhtml

qa: lint lint-ratchet fmt-check spdx release-check cov-xml docs ## Everything that should pass before a commit
	@printf "$(CG)$(CB)qa passed$(C0)\n"

clean: ## Remove caches and build output
	rm -rf build dist .pytest_cache .ruff_cache .coverage coverage.xml htmlcov site docs/_build scubas.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
