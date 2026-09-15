# fem-project developer & shipping targets.
# CMake artifacts live under build/; this Makefile owns the user-facing steps.

PYTHON := $(if $(wildcard .venv/bin/python),.venv/bin/python,python)
PIP    ?= $(PYTHON) -m pip

.PHONY: help build test benchmark lint format verify registry clean run-input

help: ## List available targets
	@printf '%s\n' 'fem-project targets:' \
	  '  make build      compile + install femcore and the fem package' \
	  '  make test       run the pytest suite (unit + integration)' \
	  '  make benchmark  run the analytical benchmark sweep' \
	  '  make lint       static checks (flake8 + black --check)' \
	  '  make format     format Python and C++ in place' \
	  '  make registry   print the registered element list' \
	  '  make verify     build + registry + tests + benchmarks (shipping gate)' \
	  '  make run-input  solve IN -> OUT (e.g. IN=tests/inputs/cantilever.inp OUT=out.vtu)' \
	  '  make clean      remove python build artifacts'

build: ## Compile the femcore C++ extension and install the package
	$(PIP) install . -v

test: ## Run the pytest suite (unit + integration, no slow nightly runs)
	$(PYTHON) -m pytest python/tests -m "not slow" -v

benchmark: ## Run the analytical benchmark sweep (docs/06 §2)
	$(PYTHON) tests/run_benchmarks.py

lint: ## Static checks: flake8 + black --check
	$(PYTHON) -m flake8 python/fem --max-line-length=100
	$(PYTHON) -m black --check --diff python/fem python/tests

format: ## Format Python (black) and C++ (clang-format) in place
	$(PYTHON) -m black python/fem python/tests
	clang-format -i cpp/include/*.hpp cpp/src/*.cpp cpp/bindings/*.cpp

registry: ## Print the registered element names (doc 05 §2.2 verification)
	$(PYTHON) -c "import femcore; print(femcore.registered_elements())"

verify: build registry test benchmark ## Full shipping gate: build -> registry -> tests -> benchmarks
	@echo "=== fem-project shipping gate: ALL GREEN ==="

run-input: ## Solve a deck: make run-input IN=path/to.deck.inp OUT=result.vtu
	$(PYTHON) -m fem.cli --in "$(IN)" --out "$(OUT)" --verbose

clean: ## Remove build artifacts
	rm -rf build dist *.egg-info python/femcore*.so python/build
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache