# Contributing to GRsuite

Bug reports, models, criteria and examples are all welcome.

GRsuite is an unofficial, personal port of airGR, maintained on personal time and
with no institutional backing. Expect the review pace of a one-person project —
and please never route a GRsuite problem to the airGR maintainers.

## The one hard rule

**Anything that changes model numerics must keep the airGR comparison tests
green.**

GRsuite is a faithful port of airGR: its value comes from producing the same
numbers as the reference implementation. The test suite compares every model
output, criterion and calibrated parameter set against reference files produced
by airGR 1.7.9 itself (see `tests/data`). If your change makes one of those
tests fail, the change is wrong — even when it looks like an improvement.

Two subtleties you must not "fix" (see [docs/FIDELITY.md](docs/FIDELITY.md)
for the details):

- the model constants (`0.9`, `0.4`, …) are evaluated in single precision
  before promotion, exactly as airGR's Fortran does;
- the automatic warm-up period carries airGR's leap-year correction.

## Development setup

```bash
git clone https://github.com/Najim33/GRsuite
cd GRsuite
pip install -e ".[dev]"
pytest
```

## Style

- Code is linted with `ruff` (see `pyproject.toml`); run `ruff check src tests examples`.
- The numerical kernels (`src/grsuite/_kernels*.py`) are translated line by line
  from airGR's Fortran. Variable names and operation order are airGR's own —
  do not rename, reorder or "simplify" them.
- The rest of the package follows the existing naming: the airGR-mirroring
  API keeps airGR's names, lowercased (`run_model_gr4j`, `calibration_michel`),
  the high-level API (`api.py`) is idiomatic Python.
- Docstrings of the translated modules are in French (they mirror airGR's
  sources); user-facing documentation and examples are in English.

## Adding a model or a criterion

1. Translate the reference implementation, keeping the exact operation order.
2. Add the corresponding test in `tests/`, against a reference file exported
   from airGR itself — a test without an airGR reference is not accepted for
   numerics.
3. Document the mapping in [docs/MIGRATING_FROM_AIRGR.md](docs/MIGRATING_FROM_AIRGR.md).
