# Changelog

All notable changes to GRsuite are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-08-15

First public release.

### Added

- The full airGR 1.7.9 model suite, translated from the Fortran and R sources
  and validated bit-for-bit against airGR: GR4J, GR5J, GR6J (daily), GR4H and
  GR5H with or without the interception store (hourly), GR2M (monthly), GR1A
  (annual).
- The CemaNeige snow module, with or without linear hysteresis, coupled to
  every daily and hourly model, with Valéry elevation-band extrapolation and
  daily temperature gradients.
- Michel's calibration algorithm (grid screening then steepest descent),
  reproducing airGR's iteration counts exactly.
- Error criteria: NSE, KGE, KGE′, RMSE, weighted composites, and the
  `sqrt`, `log`, `inv`, `sort`, `boxcox` and power transformations.
- Utilities: Oudin potential evapotranspiration, interception capacity
  estimation, time-series aggregation, semi-distributed routing
  (`run_model_lag`).
- A high-level API (`Catchment`, `calibrate`, `split_sample`, `evaluate`)
  on top of the airGR-mirroring one.
- A test suite of 159 tests comparing GRsuite to reference outputs produced by
  airGR 1.7.9, shipped as gzipped fixtures under `tests/data`.
- Validation material: campaign over 100 CAMELS-FR catchments, reported in
  [docs/VALIDATION.md](docs/VALIDATION.md), and the numerical fidelity notes
  in [docs/FIDELITY.md](docs/FIDELITY.md).
- `docs/make_validation_tables.py`: recomputes the whole airGR comparison from the
  committed reference files and writes one row per compared quantity to
  `docs/data/series_deviations.csv` — 502 quantities, 1 230 966 values, worst
  deviation 2.0 × 10⁻¹³. Section 2.1 of the validation report is its output, so
  the figures quoted there can be audited rather than trusted.
- [docs/REFERENCES.md](docs/REFERENCES.md): the published source of every model,
  snow routine, criterion, transformation and algorithm implemented here, with
  DOIs, plus the reference data used for validation. The READMEs, `NOTICE`,
  `CITATION.cff` and the module docstrings point to it.
