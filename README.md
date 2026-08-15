# GRsuite

**An unofficial, independent Python port of the airGR suite of GR rainfall-runoff
models — a personal project, checked value by value against airGR.**

> **This is not an official package, and it is not airGR.** GRsuite is one
> person's work, developed independently and on personal time. It is not
> produced, endorsed, reviewed or supported by INRAE, and it carries no
> institutional guarantee. The models, the science and the reference
> implementation belong to INRAE's HYCAR unit; what is mine is the translation
> into Python and the mistakes in it.
>
> If your work needs the reference implementation, use
> **[airGR](https://cran.r-project.org/package=airGR)** — and please report
> problems found here to [this repository's issue
> tracker](https://github.com/Najim33/GRsuite/issues), never to the airGR
> maintainers.

[![unofficial port](https://img.shields.io/badge/status-unofficial%20port%20of%20airGR-c2521f.svg)](https://github.com/Najim33/GRsuite/blob/main/NOTICE)
[![tests](https://github.com/Najim33/GRsuite/actions/workflows/ci.yml/badge.svg)](https://github.com/Najim33/GRsuite/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/grsuite.svg)](https://pypi.org/project/grsuite/)
[![Python](https://img.shields.io/pypi/pyversions/grsuite.svg)](https://pypi.org/project/grsuite/)
[![License: GPL v2](https://img.shields.io/badge/license-GPL--2.0--or--later-blue.svg)](https://github.com/Najim33/GRsuite/blob/main/LICENSE)
[![validated against airGR 1.7.9](https://img.shields.io/badge/validated-airGR%201.7.9-0a7ea4.svg)](https://github.com/Najim33/GRsuite/blob/main/docs/VALIDATION.md)

GR4J, GR5J, GR6J, GR2M, GR1A, GR4H, GR5H and the CemaNeige snow module — the models
hydrologists have been running in R for twenty years, now in Python with no R
dependency, no Fortran toolchain, and no reimplementation risk.

Every model output, every criterion and every calibrated parameter set in this
package has been checked against airGR 1.7.9 itself: **502 quantities — 1.23
million individual values — agree to 2 × 10⁻¹³, and across 500 calibrations on
100 French catchments the largest deviation is 1.9 × 10⁻¹³.**

```bash
pip install grsuite
```

---

## Sixty seconds

```python
import grsuite as gr

catchment = gr.Catchment(dates, precip=P, pot_evap=E, obs_discharge=Q)

fit = catchment.calibrate("GR4J", period=("2000-01-01", "2009-12-31"))
print(fit.describe())

check = fit.evaluate(period=("2010-01-01", "2019-12-31"))
print(check.nse(), check.kge())
```

```
GR4J calibrated on KGE = 0.8694
  X1       170.7158   production store capacity [mm]
  X2         0.5897   groundwater exchange coefficient [mm/step]
  X3        78.2571   routing store capacity [mm]
  X4         2.2568   unit hydrograph time constant [step]
  17 iterations, 202 model runs
0.799 0.776
```

A split-sample test is one call:

```python
fit, validation = catchment.split_sample(
    "GR6J",
    calibration=("2000-01-01", "2009-12-31"),
    validation=("2010-01-01", "2019-12-31"),
    criterion="KGE",
)
```

Snow, on five elevation bands, needs only temperature and a hypsometric curve:

```python
alpine = gr.Catchment(dates, precip=P, pot_evap=E, obs_discharge=Q,
                      temperature=T, hypsometry=elevation_quantiles)
fit = alpine.calibrate("CemaNeigeGR4J")
```

Every internal flux is available, under airGR's own names:

```python
sim = fit.simulate()
sim.to_dataframe()[["Qobs", "Qsim", "Prod", "Rout", "Perc", "Exch"]]
sim.plot(log=True)
```

---

## Why this exists

airGR is the reference implementation of the GR models: careful, peer-reviewed,
maintained at INRAE. Its one constraint is that it lives in R. Python users have
had to choose between calling R over a bridge, or trusting a reimplementation
nobody checked.

GRsuite removes that choice. It was written by translating airGR's Fortran
kernels and R logic line by line, and it is continuously verified against airGR's
own outputs, which ship with the test suite. If GRsuite and airGR ever disagree,
CI fails.

| | |
|---|---|
| **Same numbers** | 502 quantities compared value by value against airGR — [every internal variable](https://github.com/Najim33/GRsuite/blob/main/docs/VALIDATION.md#the-internal-state-variable-by-variable), not just discharge |
| **Same calibration** | Michel's algorithm reproduced step for step — identical iteration counts on 500 calibrations |
| **No R needed** | pure Python + NumPy, JIT-compiled with Numba |
| **Faster** | 3.2× quicker than airGR on the same 500 calibrations |
| **Actually usable** | a one-line calibration API on top of the faithful airGR-mirroring one |

---

## What is covered

| Component | Models and functions | Published as |
|---|---|---|
| Daily | GR4J, GR5J, GR6J | Perrin et al. (2003); Le Moine (2008); Pushpalatha et al. (2011) |
| Hourly | GR4H, GR5H (with or without the interception store) | Mathevet (2005); Ficchì (2017); Ficchì et al. (2019) |
| Monthly / annual | GR2M, GR1A | Mouelhi et al. (2006a, 2006b) |
| Snow | CemaNeige, with or without linear hysteresis | Valéry (2010); Valéry et al. (2014); Riboust et al. (2019) |
| Coupled | CemaNeige + GR4J / GR5J / GR6J / GR4H / GR5H | as above |
| Evapotranspiration | Oudin formula | Oudin et al. (2005) |
| Elevation bands | Valéry extrapolation, daily temperature gradients | Valéry (2010) |
| Criteria | NSE, KGE, KGE′, RMSE, weighted composites | Nash & Sutcliffe (1970); Gupta et al. (2009); Kling et al. (2012) |
| Transformations | `sqrt`, `log`, `inv`, `sort`, `boxcox`, powers | Santos et al. (2018) |
| Calibration | Michel's algorithm (grid screening then steepest descent) | Michel (1991) |
| Utilities | time-series aggregation, semi-distributed routing, interception capacity | Lobligeois (2014); Ficchì et al. (2019) |

Full bibliographic entries, with DOIs, are in **[docs/REFERENCES.md](https://github.com/Najim33/GRsuite/blob/main/docs/REFERENCES.md)**.

`CreateErrorCrit_GAPX` (de Lavenne et al., 2019) and airGR's plotting helpers are
the only pieces left out. GR6H does not exist in airGR 1.7.9, so it does not
exist here either.

---

## Validation

![Measured agreement against the 5 % requirement](https://raw.githubusercontent.com/Najim33/GRsuite/main/docs/assets/deviation_scale.png)

The brief was to stay within 5 %. The measured agreement is eleven orders of
magnitude better than that, and sits at the precision limit of the reference files
themselves.

![500 calibrations](https://raw.githubusercontent.com/Najim33/GRsuite/main/docs/assets/deviation_histogram.png)

Across 500 independent calibrations — five model configurations on 100 CAMELS-FR
catchments — 34 parameter sets came out *bit-identical* to airGR's, and the rest
differ only in the last significant digit of a double.

**[Read the full validation report →](https://github.com/Najim33/GRsuite/blob/main/docs/VALIDATION.md)**

Two implementation details decide whether a re-implementation of airGR is merely
correct or actually identical. Both are documented in
**[docs/FIDELITY.md](https://github.com/Najim33/GRsuite/blob/main/docs/FIDELITY.md)**; the short version is that airGR's Fortran
evaluates `0.9` in single precision before promoting it, and that its automatic
warm-up period carries a leap-year correction. Miss either and your discharge is
off by 2.4 × 10⁻⁷ or 3 × 10⁻⁶ respectively.

---

## Two APIs, on purpose

**The high-level one** — for getting work done.

```python
catchment = gr.Catchment(dates, precip=P, pot_evap=E, obs_discharge=Q)
fit = catchment.calibrate("GR6J", criterion="KGE", transfo="sqrt")
```

**The faithful one** — a direct translation of airGR, function by function, for
porting an existing R workflow without rethinking it.

```python
inputs = gr.InputsModel(dates, precip=P, pot_evap=E)
options = gr.RunOptions(inputs, "GR4J", ind_period_run=index)
crit = gr.InputsCrit("NSE", obs=Q[index])
result = gr.calibration_michel(inputs, options, crit, gr.CalibOptions("GR4J"))
```

The mapping from every airGR function to its GRsuite counterpart is in
**[docs/MIGRATING_FROM_AIRGR.md](https://github.com/Najim33/GRsuite/blob/main/docs/MIGRATING_FROM_AIRGR.md)**.

---

## Examples

| File | What it shows |
|---|---|
| [`examples/01_simulate.py`](https://github.com/Najim33/GRsuite/blob/main/examples/01_simulate.py) | Run GR4J with known parameters, inspect the internal fluxes |
| [`examples/02_calibrate.py`](https://github.com/Najim33/GRsuite/blob/main/examples/02_calibrate.py) | Calibrate, then validate on an independent period |
| [`examples/03_compare_models.py`](https://github.com/Najim33/GRsuite/blob/main/examples/03_compare_models.py) | GR4J vs GR5J vs GR6J on the same catchment |
| [`examples/04_snow.py`](https://github.com/Najim33/GRsuite/blob/main/examples/04_snow.py) | CemaNeige on elevation bands, snow pack dynamics |
| [`examples/05_hourly.py`](https://github.com/Najim33/GRsuite/blob/main/examples/05_hourly.py) | GR4H and GR5H with an interception store |
| [`examples/06_low_flows.py`](https://github.com/Najim33/GRsuite/blob/main/examples/06_low_flows.py) | Objective functions for low-flow performance |
| [`examples/07_semi_distributed.py`](https://github.com/Najim33/GRsuite/blob/main/examples/07_semi_distributed.py) | Routing upstream sub-catchments to an outlet |
| [`examples/08_from_csv.py`](https://github.com/Najim33/GRsuite/blob/main/examples/08_from_csv.py) | Going from a plain CSV to a calibrated model |

All of them run against the demonstration catchment shipped in `tests/data`, so
they work straight after a clone.

---

## Installation

```bash
pip install grsuite            # runtime: numpy, numba
pip install "grsuite[io]"      # adds pandas for the DataFrame helpers
```

From source:

```bash
git clone https://github.com/Najim33/GRsuite
cd GRsuite
pip install -e ".[dev]"
pytest
```

Python 3.9 to 3.13, on Linux, macOS and Windows.

---

## Performance

The kernels are JIT-compiled by Numba, so the first call to a model pays a
compilation cost of a second or two and everything afterwards runs at compiled
speed. Recompilation is cached on disk between sessions.

On the 500-calibration benchmark (5 configurations × 100 catchments, 20 years of
daily data each):

| | Wall clock |
|---|---|
| airGR 1.7.9 (R + Fortran) | 17.5 min |
| GRsuite (Python + Numba) | 5.4 min |

Same parameters, same criteria, same iteration counts.

---

## Citing

GRsuite introduces no new hydrology. The models, the calibration algorithm and
the reference implementation belong to their authors, at Cemagref / Irstea /
INRAE. If GRsuite supports published work, cite **airGR**:

> Coron, L., Thirel, G., Delaigue, O., Perrin, C., Andréassian, V. (2017).
> The suite of lumped GR hydrological models in an R package.
> *Environmental Modelling & Software*, 94, 166–171.
> doi:[10.1016/j.envsoft.2017.05.002](https://doi.org/10.1016/j.envsoft.2017.05.002)

> Coron, L., Delaigue, O., Thirel, G., Dorchies, D., Perrin, C., Michel, C.
> airGR: Suite of GR Hydrological Models for Precipitation-Runoff Modelling.
> R package version 1.7.9, INRAE, HYCAR Research Unit, Antony, France.
> doi:[10.15454/EX11NA](https://doi.org/10.15454/EX11NA)

**and the model you actually ran** — GR4J is Perrin et al. (2003), CemaNeige is
Valéry (2010), and so on. Every one of them, with its DOI, is listed in
**[docs/REFERENCES.md](https://github.com/Najim33/GRsuite/blob/main/docs/REFERENCES.md)**.

If the port itself is relevant to your methods, cite this repository too — see
[CITATION.cff](https://github.com/Najim33/GRsuite/blob/main/CITATION.cff).

---

## Contributing

Bug reports, models, criteria and examples are all welcome. The one hard rule:
**anything that changes model numerics must keep the airGR comparison tests
green.** See [CONTRIBUTING.md](https://github.com/Najim33/GRsuite/blob/main/CONTRIBUTING.md).

---

## Licence

GPL-2.0-or-later. GRsuite is a derivative work of airGR (INRAE, GPL-2) and
carries the same licence — see [NOTICE](https://github.com/Najim33/GRsuite/blob/main/NOTICE) for the full attribution.

GRsuite is an independent, unofficial project. It is not produced, endorsed or
supported by INRAE, and it is not a release of airGR.

*Ce README existe aussi [en français](https://github.com/Najim33/GRsuite/blob/main/README.fr.md).*
