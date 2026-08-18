# Migrating from airGR to GRsuite

GRsuite is a function-by-function translation of airGR 1.7.9 — an unofficial one,
written independently of INRAE. This page maps every airGR function to its GRsuite
counterpart and lists the (few) differences you need to know when porting an R
script. Porting away from airGR means leaving the reference implementation behind;
that is a decision worth making deliberately.

The short version: take your R workflow, lowercase the function names, convert
your indices from 1-based to 0-based, and you get the same numbers.

## Core workflow

| airGR (R) | GRsuite (Python) |
|---|---|
| `CreateInputsModel(FUN_MOD, DatesR, Precip, PotEvap, ...)` | `InputsModel(dates, precip, pot_evap, ...)` |
| `CreateRunOptions(FUN_MOD, InputsModel, IndPeriod_Run, ...)` | `RunOptions(inputs_model, model, ind_period_run, ...)` |
| `CreateIniStates(FUN_MOD, InputsModel, ProdStore, RoutStore, ...)` | `ini_states` / `ini_res_levels` arguments of `RunOptions` |
| `CreateInputsCrit(FUN_CRIT, Obs, ...)` | `InputsCrit(fun_crit, obs, transfo, ...)` |
| `CreateInputsCrit` (composite, `Weights`) | `InputsCritCompo([InputsCrit(..., weights=w), ...])` |
| `CreateCalibOptions(FUN_MOD, ...)` | `CalibOptions(model, fixed_param, ...)` |
| `Calibration_Michel(InputsModel, RunOptions, InputsCrit, CalibOptions, FUN_MOD)` | `calibration_michel(inputs_model, run_options, inputs_crit, calib_options)` |
| `RunModel(InputsModel, RunOptions, Param, FUN_MOD)` | `run_model(inputs_model, run_options, param, model=...)` |
| `RunModel_GR4J` … `RunModel_CemaNeigeGR5H` | `run_model_gr4j` … `run_model_cemaneige_gr5h` |
| `RunModel_Lag(InputsModelSD, RunOptions, Param, QcontribDown)` | `run_model_lag(inputs_sd, run_options, param, q_contrib_down)` |
| `InputsModelSD` helper in semi-distributed setups | `InputsModelSD(q_upstream, length_hydro, basin_areas)` |

## Criteria and utilities

| airGR (R) | GRsuite (Python) |
|---|---|
| `ErrorCrit_NSE(InputsCrit, OutputsModel)` | `error_crit_nse(inputs_crit, outputs)` |
| `ErrorCrit_KGE` / `ErrorCrit_KGE2` / `ErrorCrit_RMSE` | `error_crit_kge` / `error_crit_kge2` / `error_crit_rmse` |
| `ErrorCrit(InputsCrit, OutputsModel)` | `error_crit(inputs_crit, outputs)` |
| `TransfoParam(ParamIn, Direction, FUN_TRANSFO)` | `transfo_param(param, direction, model)` |
| `SeriesAggreg(x, Format, ConvertFun, ...)` | `series_aggreg(dates, data, fmt, convert_fun, ...)` |
| `PE_Oudin(JD, Temp, Lat, LatUnit)` | `pe_oudin(julian_day, temp, lat, lat_unit)` |
| `DataAltiExtrapolation_Valery(DatesR, Precip, Temp, ...)` | `data_alti_extrapolation_valery(dates, precip, temp, ...)` |
| `Imax(InputsModel, IndPeriod_Run, TestedValues)` | `imax_estimate(inputs_model, ind_period_run, tested_values)` |
| `CreateErrorCrit_GAPX` | not ported |
| `plot_OutputsModel`, `plot_OutputsCalib` | not ported (use `Simulation.plot()`) |

## Data assimilation (airGRdatassim)

| airGRdatassim (R) | GRsuite (Python) |
|---|---|
| `CreateInputsPert(FUN_MOD, DatesR, Precip, PotEvap, ..., NbMbr, Seed)` | `InputsPert(model, dates, precip=..., pot_evap=..., nb_mbr=..., seed=...)` |
| `RunModel_DA(InputsModel, InputsPert, Qobs, IndRun, FUN_MOD, Param, DaMethod, NbMbr, StateEnKF, StatePert, Seed)` | `run_model_da(inputs_model, ind_run, model, param, inputs_pert=..., qobs=..., da_method=..., nb_mbr=..., state_enkf=..., state_pert=..., seed=...)` |
| `plot(OutputsModelDA)` | `Assimilation.plot()` (high level: `Catchment.assimilate()`) |

One deliberate difference to know about: **airGRdatassim 0.1.4's EnKF never
updates any state.** `DA_EnKF` selects the variables to update with
`which(StateEnKF == 1)`, but `RunModel_DA` passes `StateEnKF` as a character
vector (the documented `c("Prod", "Rout", ...)`), so the selection is always
empty and the Kalman update is a silent no-op. GRsuite implements the
documented behaviour: the names in `state_enkf` select the variables to
update. An EnKF ported from R therefore gives *different* results here —
different because it actually assimilates. The particle filter and
`CreateInputsPert` are unaffected and agree with the R package value by
value. See [VALIDATION.md section 9](VALIDATION.md#92-the-one-deliberate-deviation-the-enkf-fix).

## Differences to know about

- **Indexing.** R is 1-based, GRsuite is 0-based. `IndPeriod_Run = 365:730` in R
  is `ind_period_run=np.arange(364, 730)` in GRsuite.
- **Model selection.** airGR passes `FUN_MOD = RunModel_GR4J`; GRsuite passes the
  model name as a string (`"GR4J"`) to `RunOptions`, `CalibOptions` and
  `run_model`.
- **Outputs.** `RunModel_*` returns a dictionary keyed by airGR's own output
  names (`"Qsim"`, `"Prod"`, `"Perc"`, `"Exch"`, `"StateEnd"`, …), with NumPy
  arrays instead of R vectors. `OutputsModel$DatesR` is `outputs["DatesR"]`.
- **Criteria.** Same convention as airGR: NSE, KGE and KGE′ are maximised,
  RMSE is minimised, and calibration always minimises `crit_value * multiplier`.
- **Missing values.** `NA` becomes `np.nan`; series are plain NumPy arrays.

## Porting example

```r
## airGR
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = dates,
                                 Precip = P, PotEvap = E)
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel,
                               IndPeriod_Run = seq_len(length(dates)))
InputsCrit  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                                RunOptions = RunOptions, Obs = Q)
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)
OutputsCalib <- Calibration_Michel(InputsModel, RunOptions, InputsCrit,
                                   CalibOptions, FUN_MOD = RunModel_GR4J)
```

```python
## GRsuite
inputs  = gr.InputsModel(dates, precip=P, pot_evap=E)
options = gr.RunOptions(inputs, "GR4J", ind_period_run=np.arange(len(dates)))
crit    = gr.InputsCrit("NSE", obs=Q)
result  = gr.calibration_michel(inputs, options, crit, gr.CalibOptions("GR4J"))
```

Or skip all of it and use the high-level API:

```python
catchment = gr.Catchment(dates, precip=P, pot_evap=E, obs_discharge=Q)
fit = catchment.calibrate("GR4J", criterion="NSE")
```

## References

Function names change; the science does not. The reference for each model,
criterion and algorithm — airGR's own and the underlying publications — is in
**[REFERENCES.md](REFERENCES.md)**. A paper that cited airGR's `RunModel_GR4J`
should still cite Perrin et al. (2003) and airGR (Coron et al., 2017) after
porting to `run_model_gr4j`.
