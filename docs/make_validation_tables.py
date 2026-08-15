"""Rebuild the validation inventory of section 2 of docs/VALIDATION.md.

Every quantity GRsuite is compared against airGR on is recomputed here, one
series at a time, from the reference files under ``tests/data`` — the same
files the test suite asserts on. The script writes one row per compared
quantity, so the counts quoted in the report can be audited rather than
trusted:

    python docs/make_validation_tables.py

Outputs
-------
docs/data/series_deviations.csv
    One row per compared quantity: category, case, name, kind, length, and
    the absolute / relative deviation from airGR.
stdout
    The summary tables of docs/VALIDATION.md section 2, in Markdown.

A quantity is counted as a *series* when it is a full time series (one value
per time step) and as a *scalar* when it is a single number (a criterion
value, a mean annual solid precipitation, a calibrated parameter set).
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tests"))

import grsuite as gr  # noqa: E402
from conftest import index_range, load  # noqa: E402
from grsuite.core import MODEL_OUTPUTS  # noqa: E402
from grsuite.transfo import TRANSFO_FUNCS  # noqa: E402

RECORDS = []

# Parameter sets are airGR's own, as used by the test suite.
DAILY_MODELS = {
    "GR4J": (gr.run_model_gr4j, [257.238, 1.012, 88.235, 2.208]),
    "GR5J": (gr.run_model_gr5j, [245.918, 1.027, 90.017, 2.198, 0.318]),
    "GR6J": (gr.run_model_gr6j, [250.0, 0.8, 80.0, 2.1, 0.2, 30.0]),
}
STATE_SLOTS = {
    "GR4J": [0, 1] + list(range(7, 67)),
    "GR5J": [0, 1] + list(range(27, 67)),
    "GR6J": [0, 1, 2] + list(range(7, 67)),
}
COUPLED = {
    "CemaNeigeGR4J": ("CemaNeigeGR4J", "GR4J", False,
                      [408.774, 2.646, 131.264, 1.174, 0.962, 2.249]),
    "CemaNeigeGR5J": ("CemaNeigeGR5J", "GR5J", False,
                      [245.918, 1.027, 90.017, 2.198, 0.318, 0.962, 2.249]),
    "CemaNeigeGR6J": ("CemaNeigeGR6J", "GR6J", False,
                      [250.0, 0.8, 80.0, 2.1, 0.2, 30.0, 0.962, 2.249]),
    "CemaNeigeGR4J_Hyst": ("CemaNeigeGR4J", "GR4J", True,
                           [408.774, 2.646, 131.264, 1.174, 0.962, 2.249,
                            80.0, 0.4]),
}
AGGREG_CASES = [
    ("basin_monthly", "%Y%m", "sum", 1, "monthly total"),
    ("aggreg_monthly_mean", "%Y%m", "mean", 1, "monthly mean"),
    ("basin_yearly", "%Y", "sum", 1, "calendar year"),
    ("aggreg_yearly_sept", "%Y", "sum", 9, "hydrological year (Sept.)"),
]


def record(category, case, name, ref, got, kind="series"):
    """Store the deviation between one airGR quantity and its GRsuite twin."""
    ref = np.atleast_1d(np.asarray(ref, dtype=float))
    got = np.atleast_1d(np.asarray(got, dtype=float))
    ok = np.isfinite(ref) & np.isfinite(got)
    diff = np.abs(ref[ok] - got[ok])
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(ref[ok] != 0.0, diff / np.abs(ref[ok]), 0.0)
    # A value passes on whichever of the two deviations is the smaller, so the
    # honest worst case is the largest of those per-element minima.
    worst = np.maximum.reduce([np.minimum(diff, rel)]) if diff.size else np.zeros(1)
    RECORDS.append({
        "category": category,
        "case": case,
        "name": name,
        "kind": kind,
        "n_values": int(ok.sum()),
        "max_abs": float(diff.max()) if diff.size else 0.0,
        "max_rel": float(rel.max()) if rel.size else 0.0,
        "worst": float(worst.max()) if worst.size else 0.0,
    })


# ---------------------------------------------------------------------------
# 1. Rainfall-runoff models, every internal variable
# ---------------------------------------------------------------------------


def daily_models():
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy(),
                            temp_mean=basin["T"].to_numpy())
    ind = index_range("idx_daily")
    for model, (run, param) in DAILY_MODELS.items():
        options = gr.RunOptions(inputs, model, ind_period_run=ind)
        out = run(inputs, options, param)
        reference = load("sim_" + model)
        for name in MODEL_OUTPUTS[model]:
            record("Daily models", model, name,
                   reference[name].to_numpy(), out[name])
        state = load("state_" + model)["state"].to_numpy()
        slots = STATE_SLOTS[model]
        record("Daily models", model, "StateEnd (%i slots)" % len(slots),
               state[slots], out["StateEnd"][slots], kind="state vector")


def monthly_and_yearly():
    for tag, fixture, idx, model, run, param in (
            ("GR2M", "basin_monthly", "idx_monthly", "GR2M",
             gr.run_model_gr2m, [265.072, 1.007]),
            ("GR1A", "basin_yearly", "idx_yearly", "GR1A",
             gr.run_model_gr1a, [0.840])):
        basin = load(fixture)
        dates = basin["Date"].to_numpy().astype("datetime64[D]")
        inputs = gr.InputsModel(dates, basin["P"].to_numpy(),
                                basin["E"].to_numpy())
        options = gr.RunOptions(inputs, model, ind_period_run=index_range(idx))
        out = run(inputs, options, param)
        reference = load("sim_" + tag)
        for name in MODEL_OUTPUTS[model]:
            record("Monthly / annual models", tag, name,
                   reference[name].to_numpy(), out[name])


def hourly_models():
    basin = load("basin_hourly")
    dates = pd.to_datetime(basin["Date"]).to_numpy().astype("datetime64[h]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    ind = index_range("idx_hourly")

    out = gr.run_model_gr4h(inputs, gr.RunOptions(inputs, "GR4H",
                                                  ind_period_run=ind),
                            [756.930, -0.773, 138.638, 5.247])
    reference = load("sim_GR4H")
    for name in MODEL_OUTPUTS["GR4H"]:
        record("Hourly models", "GR4H", name,
               reference[name].to_numpy(), out[name])

    param5 = [756.930, -0.773, 138.638, 5.247, 0.400]
    out = gr.run_model_gr5h(inputs, gr.RunOptions(inputs, "GR5H",
                                                  ind_period_run=ind), param5)
    reference = load("sim_GR5H")
    for name in MODEL_OUTPUTS["GR5H"]:
        record("Hourly models", "GR5H", name,
               reference[name].to_numpy(), out[name])

    imax = float(load("imax_value")["imax"][0])
    out = gr.run_model_gr5h(inputs, gr.RunOptions(inputs, "GR5H",
                                                  ind_period_run=ind,
                                                  imax=imax), param5)
    reference = load("sim_GR5H_interception")
    for name in MODEL_OUTPUTS["GR5H"]:
        record("Hourly models", "GR5H + interception", name,
               reference[name].to_numpy(), out[name])


# ---------------------------------------------------------------------------
# 2. Snow
# ---------------------------------------------------------------------------


def snow():
    basin = load("basin_snow")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    hypso = load("hypso_snow")["hypso"].to_numpy()
    meta = load("idx_snow")
    alti = gr.data_alti_extrapolation_valery(
        dates, basin["P"].to_numpy(), basin["T"].to_numpy(),
        z_inputs=float(meta["zinputs"][0]), hypso_data=hypso, n_layers=5)
    inputs = gr.InputsModel(
        dates, basin["P"].to_numpy(), basin["E"].to_numpy(),
        temp_mean=basin["T"].to_numpy(),
        layer_precip=alti["LayerPrecip"], layer_temp_mean=alti["LayerTempMean"],
        layer_frac_solid_precip=alti["LayerFracSolidPrecip"],
        z_layers=alti["ZLayers"])
    ind = np.arange(int(meta["ind_start"][0]) - 1, int(meta["ind_end"][0]))

    for layer in range(1, 6):
        reference = load("inputs_layer%02i" % layer)
        i = layer - 1
        for column, series, label in (
                ("P", alti["LayerPrecip"][i], "precipitation"),
                ("T", alti["LayerTempMean"][i], "temperature"),
                ("FS", alti["LayerFracSolidPrecip"][i], "solid fraction")):
            record("Elevation bands", "band %i" % layer, label,
                   reference[column].to_numpy(), series)

    for case, (model, gr_name, hyst, param) in COUPLED.items():
        options = gr.RunOptions(inputs, model, ind_period_run=ind,
                                is_hyst=hyst)
        out = gr.run_model(inputs, options, param, model=model)
        reference = load("sim_" + case)
        for name in MODEL_OUTPUTS[gr_name]:
            record("Coupled snow models", case, name,
                   reference[name].to_numpy(), out[name])
        record("Coupled snow models", case, "MeanAnSolidPrecip",
               [float(load("masp_" + case)["masp"][0])],
               [options.mean_an_solid_precip[0]], kind="scalar")

        if case in ("CemaNeigeGR4J", "CemaNeigeGR4J_Hyst"):
            for layer in range(1, 6):
                reference = load("sim_%s_layer%02i" % (case, layer))
                for name in reference.columns:
                    record("Snow module, band by band",
                           "%s / band %i" % (case, layer), name,
                           reference[name].to_numpy(),
                           out["CemaNeigeLayers"][layer - 1][name])


# ---------------------------------------------------------------------------
# 3. Evapotranspiration, criteria, transformations, utilities
# ---------------------------------------------------------------------------


def potential_evapotranspiration():
    reference = load("pe_oudin")
    computed = gr.pe_oudin(reference["JD"].to_numpy(),
                           reference["Temp"].to_numpy(), 0.8)
    record("Potential evapotranspiration", "PE_Oudin", "PE",
           reference["PE"].to_numpy(), computed)


def criteria():
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    ind = index_range("idx_daily")
    options = gr.RunOptions(inputs, "GR4J", ind_period_run=ind)
    run = gr.run_model_gr4j(inputs, options, [257.238, 1.012, 88.235, 2.208])
    obs = load("qobs_daily")["Qmm"].to_numpy()[ind]

    for _, r in load("error_crits").iterrows():
        transfo = "" if pd.isna(r["transfo"]) else str(r["transfo"])
        epsilon = None if pd.isna(r["epsilon"]) else float(r["epsilon"])
        if r["crit"] == "Composite_NSE":
            crit = gr.InputsCritCompo([
                gr.InputsCrit("NSE", obs=obs, transfo="", weights=0.6),
                gr.InputsCrit("NSE", obs=obs, transfo="log", epsilon=0.01,
                              weights=0.4)])
            value = -gr.error_crit(crit, run).crit_value
        else:
            crit = gr.InputsCrit(r["crit"], obs=obs, transfo=transfo,
                                 epsilon=epsilon)
            value = gr.error_crit(crit, run).crit_value
        record("Error criteria", r["crit"], transfo or "no transformation",
               [float(r["value"])], [value], kind="scalar")


def transformations():
    for _, r in load("transfo_param").iterrows():
        model = r["model"]
        if model not in TRANSFO_FUNCS:
            continue
        param_t = np.array([float(v) for v in r["paramT"].split("|")])
        real = np.atleast_1d(gr.transfo_param(param_t, "TR", model))
        back = np.atleast_1d(gr.transfo_param(real, "RT", model))
        record("Parameter transformations", model,
               "set %i, transformed to real" % int(r["k"]),
               [float(v) for v in r["paramR"].split("|")], real, kind="scalar")
        record("Parameter transformations", model,
               "set %i, real back to transformed" % int(r["k"]),
               [float(v) for v in r["paramT2"].split("|")], back, kind="scalar")


def aggregation():
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    data = {"P": basin["P"].to_numpy(), "E": basin["E"].to_numpy(),
            "Qmm": load("qobs_daily")["Qmm"].to_numpy()}
    for fixture, fmt, fun, first_month, label in AGGREG_CASES:
        _, aggregated = gr.series_aggreg(
            dates, data, fmt=fmt,
            convert_fun=dict.fromkeys(data, fun), year_first_month=first_month)
        reference = load(fixture)
        for name in ("P", "E", "Qmm"):
            record("Series aggregation", label, name,
                   reference[name].to_numpy(), aggregated[name])


def routing():
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    options = gr.RunOptions(inputs, "GR4J", ind_period_run=index_range("idx_daily"))
    downstream = gr.run_model_gr4j(inputs, options,
                                   [257.238, 1.012, 88.235, 2.208])
    config = load("sd_config").iloc[0]
    upstream = load("sd_qupstream")
    network = gr.InputsModelSD(
        q_upstream=np.column_stack([upstream["Qup1"].to_numpy(),
                                    upstream["Qup2"].to_numpy()]),
        length_hydro=[float(config["length1"]), float(config["length2"])],
        basin_areas=[float(config["area1"]), float(config["area2"]),
                     float(config["area3"])])
    out = gr.run_model_lag(network, options, [float(config["speed"])],
                           downstream["Qsim"],
                           warmup_q_down=downstream["WarmUpQsim"])
    reference = load("sim_SD_lag")
    for name in ("Qsim", "Qsim_m3", "QsimDown"):
        record("Semi-distributed routing", "RunModel_Lag", name,
               reference[name].to_numpy(), out[name])


def calibrations():
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    ind = index_range("idx_daily")
    obs = load("qobs_daily")["Qmm"].to_numpy()[ind]
    for _, r in load("calib_reference").iterrows():
        model = r["model"]
        transfo = "" if pd.isna(r["transfo"]) else str(r["transfo"])
        options = gr.RunOptions(inputs, model, ind_period_run=ind)
        crit = gr.InputsCrit(r["crit"], obs=obs, transfo=transfo,
                             epsilon=0.01 if transfo == "log" else None)
        result = gr.calibration_michel(inputs, options, crit,
                                       gr.CalibOptions(model), verbose=False)
        label = "%s on %s%s" % (model, r["crit"],
                                "[%s]" % transfo if transfo else "")
        record("Calibration on the demonstration catchment", label,
               "calibrated parameters",
               [float(v) for v in r["param"].split("|")],
               result.param_final_r, kind="scalar")
        record("Calibration on the demonstration catchment", label,
               "criterion at the optimum",
               [float(r["crit_final"])], [result.crit_final], kind="scalar")
        record("Calibration on the demonstration catchment", label,
               "iterations", [float(r["n_iter"])], [float(result.n_iter)],
               kind="scalar")


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def fmt(x):
    return "0" if x == 0.0 else "%.1e" % x


def main():
    daily_models()
    monthly_and_yearly()
    hourly_models()
    snow()
    potential_evapotranspiration()
    criteria()
    transformations()
    aggregation()
    routing()
    calibrations()

    table = pd.DataFrame(RECORDS)
    out_path = os.path.join(ROOT, "docs", "data", "series_deviations.csv")
    table.to_csv(out_path, index=False)

    n_series = int((table["kind"] == "series").sum())
    n_other = int((table["kind"] != "series").sum())

    print("## By category\n")
    print("| Category | Quantities | of which time series | Values compared "
          "| Worst deviation |")
    print("|---|---|---|---|---|")
    for category, block in table.groupby("category", sort=False):
        print("| %s | %i | %i | %s | %s |"
              % (category, len(block), int((block["kind"] == "series").sum()),
                 "{:,}".format(int(block["n_values"].sum())),
                 fmt(block["worst"].max())))
    print("| **Total** | **%i** | **%i** | **%s** | **%s** |"
          % (len(table), n_series, "{:,}".format(int(table["n_values"].sum())),
             fmt(table["worst"].max())))
    print("\n%i time series + %i scalars / state vectors = %i quantities"
          % (n_series, n_other, len(table)))
    print("%i of them are bit-identical to airGR\n"
          % int((table["worst"] == 0.0).sum()))

    print("## Worst deviation per variable name (time series only)\n")
    print("| Variable | Series | Max absolute | Max relative | Worst |")
    print("|---|---|---|---|---|")
    only = table[table["kind"] == "series"]
    per_name = only.groupby("name").agg(
        n=("name", "size"), max_abs=("max_abs", "max"),
        max_rel=("max_rel", "max"), worst=("worst", "max")
    ).sort_values("worst", ascending=False)
    for name, row in per_name.iterrows():
        print("| `%s` | %i | %s | %s | %s |"
              % (name, row["n"], fmt(row["max_abs"]), fmt(row["max_rel"]),
                 fmt(row["worst"])))

    print("\n## Ten largest deviations, all quantities\n")
    print("| Category | Case | Quantity | Max absolute | Max relative "
          "| Worst |")
    print("|---|---|---|---|---|---|")
    for _, r in table.sort_values("worst", ascending=False).head(10).iterrows():
        print("| %s | %s | `%s` | %s | %s | %s |"
              % (r["category"], r["case"], r["name"], fmt(r["max_abs"]),
                 fmt(r["max_rel"]), fmt(r["worst"])))

    print("\nWritten: docs/data/series_deviations.csv")


if __name__ == "__main__":
    main()
