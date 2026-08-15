"""Error criteria, parameter transformations, aggregation, routing, calibration."""

import numpy as np
import pandas as pd
import pytest

import grsuite as gr
from conftest import assert_matches, load
from grsuite.transfo import TRANSFO_FUNCS

CRITERIA = load("error_crits")
TRANSFOS = load("transfo_param")
CALIBRATIONS = load("calib_reference")

CALIB_MODELS = {"GR4J": gr.run_model_gr4j, "GR5J": gr.run_model_gr5j,
                "GR6J": gr.run_model_gr6j}


@pytest.fixture(scope="module")
def gr4j_run(daily):
    options = gr.RunOptions(daily["inputs"], "GR4J",
                            ind_period_run=daily["ind_run"])
    return gr.run_model_gr4j(daily["inputs"], options,
                             [257.238, 1.012, 88.235, 2.208])


# ---------------------------------------------------------------------------
# Error criteria
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("row", range(len(CRITERIA)))
def test_error_criteria(daily, gr4j_run, row):
    r = CRITERIA.iloc[row]
    obs = daily["Q"][daily["ind_run"]]
    transfo = "" if pd.isna(r["transfo"]) else str(r["transfo"])
    epsilon = None if pd.isna(r["epsilon"]) else float(r["epsilon"])
    expected = float(r["value"])

    if r["crit"] == "Composite_NSE":
        criterion = gr.InputsCritCompo([
            gr.InputsCrit("NSE", obs=obs, transfo="", weights=0.6),
            gr.InputsCrit("NSE", obs=obs, transfo="log", epsilon=0.01,
                          weights=0.4)])
        # airGR reports the composite already multiplied by its sign convention
        value = -gr.error_crit(criterion, gr4j_run).crit_value
    else:
        criterion = gr.InputsCrit(r["crit"], obs=obs, transfo=transfo,
                                  epsilon=epsilon)
        value = gr.error_crit(criterion, gr4j_run).crit_value

    assert_matches(np.array([expected]), np.array([value]),
                   "%s[%s]" % (r["crit"], transfo or "-"))


def test_perfect_simulation_scores_perfectly(gr4j_run):
    perfect = gr.InputsCrit("NSE", obs=gr4j_run["Qsim"])
    assert gr.error_crit(perfect, gr4j_run).crit_value == pytest.approx(1.0)
    perfect_kge = gr.InputsCrit("KGE", obs=gr4j_run["Qsim"])
    assert gr.error_crit(perfect_kge, gr4j_run).crit_value == pytest.approx(1.0)


def test_bool_crit_restricts_the_period(daily, gr4j_run):
    obs = daily["Q"][daily["ind_run"]]
    mask = np.zeros(obs.shape[0], dtype=bool)
    mask[: obs.shape[0] // 2] = True
    full = gr.error_crit(gr.InputsCrit("NSE", obs=obs), gr4j_run).crit_value
    half = gr.error_crit(gr.InputsCrit("NSE", obs=obs, bool_crit=mask),
                         gr4j_run).crit_value
    assert full != half
    assert np.isfinite(half)


def test_unknown_criterion_is_rejected(daily):
    with pytest.raises(ValueError, match="critere inconnu"):
        gr.InputsCrit("R2", obs=daily["Q"])


# ---------------------------------------------------------------------------
# Parameter transformations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("row", range(len(TRANSFOS)))
def test_parameter_transformation(row):
    r = TRANSFOS.iloc[row]
    model = r["model"]
    if model not in TRANSFO_FUNCS:
        pytest.skip("model %s has no transformation" % model)
    transformed = np.array([float(v) for v in r["paramT"].split("|")])
    expected_real = np.array([float(v) for v in r["paramR"].split("|")])
    expected_back = np.array([float(v) for v in r["paramT2"].split("|")])

    real = np.atleast_1d(gr.transfo_param(transformed, "TR", model))
    back = np.atleast_1d(gr.transfo_param(real, "RT", model))
    assert_matches(expected_real, real, "%s TR" % model)
    assert_matches(expected_back, back, "%s RT" % model)


@pytest.mark.parametrize("model", sorted(TRANSFO_FUNCS))
def test_transformation_round_trip(model):
    from grsuite.transfo import N_PARAM
    rng = np.random.default_rng(11)
    transformed = rng.uniform(-9.0, 9.0, N_PARAM[model])
    real = gr.transfo_param(transformed, "TR", model)
    back = np.atleast_1d(gr.transfo_param(real, "RT", model))
    assert_matches(transformed, back, "%s round trip" % model, tol=1e-10)


def test_wrong_parameter_count_is_rejected():
    with pytest.raises(ValueError, match="requiert 4 parametres"):
        gr.transfo_param([1.0, 2.0], "TR", "GR4J")


# ---------------------------------------------------------------------------
# Time-series aggregation
# ---------------------------------------------------------------------------

AGGREG_CASES = [
    ("basin_monthly", "%Y%m", {"P": "sum", "E": "sum", "Qmm": "sum"}, 1),
    ("aggreg_monthly_mean", "%Y%m", {"P": "mean", "E": "mean", "Qmm": "mean"}, 1),
    ("basin_yearly", "%Y", {"P": "sum", "E": "sum", "Qmm": "sum"}, 1),
    ("aggreg_yearly_sept", "%Y", {"P": "sum", "E": "sum", "Qmm": "sum"}, 9),
]


@pytest.mark.parametrize("reference_file,fmt,funs,first_month", AGGREG_CASES)
def test_series_aggreg(daily, reference_file, fmt, funs, first_month):
    data = {"P": daily["P"], "E": daily["E"], "Qmm": daily["Q"]}
    dates, aggregated = gr.series_aggreg(daily["dates"], data, fmt=fmt,
                                         convert_fun=funs,
                                         year_first_month=first_month)
    reference = load(reference_file)
    expected_dates = reference["Date"].to_numpy().astype("datetime64[D]")
    assert np.array_equal(dates, expected_dates), (
        "%s: %i aggregated steps, airGR produced %i"
        % (reference_file, dates.size, expected_dates.size))
    for name in ("P", "E", "Qmm"):
        assert_matches(reference[name].to_numpy(), aggregated[name],
                       "%s / %s" % (reference_file, name))


def test_unsupported_aggregation_format(daily):
    with pytest.raises(ValueError, match="format non pris en charge"):
        gr.series_aggreg(daily["dates"], {"P": daily["P"]}, fmt="%W")


# ---------------------------------------------------------------------------
# Semi-distributed routing
# ---------------------------------------------------------------------------


def test_semi_distributed_lag(daily):
    config = load("sd_config").iloc[0]
    upstream = load("sd_qupstream")
    reference = load("sim_SD_lag")

    options = gr.RunOptions(daily["inputs"], "GR4J",
                            ind_period_run=daily["ind_run"])
    downstream = gr.run_model_gr4j(daily["inputs"], options,
                                   [257.238, 1.012, 88.235, 2.208])
    network = gr.InputsModelSD(
        q_upstream=np.column_stack([upstream["Qup1"].to_numpy(),
                                    upstream["Qup2"].to_numpy()]),
        length_hydro=[float(config["length1"]), float(config["length2"])],
        basin_areas=[float(config["area1"]), float(config["area2"]),
                     float(config["area3"])])
    out = gr.run_model_lag(network, options, [float(config["speed"])],
                           downstream["Qsim"],
                           warmup_q_down=downstream["WarmUpQsim"])
    for name in ("Qsim", "Qsim_m3", "QsimDown"):
        assert_matches(reference[name].to_numpy(), out[name], "lag / " + name)


def test_lag_rejects_inconsistent_network():
    with pytest.raises(ValueError, match="basin_areas"):
        gr.InputsModelSD(q_upstream=np.zeros((10, 2)), length_hydro=[10.0, 20.0],
                         basin_areas=[100.0, 200.0])


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("row", range(len(CALIBRATIONS)))
def test_calibration_matches_airgr(daily, row):
    r = CALIBRATIONS.iloc[row]
    model = r["model"]
    transfo = "" if pd.isna(r["transfo"]) else str(r["transfo"])
    epsilon = 0.01 if transfo == "log" else None

    options = gr.RunOptions(daily["inputs"], model,
                            ind_period_run=daily["ind_run"])
    criterion = gr.InputsCrit(r["crit"], obs=daily["Q"][daily["ind_run"]],
                              transfo=transfo, epsilon=epsilon)
    result = gr.calibration_michel(daily["inputs"], options, criterion,
                                   gr.CalibOptions(model), verbose=False)

    expected_param = np.array([float(v) for v in r["param"].split("|")])
    assert_matches(expected_param, result.param_final_r,
                   "%s %s[%s] parameters" % (model, r["crit"], transfo or "-"),
                   tol=1e-8)
    assert_matches(np.array([float(r["crit_final"])]),
                   np.array([result.crit_final]),
                   "%s %s[%s] criterion" % (model, r["crit"], transfo or "-"),
                   tol=1e-8)
    assert result.n_iter == int(r["n_iter"]), (
        "calibration took %i iterations, airGR took %i"
        % (result.n_iter, int(r["n_iter"])))


def test_calibration_improves_on_the_starting_point(daily):
    options = gr.RunOptions(daily["inputs"], "GR4J",
                            ind_period_run=daily["ind_run"])
    obs = daily["Q"][daily["ind_run"]]
    criterion = gr.InputsCrit("NSE", obs=obs)
    result = gr.calibration_michel(daily["inputs"], options, criterion,
                                   gr.CalibOptions("GR4J"), verbose=False)
    assert result.crit_final > 0.5
    assert result.n_runs > result.n_iter
    assert result.hist_crit.shape[0] == result.n_iter


def test_fixed_parameters_are_not_optimised(daily):
    options = gr.RunOptions(daily["inputs"], "GR4J",
                            ind_period_run=daily["ind_run"])
    criterion = gr.InputsCrit("NSE", obs=daily["Q"][daily["ind_run"]])
    calib = gr.CalibOptions("GR4J", fixed_param=[np.nan, 0.0, np.nan, np.nan])
    result = gr.calibration_michel(daily["inputs"], options, criterion, calib,
                                   verbose=False)
    assert result.param_final_r[1] == 0.0
