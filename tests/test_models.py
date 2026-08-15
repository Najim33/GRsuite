"""Every rainfall-runoff model reproduces airGR, variable by variable."""

import numpy as np
import pytest

import grsuite as gr
from conftest import assert_matches, index_range, load
from grsuite.core import MODEL_OUTPUTS

DAILY_MODELS = {
    "GR4J": (gr.run_model_gr4j, [257.238, 1.012, 88.235, 2.208]),
    "GR5J": (gr.run_model_gr5j, [245.918, 1.027, 90.017, 2.198, 0.318]),
    "GR6J": (gr.run_model_gr6j, [250.0, 0.8, 80.0, 2.1, 0.2, 30.0]),
}

# airGR blanks out the state slots a model does not use; only the slots the
# model actually fills are meaningful to compare.
STATE_SLOTS = {
    "GR4J": [0, 1] + list(range(7, 67)),
    "GR5J": [0, 1] + list(range(27, 67)),       # GR5J has no UH1
    "GR6J": [0, 1, 2] + list(range(7, 67)),     # plus the exponential store
}


@pytest.mark.parametrize("model", sorted(DAILY_MODELS))
def test_daily_model_outputs(daily, model):
    run, param = DAILY_MODELS[model]
    options = gr.RunOptions(daily["inputs"], model,
                            ind_period_run=daily["ind_run"])
    out = run(daily["inputs"], options, param)
    reference = load("sim_" + model)

    for name in MODEL_OUTPUTS[model]:
        assert_matches(reference[name].to_numpy(), out[name],
                       "%s / %s" % (model, name))


@pytest.mark.parametrize("model", sorted(DAILY_MODELS))
def test_daily_model_final_state(daily, model):
    run, param = DAILY_MODELS[model]
    options = gr.RunOptions(daily["inputs"], model,
                            ind_period_run=daily["ind_run"])
    out = run(daily["inputs"], options, param)
    reference = load("state_" + model)["state"].to_numpy()
    slots = STATE_SLOTS[model]
    assert_matches(reference[slots], out["StateEnd"][slots],
                   "%s / StateEnd" % model)


def test_monthly_model_gr2m():
    basin = load("basin_monthly")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    options = gr.RunOptions(inputs, "GR2M",
                            ind_period_run=index_range("idx_monthly"))
    out = gr.run_model_gr2m(inputs, options, [265.072, 1.007])
    reference = load("sim_GR2M")
    for name in MODEL_OUTPUTS["GR2M"]:
        assert_matches(reference[name].to_numpy(), out[name], "GR2M / " + name)


def test_yearly_model_gr1a():
    basin = load("basin_yearly")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    inputs = gr.InputsModel(dates, basin["P"].to_numpy(), basin["E"].to_numpy())
    options = gr.RunOptions(inputs, "GR1A",
                            ind_period_run=index_range("idx_yearly"))
    out = gr.run_model_gr1a(inputs, options, [0.840])
    reference = load("sim_GR1A")
    for name in MODEL_OUTPUTS["GR1A"]:
        assert_matches(reference[name].to_numpy(), out[name], "GR1A / " + name)


def test_hourly_model_gr4h(hourly):
    options = gr.RunOptions(hourly["inputs"], "GR4H",
                            ind_period_run=hourly["ind_run"])
    out = gr.run_model_gr4h(hourly["inputs"], options,
                            [756.930, -0.773, 138.638, 5.247])
    reference = load("sim_GR4H")
    for name in MODEL_OUTPUTS["GR4H"]:
        assert_matches(reference[name].to_numpy(), out[name], "GR4H / " + name)


def test_hourly_model_gr5h(hourly):
    options = gr.RunOptions(hourly["inputs"], "GR5H",
                            ind_period_run=hourly["ind_run"])
    out = gr.run_model_gr5h(hourly["inputs"], options,
                            [756.930, -0.773, 138.638, 5.247, 0.400])
    reference = load("sim_GR5H")
    for name in MODEL_OUTPUTS["GR5H"]:
        assert_matches(reference[name].to_numpy(), out[name], "GR5H / " + name)


def test_hourly_model_gr5h_with_interception(hourly):
    imax = float(load("imax_value")["imax"][0])
    options = gr.RunOptions(hourly["inputs"], "GR5H",
                            ind_period_run=hourly["ind_run"], imax=imax)
    out = gr.run_model_gr5h(hourly["inputs"], options,
                            [756.930, -0.773, 138.638, 5.247, 0.400])
    reference = load("sim_GR5H_interception")
    for name in MODEL_OUTPUTS["GR5H"]:
        assert_matches(reference[name].to_numpy(), out[name],
                       "GR5H interception / " + name)


def test_interception_store_changes_the_simulation(hourly):
    """A sanity check that Imax is actually wired through, not ignored."""
    param = [756.930, -0.773, 138.638, 5.247, 0.400]
    plain = gr.run_model_gr5h(
        hourly["inputs"],
        gr.RunOptions(hourly["inputs"], "GR5H", ind_period_run=hourly["ind_run"]),
        param)
    intercepted = gr.run_model_gr5h(
        hourly["inputs"],
        gr.RunOptions(hourly["inputs"], "GR5H", ind_period_run=hourly["ind_run"],
                      imax=0.7),
        param)
    assert not np.allclose(plain["Qsim"], intercepted["Qsim"])
    assert np.all(intercepted["Interc"] <= 0.7 + 1e-12)


def test_pe_oudin():
    reference = load("pe_oudin")
    computed = gr.pe_oudin(reference["JD"].to_numpy(),
                           reference["Temp"].to_numpy(), 0.8)
    assert_matches(reference["PE"].to_numpy(), computed, "PE_Oudin")


def test_pe_oudin_accepts_degrees():
    jd = np.arange(1, 366, dtype=float)
    temp = 10.0 + 8.0 * np.sin(jd / 58.1)
    in_radians = gr.pe_oudin(jd, temp, 0.8, lat_unit="rad")
    in_degrees = gr.pe_oudin(jd, temp, 0.8 * 180.0 / np.pi, lat_unit="deg")
    assert_matches(in_radians, in_degrees, "PE_Oudin degrees vs radians")
