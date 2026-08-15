"""CemaNeige: elevation banding, snow module, and coupling to the GR models."""

import numpy as np
import pytest

import grsuite as gr
from conftest import assert_matches, load
from grsuite.core import MODEL_OUTPUTS

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


@pytest.mark.parametrize("layer", range(1, 6))
def test_elevation_extrapolation(snow, layer):
    """Precipitation, temperature and solid fraction per elevation band."""
    reference = load("inputs_layer%02i" % layer)
    alti = snow["alti"]
    i = layer - 1
    assert_matches(reference["P"].to_numpy(), alti["LayerPrecip"][i],
                   "layer %i precipitation" % layer)
    assert_matches(reference["T"].to_numpy(), alti["LayerTempMean"][i],
                   "layer %i temperature" % layer)
    assert_matches(reference["FS"].to_numpy(),
                   alti["LayerFracSolidPrecip"][i],
                   "layer %i solid fraction" % layer)


def test_elevation_bands_are_ordered(snow):
    """Bands run from the lowest to the highest elevation."""
    z = snow["alti"]["ZLayers"]
    assert len(z) == 5
    assert np.all(np.diff(z) > 0)


@pytest.mark.parametrize("case", sorted(COUPLED))
def test_coupled_model_outputs(snow, case):
    model, gr_name, hysteresis, param = COUPLED[case]
    options = gr.RunOptions(snow["inputs"], model,
                            ind_period_run=snow["ind_run"],
                            is_hyst=hysteresis)
    out = gr.run_model(snow["inputs"], options, param, model=model)
    reference = load("sim_" + case)
    for name in MODEL_OUTPUTS[gr_name]:
        assert_matches(reference[name].to_numpy(), out[name],
                       "%s / %s" % (case, name))


@pytest.mark.parametrize("case", sorted(COUPLED))
def test_mean_annual_solid_precipitation(snow, case):
    model, _, hysteresis, param = COUPLED[case]
    options = gr.RunOptions(snow["inputs"], model,
                            ind_period_run=snow["ind_run"],
                            is_hyst=hysteresis)
    gr.run_model(snow["inputs"], options, param, model=model)
    reference = float(load("masp_" + case)["masp"][0])
    assert_matches(np.array([reference]),
                   np.array([options.mean_an_solid_precip[0]]),
                   "%s / MeanAnSolidPrecip" % case)


@pytest.mark.parametrize("case", ["CemaNeigeGR4J", "CemaNeigeGR4J_Hyst"])
@pytest.mark.parametrize("layer", range(1, 6))
def test_snow_layer_series(snow, case, layer):
    """Snow pack, thermal state, melt and Gratio, band by band."""
    model, _, hysteresis, param = COUPLED[case]
    options = gr.RunOptions(snow["inputs"], model,
                            ind_period_run=snow["ind_run"],
                            is_hyst=hysteresis)
    out = gr.run_model(snow["inputs"], options, param, model=model)
    reference = load("sim_%s_layer%02i" % (case, layer))
    for name in reference.columns:
        assert_matches(reference[name].to_numpy(),
                       out["CemaNeigeLayers"][layer - 1][name],
                       "%s / band %i %s" % (case, layer, name))


def test_hysteresis_changes_the_snow_dynamics(snow):
    """The hysteresis variant must not silently behave like the plain one."""
    base = gr.run_model(
        snow["inputs"],
        gr.RunOptions(snow["inputs"], "CemaNeigeGR4J",
                      ind_period_run=snow["ind_run"]),
        [408.774, 2.646, 131.264, 1.174, 0.962, 2.249], model="CemaNeigeGR4J")
    hyst = gr.run_model(
        snow["inputs"],
        gr.RunOptions(snow["inputs"], "CemaNeigeGR4J",
                      ind_period_run=snow["ind_run"], is_hyst=True),
        [408.774, 2.646, 131.264, 1.174, 0.962, 2.249, 80.0, 0.4],
        model="CemaNeigeGR4J")
    assert not np.allclose(base["Qsim"], hyst["Qsim"])


def test_snow_module_alone(snow):
    """CemaNeige can run without a rainfall-runoff model attached."""
    options = gr.RunOptions(snow["inputs"], "CemaNeigeGR4J",
                            ind_period_run=snow["ind_run"])
    out = gr.run_model_cemaneige(snow["inputs"], options, [0.962, 2.249])
    assert len(out["CemaNeigeLayers"]) == 5
    assert out["PliqAndMelt"].shape == snow["ind_run"].shape
    assert np.all(out["PliqAndMelt"] >= 0)
