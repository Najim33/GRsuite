"""Data assimilation reproduces airGRdatassim, value by value.

The ``tests/data/da_*`` reference files were produced by R running the
instrumented copies of airGRdatassim 0.1.4 in ``tools/oracle/
da_instrumented.R`` (see ``tools/oracle/export_da_fixtures.R``). R's random
draws cannot be reproduced by numpy, so the exact draws are exported too and
replayed here through ``_ReplayRNG``: a failure means a real divergence from
the reference, including in the sequence of random calls (each test also
checks that the exported draws are consumed exactly).

Note that airGRdatassim 0.1.4's EnKF is a silent no-op with the documented
character ``StateEnKF``; GRsuite implements the documented behaviour (the
GRSUITE-FIX marked in both code bases), so the EnKF references are produced
with that same one-line fix applied. See docs/VALIDATION.md.
"""

import numpy as np
import pytest

import grsuite as gr
from conftest import assert_matches, load
from grsuite._rng import _ReplayRNG

DA_CONFIGS = ["enkf_GR4J", "enkf_GR5J", "enkf_GR6J", "pf_GR4J", "none_GR4J",
              "enkf_CemaNeigeGR4J", "pf_CemaNeigeGR6J"]


def _draws(name):
    df = load(name)
    return [(g["kind"].iloc[0], g["value"].to_numpy())
            for _, g in df.groupby("call", sort=True)]


def _meta(name):
    return load(name).iloc[0].to_dict()


def _split(value):
    s = str(value)
    return None if s in ("", "nan") else s.split("|")


def _inputs(basin):
    if basin == "daily":
        b = load("basin_daily")
        dates = b["Date"].to_numpy().astype("datetime64[D]")
        return gr.InputsModel(dates, b["P"].to_numpy(), b["E"].to_numpy(),
                              temp_mean=b["T"].to_numpy())
    b = load("basin_snow")
    hypso = load("hypso_snow")["hypso"].to_numpy()
    z_inputs = float(load("idx_snow")["zinputs"][0])
    dates = b["Date"].to_numpy().astype("datetime64[D]")
    alti = gr.data_alti_extrapolation_valery(
        dates, b["P"].to_numpy(), b["T"].to_numpy(),
        z_inputs=z_inputs, hypso_data=hypso, n_layers=5)
    return gr.InputsModel(
        dates, b["P"].to_numpy(), b["E"].to_numpy(),
        temp_mean=b["T"].to_numpy(),
        layer_precip=alti["LayerPrecip"],
        layer_temp_mean=alti["LayerTempMean"],
        layer_frac_solid_precip=alti["LayerFracSolidPrecip"],
        z_layers=alti["ZLayers"])


@pytest.fixture(scope="module")
def daily():
    b = load("basin_daily")
    return _inputs("daily"), b["Qmm"].to_numpy()


@pytest.fixture(scope="module")
def snow():
    b = load("basin_snow")
    return _inputs("snow"), b["Qmm"].to_numpy()


# ---------------------------------------------------------------------------
# Reference comparisons (R-produced fixtures, replayed draws)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cfg", ["GR4J", "CemaNeigeGR4J"])
def test_inputs_pert_matches_reference(cfg, daily, snow):
    inputs, _ = (daily if cfg == "GR4J" else snow)
    m = _meta("da_pert_%s_meta" % cfg)
    i0, i1 = int(m["ind_start"]) - 1, int(m["ind_end"])
    kw = {}
    if cfg != "GR4J":
        b = load("basin_snow")
        kw = {"temp_mean": b["T"].to_numpy()[i0:i1],
              "z_inputs": float(load("idx_snow")["zinputs"][0]),
              "hypso_data": load("hypso_snow")["hypso"].to_numpy(),
              "n_layers": 5}
    rng = _ReplayRNG(_draws("da_pert_%s_draws" % cfg))
    ip = gr.InputsPert(m["model"], inputs.dates[i0:i1],
                       precip=inputs.precip[i0:i1],
                       pot_evap=inputs.pot_evap[i0:i1],
                       nb_mbr=int(m["nb_mbr"]), seed=int(m["seed"]),
                       _rng=rng, **kw)
    for var, attr in (("precip", "precip"), ("potevap", "pot_evap")):
        ref = load("da_pert_%s_%s" % (cfg, var))
        got = getattr(ip, attr)
        assert got is not None, "%s: %s should be perturbed" % (cfg, attr)
        assert_matches(ref.drop(columns=["Date"]).to_numpy(), got,
                       "%s / %s" % (cfg, var))
    assert rng.exhausted(), "%s: not all exported draws were consumed" % cfg


@pytest.mark.parametrize("cfg", DA_CONFIGS)
def test_run_model_da_matches_reference(cfg, daily, snow):
    m = _meta("da_%s_meta" % cfg)
    inputs, qmm = (daily if m["basin"] == "daily" else snow)
    ind_run = np.arange(int(m["ind_start"]) - 1, int(m["ind_end"]))
    param = [float(x) for x in str(m["param"]).split("|")]
    rng = _ReplayRNG(_draws("da_%s_draws" % cfg))
    out = gr.run_model_da(
        inputs, ind_run, m["model"], param,
        qobs=qmm, da_method=m["da_method"], nb_mbr=int(m["nb_mbr"]),
        state_enkf=_split(m["state_enkf"]),
        state_pert=_split(m["state_pert"]),
        seed=int(m["seed"]), _rng=rng)

    ref_q = load("da_%s_qsimens" % cfg)
    assert_matches(ref_q.drop(columns=["Date"]).to_numpy(), out["QsimEns"],
                   cfg + " / QsimEns")
    for part, key in (("bkg", "EnsStateBkg"), ("ana", "EnsStateA")):
        ref = load("da_%s_%s" % (cfg, part))
        for k, sname in enumerate(out["StateNames"]):
            assert_matches(ref[sname].to_numpy(), out[key][:, :, k].ravel(),
                           "%s / %s / %s" % (cfg, key, sname))
    if m["da_method"] == "EnKF":
        ref_o = load("da_%s_obspert" % cfg)
        assert_matches(ref_o.drop(columns=["Date"]).to_numpy(),
                       out["ObsPert"], cfg + " / ObsPert")
    assert rng.exhausted(), "%s: not all exported draws were consumed" % cfg


def test_run_model_da_meteo_matches_reference(daily):
    cfg = "enkfmet_GR4J"
    inputs, qmm = daily
    m = _meta("da_%s_meta" % cfg)
    w0, w1 = int(m["win_start"]) - 1, int(m["win_end"])
    ind_run = np.arange(int(m["ind_start"]) - 1, int(m["ind_end"]))
    param = [float(x) for x in str(m["param"]).split("|")]
    sub = gr.InputsModel(inputs.dates[w0:w1], inputs.precip[w0:w1],
                         inputs.pot_evap[w0:w1])

    rng_p = _ReplayRNG(_draws("da_%s_pdraws" % cfg))
    ip = gr.InputsPert("GR4J", sub.dates, precip=sub.precip,
                       pot_evap=sub.pot_evap, nb_mbr=int(m["nb_mbr"]),
                       seed=int(m["pert_seed"]), _rng=rng_p)
    for var, attr in (("precip", "precip"), ("potevap", "pot_evap")):
        ref = load("da_%s_%s" % (cfg, var))
        assert_matches(ref.drop(columns=["Date"]).to_numpy(),
                       getattr(ip, attr), "%s / %s" % (cfg, var))
    assert rng_p.exhausted()

    rng = _ReplayRNG(_draws("da_%s_draws" % cfg))
    out = gr.run_model_da(
        sub, ind_run, "GR4J", param, inputs_pert=ip,
        qobs=qmm[w0:w1], da_method="EnKF", nb_mbr=int(m["nb_mbr"]),
        state_enkf=["Prod", "Rout", "UH1", "UH2"],
        state_pert=["Prod", "Rout"], seed=int(m["seed"]), _rng=rng)

    ref_q = load("da_%s_qsimens" % cfg)
    assert_matches(ref_q.drop(columns=["Date"]).to_numpy(), out["QsimEns"],
                   cfg + " / QsimEns")
    for part, key in (("bkg", "EnsStateBkg"), ("ana", "EnsStateA")):
        ref = load("da_%s_%s" % (cfg, part))
        for k, sname in enumerate(out["StateNames"]):
            assert_matches(ref[sname].to_numpy(), out[key][:, :, k].ravel(),
                           "%s / %s / %s" % (cfg, key, sname))
    ref_o = load("da_%s_obspert" % cfg)
    assert_matches(ref_o.drop(columns=["Date"]).to_numpy(), out["ObsPert"],
                   cfg + " / ObsPert")
    assert rng.exhausted()


# ---------------------------------------------------------------------------
# Self-contained behaviour checks (no R reference needed)
# ---------------------------------------------------------------------------

PARAM4J = [257.238, 1.012, 88.235, 2.208]


def test_open_loop_matches_plain_run(daily):
    """DaMethod "none" is the deterministic model, member by member."""
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 120)
    out = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                          da_method="none", nb_mbr=5, seed=1)
    options = gr.RunOptions(inputs, "GR4J", ind_period_run=ind_run,
                            ind_period_warmup=[])
    ref = gr.run_model_gr4j(inputs, options, PARAM4J)
    for m in range(out["NbMbr"]):
        assert_matches(ref["Qsim"], out["QsimEns"][:, m],
                       "open loop / member %i" % m, tol=1e-12)


def test_first_step_is_degenerate(daily):
    """As in R, all members share the same unperturbed first step."""
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 30)
    ip = gr.InputsPert("GR4J", inputs.dates, precip=inputs.precip,
                       pot_evap=inputs.pot_evap, nb_mbr=6, seed=3)
    out = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, inputs_pert=ip,
                          qobs=qmm, da_method="PF", seed=3)
    for k in range(out["NbState"]):
        assert np.ptp(out["EnsStateBkg"][0, :, k]) == 0.0
    assert np.ptp(out["QsimEns"][0, :]) == 0.0


def test_enkf_actually_updates_states(daily):
    """Regression guard: the EnKF update must not be a no-op (GRSUITE-FIX).

    Without any perturbation the ensemble stays degenerate (identical
    members), the Kalman gain is zero and the update is legitimately null;
    ``state_pert`` is what keeps the ensemble spread, as in airGRdatassim's
    documented usage.
    """
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 30)
    out = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                          da_method="EnKF", nb_mbr=6,
                          state_enkf=["Prod", "Rout", "UH1", "UH2"],
                          state_pert=["Prod", "Rout"], seed=1)
    assert np.nanmax(np.abs(out["EnsStateA"] - out["EnsStateBkg"])) > 0.0


def test_gr5j_has_no_uh1(daily):
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 10)
    out = gr.run_model_da(inputs, ind_run, "GR5J",
                          [245.918, 1.027, 90.017, 2.198, 0.318], qobs=qmm,
                          da_method="EnKF", nb_mbr=4,
                          state_enkf=["Prod", "Rout", "UH2"], seed=1)
    assert out["StateNames"] == ("Prod", "Rout", "UH2")


def test_all_nan_qobs_falls_back_to_none(daily):
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 30)
    nan_q = np.full_like(qmm, np.nan)
    with pytest.warns(UserWarning, match="none"):
        out = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=nan_q,
                              da_method="EnKF", nb_mbr=5,
                              state_enkf=["Prod", "Rout"], seed=1)
    ref = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                          da_method="none", nb_mbr=5, seed=1)
    assert_matches(ref["QsimEns"], out["QsimEns"], "fallback / QsimEns",
                   tol=1e-12)


def test_negative_qobs_is_nan(daily):
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 30)
    qmod = qmm.copy()
    qmod[ind_run[10]] = -1.0
    with pytest.warns(UserWarning, match="negative"):
        out = gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmod,
                              da_method="EnKF", nb_mbr=5,
                              state_enkf=["Prod", "Rout"], seed=1)
    # no assimilation at that step: analysis equals background
    assert_matches(out["EnsStateBkg"][10], out["EnsStateA"][10],
                   "negative qobs / step without assimilation", tol=1e-12)


def test_run_model_da_errors(daily):
    inputs, qmm = daily
    ind_run = np.arange(3653, 3653 + 10)
    with pytest.raises(ValueError, match="state_enkf"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                        da_method="EnKF")
    with pytest.raises(ValueError, match="state_pert"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                        da_method="EnKF", nb_mbr=4,
                        state_enkf=["Prod", "Rout"], state_pert=["UH1"])
    with pytest.raises(ValueError, match="nb_mbr"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                        da_method="none", nb_mbr=1)
    with pytest.raises(ValueError, match="assimilation"):
        gr.run_model_da(inputs, ind_run, "GR2M", [265.072, 1.007], qobs=qmm,
                        da_method="none", nb_mbr=4)
    with pytest.raises(ValueError, match="inconnue"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, qobs=qmm,
                        da_method="EnKF", nb_mbr=4, state_enkf=["Store"])
    ip = gr.InputsPert("GR4J", inputs.dates[:100],
                       precip=inputs.precip[:100],
                       pot_evap=inputs.pot_evap[:100], nb_mbr=4)
    with pytest.raises(ValueError, match="longueur"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J, inputs_pert=ip,
                        qobs=qmm, da_method="PF")
    with pytest.raises(ValueError, match="superieur"):
        gr.run_model_da(inputs, ind_run, "GR4J", PARAM4J,
                        inputs_pert=gr.InputsPert(
                            "GR4J", inputs.dates, precip=inputs.precip,
                            pot_evap=inputs.pot_evap, nb_mbr=4),
                        qobs=qmm, da_method="PF", nb_mbr=5)


def test_inputs_pert_errors(daily):
    inputs, _ = daily
    with pytest.raises(ValueError, match="nb_mbr"):
        gr.InputsPert("GR4J", inputs.dates, precip=inputs.precip, nb_mbr=1)
    with pytest.raises(ValueError, match="absents"):
        gr.InputsPert("GR4J", inputs.dates)
    with pytest.raises(ValueError, match="assimilation"):
        gr.InputsPert("GR2M", inputs.dates, precip=inputs.precip)
