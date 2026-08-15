"""Shared fixtures: airGR reference data and equality assertions.

Every reference file under ``tests/data`` was produced by airGR 1.7.9 itself,
run in R on the demonstration catchments that package ships (see
``docs/VALIDATION.md``). The test suite asserts that GRsuite reproduces those
numbers, so a failure means a real divergence from the reference
implementation.
"""

import gzip
import io
import os

import numpy as np
import pandas as pd
import pytest

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

#: Reference outputs are written with 17 significant digits, so agreement is
#: expected down to double-precision round-off, far below this tolerance.
TOL = 1e-9


def load(name):
    """Read a gzipped reference table shipped with the test suite."""
    path = os.path.join(DATA, name + ".csv.gz")
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return pd.read_csv(io.StringIO(fh.read()))


def assert_matches(reference, computed, label, tol=TOL):
    """Assert a computed series matches airGR, elementwise.

    Passes when either the absolute or the relative deviation is within
    ``tol``; the pair matters because several model outputs legitimately
    cross zero, where a relative test is meaningless.
    """
    ref = np.asarray(reference, dtype=float)
    got = np.asarray(computed, dtype=float)
    assert ref.shape == got.shape, (
        "%s: shape mismatch, airGR %s vs GRsuite %s"
        % (label, ref.shape, got.shape))

    nan_mismatch = int((np.isnan(ref) != np.isnan(got)).sum())
    assert nan_mismatch == 0, (
        "%s: %i time steps differ in missing-value placement"
        % (label, nan_mismatch))

    ok = ~(np.isnan(ref) | np.isnan(got))
    if not ok.any():
        return
    diff = np.abs(ref[ok] - got[ok])
    rel = diff / np.maximum(np.abs(ref[ok]), 1e-300)
    within = (diff <= tol) | (rel <= tol)
    if within.all():
        return

    bad = np.where(~within)[0]
    i = bad[int(np.argmax(diff[bad]))]
    raise AssertionError(
        "%s: %i / %i values outside tolerance %.0e\n"
        "  worst at index %i: airGR %.17g vs GRsuite %.17g "
        "(absolute %.3e, relative %.3e)"
        % (label, bad.size, ok.sum(), tol, i, ref[ok][i], got[ok][i],
           diff[i], rel[i]))


def index_range(name):
    """Convert a 1-based R index range into a 0-based numpy range."""
    idx = load(name)
    return np.arange(int(idx["ind_start"][0]) - 1, int(idx["ind_end"][0]))


# ---------------------------------------------------------------------------
# Catchment fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def daily():
    """Daily demonstration catchment (airGR L0123001)."""
    import grsuite as gr
    basin = load("basin_daily")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    return {
        "dates": dates,
        "P": basin["P"].to_numpy(),
        "E": basin["E"].to_numpy(),
        "T": basin["T"].to_numpy(),
        "Q": load("qobs_daily")["Qmm"].to_numpy(),
        "ind_run": index_range("idx_daily"),
        "inputs": gr.InputsModel(dates, basin["P"].to_numpy(),
                                 basin["E"].to_numpy(),
                                 temp_mean=basin["T"].to_numpy()),
    }


@pytest.fixture(scope="session")
def hourly():
    """Hourly demonstration catchment (airGR L0123003)."""
    import grsuite as gr
    basin = load("basin_hourly")
    dates = pd.to_datetime(basin["Date"]).to_numpy().astype("datetime64[h]")
    return {
        "dates": dates,
        "ind_run": index_range("idx_hourly"),
        "inputs": gr.InputsModel(dates, basin["P"].to_numpy(),
                                 basin["E"].to_numpy()),
    }


@pytest.fixture(scope="session")
def snow():
    """Snow-affected catchment with five elevation bands (airGR L0123002)."""
    import grsuite as gr
    basin = load("basin_snow")
    dates = basin["Date"].to_numpy().astype("datetime64[D]")
    hypso = load("hypso_snow")["hypso"].to_numpy()
    meta = load("idx_snow")
    z_inputs = float(meta["zinputs"][0])

    alti = gr.data_alti_extrapolation_valery(
        dates, basin["P"].to_numpy(), basin["T"].to_numpy(),
        z_inputs=z_inputs, hypso_data=hypso, n_layers=5)

    return {
        "dates": dates,
        "hypso": hypso,
        "z_inputs": z_inputs,
        "alti": alti,
        "ind_run": np.arange(int(meta["ind_start"][0]) - 1,
                             int(meta["ind_end"][0])),
        "inputs": gr.InputsModel(
            dates, basin["P"].to_numpy(), basin["E"].to_numpy(),
            temp_mean=basin["T"].to_numpy(),
            layer_precip=alti["LayerPrecip"],
            layer_temp_mean=alti["LayerTempMean"],
            layer_frac_solid_precip=alti["LayerFracSolidPrecip"],
            z_layers=alti["ZLayers"]),
    }
