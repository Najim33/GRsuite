"""Load the demonstration catchments bundled with the test suite.

Every example in this folder runs straight after a clone: the data comes from
``tests/data``, which ships airGR's own demonstration catchments.
"""

import gzip
import io
import os

import numpy as np
import pandas as pd

import grsuite as gr

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "tests", "data")


def _read(name):
    with gzip.open(os.path.join(DATA, name + ".csv.gz"), "rt",
                   encoding="utf-8") as fh:
        return pd.read_csv(io.StringIO(fh.read()))


def daily_catchment():
    """A temperate lowland catchment, daily time step, 29 years."""
    basin = _read("basin_daily")
    return gr.Catchment(
        basin["Date"].to_numpy().astype("datetime64[D]"),
        precip=basin["P"].to_numpy(),
        pot_evap=basin["E"].to_numpy(),
        obs_discharge=_read("qobs_daily")["Qmm"].to_numpy(),
        temperature=basin["T"].to_numpy(),
        name="Demonstration catchment (daily)")


def snow_catchment():
    """A mountainous catchment with a hypsometric curve, daily time step."""
    basin = _read("basin_snow")
    return gr.Catchment(
        basin["Date"].to_numpy().astype("datetime64[D]"),
        precip=basin["P"].to_numpy(),
        pot_evap=basin["E"].to_numpy(),
        obs_discharge=basin["Qmm"].to_numpy(),
        temperature=basin["T"].to_numpy(),
        hypsometry=_read("hypso_snow")["hypso"].to_numpy(),
        n_layers=5,
        name="Demonstration catchment (snow)")


def hourly_catchment():
    """A flashy catchment at hourly time step."""
    basin = _read("basin_hourly")
    return gr.Catchment(
        pd.to_datetime(basin["Date"]).to_numpy().astype("datetime64[h]"),
        precip=basin["P"].to_numpy(),
        pot_evap=basin["E"].to_numpy(),
        obs_discharge=basin["Qmm"].to_numpy(),
        name="Demonstration catchment (hourly)")


def upstream_discharge(catchment, areas):
    """Synthetic upstream inflows [m3/step], for the routing example."""
    q = np.nan_to_num(catchment.obs_discharge)
    return np.column_stack([q * areas[0] * 1e3 * 0.9,
                            q * areas[1] * 1e3 * 1.1])
