"""From a plain CSV file to a calibrated model.

The only format GRsuite needs is three columns — date, precipitation, potential
evapotranspiration — plus observed discharge if you want to calibrate. This
example starts from such a CSV and ends with a validated GR4J.

    python examples/08_from_csv.py

Your own file only needs to look like this (any column order, any step):

    date,P,E,Q
    2000-01-01,4.2,0.1,0.83
    2000-01-02,0.0,0.2,0.80
    ...

To keep this example runnable straight after a clone, it reads the
demonstration catchment bundled in ``tests/data`` — a gzipped CSV, which
pandas reads transparently; a plain ``.csv`` works the same way.
"""

import os

import pandas as pd

import grsuite as gr

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO_CSV = os.path.join(HERE, "..", "tests", "data", "basin_daily.csv.gz")
DEMO_QOBS = os.path.join(HERE, "..", "tests", "data", "qobs_daily.csv.gz")


def main():
    # 1. Read the CSV. With your own file: df = pd.read_csv("my_catchment.csv")
    df = pd.read_csv(DEMO_CSV)
    qobs = pd.read_csv(DEMO_QOBS)

    # 2. Build the catchment: dates as datetime64, series in mm per time step.
    catchment = gr.Catchment(
        df["Date"].to_numpy().astype("datetime64[D]"),
        precip=df["P"].to_numpy(),
        pot_evap=df["E"].to_numpy(),
        obs_discharge=qobs["Qmm"].to_numpy(),
        name="From-CSV catchment")
    print(catchment)

    # 3. Calibrate on the first ten years, validate on the next ten.
    fit, validation = catchment.split_sample(
        "GR4J",
        calibration=("1990-01-01", "1999-12-31"),
        validation=("2000-01-01", "2009-12-31"),
        criterion="KGE")

    print("\n", fit.describe(), sep="")
    print("\nValidation period:", validation)


if __name__ == "__main__":
    main()
