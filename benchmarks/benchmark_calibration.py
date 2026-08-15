"""Calibration benchmark: time Michel's algorithm on the bundled catchment.

The numbers in the README (5.4 min vs airGR's 17.5 min) come from a campaign
of 500 calibrations on CAMELS-FR catchments, which are too large to ship with
the repository. This script reproduces the *measurement* on the demonstration
catchment bundled in ``tests/data``: each daily model is calibrated
``--repeat`` times on the same 20-year period, and the wall clock per
calibration is reported.

    python benchmarks/benchmark_calibration.py [--repeat 5]

The first calibration of each model pays Numba's JIT compilation (a second or
two); it is excluded from the statistics. Compilation is cached on disk, so a
second run of this script has no compilation cost at all.
"""

import argparse
import gzip
import io
import os
import time

import pandas as pd

import grsuite as gr

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "tests", "data")
MODELS = ["GR4J", "GR5J", "GR6J"]
PERIOD = ("1990-01-01", "2009-12-31")


def demo_catchment():
    with gzip.open(os.path.join(DATA, "basin_daily.csv.gz"), "rt",
                   encoding="utf-8") as fh:
        basin = pd.read_csv(io.StringIO(fh.read()))
    with gzip.open(os.path.join(DATA, "qobs_daily.csv.gz"), "rt",
                   encoding="utf-8") as fh:
        qobs = pd.read_csv(io.StringIO(fh.read()))
    return gr.Catchment(
        basin["Date"].to_numpy().astype("datetime64[D]"),
        precip=basin["P"].to_numpy(),
        pot_evap=basin["E"].to_numpy(),
        obs_discharge=qobs["Qmm"].to_numpy())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeat", type=int, default=5,
                        help="calibrations per model (default: 5)")
    args = parser.parse_args()

    catchment = demo_catchment()
    n_steps = catchment.period_index(PERIOD).size
    print(f"{args.repeat} calibrations per model on {PERIOD[0]} to {PERIOD[1]} "
          f"({n_steps} steps)\n")
    print(f"{'model':<8} {'runs':>10} {'total [s]':>12} {'per cal [s]':>12}")
    print("-" * 46)

    for model in MODELS:
        # Warm-up: pays the one-off Numba compilation, excluded from timings.
        catchment.calibrate(model, period=PERIOD)

        tic = time.perf_counter()
        for _ in range(args.repeat):
            catchment.calibrate(model, period=PERIOD)
        total = time.perf_counter() - tic

        print(f"{model:<8} {args.repeat:>10} {total:>12.2f} "
              f"{total / args.repeat:>12.3f}")


if __name__ == "__main__":
    main()
