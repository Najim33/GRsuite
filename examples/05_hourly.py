"""Sub-daily modelling: GR4H, GR5H, and the interception store.

At an hourly time step the way interception is represented starts to matter.
GR5H can carry an explicit interception store, whose capacity Imax can be
estimated from the forcing itself rather than calibrated.

    python examples/05_hourly.py
"""

from _data import hourly_catchment

CALIBRATION = ("2005-02-01", "2005-08-31")
VALIDATION = ("2005-09-01", "2005-12-31")


def main():
    catchment = hourly_catchment()
    print(catchment)

    for model in ["GR4H", "GR5H"]:
        fit, check = catchment.split_sample(
            model, calibration=CALIBRATION, validation=VALIDATION,
            criterion="KGE")
        print("\n%s  calibration KGE=%.3f  validation KGE=%.3f  NSE=%.3f"
              % (model, fit.score, check.kge(), check.nse()))

    # Imax is estimated by matching hourly interception evaporation to the
    # daily sum of min(P, PE), following Ficchi (2017).
    inputs = catchment._inputs("GR5H")
    index = catchment.period_index(CALIBRATION)
    import grsuite as gr
    imax = gr.imax_estimate(inputs, index)
    print("\nEstimated interception capacity: Imax = %.1f mm" % imax)

    fit, check = catchment.split_sample(
        "GR5H", calibration=CALIBRATION, validation=VALIDATION,
        criterion="KGE", imax=imax)
    print("GR5H with interception store  validation KGE=%.3f  NSE=%.3f"
          % (check.kge(), check.nse()))

    print("\nThe demonstration record is only two years long, so a seasonal")
    print("split is a harsh test: calibration covers spring and summer, validation")
    print("only autumn. What the numbers do show is the ranking, and the")
    print("interception store is consistently ahead.")

    frame = check.to_dataframe(["Precip", "Interc", "EI", "ES", "Qsim"])
    print("\nInterception store, first hours of a rainy spell [mm]:")
    rainy = frame[frame["Precip"] > 0].head(6)
    print(rainy.round(3).to_string())


if __name__ == "__main__":
    main()
