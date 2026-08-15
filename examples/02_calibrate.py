"""Calibrate a model, then check it on a period it has never seen.

The split-sample test of Klemes (1986) is the standard way to find out whether a
calibrated model transfers, rather than memorises.

    python examples/02_calibrate.py
"""

from _data import daily_catchment


def main():
    catchment = daily_catchment()

    fit, validation = catchment.split_sample(
        "GR4J",
        calibration=("1990-01-01", "1999-12-31"),
        validation=("2000-01-01", "2009-12-31"),
        criterion="KGE")

    print(fit.describe())

    print("\nCalibration period")
    print("  ", fit.simulate(period=("1990-01-01", "1999-12-31")))
    print("Validation period (never seen during calibration)")
    print("  ", validation)

    # Parameters are available by name, not just by position.
    print("\nParameters:", {k: round(v, 3)
                            for k, v in fit.named_params().items()})

    # A different objective function moves the parameters somewhere else.
    print("\nSame model, different objective:")
    for criterion, transfo in [("NSE", ""), ("KGE", ""), ("NSE", "log"),
                               ("NSE", "sqrt")]:
        other = catchment.calibrate("GR4J", criterion=criterion,
                                    transfo=transfo,
                                    period=("1990-01-01", "1999-12-31"))
        check = other.evaluate(period=("2000-01-01", "2009-12-31"))
        label = "%s[%s]" % (criterion, transfo or "-")
        print("  %-10s  X1=%7.1f  X3=%6.1f   validation NSE=%.3f  KGE=%.3f"
              % (label, other.params[0], other.params[2],
                 check.nse(), check.kge()))


if __name__ == "__main__":
    main()
