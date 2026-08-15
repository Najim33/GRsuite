"""Which GR model is worth the extra parameters on this catchment?

Calibrates every daily model on the same period with the same objective, then
scores them on an independent one.

    python examples/03_compare_models.py
"""

from _data import daily_catchment

CALIBRATION = ("1990-01-01", "1999-12-31")
VALIDATION = ("2000-01-01", "2009-12-31")


def main():
    catchment = daily_catchment()
    rows = []

    for model in ["GR4J", "GR5J", "GR6J"]:
        fit, check = catchment.split_sample(
            model, calibration=CALIBRATION, validation=VALIDATION,
            criterion="KGE")
        rows.append((model, len(fit.params), fit.score, check.nse(),
                     check.kge(), check.nse("log"), check.bias()))

    header = ("model", "params", "cal KGE", "val NSE", "val KGE",
              "val NSE_log", "val bias")
    print("%-7s %7s %9s %9s %9s %12s %10s" % header)
    print("-" * 68)
    for r in rows:
        print("%-7s %7d %9.3f %9.3f %9.3f %12.3f %9.1f%%"
              % (r[0], r[1], r[2], r[3], r[4], r[5], 100 * r[6]))

    best = max(rows, key=lambda r: r[4])
    print("\nBest validation KGE: %s (%.3f)" % (best[0], best[4]))
    print("NSE_log rewards low-flow performance, where GR6J's exponential store")
    print("usually earns its keep.")


if __name__ == "__main__":
    main()
