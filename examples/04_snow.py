"""Snow: CemaNeige on elevation bands.

CemaNeige splits the catchment into elevation bands, extrapolates precipitation
and temperature to each one, and runs an independent snow pack per band. The
liquid water it releases feeds the rainfall-runoff model.

    python examples/04_snow.py
"""

import numpy as np

from _data import snow_catchment


def main():
    catchment = snow_catchment()
    print(catchment)
    print("Elevation bands from the hypsometric curve: %s m"
          % np.round(catchment._inputs("CemaNeigeGR4J").z_layers).astype(int))

    plain = catchment.split_sample(
        "GR4J", calibration=("1990-01-01", "1999-12-31"),
        validation=("2000-01-01", "2009-12-31"), criterion="KGE")[1]
    fit, with_snow = catchment.split_sample(
        "CemaNeigeGR4J", calibration=("1990-01-01", "1999-12-31"),
        validation=("2000-01-01", "2009-12-31"), criterion="KGE")

    print("\nDoes the snow module earn its two extra parameters?")
    print("  GR4J alone       validation KGE = %.3f" % plain.kge())
    print("  CemaNeige-GR4J   validation KGE = %.3f" % with_snow.kge())

    print("\n" + fit.describe())

    print("\nSnow pack per elevation band [mm water equivalent]")
    print("  %-8s %10s %10s %12s" % ("band", "mean", "max", "days with snow"))
    for i, layer in enumerate(with_snow["CemaNeigeLayers"], start=1):
        pack = layer["SnowPack"]
        print("  band %-3d %10.1f %10.1f %12d"
              % (i, pack.mean(), pack.max(), int((pack > 1).sum())))

    melt = sum(layer["Melt"].sum() for layer in with_snow["CemaNeigeLayers"])
    solid = sum(layer["Psol"].sum() for layer in with_snow["CemaNeigeLayers"])
    n_bands = len(with_snow["CemaNeigeLayers"])
    print("\nOver the validation period, averaged over the catchment:")
    print("  solid precipitation %.0f mm, melt %.0f mm"
          % (solid / n_bands, melt / n_bands))

    # Hysteresis adds two parameters describing how snow cover builds and decays.
    hyst = catchment.calibrate("CemaNeigeGR4J", hysteresis=True,
                               criterion="KGE",
                               period=("1990-01-01", "1999-12-31"))
    print("\nWith linear hysteresis (2 more parameters):")
    print("  validation KGE = %.3f"
          % hyst.evaluate(period=("2000-01-01", "2009-12-31")).kge())


if __name__ == "__main__":
    main()
