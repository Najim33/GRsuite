"""Run GR4J with a known parameter set and look inside the model.

    python examples/01_simulate.py
"""

from _data import daily_catchment


def main():
    catchment = daily_catchment()
    print(catchment)

    # X1 production store [mm], X2 exchange [mm/d], X3 routing store [mm],
    # X4 unit hydrograph time constant [d]
    simulation = catchment.simulate("GR4J", [257.238, 1.012, 88.235, 2.208],
                                    period=("1995-01-01", "1999-12-31"))

    print("\nScores against observed discharge")
    for name, value in simulation.summary().items():
        print("  %-10s %s" % (name, round(value, 4)
                              if isinstance(value, float) else value))

    print("\nEvery internal variable the model exposes:")
    print(" ", ", ".join(simulation.variables()))

    print("\nFirst days of the water balance [mm/d]:")
    frame = simulation.to_dataframe()
    print(frame[["Precip", "PotEvap", "AE", "Prod", "Perc", "Rout",
                 "Exch", "Qsim", "Qobs"]].head(8).round(3).to_string())

    # The store levels are what a modeller usually wants to plot next.
    print("\nProduction store: %.1f to %.1f mm, routing store: %.1f to %.1f mm"
          % (frame["Prod"].min(), frame["Prod"].max(),
             frame["Rout"].min(), frame["Rout"].max()))


if __name__ == "__main__":
    main()
