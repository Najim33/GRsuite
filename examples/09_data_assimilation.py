"""Assimilate observed discharge into a calibrated model.

Data assimilation corrects the model's internal states step by step, every
time an observation is available: the ensemble Kalman filter (EnKF) or a
particle filter (PF) pulls the production store, the routing store and the
unit hydrograph levels towards the observed discharge. It is the standard
way to put a rainfall-runoff model back on track before a forecast.

    python examples/09_data_assimilation.py
"""

from _data import daily_catchment


def main():
    catchment = daily_catchment()

    fit = catchment.calibrate("GR4J", criterion="KGE",
                              period=("1990-01-01", "1999-12-31"))
    print(fit.describe())

    period = ("2000-01-01", "2005-12-31")

    # Open loop: the calibrated model on its own, no correction.
    check = fit.evaluate(period=period)
    print("\nOpen loop over %s to %s" % period)
    print("   NSE = %.3f" % check.nse())

    # EnKF: states updated at every observed step. The ensemble is kept
    # spread by perturbing the production and routing stores after each
    # update (without it the members would stay identical and the Kalman
    # gain would be zero).
    enkf = catchment.assimilate("GR4J", fit.params, method="EnKF",
                                period=period, nb_mbr=50,
                                state_pert=["Prod", "Rout"], seed=42)
    print("EnKF, ensemble median over the same period")
    print("   NSE = %.3f   (%s)" % (_nse(enkf.qobs, enkf.qsim_median()), enkf))

    # Particle filter on the same setup.
    pf = catchment.assimilate("GR4J", fit.params, method="PF",
                              period=period, nb_mbr=50,
                              state_pert=["Prod", "Rout"], seed=42)
    print("PF, ensemble median over the same period")
    print("   NSE = %.3f   (%s)" % (_nse(pf.qobs, pf.qsim_median()), pf))

    # State trajectories are available member by member.
    prod = enkf.states("Prod")
    print("\nProduction store level [mm], first and last day (EnKF)")
    print("   day 1 :", prod[0].round(1))
    print("   day %i: %s ... %s"
          % (prod.shape[0], prod[-1].min().round(1), prod[-1].max().round(1)))


def _nse(obs, sim):
    ok = obs >= 0
    return 1.0 - float(((sim[ok] - obs[ok]) ** 2).sum()
                       / ((obs[ok] - obs[ok].mean()) ** 2).sum())


if __name__ == "__main__":
    main()
