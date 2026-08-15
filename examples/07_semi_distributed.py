"""Routing upstream sub-catchments down to a gauged outlet.

A lumped model on the downstream sub-catchment, plus a pure translation of the
upstream inflows, gives a simple semi-distributed scheme. The single parameter
is the propagation velocity.

    python examples/07_semi_distributed.py
"""

import numpy as np

import grsuite as gr
from _data import daily_catchment, upstream_discharge

AREAS = [180.0, 240.0, 360.0]     # km2: upstream 1, upstream 2, downstream
LENGTHS = [30.0, 55.0]            # km of river between each inflow and the outlet


def main():
    catchment = daily_catchment()
    period = ("1995-01-01", "1999-12-31")

    # Local runoff of the downstream sub-catchment.
    inputs = catchment._inputs("GR4J")
    options = gr.RunOptions(inputs, "GR4J",
                            ind_period_run=catchment.period_index(period))
    local = gr.run_model_gr4j(inputs, options, [257.238, 1.012, 88.235, 2.208])

    network = gr.InputsModelSD(
        q_upstream=upstream_discharge(catchment, AREAS),
        length_hydro=LENGTHS,
        basin_areas=AREAS)

    print("A pure translation conserves volume and only moves water in time.\n")
    print("%-11s %11s %11s %15s" % ("velocity", "mean Qsim", "peak Qsim",
                                    "travel time"))
    print("-" * 52)
    for speed in [0.5, 1.0, 1.5, 2.5]:
        routed = gr.run_model_lag(network, options, [speed], local["Qsim"],
                                  warmup_q_down=local["WarmUpQsim"])
        lag_days = np.array(LENGTHS) * 1e3 / speed / 86400.0
        print("%7.1f m/s %11.3f %11.3f %14s d"
              % (speed, routed["Qsim"].mean(), routed["Qsim"].max(),
                 np.round(lag_days, 2).tolist()))

    print("\nMean discharge is identical across velocities: routing moves the")
    print("hydrograph without creating or destroying water. The peak does move,")
    print("because a slower velocity spreads the two inflows further apart and")
    print("they stop coinciding with the local flood.")

    routed = gr.run_model_lag(network, options, [1.2], local["Qsim"],
                              warmup_q_down=local["WarmUpQsim"])
    print("\nAt the outlet the catchment is %.0f km2 instead of %.0f km2,"
          % (sum(AREAS), AREAS[-1]))
    print("so specific discharge is the weighted mix of local runoff and inflows:")
    print("  local sub-catchment only  %.3f mm/d" % local["Qsim"].mean())
    print("  at the outlet             %.3f mm/d" % routed["Qsim"].mean())
    print("  total volume              %.0f x 10^6 m3"
          % (routed["Qsim_m3"].sum() / 1e6))


if __name__ == "__main__":
    main()
