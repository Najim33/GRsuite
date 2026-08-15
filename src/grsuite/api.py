"""High-level interface: get from a CSV of rainfall and streamflow to a
calibrated model in a few lines.

The rest of the package mirrors airGR one function at a time, which is what
you want when porting an R workflow. This module is the shortcut for
everything else::

    import grsuite as gr

    catchment = gr.Catchment(dates, precip=P, pot_evap=E, obs_discharge=Q)
    fit = catchment.calibrate("GR4J", period=("2000-01-01", "2009-12-31"))
    check = fit.evaluate(period=("2010-01-01", "2019-12-31"))

    print(fit.params, check.nse(), check.kge())

Everything here is a thin, explicit wrapper: no hidden defaults beyond the
ones airGR itself applies.
"""

import numpy as np

from .calib import calibration_michel
from .core import CalibOptions, InputsModel, RunOptions
from .crit import CRIT_FUNCS, InputsCrit, InputsCritCompo, error_crit
from .models import MODEL_FUNCS, data_alti_extrapolation_valery, run_model

__all__ = ["Catchment", "Simulation", "CalibratedModel", "list_models"]

SNOW_MODELS = {"CemaNeigeGR4J", "CemaNeigeGR5J", "CemaNeigeGR6J",
               "CemaNeigeGR4H", "CemaNeigeGR5H"}

#: Parameter names, in the order each model expects them.
PARAM_NAMES = {
    "GR1A": ["X1"],
    "GR2M": ["X1", "X2"],
    "GR4J": ["X1", "X2", "X3", "X4"],
    "GR5J": ["X1", "X2", "X3", "X4", "X5"],
    "GR6J": ["X1", "X2", "X3", "X4", "X5", "X6"],
    "GR4H": ["X1", "X2", "X3", "X4"],
    "GR5H": ["X1", "X2", "X3", "X4", "X5"],
}
SNOW_PARAM_NAMES = ["CN1", "CN2"]
HYST_PARAM_NAMES = ["CN3", "CN4"]

PARAM_UNITS = {
    "X1": "production store capacity [mm]",
    "X2": "groundwater exchange coefficient [mm/step]",
    "X3": "routing store capacity [mm]",
    "X4": "unit hydrograph time constant [step]",
    "X5": "exchange threshold [-]",
    "X6": "exponential store depletion coefficient [mm]",
    "CN1": "snow pack thermal state weighting [-]",
    "CN2": "degree-day melt factor [mm/degC/step]",
    "CN3": "hysteresis accumulation threshold [mm]",
    "CN4": "hysteresis melt coefficient [-]",
}


def list_models():
    """Return the model names accepted by :class:`Catchment`."""
    return sorted(MODEL_FUNCS)


def param_names(model, hysteresis=False):
    """Parameter names for a model, snow parameters included."""
    base = model.replace("CemaNeige", "") if model in SNOW_MODELS else model
    names = list(PARAM_NAMES[base])
    if model in SNOW_MODELS:
        names += SNOW_PARAM_NAMES
        if hysteresis:
            names += HYST_PARAM_NAMES
    return names


#: Criteria airGR maximises (their ``multiplier`` is -1).
MAXIMISED_CRITERIA = ("NSE", "KGE", "KGE2")


def _to_date(value, dtype):
    return np.asarray(np.datetime64(value)).astype(dtype)


def _oriented(inputs_crit, value):
    """Give a composite criterion the orientation of a simple one.

    ``error_crit`` follows airGR and returns a composite already multiplied by
    its sign convention, so a composite of two NSE would come out negative.
    The high-level API flips it back and leaves simple criteria untouched.
    """
    if not isinstance(inputs_crit, InputsCritCompo):
        return value
    if all(item.fun_crit in MAXIMISED_CRITERIA for item in inputs_crit.items):
        return -value
    return value


class Catchment:
    """A catchment: forcing series, observed discharge, optional elevation.

    Parameters
    ----------
    dates : array of numpy.datetime64
        Regular, sorted time steps. The time unit drives the model family:
        daily arrays feed GR4J/GR5J/GR6J, hourly arrays feed GR4H/GR5H.
    precip : array
        Total precipitation [mm per time step].
    pot_evap : array
        Potential evapotranspiration [mm per time step]. If you only have
        temperature, build it first with :func:`grsuite.pe_oudin`.
    obs_discharge : array, optional
        Observed discharge [mm per time step]. Required to calibrate or to
        score a simulation; missing values are allowed and simply skipped.
    temperature : array, optional
        Mean air temperature [degC]. Required by the snow models.
    hypsometry : array of 101 floats, optional
        Elevation quantiles of the catchment [m], from minimum to maximum.
        Required by the snow models.
    elevation : float, optional
        Elevation the forcing series refer to [m]. Defaults to the median
        of ``hypsometry``.
    n_layers : int
        Number of elevation bands used by CemaNeige.
    name : str, optional
        Free-form label, used in reports and plots.
    """

    def __init__(self, dates, precip, pot_evap, obs_discharge=None,
                 temperature=None, hypsometry=None, elevation=None,
                 n_layers=5, name=None):
        self.dates = np.asarray(dates)
        self.precip = np.asarray(precip, dtype=float)
        self.pot_evap = np.asarray(pot_evap, dtype=float)
        self.obs_discharge = (None if obs_discharge is None
                              else np.asarray(obs_discharge, dtype=float))
        self.temperature = (None if temperature is None
                            else np.asarray(temperature, dtype=float))
        self.hypsometry = (None if hypsometry is None
                           else np.asarray(hypsometry, dtype=float))
        self.n_layers = int(n_layers)
        self.name = name

        if self.hypsometry is not None:
            self.elevation = (float(np.median(self.hypsometry))
                              if elevation is None else float(elevation))
        else:
            self.elevation = elevation

        n = self.precip.shape[0]
        for label, series in (("pot_evap", self.pot_evap),
                              ("obs_discharge", self.obs_discharge),
                              ("temperature", self.temperature)):
            if series is not None and series.shape[0] != n:
                raise ValueError(
                    "'%s' has %i values but 'precip' has %i"
                    % (label, series.shape[0], n))
        if self.dates.shape[0] != n:
            raise ValueError("'dates' and 'precip' must have the same length")

        self._plain = None
        self._snow = None

    def __len__(self):
        return self.precip.shape[0]

    def __repr__(self):
        return ("Catchment(%s%i steps, %s to %s)"
                % ("%s, " % self.name if self.name else "", len(self),
                   self.dates[0], self.dates[-1]))

    # -- internals ---------------------------------------------------------

    def _inputs(self, model):
        if model in SNOW_MODELS:
            if self._snow is None:
                if self.temperature is None or self.hypsometry is None:
                    raise ValueError(
                        "%s needs 'temperature' and 'hypsometry'; pass them to "
                        "Catchment(), or use a model without snow (%s)"
                        % (model, ", ".join(sorted(PARAM_NAMES))))
                alti = data_alti_extrapolation_valery(
                    self.dates, self.precip, self.temperature,
                    z_inputs=self.elevation, hypso_data=self.hypsometry,
                    n_layers=self.n_layers)
                self._snow = InputsModel(
                    self.dates, self.precip, self.pot_evap,
                    temp_mean=self.temperature,
                    layer_precip=alti["LayerPrecip"],
                    layer_temp_mean=alti["LayerTempMean"],
                    layer_frac_solid_precip=alti["LayerFracSolidPrecip"],
                    z_layers=alti["ZLayers"])
            return self._snow
        if self._plain is None:
            self._plain = InputsModel(self.dates, self.precip, self.pot_evap,
                                      temp_mean=self.temperature)
        return self._plain

    def period_index(self, period=None):
        """Turn a ``(start, end)`` pair of dates into an index array.

        ``None`` means the whole series minus a one-year warm-up, which is
        the sensible default for a first look at a catchment.
        """
        if period is None:
            step = self.dates[1] - self.dates[0]
            per_year = int(np.timedelta64(365, "D") / step)
            start = min(per_year, len(self) // 2)
            return np.arange(start, len(self))
        start, end = period
        dtype = self.dates.dtype
        i0 = int(np.searchsorted(self.dates, _to_date(start, dtype)))
        i1 = int(np.searchsorted(self.dates, _to_date(end, dtype), side="right"))
        if i1 <= i0:
            raise ValueError("empty period: %s to %s" % (start, end))
        return np.arange(i0, i1)

    def _options(self, model, period, hysteresis, imax, warmup):
        index = self.period_index(period)
        warmup_index = None
        if warmup is not None:
            warmup_index = self.period_index(warmup)
        return RunOptions(self._inputs(model), model, ind_period_run=index,
                          ind_period_warmup=warmup_index, is_hyst=hysteresis,
                          imax=imax)

    # -- public API --------------------------------------------------------

    def simulate(self, model, params, period=None, hysteresis=False,
                 imax=None, warmup=None):
        """Run a model with a given parameter set.

        Parameters
        ----------
        model : str
            One of :func:`list_models`.
        params : sequence of float
            Parameters in the model's own order, see :func:`param_names`.
        period : (str, str), optional
            Simulation period as ``("YYYY-MM-DD", "YYYY-MM-DD")``.
        hysteresis : bool
            Use CemaNeige with linear hysteresis (two extra parameters).
        imax : float, optional
            Interception store capacity [mm], GR5H only.
        warmup : (str, str), optional
            Explicit warm-up period. Left out, airGR's own default applies:
            the year preceding the simulation, if available.

        Returns
        -------
        Simulation
        """
        if model not in MODEL_FUNCS:
            raise ValueError("unknown model %r; available: %s"
                             % (model, ", ".join(list_models())))
        options = self._options(model, period, hysteresis, imax, warmup)
        outputs = run_model(self._inputs(model), options, params, model=model)
        return Simulation(self, model, np.asarray(params, dtype=float),
                          options, outputs, hysteresis=hysteresis)

    def calibrate(self, model, criterion="KGE", transfo="", period=None,
                  hysteresis=False, imax=None, warmup=None, epsilon=None,
                  fixed_params=None, verbose=False):
        """Calibrate a model with Michel's algorithm, exactly as airGR does.

        Parameters
        ----------
        model : str
        criterion : str or sequence
            ``"NSE"``, ``"KGE"``, ``"KGE2"``, ``"RMSE"``, or a sequence of
            ``(criterion, transfo, weight)`` tuples for a composite objective.
        transfo : str
            Discharge transformation applied before scoring: ``""``,
            ``"sqrt"``, ``"log"``, ``"inv"``, ``"sort"``, ``"boxcox"`` or a
            power such as ``"^0.5"``. ``"log"`` and ``"inv"`` default to an
            epsilon of one hundredth of the mean observed discharge.
        period : (str, str), optional
            Calibration period.
        fixed_params : sequence, optional
            Values to hold fixed; use ``None`` or ``nan`` for the parameters
            that should still be optimised.

        Returns
        -------
        CalibratedModel
        """
        if self.obs_discharge is None:
            raise ValueError("calibration needs 'obs_discharge'")
        if model not in MODEL_FUNCS:
            raise ValueError("unknown model %r; available: %s"
                             % (model, ", ".join(list_models())))

        options = self._options(model, period, hysteresis, imax, warmup)
        obs = self.obs_discharge[options.ind_period_run]
        inputs_crit = self._build_criterion(criterion, transfo, obs, epsilon)

        calib_model = model + ("Hyst" if hysteresis and model in SNOW_MODELS
                               else "")
        fixed = None
        if fixed_params is not None:
            fixed = np.array([np.nan if v is None else float(v)
                              for v in fixed_params])
        calib_options = CalibOptions(calib_model, fixed_param=fixed,
                                     is_int_store=imax is not None)

        result = calibration_michel(self._inputs(model), options, inputs_crit,
                                    calib_options, model=model, verbose=verbose)
        return CalibratedModel(self, model, result, options,
                               hysteresis=hysteresis, imax=imax,
                               criterion=criterion, transfo=transfo,
                               inputs_crit=inputs_crit)

    def _build_criterion(self, criterion, transfo, obs, epsilon):
        default_eps = None
        if epsilon is None:
            finite = obs[np.isfinite(obs)]
            if finite.size:
                default_eps = 0.01 * float(np.mean(finite))

        def one(name, tf, weight=None):
            name = name.upper()
            if name not in CRIT_FUNCS:
                raise ValueError("unknown criterion %r; available: %s"
                                 % (name, ", ".join(sorted(CRIT_FUNCS))))
            eps = epsilon
            if eps is None and tf in ("log", "inv"):
                eps = default_eps
            return InputsCrit(name, obs=obs, transfo=tf, epsilon=eps,
                              weights=weight)

        if isinstance(criterion, str):
            return one(criterion, transfo)
        parts = []
        for item in criterion:
            name, tf, weight = (list(item) + ["", 1.0])[:3]
            parts.append(one(name, tf, weight))
        return InputsCritCompo(parts)

    def split_sample(self, model, calibration, validation, **kwargs):
        """Klemes split-sample test: calibrate on one period, score on another.

        Returns
        -------
        (CalibratedModel, Simulation)
            The fit and its evaluation over the validation period.
        """
        fit = self.calibrate(model, period=calibration, **kwargs)
        return fit, fit.evaluate(period=validation)


class Simulation:
    """The result of one model run, with the usual scores attached."""

    def __init__(self, catchment, model, params, options, outputs,
                 hysteresis=False):
        self.catchment = catchment
        self.model = model
        self.params = params
        self.options = options
        self.outputs = outputs
        self.hysteresis = hysteresis

    # -- data --------------------------------------------------------------

    @property
    def dates(self):
        return self.outputs["DatesR"]

    @property
    def qsim(self):
        """Simulated discharge [mm per time step]."""
        return self.outputs["Qsim"]

    @property
    def qobs(self):
        """Observed discharge over the simulated period, or None."""
        if self.catchment.obs_discharge is None:
            return None
        return self.catchment.obs_discharge[self.options.ind_period_run]

    def __getitem__(self, key):
        return self.outputs[key]

    def variables(self):
        """Names of every internal variable the model exposes.

        Dates, parameters and the final state vector are metadata, not
        variables, so they are left out.
        """
        skip = ("DatesR", "Param", "StateEnd", "WarmUpQsim", "CemaNeigeLayers")
        return [k for k, v in self.outputs.items()
                if k not in skip and isinstance(v, np.ndarray)
                and v.shape == self.qsim.shape
                and np.issubdtype(v.dtype, np.number)]

    def to_dataframe(self, variables=None):
        """Return the simulation as a pandas DataFrame indexed by date."""
        import pandas as pd
        names = variables or self.variables()
        frame = pd.DataFrame({n: self.outputs[n] for n in names},
                             index=pd.DatetimeIndex(self.dates, name="date"))
        if self.qobs is not None:
            frame.insert(0, "Qobs", self.qobs)
        return frame

    # -- scores ------------------------------------------------------------

    def score(self, criterion="KGE", transfo="", epsilon=None):
        """Score the simulation against the observed discharge.

        A composite criterion comes back with the same orientation as a simple
        one: higher is better for NSE, KGE and KGE2, lower is better for RMSE.
        (:func:`grsuite.error_crit` keeps airGR's own raw sign convention.)
        """
        if self.qobs is None:
            raise ValueError("scoring needs 'obs_discharge' on the catchment")
        crit = self.catchment._build_criterion(criterion, transfo, self.qobs,
                                               epsilon)
        return _oriented(crit, error_crit(crit, self.outputs).crit_value)

    def nse(self, transfo=""):
        """Nash-Sutcliffe efficiency (1 is perfect)."""
        return self.score("NSE", transfo)

    def kge(self, transfo=""):
        """Kling-Gupta efficiency (1 is perfect)."""
        return self.score("KGE", transfo)

    def rmse(self, transfo=""):
        """Root mean square error, in mm per time step (0 is perfect)."""
        return self.score("RMSE", transfo)

    def bias(self):
        """Relative volume bias, simulated over observed minus one."""
        obs, sim = self.qobs, self.qsim
        ok = np.isfinite(obs) & np.isfinite(sim)
        return float(np.mean(sim[ok]) / np.mean(obs[ok]) - 1.0)

    def summary(self):
        """A small dictionary of the scores people usually report."""
        out = {"model": self.model, "n_steps": int(self.qsim.shape[0]),
               "mean_qsim": float(np.nanmean(self.qsim))}
        if self.qobs is not None:
            out.update({
                "NSE": self.nse(), "KGE": self.kge(),
                "NSE_log": self.nse("log"), "RMSE": self.rmse(),
                "bias": self.bias(),
                "mean_qobs": float(np.nanmean(self.qobs)),
            })
        return out

    def __repr__(self):
        if self.qobs is None:
            return "Simulation(%s, %i steps)" % (self.model, len(self.qsim))
        return ("Simulation(%s, %i steps, NSE=%.3f, KGE=%.3f)"
                % (self.model, len(self.qsim), self.nse(), self.kge()))

    # -- plotting ----------------------------------------------------------

    def plot(self, ax=None, log=False, title=None):
        """Plot observed and simulated hydrographs (needs matplotlib)."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(11, 4))
        if self.qobs is not None:
            ax.plot(self.dates, self.qobs, lw=1.1, color="#3c4c54",
                    label="observed")
        ax.plot(self.dates, self.qsim, lw=1.1, color="#0a7ea4",
                label="%s simulated" % self.model)
        if log:
            ax.set_yscale("log")
        ax.set_ylabel("discharge [mm/step]")
        ax.legend(frameon=False, ncol=2)
        ax.set_title(title or (self.catchment.name or "") or None)
        ax.spines[["top", "right"]].set_visible(False)
        return ax


class CalibratedModel:
    """A calibrated parameter set, plus what it took to get there."""

    def __init__(self, catchment, model, result, options, hysteresis=False,
                 imax=None, criterion="KGE", transfo="", inputs_crit=None):
        self.catchment = catchment
        self.model = model
        self.result = result
        self.options = options
        self.hysteresis = hysteresis
        self.imax = imax
        self.criterion = criterion
        self.transfo = transfo
        self._inputs_crit = inputs_crit

    @property
    def params(self):
        """Calibrated parameters, in the model's own order."""
        return self.result.param_final_r

    @property
    def score(self):
        """Value of the objective function at the optimum.

        Oriented as the criterion itself is: higher is better for NSE, KGE
        and KGE2, lower is better for RMSE.
        """
        return _oriented(self._inputs_crit, self.result.crit_final)

    def named_params(self):
        """Parameters as a dictionary, keyed by their conventional names."""
        names = param_names(self.model, self.hysteresis)
        return dict(zip(names, (float(v) for v in self.params)))

    def describe(self):
        """A readable, unit-annotated parameter listing."""
        if isinstance(self.criterion, str):
            label = "%s%s" % (self.criterion,
                              "[%s]" % self.transfo if self.transfo else "")
        else:
            label = "composite(%s)" % ", ".join(
                "%s%s" % (str(item[0]).upper(),
                          "[%s]" % item[1] if len(item) > 1 and item[1] else "")
                for item in self.criterion)
        lines = ["%s calibrated on %s = %.4f" % (self.model, label, self.score)]
        for name, value in self.named_params().items():
            lines.append("  %-4s %12.4f   %s"
                         % (name, value, PARAM_UNITS.get(name, "")))
        lines.append("  %i iterations, %i model runs"
                     % (self.result.n_iter, self.result.n_runs))
        return "\n".join(lines)

    def simulate(self, period=None, **kwargs):
        """Re-run the calibrated model, over any period."""
        return self.catchment.simulate(
            self.model, self.params, period=period,
            hysteresis=self.hysteresis, imax=self.imax, **kwargs)

    def evaluate(self, period=None, **kwargs):
        """Alias of :meth:`simulate`, for use on a validation period."""
        return self.simulate(period=period, **kwargs)

    def __repr__(self):
        name = self.criterion if isinstance(self.criterion, str) else "composite"
        return ("CalibratedModel(%s, %s=%.4f, params=%s)"
                % (self.model, name, self.score, np.round(self.params, 3)))
