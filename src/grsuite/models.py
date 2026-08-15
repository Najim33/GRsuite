"""Modeles pluie-debit de GRsuite.

Traduction de airGR/R/RunModel_*.R et DataAltiExtrapolation_Valery.R
(v1.7.9). Chaque fonction `run_model_*` renvoie un dictionnaire dont les
cles reprennent exactement les noms de sortie d'airGR.

References
----------
Perrin, C., Michel, C., Andreassian, V. (2003). Improvement of a
    parsimonious model for streamflow simulation. Journal of Hydrology
    279, 275-289. doi:10.1016/S0022-1694(03)00225-7               [GR4J]
Le Moine, N. (2008). Le bassin versant de surface vu par le souterrain.
    These de doctorat, UPMC / Cemagref, Antony.                   [GR5J]
Pushpalatha, R., Perrin, C., Le Moine, N., Mathevet, T., Andreassian, V.
    (2011). A downward structural sensitivity analysis of hydrological
    models to improve low-flow simulation. Journal of Hydrology 411,
    66-76. doi:10.1016/j.jhydrol.2011.09.034                [GR6J, GR5J]
Mathevet, T. (2005). Quels modeles pluie-debit globaux au pas de temps
    horaire ? These de doctorat, ENGREF / Cemagref, Antony.       [GR4H]
Ficchi, A. (2017). An adaptive hydrological model for multiple time
    steps. These de doctorat, UPMC / Irstea, Antony.              [GR5H]
Ficchi, A., Perrin, C., Andreassian, V. (2019). Hydrological modelling
    at multiple sub-daily time steps: model improvement via
    flux-matching. Journal of Hydrology 575, 1308-1327.
    doi:10.1016/j.jhydrol.2019.05.084          [GR5H, interception, Imax]
Mouelhi, S., Michel, C., Perrin, C., Andreassian, V. (2006a). Stepwise
    development of a two-parameter monthly water balance model. Journal
    of Hydrology 318, 200-214. doi:10.1016/j.jhydrol.2005.06.014  [GR2M]
Mouelhi, S., Michel, C., Perrin, C., Andreassian, V. (2006b). Linking
    stream flow to rainfall at the annual time step: the Manabe bucket
    model revisited. Journal of Hydrology 328, 283-296.
    doi:10.1016/j.jhydrol.2005.12.022                             [GR1A]
Valery, A. (2010). Modelisation precipitations-debit sous influence
    nivale. These de doctorat, AgroParisTech / Cemagref, Antony.
                                       [CemaNeige, bandes d'altitude]
Riboust, P., Thirel, G., Le Moine, N., Ribstein, P. (2019). Revisiting a
    simple degree-day model for integrating satellite data:
    implementation of SWE-SCA hystereses. Journal of Hydrology and
    Hydromechanics 67, 70-81. doi:10.2478/johh-2018-0004     [hysteresis]
Oudin, L., Hervieu, F., Michel, C., Perrin, C., Andreassian, V., Anctil,
    F., Loumagne, C. (2005). Which potential evapotranspiration input
    for a lumped rainfall-runoff model? Part 2. Journal of Hydrology
    303, 290-306. doi:10.1016/j.jhydrol.2004.08.026                [ETP]
"""

import os

import numpy as np

from . import _kernels as K
from . import _kernels_h as KH
from .core import CEMANEIGE_OUTPUTS, MODEL_OUTPUTS, N_STATES_BASE, N_STATES_BASE_H

__all__ = [
    "run_model_gr4j", "run_model_gr5j", "run_model_gr6j",
    "run_model_gr2m", "run_model_gr1a", "run_model_gr4h", "run_model_gr5h",
    "run_model_cemaneige", "run_model_cemaneige_gr4j",
    "run_model_cemaneige_gr5j",
    "run_model_cemaneige_gr6j", "run_model_cemaneige_gr4h",
    "run_model_cemaneige_gr5h",
    "run_model", "pe_oudin", "data_alti_extrapolation_valery", "imax_estimate",
    "MODEL_FUNCS",
]

_GRADT_PATH = os.path.join(os.path.dirname(__file__), "data",
                           "gradT_valery2010.csv")
_GRADT = None


def _grad_t():
    """Table journaliere des gradients de temperature (Valery, 2010)."""
    global _GRADT
    if _GRADT is None:
        raw = np.genfromtxt(_GRADT_PATH, delimiter=",", names=True)
        _GRADT = {
            "day": raw["day"].astype(int),
            "month": raw["month"].astype(int),
            "grad_Tmean": raw["grad_Tmean"],
            "grad_Tmin": raw["grad_Tmin"],
            "grad_Tmax": raw["grad_Tmax"],
        }
    return _GRADT


# ---------------------------------------------------------------------------
# Assemblage des sorties
# ---------------------------------------------------------------------------


def _assemble(model, outputs, state_end, run_options, inputs_model,
              param, cemaneige_layers=None):
    n_warmup = run_options.ind_period_warmup.shape[0]
    ind2 = slice(n_warmup, outputs.shape[0])
    names = MODEL_OUTPUTS[model]

    res = {"DatesR": inputs_model.dates[run_options.ind_period_run]}
    for i, name in enumerate(names):
        col = outputs[ind2, i].copy()
        col[col <= -99e8] = np.nan
        res[name] = col
    if n_warmup > 0:
        qcol = names.index("Qsim")
        res["WarmUpQsim"] = outputs[:n_warmup, qcol].copy()
    res["StateEnd"] = state_end
    res["Param"] = np.asarray(param, dtype=float)
    if cemaneige_layers is not None:
        res["CemaNeigeLayers"] = cemaneige_layers
    return res


def _ind_period1(run_options):
    return np.concatenate([run_options.ind_period_warmup,
                           run_options.ind_period_run])


def _check_param(param, n_expected, model):
    p = np.asarray(param, dtype=float).copy()
    if p.shape[0] != n_expected:
        raise ValueError("le modele %s requiert %i parametres (recu %i)"
                         % (model, n_expected, p.shape[0]))
    return p


def _apply_thresholds_gr4j_like(p):
    """Seuils de securite appliques par airGR sur X1, X3 et X4."""
    if p[0] < 1e-2:
        p[0] = 1e-2
    if p[2] < 1e-2:
        p[2] = 1e-2
    if p[3] < 0.5:
        p[3] = 0.5
    return p


# ---------------------------------------------------------------------------
# Modeles GR sans neige
# ---------------------------------------------------------------------------


def run_model_gr4j(inputs_model, run_options, param):
    """GR4J (Perrin et al., 2003), pas de temps journalier."""
    p = _apply_thresholds_gr4j_like(_check_param(param, 4, "GR4J"))
    ind1 = _ind_period1(run_options)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = lv[1] * p[2]

    outputs, state_end = K.run_gr4j(inputs_model.precip[ind1],
                                    inputs_model.pot_evap[ind1], p,
                                    ini[:N_STATES_BASE])
    return _assemble("GR4J", outputs, state_end, run_options, inputs_model, p)


def run_model_gr5j(inputs_model, run_options, param):
    """GR5J (Le Moine, 2008), pas de temps journalier."""
    p = _apply_thresholds_gr4j_like(_check_param(param, 5, "GR5J"))
    ind1 = _ind_period1(run_options)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = lv[1] * p[2]

    outputs, state_end = K.run_gr5j(inputs_model.precip[ind1],
                                    inputs_model.pot_evap[ind1], p,
                                    ini[:N_STATES_BASE])
    return _assemble("GR5J", outputs, state_end, run_options, inputs_model, p)


def run_model_gr6j(inputs_model, run_options, param):
    """GR6J (Pushpalatha et al., 2011), pas de temps journalier."""
    p = _check_param(param, 6, "GR6J")
    if p[0] < 1e-2:
        p[0] = 1e-2
    if p[2] < 1e-2:
        p[2] = 1e-2
    if p[3] < 0.5:
        p[3] = 0.5
    if p[5] < 1e-2:
        p[5] = 1e-2
    ind1 = _ind_period1(run_options)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = lv[1] * p[2]
        ini[2] = lv[2]

    outputs, state_end = K.run_gr6j(inputs_model.precip[ind1],
                                    inputs_model.pot_evap[ind1], p,
                                    ini[:N_STATES_BASE])
    return _assemble("GR6J", outputs, state_end, run_options, inputs_model, p)


def run_model_gr2m(inputs_model, run_options, param):
    """GR2M (Mouelhi et al., 2006a), pas de temps mensuel."""
    p = _check_param(param, 2, "GR2M")
    if p[0] < 1e-2:
        p[0] = 1e-2
    if p[1] < 1e-2:
        p[1] = 1e-2
    ind1 = _ind_period1(run_options)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = 60.0 * lv[1] * p[1]

    outputs, state_end = K.run_gr2m(inputs_model.precip[ind1],
                                    inputs_model.pot_evap[ind1], p,
                                    ini[:N_STATES_BASE])
    return _assemble("GR2M", outputs, state_end, run_options, inputs_model, p)


def run_model_gr1a(inputs_model, run_options, param):
    """GR1A (Mouelhi et al., 2006b), pas de temps annuel."""
    p = _check_param(param, 1, "GR1A")
    ind1 = _ind_period1(run_options)
    outputs = K.run_gr1a(inputs_model.precip[ind1],
                         inputs_model.pot_evap[ind1], p)
    state_end = np.zeros(N_STATES_BASE)
    return _assemble("GR1A", outputs, state_end, run_options, inputs_model, p)


# ---------------------------------------------------------------------------
# Modeles horaires
# ---------------------------------------------------------------------------


def run_model_gr4h(inputs_model, run_options, param):
    """GR4H (Mathevet, 2005), pas de temps horaire."""
    p = _apply_thresholds_gr4j_like(_check_param(param, 4, "GR4H"))
    ind1 = _ind_period1(run_options)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = lv[1] * p[2]

    outputs, state_end = KH.run_gr4h(inputs_model.precip[ind1],
                                     inputs_model.pot_evap[ind1], p,
                                     ini[:N_STATES_BASE_H])
    return _assemble("GR4H", outputs, state_end, run_options, inputs_model, p)


def run_model_gr5h(inputs_model, run_options, param):
    """GR5H (Ficchi, 2017 ; Ficchi et al., 2019), pas de temps horaire.

    Le reservoir d'interception est active en fournissant `imax` dans
    les options de simulation.
    """
    p = _apply_thresholds_gr4j_like(_check_param(param, 5, "GR5H"))
    ind1 = _ind_period1(run_options)
    imax = run_options.imax
    imax_val = -1.0 if imax is None else float(imax)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * p[0]
        ini[1] = lv[1] * p[2]
        if imax is not None and lv[3] is not None:
            ini[3] = lv[3] * imax_val

    outputs, state_end = KH.run_gr5h(inputs_model.precip[ind1],
                                     inputs_model.pot_evap[ind1], p,
                                     ini[:N_STATES_BASE_H], imax_val)
    return _assemble("GR5H", outputs, state_end, run_options, inputs_model, p)


def imax_estimate(inputs_model, ind_period_run,
                  tested_values=None):
    """Estime la capacite du reservoir d'interception (airGR/Imax.R).

    Retient la valeur qui egalise l'evaporation d'interception cumulee au
    pas horaire et la somme journaliere de min(P, ETP). Methode de
    Ficchi et al. (2019).
    """
    if tested_values is None:
        tested_values = np.arange(0.1, 3.0 + 1e-9, 0.1)
    ind = np.asarray(ind_period_run, dtype=np.int64)
    P = inputs_model.precip[ind]
    E = inputs_model.pot_evap[ind]

    dates = inputs_model.dates[ind].astype("datetime64[D]")
    _, inv = np.unique(dates, return_inverse=True)
    n_days = inv.max() + 1
    p_day = np.bincount(inv, weights=P, minlength=n_days)
    e_day = np.bincount(inv, weights=E, minlength=n_days)
    cum_daily = float(np.sum(np.minimum(p_day, e_day)))

    best_val, best_diff = None, np.inf
    for imax in tested_values:
        c0 = 0.0
        cum_hourly = 0.0
        for i in range(P.shape[0]):
            ec = min(E[i], P[i] + c0)
            pth = max(0.0, P[i] - (imax - c0) - ec)
            c0 = c0 + P[i] - ec - pth
            cum_hourly += ec
        diff = abs(cum_hourly - cum_daily)
        if diff < best_diff:
            best_diff, best_val = diff, imax
    return best_val


# ---------------------------------------------------------------------------
# Extrapolation altitudinale (Valery, 2010)
# ---------------------------------------------------------------------------


def data_alti_extrapolation_valery(dates, precip, temp_mean, z_inputs,
                                   hypso_data, n_layers, precip_scale=True,
                                   temp_min=None, temp_max=None):
    """Repartit P et T sur `n_layers` couches d'altitude.

    Reproduit DataAltiExtrapolation_Valery : gradient de precipitation
    de 0.00041 m-1 plafonne a 4000 m, gradients de temperature
    journaliers de Valery (2010), fraction solide USACE ou Hydrotel.
    """
    grad_p = 0.00041
    hypso_data = np.asarray(hypso_data, dtype=float)
    precip = np.asarray(precip, dtype=float)
    temp_mean = np.asarray(temp_mean, dtype=float)

    z_layers = np.full(n_layers, np.nan)
    if not np.all(np.isnan(hypso_data)):
        nmoy = 100 // n_layers
        nreste = 100 % n_layers
        ncont = 0
        for i in range(n_layers):
            if nreste > 0:
                nn = nmoy + 1
                nreste -= 1
            else:
                nn = nmoy
            if nn == 1:
                z_layers[i] = hypso_data[ncont]
            elif nn == 2:
                z_layers[i] = 0.5 * (hypso_data[ncont] + hypso_data[ncont + 1])
            else:
                z_layers[i] = hypso_data[ncont + nn // 2]
            ncont += nn

    single_layer = (n_layers == 1 and z_inputs == hypso_data[50])

    if single_layer:
        layer_precip = [precip.astype(float).copy()]
    else:
        z_threshold = 4000.0
        cols = []
        for i in range(n_layers):
            if z_layers[i] <= z_threshold:
                cols.append(precip * np.exp(grad_p * (z_layers[i] - z_inputs)))
            elif z_inputs <= z_threshold:
                cols.append(precip * np.exp(grad_p * (z_threshold - z_inputs)))
            else:
                cols.append(precip.copy())
        mat = np.column_stack(cols)
        if precip_scale:
            row_means = mat.mean(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                mat = mat / row_means[:, None] * precip[:, None]
            mat[np.isnan(mat)] = 0.0
        layer_precip = [mat[:, i].copy() for i in range(n_layers)]

    layer_temp_mean = []
    layer_temp_min = []
    layer_temp_max = []
    if single_layer:
        layer_temp_mean.append(temp_mean.astype(float).copy())
        if temp_min is not None and temp_max is not None:
            layer_temp_min.append(np.asarray(temp_min, float).copy())
            layer_temp_max.append(np.asarray(temp_max, float).copy())
    else:
        g = _grad_t()
        key_table = g["day"] * 100 + g["month"]
        d = dates.astype("datetime64[D]")
        days = (d - d.astype("datetime64[M]")).astype(int) + 1
        months = d.astype("datetime64[M]").astype(int) % 12 + 1
        keys = days * 100 + months
        idx = np.searchsorted(np.sort(key_table), keys)
        order = np.argsort(key_table)
        idx = order[np.clip(idx, 0, len(order) - 1)]
        gtm = g["grad_Tmean"][idx]
        gtn = g["grad_Tmin"][idx]
        gtx = g["grad_Tmax"][idx]
        for i in range(n_layers):
            dz = (z_inputs - z_layers[i])
            layer_temp_mean.append(temp_mean + dz * np.abs(gtm) / 100.0)
            if temp_min is not None and temp_max is not None:
                layer_temp_min.append(np.asarray(temp_min, float)
                                      + dz * np.abs(gtn) / 100.0)
                layer_temp_max.append(np.asarray(temp_max, float)
                                      + dz * np.abs(gtx) / 100.0)

    layer_frac_solid = []
    option = "USACE"
    if not np.isnan(z_inputs):
        if z_inputs < 1500.0 and temp_min is not None and temp_max is not None:
            option = "Hydrotel"
    for i in range(n_layers):
        if option == "Hydrotel":
            tmin = layer_temp_min[i]
            tmax = layer_temp_max[i]
            with np.errstate(invalid="ignore", divide="ignore"):
                frac = 1.0 - tmax / (tmax - tmin)
            frac[tmin >= 0.0] = 0.0
            frac[tmax <= 0.0] = 1.0
        else:
            tmean = layer_temp_mean[i]
            frac = 1.0 - (tmean - (-1.0)) / (3.0 - (-1.0))
            frac[tmean > 3.0] = 0.0
            frac[tmean < -1.0] = 1.0
        layer_frac_solid.append(frac)

    return {
        "LayerPrecip": layer_precip,
        "LayerTempMean": layer_temp_mean,
        "LayerTempMin": layer_temp_min,
        "LayerTempMax": layer_temp_max,
        "LayerFracSolidPrecip": layer_frac_solid,
        "ZLayers": z_layers,
    }


def mean_an_solid_precip(inputs_model, factor=365.25):
    """Precipitation solide annuelle moyenne, par couche (valeur commune)."""
    n_layers = inputs_model.n_layers
    total = None
    for i in range(n_layers):
        contrib = (inputs_model.layer_frac_solid_precip[i]
                   * inputs_model.layer_precip[i] / n_layers)
        total = contrib if total is None else total + contrib
    return np.full(n_layers, float(np.mean(total)) * factor)


# ---------------------------------------------------------------------------
# Modeles couples CemaNeige + GR
# ---------------------------------------------------------------------------


def _run_cemaneige_layers(inputs_model, run_options, param_cn, ind1, is_hyst,
                          n_base=N_STATES_BASE, factor=365.25):
    n_layers = inputs_model.n_layers
    masp = getattr(run_options, "mean_an_solid_precip", None)
    if masp is None:
        masp = mean_an_solid_precip(inputs_model, factor=factor)
        run_options.mean_an_solid_precip = masp

    n_warmup = run_options.ind_period_warmup.shape[0]
    layers = []
    catch_melt_pliq = None
    state_end_cn = []

    for i in range(n_layers):
        base = n_base
        st = np.zeros(4)
        st[0] = run_options.ini_states[base + i]
        st[1] = run_options.ini_states[base + n_layers + i]
        if is_hyst:
            st[2] = run_options.ini_states[base + 2 * n_layers + i]
            st[3] = run_options.ini_states[base + 3 * n_layers + i]

        out, se = K.run_cemaneige(
            inputs_model.layer_precip[i][ind1],
            inputs_model.layer_frac_solid_precip[i][ind1],
            inputs_model.layer_temp_mean[i][ind1],
            float(masp[i]), param_cn, st, is_hyst)

        layer = {name: out[n_warmup:, j].copy()
                 for j, name in enumerate(CEMANEIGE_OUTPUTS)}
        layers.append(layer)
        pm = out[:, CEMANEIGE_OUTPUTS.index("PliqAndMelt")]
        catch_melt_pliq = pm / n_layers if i == 0 else catch_melt_pliq + pm / n_layers
        state_end_cn.append(se)

    return layers, catch_melt_pliq, state_end_cn


def _run_cemaneige_gr(inputs_model, run_options, param, gr_name, n_gr,
                      kernel, thresholds, hourly=False):
    is_hyst = run_options.is_hyst
    n_cn = 4 if is_hyst else 2
    p = _check_param(param, n_gr + n_cn, "CemaNeige" + gr_name)
    p = thresholds(p)
    param_mod = p[:n_gr]
    param_cn = p[n_gr:]

    n_base = N_STATES_BASE_H if hourly else N_STATES_BASE
    factor = 365.25 * 24 if hourly else 365.25

    ind1 = _ind_period1(run_options)
    layers, catch_melt_pliq, state_end_cn = _run_cemaneige_layers(
        inputs_model, run_options, param_cn, ind1, is_hyst, n_base, factor)

    ini = run_options.ini_states.copy()
    lv = run_options.ini_res_levels
    if lv is not None:
        ini[0] = lv[0] * param_mod[0]
        ini[1] = lv[1] * param_mod[2]
        if gr_name == "GR6J":
            ini[2] = lv[2]

    if gr_name == "GR5H":
        imax = run_options.imax
        imax_val = -1.0 if imax is None else float(imax)
        if lv is not None and imax is not None and lv[3] is not None:
            ini[3] = lv[3] * imax_val
        outputs, state_end = kernel(catch_melt_pliq,
                                    inputs_model.pot_evap[ind1], param_mod,
                                    ini[:n_base], imax_val)
    else:
        outputs, state_end = kernel(catch_melt_pliq,
                                    inputs_model.pot_evap[ind1], param_mod,
                                    ini[:n_base])

    # airGR restitue la precipitation observee du bassin, pas celle des couches
    outputs[:, MODEL_OUTPUTS[gr_name].index("Precip")] = \
        inputs_model.precip[ind1]

    n_layers = inputs_model.n_layers
    full_state_end = np.zeros(n_base + 4 * n_layers)
    full_state_end[:n_base] = state_end
    for i, se in enumerate(state_end_cn):
        full_state_end[n_base + i] = se[0]
        full_state_end[n_base + n_layers + i] = se[1]
        full_state_end[n_base + 2 * n_layers + i] = se[2]
        full_state_end[n_base + 3 * n_layers + i] = se[3]

    return _assemble(gr_name, outputs, full_state_end, run_options,
                     inputs_model, p, cemaneige_layers=layers)


def run_model_cemaneige(inputs_model, run_options, param, hourly=False):
    """Module de neige CemaNeige seul (sans modele pluie-debit).

    Renvoie les series par couche d'altitude ainsi que la lame
    `PliqAndMelt` moyenne du bassin, qui alimente le modele GR.

    CemaNeige : Valery (2010) ; hysteresis lineaire : Riboust et al.
    (2019).
    """
    is_hyst = run_options.is_hyst
    n_cn = 4 if is_hyst else 2
    p = _check_param(param, n_cn, "CemaNeige")
    n_base = N_STATES_BASE_H if hourly else N_STATES_BASE
    factor = 365.25 * 24 if hourly else 365.25

    ind1 = _ind_period1(run_options)
    layers, catch_melt_pliq, state_end_cn = _run_cemaneige_layers(
        inputs_model, run_options, p, ind1, is_hyst, n_base, factor)

    n_warmup = run_options.ind_period_warmup.shape[0]
    return {
        "DatesR": inputs_model.dates[run_options.ind_period_run],
        "CemaNeigeLayers": layers,
        "PliqAndMelt": catch_melt_pliq[n_warmup:],
        "StateEnd": np.concatenate(state_end_cn) if state_end_cn else None,
        "Param": p,
    }


def run_model_cemaneige_gr4j(inputs_model, run_options, param):
    """CemaNeige couple a GR4J."""
    return _run_cemaneige_gr(inputs_model, run_options, param, "GR4J", 4,
                             K.run_gr4j, _apply_thresholds_gr4j_like)


def run_model_cemaneige_gr5j(inputs_model, run_options, param):
    """CemaNeige couple a GR5J."""
    return _run_cemaneige_gr(inputs_model, run_options, param, "GR5J", 5,
                             K.run_gr5j, _apply_thresholds_gr4j_like)


def _thresholds_gr6j(p):
    if p[0] < 1e-2:
        p[0] = 1e-2
    if p[2] < 1e-2:
        p[2] = 1e-2
    if p[3] < 0.5:
        p[3] = 0.5
    if p[5] < 1e-2:
        p[5] = 1e-2
    return p


def run_model_cemaneige_gr6j(inputs_model, run_options, param):
    """CemaNeige couple a GR6J."""
    return _run_cemaneige_gr(inputs_model, run_options, param, "GR6J", 6,
                             K.run_gr6j, _thresholds_gr6j)


def run_model_cemaneige_gr4h(inputs_model, run_options, param):
    """CemaNeige couple a GR4H."""
    return _run_cemaneige_gr(inputs_model, run_options, param, "GR4H", 4,
                             KH.run_gr4h, _apply_thresholds_gr4j_like,
                             hourly=True)


def run_model_cemaneige_gr5h(inputs_model, run_options, param):
    """CemaNeige couple a GR5H."""
    return _run_cemaneige_gr(inputs_model, run_options, param, "GR5H", 5,
                             KH.run_gr5h, _apply_thresholds_gr4j_like,
                             hourly=True)


# ---------------------------------------------------------------------------
# ETP d'Oudin
# ---------------------------------------------------------------------------


def pe_oudin(julian_day, temp, lat, lat_unit="rad"):
    """ETP journaliere d'Oudin et al. (2005).

    `julian_day` est le quantieme (1 a 366), `lat` un scalaire ou un
    vecteur de meme longueur que `temp`.
    """
    temp = np.asarray(temp, dtype=float)
    jd = np.asarray(julian_day, dtype=float)
    if lat_unit == "deg":
        lat_rad = np.asarray(lat, dtype=float) / (180.0 / np.pi)
    elif lat_unit == "rad":
        lat_rad = np.asarray(lat, dtype=float)
    else:
        raise ValueError("lat_unit doit valoir 'rad' ou 'deg'")
    if lat_rad.ndim == 0:
        lat_rad = np.full(temp.shape[0], float(lat_rad))
    return K.run_pe_oudin(lat_rad, temp, jd)


MODEL_FUNCS = {
    "GR1A": run_model_gr1a,
    "GR2M": run_model_gr2m,
    "GR4J": run_model_gr4j,
    "GR5J": run_model_gr5j,
    "GR6J": run_model_gr6j,
    "GR4H": run_model_gr4h,
    "GR5H": run_model_gr5h,
    "CemaNeigeGR4J": run_model_cemaneige_gr4j,
    "CemaNeigeGR5J": run_model_cemaneige_gr5j,
    "CemaNeigeGR6J": run_model_cemaneige_gr6j,
    "CemaNeigeGR4H": run_model_cemaneige_gr4h,
    "CemaNeigeGR5H": run_model_cemaneige_gr5h,
}


def run_model(inputs_model, run_options, param, model=None):
    """Point d'entree generique, equivalent de RunModel d'airGR."""
    name = model or run_options.model
    if name not in MODEL_FUNCS:
        raise ValueError("modele inconnu : %s" % name)
    return MODEL_FUNCS[name](inputs_model, run_options, param)
