"""Assimilation de donnees : EnKF et filtre particulaire pour les modeles GR.

Traduction de airGRdatassim 0.1.4 R/RunModel_DA.R, R/DA_EnKF.R et
R/DA_PF.R (GPL-2, INRAE). Les debits observes sont assimiles pas a pas :
a chaque pas observe, les etats de l'ensemble (reservoirs de production et
de routage, premiers etats des hydrogrammes unitaires) sont corriges, puis
le modele repart de ces etats.

Comme dans le code R, l'InputsModel est d'abord restreint a la periode de
simulation : il n'y a donc jamais de periode de chauffe, et le premier pas
de temps tourne pour tous les membres avec les entrees non perturbees et
les niveaux initiaux par defaut. Deux ecritures du code R sont des
affectations mortes (EnsStateBkg au pas suivant, systematiquement re-ecrit
par la boucle des membres) : elles sont signalees mais non reproduites.

GRSUITE-FIX : dans airGRdatassim 0.1.4, DA_EnKF calcule
`IndDa <- which(StateEnKF == 1)` alors que RunModel_DA passe un StateEnKF
caractere (l'usage documente) : la selection est toujours vide et la mise a
jour de Kalman ne s'execute jamais - l'EnKF est un no-op silencieux
(verifie empiriquement : EnsStateA == EnsStateBkg de bout en bout). Ce
port implemente le comportement documente : les noms de `state_enkf`
selectionnent les variables mises a jour. Voir tools/oracle/da_instrumented.R
pour la meme correction cote reference R, et docs/VALIDATION.md.

References
----------
Piazzi, G., Thirel, G., Perrin, C., Delaigue, O. (2021). Sequential data
assimilation for streamflow forecasting: assessing the sensitivity to
uncertainties and updated variables of a conceptual hydrological model at
basin scale. Water Resources Research 57, e2020WR028390.
doi:10.1029/2020WR028390
Piazzi, G., Delaigue, O. (2025). airGRdatassim: Ensemble-Based Data
Assimilation with GR Hydrological Models. R package version 0.1.4,
INRAE, HYCAR Research Unit. doi:10.32614/CRAN.package.airGRdatassim
"""

import math
import warnings

import numpy as np

from ._rng import _NumpyRNG
from .core import InputsModel, RunOptions
from .models import MODEL_FUNCS
from .perturb import DA_MODELS

__all__ = ["run_model_da"]

#: Rang des variables d'etat dans le vecteur d'etat (voir core.py)
_SLOT = {"Prod": 0, "Rout": 1, "UH1": 7, "UH2": 27}

#: 1/sqrt(2*pi), constante M_1_SQRT_2PI de R (dnorm)
_M_1_SQRT_2PI = 0.3989422804014327


def _state_names(model):
    """Variables d'etat assimilees ; GR5J n'a pas d'hydrogramme UH1."""
    if model in ("GR5J", "CemaNeigeGR5J"):
        return ("Prod", "Rout", "UH2")
    return ("Prod", "Rout", "UH1", "UH2")


def _slice_inputs(inputs_model, ind):
    """Restreint un InputsModel a une fenetre, comme InputsModel[IndRun]."""
    n_layers = inputs_model.n_layers
    return InputsModel(
        inputs_model.dates[ind],
        inputs_model.precip[ind],
        inputs_model.pot_evap[ind],
        temp_mean=(None if inputs_model.temp_mean is None
                   else inputs_model.temp_mean[ind]),
        layer_precip=(None if n_layers == 0 else
                      [a[ind] for a in inputs_model.layer_precip]),
        layer_temp_mean=(None if n_layers == 0 else
                         [a[ind] for a in inputs_model.layer_temp_mean]),
        layer_frac_solid_precip=(
            None if n_layers == 0 else
            [a[ind] for a in inputs_model.layer_frac_solid_precip]),
        z_layers=inputs_model.z_layers)


def _check_state_names(values, state_names, label):
    """Equivalent de match.arg(..., several.ok = TRUE) : noms valides."""
    for v in values:
        if v not in state_names:
            raise ValueError(
                "'%s' contient une variable inconnue : %s (attendu parmi : "
                "%s)" % (label, v, ", ".join(state_names)))


def _enkf_constraints(ens, state_names, param):
    """Bornes appliquees aux etats analyses/perturbes (DA_EnKF.R)."""
    r = ens[state_names.index("Prod")]
    r[r < 0.05 * param[0]] = 0.05 * param[0]
    r = ens[state_names.index("Rout")]
    r[r <= 0.0] = 1e-3
    r = ens[state_names.index("UH2")]
    r[r < 0.0] = 1e-3
    if "UH1" in state_names:
        r = ens[state_names.index("UH1")]
        r[r < 0.0] = 1e-3
    r = ens[state_names.index("Prod")]
    r[r > param[0]] = param[0]
    r = ens[state_names.index("Rout")]
    r[r > param[2]] = param[2]


def _da_enkf(obs, qsim, ens_state, param, state_names, nb_mbr, state_enkf,
             state_pert, var_thr, rng):
    """Mise a jour par filtre de Kalman d'ensemble (traduction de DA_EnKF.R).

    GRSUITE-FIX (voir l'en-tete du module) : `state_enkf` selectionne par
    nom les lignes de l'etat a mettre a jour.
    """
    ind_da = np.array([k for k, n in enumerate(state_names)
                       if n in state_enkf])
    state_bkg = ens_state[ind_da, :]

    # covariance d'erreur d'observation et perturbation des observations
    var_obs = max(var_thr**2, (0.1 * obs)**2)
    pert = rng.rnorm(nb_mbr, 0.0, math.sqrt(var_obs))
    obs_pert = obs + pert
    obs_pert[obs_pert < 0.0] = 0.0
    obs_err = np.var(pert, ddof=1)

    # innovations
    innov = obs_pert - qsim

    # gain de Kalman a partir des anomalies d'ensemble
    ens_mean_bkg = state_bkg.mean(axis=1)
    ens_mean_q = qsim.mean()
    anom = state_bkg - ens_mean_bkg[:, None]
    anom_q = qsim - ens_mean_q
    bht = (anom * anom_q[None, :]).sum(axis=1) / (nb_mbr - 1)
    hbht = (anom_q * anom_q).sum() / (nb_mbr - 1)
    k = bht * (1.0 / (hbht + obs_err))

    # analyse
    state_a = state_bkg + k[:, None] * innov[None, :]
    ens_state_enkf = ens_state.copy()
    ens_state_enkf[ind_da, :] = state_a
    _enkf_constraints(ens_state_enkf, state_names, param)

    ans = {"ens_state_enkf": ens_state_enkf, "obs_pert": obs_pert}

    # perturbation d'etat (incertitude sur les variables d'etat)
    if state_pert is not None:
        ind_pert = np.array([n in state_pert for n in state_names])
        sd0 = ens_state_enkf.std(axis=1, ddof=1)
        sd_state = np.minimum(3.0, np.maximum(1.2, sd0))
        nb_state = len(state_names)
        tao = rng.rnorm(nb_state * nb_mbr, np.zeros(nb_state), sd_state)
        tao = tao.reshape(nb_state, nb_mbr)      # byrow = TRUE
        tao[~ind_pert, :] = 0.0
        ens_state_pert = ens_state_enkf + tao
        _enkf_constraints(ens_state_pert, state_names, param)
        ans["ens_state_pert"] = ens_state_pert

    return ans


def _da_pf(obs, qsim, states, param, state_names, nb_mbr, state_pert,
           var_thr, rng):
    """Pesee et reechantillonnage particulaires (traduction de DA_PF.R)."""
    nb_state = len(state_names)
    i_prod = state_names.index("Prod")
    i_rout = state_names.index("Rout")

    # pesee : vraisemblance gaussienne des innovations
    var_obs = max(var_thr**2, (0.1 * obs)**2)
    innov = obs - qsim
    sd_obs = math.sqrt(var_obs)
    z = innov / sd_obs
    with np.errstate(divide="ignore", invalid="ignore", under="ignore"):
        weights = (_M_1_SQRT_2PI * np.exp(-0.5 * z * z)) / sd_obs
        weights = weights / weights.sum()
    if not np.all(np.isfinite(weights)):
        # ensemble ecrase : poids uniformes
        weights = np.full(nb_mbr, 1.0 / nb_mbr)

    # reechantillonnage systematique sur la fonction de repartition
    cdf = np.cumsum(weights)
    a = cdf[0]
    b = min(1.0, a + (1.0 / (nb_mbr + 1)))
    urand0 = a + (b - a) * rng.runif(1)[0]
    step = (1.0 - urand0) / nb_mbr
    urand = urand0 + np.arange(nb_mbr) * step
    indices = np.searchsorted(cdf, urand, side="right")
    indices[indices == nb_mbr] = nb_mbr - 1   # rightmost.closed = TRUE
    indices = indices + 1

    if not np.all(np.isfinite(indices)):
        # branche de repli du code R (inaccessible avec des poids finis,
        # reproduite pour la tracabilite ligne a ligne)
        weights = np.full(nb_mbr, 1.0 / nb_mbr)
        cdf = np.cumsum(weights)
        a = cdf[0]
        b = min(1.0, a + (1.0 / (nb_mbr + 1)))
        urand0 = a + (b - a) * rng.runif(1)[0]
        step = (1.0 - urand0) / nb_mbr
        urand = urand0 + np.arange(nb_mbr) * step
        indices = np.searchsorted(cdf, urand, side="right")
        indices[indices == nb_mbr] = nb_mbr - 1
        indices = indices + 1

    uniq, counts = np.unique(indices, return_counts=True)   # table(Indices)

    if state_pert is not None:
        ind_pert = np.array([n in state_pert for n in state_names])
        ens = np.empty((nb_state, nb_mbr))
        for k, name in enumerate(state_names):
            ens[k] = np.array([s[_SLOT[name]] for s in states])
        sel = ens[:, uniq - 1]
        if sel.shape[1] > 1:
            sd0 = sel.std(axis=1, ddof=1)
        else:
            # R sd() d'une seule valeur vaut NA ; pmax(1.2, NA, na.rm=TRUE)
            sd0 = np.full(nb_state, np.nan)
        sd_state = np.minimum(3.0, np.fmax(1.2, sd0))

    ens_pf = []
    ens_pert = [] if state_pert is not None else None
    for idx, rep in zip(uniq, counts):
        base = states[idx - 1]
        if state_pert is not None and rep > 1:
            # les replications d'une particule selectionnee sont perturbees
            state_rep = np.empty((nb_state, rep))
            for k, name in enumerate(state_names):
                state_rep[k] = base[_SLOT[name]]
            noise = rng.rnorm(rep * nb_state, np.zeros(nb_state), sd_state)
            noise = noise.reshape(nb_state, rep, order="F")   # byrow = FALSE
            noise[~ind_pert, :] = 0.0
            state_rep_pert = state_rep + noise
            r = state_rep_pert[i_prod]
            r[r > param[0]] = param[0]
            r = state_rep_pert[i_rout]
            r[r > param[2]] = param[2]
            for j in range(rep):
                st = base.copy()
                for k, name in enumerate(state_names):
                    st[_SLOT[name]] = state_rep_pert[k, j]
                ens_pert.append(st)
        elif state_pert is not None:
            # particule selectionnee une seule fois : pas de perturbation
            ens_pert.extend([base] * rep)
        ens_pf.extend([base] * rep)

    ans = {"ens_state_pf": ens_pf}
    if state_pert is not None:
        ans["ens_state_pert"] = ens_pert
    return ans


def run_model_da(inputs_model, ind_run, model, param, inputs_pert=None,
                 qobs=None, da_method="EnKF", nb_mbr=None, state_enkf=None,
                 state_pert=None, seed=None, _rng=None):
    """Assimilation de debits observes (traduction de RunModel_DA.R).

    Parameters
    ----------
    inputs_model : InputsModel
        Series d'entree sur toute la periode ; elles sont restreintes a
        ``ind_run`` en interne, comme dans le code R.
    ind_run : array of int
        Indices (base 0) de la periode d'assimilation.
    model : str
        "GR4J", "GR5J", "GR6J" ou une variante CemaNeige.
    param : sequence of float
        Parametres du modele, dans son ordre habituel.
    inputs_pert : InputsPert or None
        Ensembles de forcages perturbes ; sans eux l'ensemble n'est
        diversifie que par l'assimilation elle-meme.
    qobs : array or None
        Debits observes [mm/pas de temps], meme longueur que
        ``inputs_model`` ; les valeurs negatives sont mises a NaN. Si tout
        est manquant ou negatif, ``da_method`` devient "none".
    da_method : str
        "EnKF", "PF" ou "none" (boucle ouverte).
    nb_mbr : int or None
        Nombre de membres (>= 2). Defaut : 50, ou ``inputs_pert.nb_mbr``.
    state_enkf : sequence of str or None
        Variables mises a jour par l'EnKF (requis si ``da_method="EnKF"``),
        parmi "Prod", "Rout", "UH1", "UH2" ("UH1" absent de GR5J).
    state_pert : sequence of str or None
        Variables perturbees apres assimilation ; pour l'EnKF, un
        sous-ensemble de ``state_enkf``.
    seed : int or None
        Graine aleatoire ; comme dans le code R, le generateur est
        reinitialise a chaque pas de temps avec ``seed + i`` (i base 1).

    Returns
    -------
    dict
        "DatesR", "QsimEns" (nb_time, nb_mbr), "EnsStateBkg" et "EnsStateA"
        (nb_time, nb_mbr, nb_state), "ObsPert" (nb_time, nb_mbr),
        "StateNames", "NbTime", "NbMbr", "NbState".
    """
    # ---- controles --------------------------------------------------------
    if model not in DA_MODELS:
        raise ValueError(
            "modele non pris en charge pour l'assimilation : %s "
            "(attendu : %s)" % (model, ", ".join(DA_MODELS)))
    state_names = _state_names(model)

    if da_method not in ("EnKF", "PF", "none"):
        raise ValueError("'da_method' doit valoir 'EnKF', 'PF' ou 'none'")
    if da_method == "none" and (state_enkf is not None
                                or state_pert is not None):
        warnings.warn("'state_enkf' et/ou 'state_pert' ne sont pas pris en "
                      "compte quand da_method='none'", stacklevel=2)
    if da_method == "PF" and state_enkf is not None:
        warnings.warn("'state_enkf' n'est pas pris en compte quand "
                      "da_method='PF'", stacklevel=2)
    if da_method == "EnKF" and state_enkf is None:
        raise ValueError("'state_enkf' doit etre defini quand "
                         "da_method='EnKF'")
    if da_method != "none":
        if state_enkf is not None:
            _check_state_names(state_enkf, state_names, "state_enkf")
        if state_pert is not None:
            _check_state_names(state_pert, state_names, "state_pert")
    if (da_method == "EnKF" and state_pert is not None
            and not set(state_pert) <= set(state_enkf)):
        raise ValueError(
            "la perturbation n'est permise que pour les variables mises a "
            "jour via EnKF (%s) : verifier la coherence entre 'state_pert' "
            "et 'state_enkf'" % ", ".join(state_enkf))

    # InputsPert
    if inputs_pert is None:
        is_meteo = False
        if nb_mbr is None:
            nb_mbr = 50
    else:
        if len(inputs_pert) != len(inputs_model):
            raise ValueError("les elements de 'inputs_pert' doivent avoir "
                             "la meme longueur que ceux d''inputs_model'")
        if inputs_pert.model != model:
            raise ValueError("'inputs_pert' et 'model' ne sont pas "
                             "coherents (%s vs %s)"
                             % (inputs_pert.model, model))
        is_meteo = True
        if nb_mbr is None:
            nb_mbr = inputs_pert.nb_mbr
    if not (isinstance(nb_mbr, (int, np.integer)) and nb_mbr >= 2):
        raise ValueError("'nb_mbr' doit etre un entier >= 2")
    nb_mbr = int(nb_mbr)
    if is_meteo:
        avail = inputs_pert.precip if inputs_pert.precip is not None \
            else inputs_pert.pot_evap
        nb_mbr_meteo = avail.shape[1]
        if nb_mbr > nb_mbr_meteo:
            raise ValueError(
                "impossible de prendre un nombre de membres (%i) superieur "
                "au nombre disponible pour les variables meteo perturbees "
                "(%i)" % (nb_mbr, nb_mbr_meteo))
        if nb_mbr < nb_mbr_meteo:
            warnings.warn(
                "seuls %i membres sont pris, alors que le nombre disponible "
                "pour les variables meteo perturbees est %i"
                % (nb_mbr, nb_mbr_meteo), stacklevel=2)

    # Qobs
    if qobs is None:
        qa = None
    else:
        qa = np.asarray(qobs, dtype=float)
    if qa is None or np.isnan(qa).all() or np.all(qa[~np.isnan(qa)] < 0.0):
        da_method = "none"
        warnings.warn("'da_method' est automatiquement mis a 'none' : "
                      "'qobs' est NULL, tout NaN ou tout negatif",
                      stacklevel=2)
    else:
        if is_meteo and qa.shape[0] != len(inputs_model):
            # le code R ne fait ce controle que si InputsPert est fourni
            raise ValueError("'qobs' doit avoir la meme longueur que les "
                             "elements de 'inputs_model'")
        if np.any(qa[~np.isnan(qa)] < 0.0):
            warnings.warn("les valeurs negatives de 'qobs' sont "
                          "automatiquement mises a NaN", stacklevel=2)

    # ---- initialisations --------------------------------------------------
    is_da = da_method != "none"
    ind_run = np.asarray(ind_run, dtype=np.int64)
    nb_time = ind_run.shape[0]
    nb_state = len(state_names)
    p = np.asarray(param, dtype=float)

    if qa is not None:
        qa = qa.copy()
        qa[qa < 0.0] = np.nan           # avant la restriction, comme dans R
        # inutile (et NaN) en boucle ouverte ; nanquantile sur une serie
        # tout-NaN emettrait un RuntimeWarning
        var_thr = float(np.nanquantile(qa, 0.1)) if is_da else np.nan
        qobs_run = qa[ind_run]
    else:
        var_thr = np.nan
        qobs_run = np.full(nb_time, np.nan)

    sub = _slice_inputs(inputs_model, ind_run)
    rng = _rng if _rng is not None else _NumpyRNG(seed)

    if is_meteo:
        pert_p = (np.tile(sub.precip[:, None], (1, nb_mbr))
                  if inputs_pert.precip is None
                  else inputs_pert.precip[ind_run, :nb_mbr])
        pert_e = (np.tile(sub.pot_evap[:, None], (1, nb_mbr))
                  if inputs_pert.pot_evap is None
                  else inputs_pert.pot_evap[ind_run, :nb_mbr])
        member_inputs = [
            InputsModel(sub.dates, pert_p[:, m], pert_e[:, m],
                        temp_mean=sub.temp_mean,
                        layer_precip=sub.layer_precip,
                        layer_temp_mean=sub.layer_temp_mean,
                        layer_frac_solid_precip=sub.layer_frac_solid_precip,
                        z_layers=sub.z_layers)
            for m in range(nb_mbr)]

    run_fun = MODEL_FUNCS[model]
    # a l'image des RunOptionsIni / RunOptionsIter du code R, deux objets
    # crees une fois et mutes a chaque pas ; jamais de periode de chauffe
    ro_ini = RunOptions(sub, model, ind_period_run=[0], ind_period_warmup=[])
    ro_iter = RunOptions(sub, model, ind_period_run=[0], ind_period_warmup=[])
    ro_iter.ini_res_levels = None

    obs_pert = np.full((nb_mbr, nb_time), np.nan)
    qsim_ens = np.full((nb_mbr, nb_time), np.nan)
    ens_bkg = np.full((nb_state, nb_mbr, nb_time), np.nan)
    ens_a = np.full((nb_state, nb_mbr, nb_time), np.nan)
    ini_states_ens = [None] * nb_mbr

    # ---- boucle temporelle ------------------------------------------------
    for t in range(nb_time):
        if seed is not None:
            rng.set_seed(seed + t + 1)      # iTime est base 1 dans le code R
        for m in range(nb_mbr):
            if t == 0:
                # pas de chauffe et niveaux par defaut, entrees non perturbees
                ro_ini.ind_period_run = np.array([0])
                out = run_fun(sub, ro_ini, p)
            else:
                st = ini_states_ens[m]
                st[4:7] = 0.0               # Store$Rest <- 0 dans le code R
                st[np.isnan(st)] = 0.0
                ro_iter.ini_states = st
                ro_iter.ind_period_run = np.array([t])
                src = member_inputs[m] if is_meteo else sub
                out = run_fun(src, ro_iter, p)

            st = out["StateEnd"]
            ini_states_ens[m] = st
            for k, name in enumerate(state_names):
                if name == "Prod":
                    ens_bkg[k, m, t] = out["Prod"][-1]
                elif name == "Rout":
                    ens_bkg[k, m, t] = out["Rout"][-1]
                else:
                    ens_bkg[k, m, t] = st[_SLOT[name]]
            qsim_ens[m, t] = out["Qsim"][-1]

        # ---- assimilation, si une observation est disponible --------------
        if is_da and np.isfinite(qobs_run[t]):
            if da_method == "EnKF":
                ans = _da_enkf(qobs_run[t], qsim_ens[:, t], ens_bkg[:, :, t],
                               p, state_names, nb_mbr, state_enkf,
                               state_pert, var_thr, rng)
                enkf = ans["ens_state_enkf"]
                for k, name in enumerate(state_names):
                    slot = _SLOT[name]
                    for m in range(nb_mbr):
                        ini_states_ens[m][slot] = enkf[k, m]
                if state_pert is not None:
                    pert = ans["ens_state_pert"]
                    for k, name in enumerate(state_names):
                        slot = _SLOT[name]
                        for m in range(nb_mbr):
                            ini_states_ens[m][slot] = pert[k, m]
                ens_a[:, :, t] = enkf
                obs_pert[:, t] = ans["obs_pert"]
            else:   # da_method == "PF"
                ans = _da_pf(qobs_run[t], qsim_ens[:, t], ini_states_ens,
                             p, state_names, nb_mbr, state_pert,
                             var_thr, rng)
                ens_pf = ans["ens_state_pf"]
                ini_states_ens = (ans["ens_state_pert"]
                                  if state_pert is not None else ens_pf)
                for k, name in enumerate(state_names):
                    slot = _SLOT[name]
                    ens_a[k, :, t] = np.array([s[slot] for s in ens_pf])
            # le code R ecrit aussi EnsStateBkg[, , t+1] ici : affectation
            # morte (la boucle des membres la re-ecrit au pas t+1), non
            # reproduite
        else:
            ens_a[:, :, t] = ens_bkg[:, :, t]
            if da_method == "EnKF":
                obs_pert[:, t] = qobs_run[t]

    return {
        "DatesR": sub.dates,
        "QsimEns": qsim_ens.T,
        "EnsStateBkg": np.transpose(ens_bkg, (2, 1, 0)),
        "EnsStateA": np.transpose(ens_a, (2, 1, 0)),
        "ObsPert": obs_pert.T,
        "StateNames": state_names,
        "NbTime": nb_time,
        "NbMbr": nb_mbr,
        "NbState": nb_state,
    }
