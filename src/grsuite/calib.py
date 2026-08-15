"""Calage automatique par la methode de Michel (algorithme HBAN).

Traduction de airGR/R/Calibration_Michel.R (v1.7.9).

Deux etapes :
  1. criblage d'une grille de jeux de parametres de depart ;
  2. recherche locale par plus forte pente en espace transforme, avec
     pas adaptatif et progression diagonale.

L'ordre d'evaluation des candidats est conserve a l'identique afin que
les egalites soient tranchees comme dans airGR.

References
----------
Michel, C. (1991). Hydrologie appliquee aux petits bassins ruraux.
    Cemagref, Antony.
Coron, L., Thirel, G., Delaigue, O., Perrin, C., Andreassian, V. (2017).
    The suite of lumped GR hydrological models in an R package.
    Environmental Modelling & Software 94, 166-171.
    doi:10.1016/j.envsoft.2017.05.002
"""

import numpy as np

from .crit import error_crit
from .models import MODEL_FUNCS
from .transfo import transfo_param

__all__ = ["OutputsCalib", "calibration_michel"]


class OutputsCalib:
    """Resultat de calage."""

    def __init__(self, param_final_r, crit_final, n_iter, n_runs,
                 hist_param_r, hist_crit, crit_name, crit_best_value):
        self.param_final_r = param_final_r
        self.crit_final = crit_final
        self.n_iter = n_iter
        self.n_runs = n_runs
        self.hist_param_r = hist_param_r
        self.hist_crit = hist_crit
        self.crit_name = crit_name
        self.crit_best_value = crit_best_value

    def __repr__(self):
        return ("OutputsCalib(param=%s, crit=%.4f, n_runs=%i)"
                % (np.round(self.param_final_r, 3), self.crit_final,
                   self.n_runs))


def _expand_grid(columns):
    """Equivalent de expand.grid : la premiere colonne varie le plus vite."""
    uniques = []
    for col in columns:
        vals = []
        for v in col:
            if not any((v == w) or (np.isnan(v) and np.isnan(w)) for w in vals):
                vals.append(v)
        uniques.append(np.asarray(vals, dtype=float))
    n_total = int(np.prod([len(u) for u in uniques]))
    n_col = len(uniques)
    out = np.empty((n_total, n_col))
    rep_each = 1
    for j in range(n_col):
        u = uniques[j]
        n_u = len(u)
        pattern = np.repeat(u, rep_each)
        n_tile = n_total // len(pattern)
        out[:, j] = np.tile(pattern, n_tile)
        rep_each *= n_u
    return out


def _propose_candidates_loc(new_t, old_t, ranges_t, optim_param, pace):
    """Genere jusqu'a 2*NParam candidats autour du point courant."""
    n_param = new_t.shape[0]
    candidates = []
    for i in range(n_param):
        if not optim_param[i]:
            continue
        for j in (1, 2):
            sign = 2 * j - 3  # -1 puis +1
            add = True
            cand = new_t.copy()
            cand[i] = new_t[i] + sign * pace
            if cand[i] < ranges_t[0, i]:
                cand[i] = ranges_t[0, i]
            if cand[i] > ranges_t[1, i]:
                cand[i] = ranges_t[1, i]
            if new_t[i] == ranges_t[0, i] and sign < 0:
                add = False
            if new_t[i] == ranges_t[1, i] and sign > 0:
                add = False
            if old_t is not None and np.array_equal(cand, old_t):
                add = False
            if add:
                candidates.append(cand)
    if not candidates:
        return np.empty((0, n_param))
    return np.vstack(candidates)


def calibration_michel(inputs_model, run_options, inputs_crit, calib_options,
                       model=None, verbose=True):
    """Cale un modele GR par la methode de Michel.

    Parameters
    ----------
    inputs_model : InputsModel
    run_options : RunOptions
    inputs_crit : InputsCrit or InputsCritCompo
    calib_options : CalibOptions
    model : str or None
        Nom du modele ; par defaut `run_options.model`.
    verbose : bool

    Returns
    -------
    OutputsCalib
    """
    model_name = model or run_options.model
    if model_name not in MODEL_FUNCS:
        raise ValueError("modele inconnu : %s" % model_name)
    fun_mod = MODEL_FUNCS[model_name]
    transfo_model = calib_options.model

    n_param = calib_options.n_param
    if n_param > 20:
        raise ValueError("le calage de Michel gere au maximum 20 parametres")

    optim_param = np.isnan(calib_options.fixed_param)
    crit_optim = 1e100
    crit_name = None
    crit_best_value = None
    multiplier = None
    n_runs = 0

    hist_param_r = np.full((500 * n_param, n_param), np.nan)
    hist_crit = np.full(500 * n_param, np.nan)

    def evaluate(param_r):
        outputs_model = fun_mod(inputs_model, run_options, param_r)
        return error_crit(inputs_crit, outputs_model)

    # ---------------------------------------------------------------
    # 1. Criblage des jeux de depart
    # ---------------------------------------------------------------
    if calib_options.start_param_list is not None:
        candidates_r = calib_options.start_param_list.copy()
    else:
        distrib = calib_options.start_param_distrib.copy()
        distrib[:, ~optim_param] = np.nan
        candidates_r = _expand_grid([distrib[:, j] for j in range(n_param)])

    candidates_r = np.atleast_2d(candidates_r).copy()
    candidates_r[:, ~optim_param] = calib_options.fixed_param[~optim_param]

    i_new_optim = -1
    n_candidates = candidates_r.shape[0]
    if verbose and n_candidates > 1:
        print("Criblage en cours (%i jeux)" % n_candidates)

    for i_new in range(n_candidates):
        out_crit = evaluate(candidates_r[i_new])
        if not np.isnan(out_crit.crit_value):
            if out_crit.crit_value * out_crit.multiplier < crit_optim:
                crit_optim = out_crit.crit_value * out_crit.multiplier
                i_new_optim = i_new
        if crit_name is None:
            crit_name = out_crit.crit_name
            crit_best_value = out_crit.crit_best_value
            multiplier = out_crit.multiplier

    if i_new_optim < 0:
        raise RuntimeError("aucun jeu de parametres de depart evaluable")

    param_start_r = candidates_r[i_new_optim].copy()
    param_start_t = transfo_param(param_start_r, "RT", transfo_model)
    crit_start = crit_optim
    n_runs += n_candidates

    if verbose:
        print("\t Criblage termine (%i simulations)" % n_runs)
        print("\t     Param = " + ", ".join("%8.3f" % v for v in param_start_r))
        print("\t     Crit. %-12s = %.4f" % (crit_name, crit_start * multiplier))

    hist_param_r[0] = param_start_r
    hist_crit[0] = crit_start

    # ---------------------------------------------------------------
    # 2. Recherche locale par plus forte pente
    # ---------------------------------------------------------------
    if verbose:
        print("Recherche locale par plus forte pente en cours")

    pace = 0.64
    pace_diag = np.zeros(n_param)
    clg = 0.7 ** (1.0 / n_param)
    compt = 0
    crit_optim = crit_start
    ranges_t = transfo_param(calib_options.search_ranges, "RT", transfo_model)
    new_param_optim_t = param_start_t.copy()
    old_param_optim_t = param_start_t.copy()
    new_param_optim_r = param_start_r.copy()

    iter_done = 0
    for it in range(1, 100 * n_param + 1):
        iter_done = it
        if pace < 0.01:
            break

        candidates_t = _propose_candidates_loc(
            new_param_optim_t, old_param_optim_t, ranges_t, optim_param, pace)
        if candidates_t.shape[0] > 0:
            candidates_r = np.atleast_2d(
                transfo_param(candidates_t, "TR", transfo_model)).copy()
            candidates_r[:, ~optim_param] = calib_options.fixed_param[~optim_param]

            i_new_optim = -1
            for i_new in range(candidates_r.shape[0]):
                out_crit = evaluate(candidates_r[i_new])
                if not np.isnan(out_crit.crit_value):
                    if out_crit.crit_value * out_crit.multiplier < crit_optim:
                        crit_optim = out_crit.crit_value * out_crit.multiplier
                        i_new_optim = i_new
            n_runs += candidates_r.shape[0]
        else:
            i_new_optim = -1

        if i_new_optim >= 0:
            old_param_optim_t = new_param_optim_t.copy()
            new_param_optim_t = candidates_t[i_new_optim].copy()
            compt += 1
            if compt > 2 * n_param:
                pace *= 2.0
                compt = 0
            vect_pace = new_param_optim_t - old_param_optim_t
            for i_c in range(n_param):
                if optim_param[i_c]:
                    pace_diag[i_c] = (clg * pace_diag[i_c]
                                      + (1.0 - clg) * vect_pace[i_c])
        else:
            pace /= 2.0
            compt = 0

        # --- Candidat supplementaire en progression diagonale
        if it > 4 * n_param:
            n_runs += 1
            cand_t = new_param_optim_t + pace_diag
            for i_c in range(n_param):
                if optim_param[i_c]:
                    if cand_t[i_c] < ranges_t[0, i_c]:
                        cand_t[i_c] = ranges_t[0, i_c]
                    if cand_t[i_c] > ranges_t[1, i_c]:
                        cand_t[i_c] = ranges_t[1, i_c]
            cand_r = np.asarray(
                transfo_param(cand_t, "TR", transfo_model), dtype=float).copy()
            cand_r[~optim_param] = calib_options.fixed_param[~optim_param]
            out_crit = evaluate(cand_r)
            if (not np.isnan(out_crit.crit_value)
                    and out_crit.crit_value * out_crit.multiplier < crit_optim):
                crit_optim = out_crit.crit_value * out_crit.multiplier
                old_param_optim_t = new_param_optim_t.copy()
                new_param_optim_t = cand_t.copy()

        new_param_optim_r = np.asarray(
            transfo_param(new_param_optim_t, "TR", transfo_model), dtype=float)
        new_param_optim_r[~optim_param] = calib_options.fixed_param[~optim_param]
        hist_param_r[it] = new_param_optim_r
        hist_crit[it] = crit_optim

    n_iter = iter_done
    if verbose and crit_optim == crit_start:
        print("\t Aucun progres realise")

    param_final_r = new_param_optim_r
    crit_final = crit_optim

    if verbose:
        print("\t Calage termine (%i iterations, %i simulations)"
              % (n_iter, n_runs))
        print("\t     Param = " + ", ".join("%8.3f" % v for v in param_final_r))
        print("\t     Crit. %-12s = %.4f" % (crit_name, crit_final * multiplier))

    return OutputsCalib(
        param_final_r=np.asarray(param_final_r, dtype=float),
        crit_final=crit_final * multiplier,
        n_iter=n_iter,
        n_runs=n_runs,
        hist_param_r=hist_param_r[:n_iter],
        hist_crit=hist_crit[:n_iter] * multiplier,
        crit_name=crit_name,
        crit_best_value=crit_best_value,
    )
