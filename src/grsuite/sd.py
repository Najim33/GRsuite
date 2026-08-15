"""Modelisation semi-distribuee : propagation des apports amont.

Traduction de airGR/R/RunModel_Lag.R (v1.7.9).

Les debits des sous-bassins amont sont propages vers l'exutoire par un
simple decalage temporel, la vitesse de propagation etant le parametre
du modele. Le decalage non entier est reparti lineairement entre deux
pas de temps consecutifs.

References
----------
Lobligeois, F. (2014). Mieux connaitre la distribution spatiale des
    pluies ameliore-t-il la modelisation des crues ? Diagnostic sur 181
    bassins versants francais. These de doctorat, AgroParisTech / Irstea.
de Lavenne, A., Thirel, G., Andreassian, V., Perrin, C., Ramos, M.-H.
    (2016). Spatial variability of the parameters of a semi-distributed
    hydrological model. Proceedings of the IAHS 373, 87-94.
    doi:10.5194/piahs-373-87-2016
"""

import numpy as np

__all__ = ["InputsModelSD", "run_model_lag"]


class InputsModelSD:
    """Entrees d'un schema semi-distribue.

    Parameters
    ----------
    q_upstream : array (n_steps, n_up)
        Debits des sous-bassins amont [m3 par pas de temps].
    length_hydro : array (n_up,)
        Longueur hydraulique de chaque troncon amont [km].
    basin_areas : array (n_up + 1,)
        Surfaces des sous-bassins amont puis du sous-bassin aval [km2].
    """

    def __init__(self, q_upstream, length_hydro, basin_areas):
        self.q_upstream = np.atleast_2d(np.asarray(q_upstream, dtype=float))
        if self.q_upstream.shape[0] == 1 and len(length_hydro) != 1:
            self.q_upstream = self.q_upstream.T
        self.length_hydro = np.asarray(length_hydro, dtype=float)
        self.basin_areas = np.asarray(basin_areas, dtype=float)
        if self.q_upstream.shape[1] != self.length_hydro.shape[0]:
            raise ValueError("'q_upstream' doit avoir une colonne par "
                             "sous-bassin amont")
        if self.basin_areas.shape[0] != self.length_hydro.shape[0] + 1:
            raise ValueError("'basin_areas' doit contenir les surfaces amont "
                             "puis celle du sous-bassin aval")

    @property
    def n_upstream(self):
        return self.length_hydro.shape[0]


def run_model_lag(inputs_sd, run_options, param, q_contrib_down,
                  time_step_seconds=86400.0, warmup_q_down=None):
    """Propage les apports amont et les agrege au debit local aval.

    Parameters
    ----------
    inputs_sd : InputsModelSD
    run_options : RunOptions
        Fournit `ind_period_run` et `ind_period_warmup`.
    param : sequence of one float
        Vitesse de propagation [m/s].
    q_contrib_down : array
        Debit du sous-bassin aval sur la periode de simulation [mm].
    time_step_seconds : float
        Duree du pas de temps [s] : 86400 en journalier, 3600 en horaire.
    warmup_q_down : array or None
        Debit aval sur la periode de chauffe [mm].

    Returns
    -------
    dict avec les cles "Qsim" [mm], "Qsim_m3" [m3] et "QsimDown" [mm].
    """
    speed = float(np.asarray(param, dtype=float).ravel()[0])
    ind_warmup = run_options.ind_period_warmup
    ind_run = run_options.ind_period_run
    ind1 = np.concatenate([ind_warmup, ind_run])
    n_total = ind1.shape[0]
    n_warmup = ind_warmup.shape[0]

    q_down = np.asarray(q_contrib_down, dtype=float)
    if q_down.shape[0] != ind_run.shape[0]:
        raise ValueError("'q_contrib_down' doit couvrir la periode de simulation")
    if warmup_q_down is None:
        head = np.full(n_warmup, np.nan)
    else:
        head = np.asarray(warmup_q_down, dtype=float)
    q_sim_down = np.concatenate([head, q_down])

    # --- temps de propagation, en nombre de pas de temps
    pt = inputs_sd.length_hydro * 1e3 / speed / time_step_seconds
    frac = pt - np.floor(pt)
    hu_trans = np.vstack([1.0 - frac, frac])

    area_down = inputs_sd.basin_areas[-1]
    q_sim_m3 = q_sim_down * area_down * 1e3

    for j in range(inputs_sd.n_upstream):
        n_lag = int(np.floor(pt[j]))
        i_start = max(0, int(ind1[0]) - n_lag - 1)
        i_stop = max(0, int(ind1[0]) - 1)
        ini = inputs_sd.q_upstream[i_start:i_stop + 1, j]
        need = n_lag + 1
        if ini.shape[0] < need:
            pad = np.full(need - ini.shape[0],
                          ini[0] if ini.size else
                          inputs_sd.q_upstream[0, j])
            ini = np.concatenate([pad, ini])
        q_up = np.concatenate([ini, inputs_sd.q_upstream[ind1, j]])
        q_sim_m3 = (q_sim_m3
                    + q_up[1:1 + n_total] * hu_trans[0, j]
                    + q_up[0:n_total] * hu_trans[1, j])

    sl = slice(n_warmup, n_total)
    q_out_m3 = q_sim_m3[sl]
    q_out = q_out_m3 / np.nansum(inputs_sd.basin_areas) / 1e3
    q_out = np.where(q_out < 0, 0.0, q_out)

    return {
        "Qsim": q_out,
        "Qsim_m3": q_out_m3,
        "QsimDown": q_sim_down[sl],
    }
