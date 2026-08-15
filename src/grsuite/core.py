"""Objets d'entree de GRsuite : InputsModel, RunOptions, CalibOptions.

Traduction de airGR/R/CreateInputsModel.R, CreateRunOptions.R,
CreateIniStates.R, CreateCalibOptions.R (v1.7.9).

Convention d'indexation : les indices de periode sont donnes en base 0
(Python), la conversion depuis les indices R en base 1 est faite par
l'appelant.

Vecteur d'etat (identique a airGR) :
    0        niveau du reservoir de production      [mm]
    1        niveau du reservoir de routage         [mm]
    2        niveau du reservoir exponentiel        [mm]
    3        niveau du reservoir d'interception     [mm]
    4 - 6    reserve
    7 - 26   etats de HU1                           [mm]
    27 - 66  etats de HU2                           [mm]
    67 +     couches CemaNeige : G, eTG, Gthr, Glocmax
"""

import numpy as np

from .transfo import N_PARAM, transfo_param

__all__ = ["InputsModel", "RunOptions", "CalibOptions", "MODEL_OUTPUTS",
           "CEMANEIGE_OUTPUTS", "N_STATES_BASE", "N_STATES_BASE_H",
           "n_states_base", "is_hourly"]

N_STATES_BASE = 7 + 20 + 40      # 67, pas de temps journalier / mensuel
N_STATES_BASE_H = 7 + 480 + 960  # 1447, pas de temps horaire


def is_hourly(model):
    """Vrai pour les modeles au pas de temps horaire (GR4H, GR5H)."""
    return model.endswith("H") or model.endswith("HHyst")


def n_states_base(model):
    """Longueur de la partie non-CemaNeige du vecteur d'etat."""
    return N_STATES_BASE_H if is_hourly(model) else N_STATES_BASE

MODEL_OUTPUTS = {
    "GR1A": ["PotEvap", "Precip", "Qsim"],
    "GR2M": ["PotEvap", "Precip", "Prod", "Pn", "Ps", "AE", "Perc", "PR",
             "Rout", "AExch", "Qsim"],
    "GR4J": ["PotEvap", "Precip", "Prod", "Pn", "Ps", "AE", "Perc", "PR",
             "Q9", "Q1", "Rout", "Exch", "AExch1", "AExch2", "AExch", "QR",
             "QD", "Qsim"],
    "GR5J": ["PotEvap", "Precip", "Prod", "Pn", "Ps", "AE", "Perc", "PR",
             "Q9", "Q1", "Rout", "Exch", "AExch1", "AExch2", "AExch", "QR",
             "QD", "Qsim"],
    "GR6J": ["PotEvap", "Precip", "Prod", "Pn", "Ps", "AE", "Perc", "PR",
             "Q9", "Q1", "Rout", "Exch", "AExch1", "AExch2", "AExch", "QR",
             "QRExp", "Exp", "QD", "Qsim"],
    "GR4H": ["PotEvap", "Precip", "Prod", "Pn", "Ps", "AE", "Perc", "PR",
             "Q9", "Q1", "Rout", "Exch", "AExch1", "AExch2", "AExch", "QR",
             "QD", "Qsim"],
    "GR5H": ["PotEvap", "Precip", "Interc", "Prod", "Pn", "Ps", "AE", "EI",
             "ES", "Perc", "PR", "Q9", "Q1", "Rout", "Exch", "AExch1",
             "AExch2", "AExch", "QR", "QD", "Qsim"],
}

CEMANEIGE_OUTPUTS = ["Pliq", "Psol", "SnowPack", "ThermalState", "Gratio",
                     "PotMelt", "Melt", "PliqAndMelt", "Temp", "Gthreshold",
                     "Glocalmax"]

# Points de depart du criblage (CreateCalibOptions.R), en espace transforme
START_PARAM_DISTRIB_T = {
    "GR4J": [[+5.13, -1.60, +3.03, -9.05],
             [+5.51, -0.61, +3.74, -8.51],
             [+6.07, -0.02, +4.42, -8.06]],
    "GR5J": [[+5.17, -1.13, +3.08, -9.37, -7.45],
             [+5.55, -0.46, +3.75, -9.09, -4.69],
             [+6.10, -0.11, +4.43, -8.60, -0.66]],
    "GR6J": [[+3.60, -1.00, +3.30, -9.10, -0.90, +3.00],
             [+3.90, -0.50, +4.10, -8.70, +0.10, +4.00],
             [+4.50, +0.50, +5.00, -8.10, +1.10, +5.00]],
    "GR4H": [[+5.12, -1.18, +4.34, -9.69],
             [+5.58, -0.85, +4.74, -9.47],
             [+6.01, -0.50, +5.14, -8.87]],
    # GR5H sans reservoir d'interception
    "GR5H": [[+3.28, -0.39, +4.14, -9.54, -7.49],
             [+3.62, -0.19, +4.80, -9.00, -6.31],
             [+4.01, -0.04, +5.43, -7.53, -5.33]],
    # GR5H avec reservoir d'interception
    "GR5H_int": [[+3.46, -1.25, +4.04, -9.53, -9.34],
                 [+3.74, -0.41, +4.78, -8.94, -3.33],
                 [+4.29, +0.16, +5.39, -7.39, +3.33]],
    "GR2M": [[+5.03, -7.15], [+5.22, -6.74], [+5.85, -6.37]],
    "GR1A": [[-1.69], [-0.38], [+1.39]],
    "CemaNeige": [[-9.96, +6.63], [-9.14, +6.90], [+4.10, +7.21]],
}

_SNOW_T = [[-9.96, +6.63], [-9.14, +6.90], [+4.10, +7.21]]
_HYST_T = [[-7.00, -7.00], [-0.00, -0.00], [+7.00, +7.00]]

for _gr in ("GR4J", "GR5J", "GR6J", "GR4H", "GR5H"):
    START_PARAM_DISTRIB_T["CemaNeige" + _gr] = [
        START_PARAM_DISTRIB_T[_gr][i] + _SNOW_T[i] for i in range(3)
    ]
    START_PARAM_DISTRIB_T["CemaNeige" + _gr + "Hyst"] = [
        START_PARAM_DISTRIB_T[_gr][i] + _SNOW_T[i] + _HYST_T[i]
        for i in range(3)
    ]


class InputsModel:
    """Series d'entree du modele.

    Parameters
    ----------
    dates : array of numpy.datetime64
        Dates du pas de temps, continues et triees.
    precip : array
        Precipitation totale [mm/pas de temps].
    pot_evap : array
        Evapotranspiration potentielle [mm/pas de temps].
    temp_mean : array or None
        Temperature moyenne de l'air [degC], requise pour CemaNeige.
    layer_precip, layer_temp_mean, layer_frac_solid_precip : list of array
        Series par couche d'altitude (CemaNeige).
    z_layers : array or None
        Altitude mediane de chaque couche [m].
    """

    def __init__(self, dates, precip, pot_evap, temp_mean=None,
                 layer_precip=None, layer_temp_mean=None,
                 layer_frac_solid_precip=None, z_layers=None):
        self.dates = np.asarray(dates)
        self.precip = np.asarray(precip, dtype=float)
        self.pot_evap = np.asarray(pot_evap, dtype=float)
        self.temp_mean = None if temp_mean is None else np.asarray(temp_mean, float)
        self.layer_precip = layer_precip
        self.layer_temp_mean = layer_temp_mean
        self.layer_frac_solid_precip = layer_frac_solid_precip
        self.z_layers = None if z_layers is None else np.asarray(z_layers, float)

        n = self.precip.shape[0]
        if self.pot_evap.shape[0] != n or self.dates.shape[0] != n:
            raise ValueError("'dates', 'precip' et 'pot_evap' doivent avoir "
                             "la meme longueur")

    @property
    def n_layers(self):
        return 0 if self.layer_precip is None else len(self.layer_precip)

    def __len__(self):
        return self.precip.shape[0]


class RunOptions:
    """Options de simulation : periodes, etats initiaux, sorties.

    Parameters
    ----------
    inputs_model : InputsModel
    model : str
        "GR4J", "GR5J", "GR6J", "GR2M", "GR1A", "CemaNeigeGR4J", ...
    ind_period_run : array of int
        Indices (base 0) de la periode de simulation, continus.
    ind_period_warmup : array of int or None
        Indices de la periode de chauffe. None declenche le comportement
        par defaut d'airGR (au plus une annee precedant la simulation).
    ini_states : array or None
        Vecteur d'etat initial complet.
    ini_res_levels : sequence or None
        Taux de remplissage initiaux (production, routage, exponentiel,
        interception). Defaut : (0.3, 0.5, 0, None).
    is_hyst : bool
        Active l'hysteresis lineaire de CemaNeige.
    """

    def __init__(self, inputs_model, model, ind_period_run,
                 ind_period_warmup=None, ini_states=None,
                 ini_res_levels=None, is_hyst=False, imax=None):
        self.model = model
        self.is_hyst = is_hyst
        self.inputs_model = inputs_model
        self.imax = imax

        ind_run = np.asarray(ind_period_run, dtype=np.int64)
        if ind_run.size == 0:
            raise ValueError("'ind_period_run' ne peut pas etre vide")
        if not np.all(np.diff(ind_run) == 1):
            raise ValueError("'ind_period_run' doit etre une sequence continue")
        self.ind_period_run = ind_run

        if ind_period_warmup is None:
            self.ind_period_warmup = self._default_warmup(inputs_model, ind_run)
        else:
            self.ind_period_warmup = np.asarray(ind_period_warmup, dtype=np.int64)

        n_layers = inputs_model.n_layers
        n_states = n_states_base(model) + 4 * n_layers
        if ini_states is None:
            self.ini_states = np.zeros(n_states)
        else:
            self.ini_states = np.asarray(ini_states, dtype=float).copy()
            if self.ini_states.shape[0] != n_states:
                raise ValueError("'ini_states' doit contenir %i valeurs" % n_states)

        if ini_res_levels is None:
            if model.endswith("GR6J") or model == "GR6J":
                self.ini_res_levels = (0.3, 0.5, 0.0, None)
            elif model.endswith("GR5H") or model == "GR5H":
                self.ini_res_levels = (0.3, 0.5, None,
                                       0.0 if imax is not None else None)
            else:
                self.ini_res_levels = (0.3, 0.5, None, None)
        else:
            self.ini_res_levels = tuple(ini_res_levels)

    @staticmethod
    def _day_of_month(d):
        day = d.astype("datetime64[D]")
        return int((day - day.astype("datetime64[M]")).astype(int)) + 1

    @staticmethod
    def _default_warmup(inputs_model, ind_run):
        """Reproduit le choix par defaut de CreateRunOptions.

        airGR remonte d'une annee avant le debut de la simulation. Le
        retrait de 365 jours fixes decale la date d'un jour lorsqu'une
        annee bissextile est traversee ; airGR detecte ce decalage en
        comparant le quantieme du mois et retire un jour de plus. Cette
        correction est reproduite ici a l'identique.
        """
        if ind_run[0] == 0:
            return np.array([], dtype=np.int64)
        d0 = inputs_model.dates[ind_run[0]]
        target = d0 - np.timedelta64(365, "D")
        if RunOptions._day_of_month(target) != RunOptions._day_of_month(d0):
            target = target - np.timedelta64(1, "D")
        first_available = inputs_model.dates[0]
        start_date = max(first_available, target)
        matches = np.where(inputs_model.dates == start_date)[0]
        start = int(matches[0]) if matches.size else 0
        return np.arange(start, ind_run[0], dtype=np.int64)


class CalibOptions:
    """Options de calage (algorithme HBAN de Michel).

    Parameters
    ----------
    model : str
    fixed_param : sequence or None
        Valeurs imposees (NaN pour les parametres a optimiser).
    search_ranges : array (2, n) or None
        Bornes en espace reel. Defaut : transformation de [-9.99, +9.99].
    start_param_distrib : array or None
        Grille de depart (espace reel). Defaut : valeurs d'airGR.
    start_param_list : array or None
        Liste de jeux de depart (espace reel), alternative a la grille.
    """

    def __init__(self, model, fixed_param=None, search_ranges=None,
                 start_param_distrib=None, start_param_list=None,
                 is_int_store=False):
        if model not in N_PARAM:
            raise ValueError("modele inconnu : %s" % model)
        self.model = model
        self.is_int_store = is_int_store
        n = N_PARAM[model]
        self.n_param = n

        if fixed_param is None:
            self.fixed_param = np.full(n, np.nan)
        else:
            self.fixed_param = np.asarray(fixed_param, dtype=float)
            if self.fixed_param.shape[0] != n:
                raise ValueError("'fixed_param' doit contenir %i valeurs" % n)

        if search_ranges is None:
            param_t = np.array([[-9.99] * n, [+9.99] * n])
            self.search_ranges = transfo_param(param_t, "TR", model)
        else:
            self.search_ranges = np.asarray(search_ranges, dtype=float)
            if self.search_ranges.shape != (2, n):
                raise ValueError("'search_ranges' doit etre de forme (2, %i)" % n)

        self.start_param_list = (None if start_param_list is None
                                 else np.asarray(start_param_list, dtype=float))
        if start_param_distrib is not None:
            self.start_param_distrib = np.asarray(start_param_distrib, float)
        elif start_param_list is None:
            key = model
            if is_int_store and model == "GR5H":
                key = "GR5H_int"   # jeux de depart specifiques a l'interception
            if key not in START_PARAM_DISTRIB_T:
                raise ValueError("pas de grille par defaut pour %s" % model)
            param_t = np.asarray(START_PARAM_DISTRIB_T[key], dtype=float)
            self.start_param_distrib = transfo_param(param_t, "TR", model)
        else:
            self.start_param_distrib = None
