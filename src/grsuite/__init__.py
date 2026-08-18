"""GRsuite - suite de modeles hydrologiques GR en Python.

Projet personnel et NON OFFICIEL : GRsuite n'est ni produit, ni approuve, ni
maintenu par l'INRAE, et n'est pas une version d'airGR. Les modeles et
l'implementation de reference sont ceux de l'unite HYCAR de l'INRAE (voir
NOTICE et docs/REFERENCES.md) ; seule la traduction en Python est de l'auteur.

Reimplementation complete du package R airGR 1.7.9 (INRAE, HYCAR) :
modeles GR1A, GR2M, GR4J, GR5J, GR6J, module de neige CemaNeige (avec
ou sans hysteresis), ETP d'Oudin, criteres d'erreur et calage par la
methode de Michel.

Exemple minimal
---------------
>>> import numpy as np
>>> from grsuite import InputsModel, RunOptions, run_model_gr4j
>>> dates = np.arange("2000-01-01", "2001-01-01", dtype="datetime64[D]")
>>> im = InputsModel(dates, precip=np.zeros(366), pot_evap=np.ones(366))
>>> ro = RunOptions(im, "GR4J", ind_period_run=np.arange(366))
>>> out = run_model_gr4j(im, ro, [350.0, 0.0, 90.0, 1.7])
>>> out["Qsim"].shape
(366,)
"""

__version__ = "1.1.0"

from .aggreg import series_aggreg
from .api import (
    Assimilation,
    CalibratedModel,
    Catchment,
    Simulation,
    list_models,
    param_names,
)
from .calib import OutputsCalib, calibration_michel
from .core import (
    CEMANEIGE_OUTPUTS,
    MODEL_OUTPUTS,
    CalibOptions,
    InputsModel,
    RunOptions,
)
from .crit import (
    InputsCrit,
    InputsCritCompo,
    error_crit,
    error_crit_kge,
    error_crit_kge2,
    error_crit_nse,
    error_crit_rmse,
)
from .da import run_model_da
from .models import (
    MODEL_FUNCS,
    data_alti_extrapolation_valery,
    imax_estimate,
    mean_an_solid_precip,
    pe_oudin,
    run_model,
    run_model_cemaneige,
    run_model_cemaneige_gr4h,
    run_model_cemaneige_gr4j,
    run_model_cemaneige_gr5h,
    run_model_cemaneige_gr5j,
    run_model_cemaneige_gr6j,
    run_model_gr1a,
    run_model_gr2m,
    run_model_gr4h,
    run_model_gr4j,
    run_model_gr5h,
    run_model_gr5j,
    run_model_gr6j,
)
from .perturb import InputsPert
from .sd import InputsModelSD, run_model_lag
from .transfo import transfo_param

__all__ = [
    # interface haut niveau
    "Catchment", "Simulation", "CalibratedModel", "Assimilation",
    "list_models", "param_names",
    # assimilation de donnees
    "InputsPert", "run_model_da",
    # interface fidele a airGR
    "InputsModel", "RunOptions", "CalibOptions",
    "InputsCrit", "InputsCritCompo", "error_crit",
    "error_crit_nse", "error_crit_kge", "error_crit_kge2", "error_crit_rmse",
    "calibration_michel", "OutputsCalib",
    "run_model", "run_model_gr1a", "run_model_gr2m", "run_model_gr4j",
    "run_model_gr5j", "run_model_gr6j", "run_model_gr4h", "run_model_gr5h",
    "run_model_cemaneige", "run_model_cemaneige_gr4j",
    "run_model_cemaneige_gr5j",
    "run_model_cemaneige_gr6j", "run_model_cemaneige_gr4h",
    "run_model_cemaneige_gr5h",
    "pe_oudin", "data_alti_extrapolation_valery", "mean_an_solid_precip",
    "imax_estimate", "series_aggreg", "InputsModelSD", "run_model_lag",
    "transfo_param", "MODEL_FUNCS", "MODEL_OUTPUTS",
    "CEMANEIGE_OUTPUTS", "__version__",
]
