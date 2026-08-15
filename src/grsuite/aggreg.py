"""Agregation temporelle des series (SeriesAggreg).

Traduction de airGR/R/SeriesAggreg.data.frame.R (v1.7.9).

Le comportement d'airGR est reproduit, y compris le remplissage
prealable de la serie par des pas de temps manquants : une periode
incomplete en debut ou fin de serie produit une valeur manquante, et
les groupes dont la premiere date n'existait pas dans la serie
d'origine sont ecartes.
"""

import numpy as np

__all__ = ["series_aggreg"]

_SUPPORTED = ("%Y%m%d", "%Y%m", "%Y")


def _infer_step(dates):
    """Pas de temps de la serie, deduit des deux premieres dates."""
    if dates.shape[0] < 2:
        raise ValueError("au moins deux dates sont necessaires")
    return dates[1] - dates[0]


def _group_keys(dates, fmt, year_first_month):
    """Identifiant de groupe et marqueur de debut de groupe."""
    d = dates.astype("datetime64[D]")
    months = d.astype("datetime64[M]")

    if fmt == "%Y%m%d":
        key = d.astype("datetime64[D]").astype("int64")
    elif fmt == "%Y%m":
        key = months.astype("int64")
    elif fmt == "%Y":
        key = months.astype("int64")
    else:
        raise ValueError("format non pris en charge : %s" % fmt)

    first = np.empty(key.shape[0], dtype=bool)
    first[0] = True
    first[1:] = key[1:] != key[:-1]

    if fmt == "%Y":
        month_of = (months.astype("int64") % 12) + 1
        first = first & (month_of == year_first_month)

    # Comme dans airGR, le groupe 0 rassemble les pas de temps qui precedent
    # le premier debut de periode ; il est ecarte a la fin.
    fac = np.cumsum(first)
    return fac, first


def series_aggreg(dates, data, fmt="%Y%m", convert_fun="sum",
                  year_first_month=1):
    """Agrege des series temporelles a un pas de temps plus large.

    Parameters
    ----------
    dates : array of numpy.datetime64
        Dates de la serie d'origine, regulieres et triees.
    data : dict of {str: array}
        Series a agreger.
    fmt : str
        "%Y%m%d" (journalier), "%Y%m" (mensuel) ou "%Y" (annuel).
    convert_fun : str or dict
        "sum" ou "mean", globalement ou par variable.
    year_first_month : int
        Premier mois de l'annee hydrologique (utilise si fmt vaut "%Y").

    Returns
    -------
    (dates_agg, data_agg) : tuple
    """
    if fmt not in _SUPPORTED:
        raise ValueError("format non pris en charge : %s" % fmt)
    dates = np.asarray(dates)
    names = list(data.keys())
    if isinstance(convert_fun, str):
        convert_fun = dict.fromkeys(names, convert_fun)

    step = _infer_step(dates)

    # --- remplissage par une serie reguliere complete, comme dans airGR
    y0 = dates[0].astype("datetime64[Y]").astype(int) + 1970 - 1
    y1 = dates[-1].astype("datetime64[Y]").astype(int) + 1970 + 1
    start = np.datetime64("%04d-01-01" % y0).astype(dates.dtype)
    stop = np.datetime64("%04d-12-31" % y1).astype(dates.dtype)
    full = np.arange(start, stop + step, step, dtype=dates.dtype)

    pos = np.searchsorted(full, dates)
    in_original = np.zeros(full.shape[0], dtype=bool)
    valid = (pos < full.shape[0])
    valid &= full[np.clip(pos, 0, full.shape[0] - 1)] == dates
    in_original[pos[valid]] = True

    padded = {}
    for k in names:
        col = np.full(full.shape[0], np.nan)
        col[pos[valid]] = np.asarray(data[k], dtype=float)[valid]
        padded[k] = col

    fac, first = _group_keys(full, fmt, year_first_month)
    n_groups = int(fac.max()) + 1

    out = {}
    for k in names:
        col = padded[k]
        sums = np.bincount(fac, weights=np.nan_to_num(col), minlength=n_groups)
        n_nan = np.bincount(fac, weights=np.isnan(col).astype(float),
                            minlength=n_groups)
        counts = np.bincount(fac, minlength=n_groups)
        if convert_fun[k] == "sum":
            agg = sums
        elif convert_fun[k] == "mean":
            with np.errstate(invalid="ignore", divide="ignore"):
                agg = sums / counts
        else:
            raise ValueError("fonction d'agregation inconnue : %s"
                             % convert_fun[k])
        agg[n_nan > 0] = np.nan   # na.action = na.pass : NA propage
        out[k] = agg

    # --- on ne garde que les groupes dont la premiere date existait
    keep_idx = np.where(first & in_original)[0]
    keep_groups = fac[keep_idx]
    dates_agg = full[keep_idx]
    data_agg = {k: out[k][keep_groups] for k in names}
    return dates_agg, data_agg
