"""Transformations de parametres (espace reel <-> espace transforme).

Traduction de airGR/R/TransfoParam_*.R (v1.7.9).

Direction "RT" : reel -> transforme, "TR" : transforme -> reel.
Les parametres transformes sont bornes par convention dans [-9.99, +9.99].
"""

import numpy as np

__all__ = [
    "transfo_param",
    "TRANSFO_FUNCS",
    "N_PARAM",
]


def _as_matrix(param_in):
    arr = np.asarray(param_in, dtype=float)
    is_vec = arr.ndim == 1
    if is_vec:
        arr = arr.reshape(1, -1)
    return arr.copy(), is_vec


def _finish(out, is_vec):
    return out.ravel() if is_vec else out


def transfo_param_gr4j(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = np.sinh(p[:, 1])
        o[:, 2] = np.exp(p[:, 2])
        o[:, 3] = 20.0 + 19.5 * (p[:, 3] - 9.99) / 19.98
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = np.arcsinh(p[:, 1])
        o[:, 2] = np.log(p[:, 2])
        o[:, 3] = 9.99 + 19.98 * (p[:, 3] - 20.0) / 19.5
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr5j(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = np.sinh(p[:, 1])
        o[:, 2] = np.exp(p[:, 2])
        o[:, 3] = 20.0 + 19.5 * (p[:, 3] - 9.99) / 19.98
        o[:, 4] = (p[:, 4] + 9.99) / 19.98
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = np.arcsinh(p[:, 1])
        o[:, 2] = np.log(p[:, 2])
        o[:, 3] = 9.99 + 19.98 * (p[:, 3] - 20.0) / 19.5
        o[:, 4] = p[:, 4] * 19.98 - 9.99
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr6j(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = np.sinh(p[:, 1])
        o[:, 2] = np.exp(p[:, 2])
        o[:, 3] = 20.0 + 19.5 * (p[:, 3] - 9.99) / 19.98
        o[:, 4] = p[:, 4] / 5.0
        o[:, 5] = np.exp(p[:, 5])
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = np.arcsinh(p[:, 1])
        o[:, 2] = np.log(p[:, 2])
        o[:, 3] = 9.99 + 19.98 * (p[:, 3] - 20.0) / 19.5
        o[:, 4] = p[:, 4] * 5.0
        o[:, 5] = np.log(p[:, 5])
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr2m(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = p[:, 1] / 4.0 + 2.5
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = (p[:, 1] - 2.5) * 4.0
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr1a(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    if direction == "TR":
        o = (p + 10.0) / 8.0
    elif direction == "RT":
        o = p * 8.0 - 10.0
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr4h(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = np.sinh(p[:, 1] / 3.0)
        o[:, 2] = np.exp(p[:, 2])
        o[:, 3] = 480.0 + (480.0 - 0.5) * (p[:, 3] - 9.99) / 19.98
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = 3.0 * np.arcsinh(p[:, 1])
        o[:, 2] = np.log(p[:, 2])
        o[:, 3] = (p[:, 3] - 480.0) * 19.98 / (480.0 - 0.5) + 9.99
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_gr5h(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = np.exp(p[:, 0])
        o[:, 1] = np.sinh(p[:, 1])
        o[:, 2] = np.exp(p[:, 2])
        o[:, 3] = 480.0 + (480.0 - 0.01) * (p[:, 3] - 10.0) / 20.0
        o[:, 4] = (p[:, 4] + 10.0) / 20.0
    elif direction == "RT":
        o[:, 0] = np.log(p[:, 0])
        o[:, 1] = np.arcsinh(p[:, 1])
        o[:, 2] = np.log(p[:, 2])
        o[:, 3] = (p[:, 3] - 480.0) * 20.0 / (480.0 - 0.01) + 10.0
        o[:, 4] = p[:, 4] * 20.0 - 10.0
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_cemaneige(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = (p[:, 0] + 9.99) / 19.98
        o[:, 1] = np.exp(p[:, 1]) / 200.0
    elif direction == "RT":
        o[:, 0] = p[:, 0] * 19.98 - 9.99
        o[:, 1] = np.log(p[:, 1] * 200.0)
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def transfo_param_cemaneige_hyst(param_in, direction):
    p, is_vec = _as_matrix(param_in)
    o = p.copy()
    if direction == "TR":
        o[:, 0] = (p[:, 0] + 9.99) / 19.98
        o[:, 1] = np.exp(p[:, 1]) / 200.0
        o[:, 2] = p[:, 2] * 5.0 + 50.0
        o[:, 3] = p[:, 3] / 19.98 + 0.5
    elif direction == "RT":
        o[:, 0] = p[:, 0] * 19.98 - 9.99
        o[:, 1] = np.log(p[:, 1] * 200.0)
        o[:, 2] = (p[:, 2] - 50.0) / 5.0
        o[:, 3] = (p[:, 3] - 0.5) * 19.98
    else:
        raise ValueError("direction doit valoir 'RT' ou 'TR'")
    return _finish(o, is_vec)


def _make_composite(gr_func, n_gr, snow_func, n_snow):
    """Compose la transfo d'un modele GR et celle de CemaNeige."""

    def _composite(param_in, direction):
        p, is_vec = _as_matrix(param_in)
        o = np.empty_like(p)
        o[:, :n_gr] = np.atleast_2d(gr_func(p[:, :n_gr], direction))
        o[:, n_gr:n_gr + n_snow] = np.atleast_2d(
            snow_func(p[:, n_gr:n_gr + n_snow], direction)
        )
        return _finish(o, is_vec)

    return _composite


TRANSFO_FUNCS = {
    "GR4J": transfo_param_gr4j,
    "GR5J": transfo_param_gr5j,
    "GR6J": transfo_param_gr6j,
    "GR4H": transfo_param_gr4h,
    "GR5H": transfo_param_gr5h,
    "GR2M": transfo_param_gr2m,
    "GR1A": transfo_param_gr1a,
    "CemaNeige": transfo_param_cemaneige,
    "CemaNeigeHyst": transfo_param_cemaneige_hyst,
    "CemaNeigeGR4J": _make_composite(transfo_param_gr4j, 4,
                                     transfo_param_cemaneige, 2),
    "CemaNeigeGR5J": _make_composite(transfo_param_gr5j, 5,
                                     transfo_param_cemaneige, 2),
    "CemaNeigeGR6J": _make_composite(transfo_param_gr6j, 6,
                                     transfo_param_cemaneige, 2),
    "CemaNeigeGR4JHyst": _make_composite(transfo_param_gr4j, 4,
                                         transfo_param_cemaneige_hyst, 4),
    "CemaNeigeGR5JHyst": _make_composite(transfo_param_gr5j, 5,
                                         transfo_param_cemaneige_hyst, 4),
    "CemaNeigeGR6JHyst": _make_composite(transfo_param_gr6j, 6,
                                         transfo_param_cemaneige_hyst, 4),
    "CemaNeigeGR4H": _make_composite(transfo_param_gr4h, 4,
                                     transfo_param_cemaneige, 2),
    "CemaNeigeGR5H": _make_composite(transfo_param_gr5h, 5,
                                     transfo_param_cemaneige, 2),
    "CemaNeigeGR4HHyst": _make_composite(transfo_param_gr4h, 4,
                                         transfo_param_cemaneige_hyst, 4),
    "CemaNeigeGR5HHyst": _make_composite(transfo_param_gr5h, 5,
                                         transfo_param_cemaneige_hyst, 4),
}

N_PARAM = {
    "GR1A": 1,
    "GR2M": 2,
    "GR4J": 4,
    "GR5J": 5,
    "GR6J": 6,
    "GR4H": 4,
    "GR5H": 5,
    "CemaNeige": 2,
    "CemaNeigeHyst": 4,
    "CemaNeigeGR4J": 6,
    "CemaNeigeGR5J": 7,
    "CemaNeigeGR6J": 8,
    "CemaNeigeGR4JHyst": 8,
    "CemaNeigeGR5JHyst": 9,
    "CemaNeigeGR6JHyst": 10,
    "CemaNeigeGR4H": 6,
    "CemaNeigeGR5H": 7,
    "CemaNeigeGR4HHyst": 8,
    "CemaNeigeGR5HHyst": 9,
}


def transfo_param(param_in, direction, model):
    """Applique la transformation associee au modele demande."""
    if model not in TRANSFO_FUNCS:
        raise ValueError("modele inconnu : %s" % model)
    arr = np.asarray(param_in, dtype=float)
    n_expected = N_PARAM[model]
    n_cols = arr.shape[-1] if arr.ndim > 1 else arr.shape[0]
    if n_cols != n_expected:
        raise ValueError(
            "le modele %s requiert %i parametres (recu %i)"
            % (model, n_expected, n_cols)
        )
    return TRANSFO_FUNCS[model](param_in, direction)
