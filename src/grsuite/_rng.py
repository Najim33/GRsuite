"""Generateurs de tirages pour l'assimilation de donnees.

Deux implementations du meme protocole (``rnorm`` / ``runif`` / ``set_seed``) :

- ``_NumpyRNG`` : production, s'appuie sur ``numpy.random.Generator``. L'ordre
  de consommation et la reinitialisation par pas de temps suivent le code R
  d'airGRdatassim ; le recyclage des vecteurs mean/sd de R est reproduit avec
  ``np.resize``.
- ``_ReplayRNG`` : tests uniquement. Rejoue les tirages exacts exportes par
  tools/oracle/export_da_fixtures.R et verifie que la sequence d'appels du
  port est identique a celle du code R (type et taille de chaque tirage).
"""

import numpy as np

__all__ = ["_NumpyRNG", "_ReplayRNG"]


class _NumpyRNG:
    """Tirages de production (numpy), protocole identique au code R."""

    def __init__(self, seed=None):
        self._g = np.random.default_rng(seed)

    def set_seed(self, seed):
        self._g = np.random.default_rng(seed)

    def rnorm(self, n, mean, sd):
        # R recycle mean/sd sur les n tirages ; np.resize reproduit ce cycle
        mean = np.resize(np.asarray(mean, dtype=float), n)
        sd = np.resize(np.asarray(sd, dtype=float), n)
        return mean + sd * self._g.standard_normal(n)

    def runif(self, n, low=0.0, high=1.0):
        low = np.resize(np.asarray(low, dtype=float), n)
        high = np.resize(np.asarray(high, dtype=float), n)
        return low + (high - low) * self._g.random(n)


class _ReplayRNG:
    """Rejoue une sequence de tirages exportee du code R.

    Parameters
    ----------
    calls : iterable of (str, array)
        Tirages dans l'ordre exact des appels du code R ; le type est
        "rnorm" ou "runif".
    """

    def __init__(self, calls):
        self._calls = [(k, np.asarray(v, dtype=float)) for k, v in calls]
        self._pos = 0

    def set_seed(self, seed):
        # le code R reinitialise sa graine a chaque pas de temps ; en replay
        # les tirages sont deja fixes, il n'y a rien a faire
        pass

    def _pop(self, kind, n):
        if self._pos >= len(self._calls):
            raise ValueError(
                "sequence de replay epuisee a l'appel %i (%s, n=%i)"
                % (self._pos + 1, kind, n))
        k, values = self._calls[self._pos]
        if k != kind or values.shape[0] != n:
            raise ValueError(
                "la sequence d'appels du port differe du code R : a l'appel "
                "%i, R a fait %s(n=%i) mais le port demande %s(n=%i)"
                % (self._pos + 1, k, values.shape[0], kind, n))
        self._pos += 1
        # copie : l'appelant peut muter le tableau (mise a zero des lignes
        # non perturbees) et la sequence stockee doit rester intacte
        return values.copy()

    def rnorm(self, n, mean, sd):
        return self._pop("rnorm", n)

    def runif(self, n, low=0.0, high=1.0):
        return self._pop("runif", n)

    def exhausted(self):
        """Vrai si tous les tirages exportes ont ete consommes."""
        return self._pos == len(self._calls)
