"""Noyaux numeriques horaires de GRsuite (GR4H, GR5H).

Traduction fidele de airGR/src/utils_H.f90, frun_GR4H.f90 et
frun_GR5H.f90 (v1.7.9). Comme au pas de temps journalier, les litteraux
reels du Fortran sont promus depuis la simple precision (voir _kernels).
"""

import numpy as np
from numba import njit

from ._kernels import _F32_09, ss1, ss2

NH_H = 480
_STORED_VAL_H = 759.69140625  # (21/4)**4, exactement representable


# ---------------------------------------------------------------------------
# Hydrogrammes unitaires horaires (utils_H.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def uh1_h(C, D):
    """Ordonnees de HU1 horaire (480 pas)."""
    ord_uh1 = np.zeros(NH_H)
    for i in range(1, NH_H + 1):
        ord_uh1[i - 1] = ss1(i, C, D) - ss1(i - 1, C, D)
    return ord_uh1


@njit(cache=True)
def uh2_h(C, D):
    """Ordonnees de HU2 horaire (960 pas)."""
    ord_uh2 = np.zeros(2 * NH_H)
    for i in range(1, 2 * NH_H + 1):
        ord_uh2[i - 1] = ss2(i, C, D) - ss2(i - 1, C, D)
    return ord_uh2


# ---------------------------------------------------------------------------
# GR4H (frun_GR4H.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr4h(precip, pe, param, state_start):
    """GR4H : 18 sorties MISC, pas de temps horaire."""
    n = precip.shape[0]
    outputs = np.zeros((n, 18))
    state_end = np.zeros(state_start.shape[0])

    st = np.zeros(2)
    st[0] = state_start[0]
    st[1] = state_start[1]
    st_uh1 = np.zeros(NH_H)
    st_uh2 = np.zeros(2 * NH_H)
    for i in range(NH_H):
        st_uh1[i] = state_start[7 + i]
    for i in range(2 * NH_H):
        st_uh2[i] = state_start[7 + NH_H + i]

    D = 1.25
    ord_uh1 = uh1_h(param[3], D)
    ord_uh2 = uh2_h(param[3], D)

    A = param[0]
    B = _F32_09
    nuh1 = max(1, min(NH_H - 1, int(param[3] + 1.0)))
    nuh2 = max(1, min(2 * NH_H - 1, 2 * int(param[3] + 1.0)))

    for k in range(n):
        P1 = precip[k]
        E = pe[k]

        if P1 <= E:
            EN = E - P1
            PN = 0.0
            WS = EN / A
            if WS > 13.0:
                WS = 13.0
            expWS = np.exp(2.0 * WS)
            TWS = (expWS - 1.0) / (expWS + 1.0)
            Sr = st[0] / A
            ER = st[0] * (2.0 - Sr) * TWS / (1.0 + (1.0 - Sr) * TWS)
            AE = ER + P1
            st[0] = st[0] - ER
            PR = 0.0
            PS = 0.0
        else:
            EN = 0.0
            AE = E
            PN = P1 - E
            WS = PN / A
            if WS > 13.0:
                WS = 13.0
            expWS = np.exp(2.0 * WS)
            TWS = (expWS - 1.0) / (expWS + 1.0)
            Sr = st[0] / A
            PS = A * (1.0 - Sr * Sr) * TWS / (1.0 + Sr * TWS)
            PR = PN - PS
            st[0] = st[0] + PS

        if st[0] < 0.0:
            st[0] = 0.0
        Sr = st[0] / param[0]
        Sr = Sr * Sr
        Sr = Sr * Sr
        PERC = st[0] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Sr / _STORED_VAL_H)))
        st[0] = st[0] - PERC
        PR = PR + PERC

        PRHU1 = PR * B
        PRHU2 = PR * (1.0 - B)

        for i in range(nuh1):
            st_uh1[i] = st_uh1[i + 1] + ord_uh1[i] * PRHU1
        st_uh1[NH_H - 1] = ord_uh1[NH_H - 1] * PRHU1

        for i in range(nuh2):
            st_uh2[i] = st_uh2[i + 1] + ord_uh2[i] * PRHU2
        st_uh2[2 * NH_H - 1] = ord_uh2[2 * NH_H - 1] * PRHU2

        Rr = st[1] / param[2]
        EXCH = param[1] * Rr * Rr * Rr * np.sqrt(Rr)

        AEXCH1 = EXCH
        if (st[1] + st_uh1[0] + EXCH) < 0.0:
            AEXCH1 = -st[1] - st_uh1[0]
        st[1] = st[1] + st_uh1[0] + EXCH
        if st[1] < 0.0:
            st[1] = 0.0
        Rr = st[1] / param[2]
        Rr = Rr * Rr
        Rr = Rr * Rr
        QR = st[1] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Rr)))
        st[1] = st[1] - QR

        AEXCH2 = EXCH
        if (st_uh2[0] + EXCH) < 0.0:
            AEXCH2 = -st_uh2[0]
        QD = max(0.0, st_uh2[0] + EXCH)

        Q = QR + QD
        if Q < 0.0:
            Q = 0.0

        outputs[k, 0] = E
        outputs[k, 1] = P1
        outputs[k, 2] = st[0]
        outputs[k, 3] = PN
        outputs[k, 4] = PS
        outputs[k, 5] = AE
        outputs[k, 6] = PERC
        outputs[k, 7] = PR
        outputs[k, 8] = st_uh1[0]
        outputs[k, 9] = st_uh2[0]
        outputs[k, 10] = st[1]
        outputs[k, 11] = EXCH
        outputs[k, 12] = AEXCH1
        outputs[k, 13] = AEXCH2
        outputs[k, 14] = AEXCH1 + AEXCH2
        outputs[k, 15] = QR
        outputs[k, 16] = QD
        outputs[k, 17] = Q

    state_end[0] = st[0]
    state_end[1] = st[1]
    for i in range(NH_H):
        state_end[7 + i] = st_uh1[i]
    for i in range(2 * NH_H):
        state_end[7 + NH_H + i] = st_uh2[i]
    return outputs, state_end


# ---------------------------------------------------------------------------
# GR5H (frun_GR5H.f90), avec reservoir d'interception optionnel
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr5h(precip, pe, param, state_start, imax):
    """GR5H : 21 sorties MISC. Un `imax` negatif desactive l'interception."""
    n = precip.shape[0]
    outputs = np.zeros((n, 21))
    state_end = np.zeros(state_start.shape[0])

    is_int_store = imax >= 0.0

    st = np.zeros(3)
    st[0] = state_start[0]
    st[1] = state_start[1]
    if is_int_store:
        st[2] = state_start[3]
    st_uh2 = np.zeros(2 * NH_H)
    for i in range(2 * NH_H):
        st_uh2[i] = state_start[7 + NH_H + i]

    D = 1.25
    ord_uh2 = uh2_h(param[3], D)

    A = param[0]
    B = _F32_09
    nuh2 = max(1, min(2 * NH_H - 1, 2 * int(param[3] + 1.0)))

    for k in range(n):
        P1 = precip[k]
        E = pe[k]

        if is_int_store:
            # --- Reservoir d'interception (Ficchi, 2017) :
            #     evaporation prioritaire, puis egouttement
            EI = min(E, P1 + st[2])
            PN = max(0.0, P1 - (imax - st[2]) - EI)
            st[2] = st[2] + P1 - EI - PN
            EN = max(0.0, E - EI)

            if EN > 0.0:
                WS = EN / A
                if WS > 13.0:
                    WS = 13.0
                expWS = np.exp(2.0 * WS)
                TWS = (expWS - 1.0) / (expWS + 1.0)
                Sr = st[0] / A
                ES = st[0] * (2.0 - Sr) * TWS / (1.0 + (1.0 - Sr) * TWS)
                st[0] = st[0] - ES
                AE = ES + EI
            else:
                AE = EI
                ES = 0.0

            if PN > 0.0:
                WS = PN / A
                if WS > 13.0:
                    WS = 13.0
                expWS = np.exp(2.0 * WS)
                TWS = (expWS - 1.0) / (expWS + 1.0)
                Sr = st[0] / A
                PS = A * (1.0 - Sr * Sr) * TWS / (1.0 + Sr * TWS)
                PR = PN - PS
                st[0] = st[0] + PS
            else:
                PS = 0.0
                PR = 0.0
        else:
            if P1 <= E:
                EN = E - P1
                PN = 0.0
                WS = EN / A
                if WS > 13.0:
                    WS = 13.0
                expWS = np.exp(2.0 * WS)
                TWS = (expWS - 1.0) / (expWS + 1.0)
                Sr = st[0] / A
                ES = st[0] * (2.0 - Sr) * TWS / (1.0 + (1.0 - Sr) * TWS)
                AE = ES + P1
                EI = P1
                st[0] = st[0] - ES
                PS = 0.0
                PR = 0.0
            else:
                EN = 0.0
                ES = 0.0
                AE = E
                EI = E
                PN = P1 - E
                WS = PN / A
                if WS > 13.0:
                    WS = 13.0
                expWS = np.exp(2.0 * WS)
                TWS = (expWS - 1.0) / (expWS + 1.0)
                Sr = st[0] / A
                PS = A * (1.0 - Sr * Sr) * TWS / (1.0 + Sr * TWS)
                PR = PN - PS
                st[0] = st[0] + PS

        if st[0] < 0.0:
            st[0] = 0.0
        Sr = st[0] / param[0]
        Sr = Sr * Sr
        Sr = Sr * Sr
        PERC = st[0] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Sr / _STORED_VAL_H)))
        st[0] = st[0] - PERC
        PR = PR + PERC

        for i in range(nuh2):
            st_uh2[i] = st_uh2[i + 1] + ord_uh2[i] * PR
        st_uh2[2 * NH_H - 1] = ord_uh2[2 * NH_H - 1] * PR

        Q9 = st_uh2[0] * B
        Q1 = st_uh2[0] * (1.0 - B)

        EXCH = param[1] * (st[1] / param[2] - param[4])

        AEXCH1 = EXCH
        if (st[1] + Q9 + EXCH) < 0.0:
            AEXCH1 = -st[1] - Q9
        st[1] = st[1] + Q9 + EXCH
        if st[1] < 0.0:
            st[1] = 0.0
        Rr = st[1] / param[2]
        Rr = Rr * Rr
        Rr = Rr * Rr
        QR = st[1] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Rr)))
        st[1] = st[1] - QR

        AEXCH2 = EXCH
        if (Q1 + EXCH) < 0.0:
            AEXCH2 = -Q1
        QD = max(0.0, Q1 + EXCH)

        Q = QR + QD
        if Q < 0.0:
            Q = 0.0

        outputs[k, 0] = E
        outputs[k, 1] = P1
        outputs[k, 2] = st[2]
        outputs[k, 3] = st[0]
        outputs[k, 4] = PN
        outputs[k, 5] = PS
        outputs[k, 6] = AE
        outputs[k, 7] = EI
        outputs[k, 8] = ES
        outputs[k, 9] = PERC
        outputs[k, 10] = PR
        outputs[k, 11] = Q9
        outputs[k, 12] = Q1
        outputs[k, 13] = st[1]
        outputs[k, 14] = EXCH
        outputs[k, 15] = AEXCH1
        outputs[k, 16] = AEXCH2
        outputs[k, 17] = AEXCH1 + AEXCH2
        outputs[k, 18] = QR
        outputs[k, 19] = QD
        outputs[k, 20] = Q

    state_end[0] = st[0]
    state_end[1] = st[1]
    state_end[3] = st[2]
    for i in range(2 * NH_H):
        state_end[7 + NH_H + i] = st_uh2[i]
    return outputs, state_end
