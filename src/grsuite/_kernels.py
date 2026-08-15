"""Noyaux numeriques de GRsuite.

Traduction fidele des sources Fortran d'airGR 1.7.9 (INRAE, GPL-2) :
src/utils_D.f90, frun_GR4J.f90, frun_GR5J.f90, frun_GR6J.f90,
frun_GR2M.f90, frun_GR1A.f90, frun_CEMANEIGE.f90, frun_PE.f90.

Toutes les operations sont menees en double precision et dans le meme
ordre que le Fortran d'origine, y compris les "speed-up" algebriques,
afin d'obtenir une egalite au bit pres.
"""

import numpy as np
from numba import njit

NH = 20
NMISC = 30

_STORED_VAL = 25.62890625  # (9/4)**4, exactement representable

# ---------------------------------------------------------------------------
# Litteraux reels du Fortran d'airGR
#
# Dans le code source d'airGR, les constantes comme `0.9` ou `0.4` sont
# ecrites sans suffixe de precision : le compilateur les evalue donc en
# simple precision (REAL(4)) avant de les promouvoir en double precision,
# y compris lorsqu'elles initialisent un `doubleprecision, parameter`.
# Reproduire cette promotion est indispensable pour retrouver les sorties
# d'airGR au bit pres ; l'ignorer introduit un ecart relatif de l'ordre
# de 2.4e-07 (epsilon de la simple precision).
# ---------------------------------------------------------------------------
_F32_09 = float(np.float32(0.9))        # B (GR4J, GR5J, GR6J) et 0.9 * MASP
_F32_04 = float(np.float32(0.4))        # C (GR6J)
_F32_07 = float(np.float32(0.7))        # GR1A
_F32_03 = float(np.float32(0.3))        # GR1A
_F32_01 = float(np.float32(0.1))        # MinSpeed (CemaNeige)
_F32_THIRD = float(np.float32(1.0) / np.float32(3.0))   # 1./3. (GR2M)
_F32_M999 = float(np.float32(-999.999))                  # valeur sentinelle


# ---------------------------------------------------------------------------
# Hydrogrammes unitaires journaliers (utils_D.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def ss1(i, C, D):
    """Courbe en S cumulee de HU1 (SS1)."""
    fi = float(i)
    if fi <= 0.0:
        return 0.0
    if fi < C:
        return (fi / C) ** D
    return 1.0


@njit(cache=True)
def ss2(i, C, D):
    """Courbe en S cumulee de HU2 (SS2)."""
    fi = float(i)
    if fi <= 0.0:
        return 0.0
    if fi <= C:
        return 0.5 * (fi / C) ** D
    if fi < 2.0 * C:
        return 1.0 - 0.5 * (2.0 - fi / C) ** D
    return 1.0


@njit(cache=True)
def uh1(C, D):
    """Ordonnees de HU1 par differences successives sur SS1."""
    ord_uh1 = np.zeros(NH)
    for i in range(1, NH + 1):
        ord_uh1[i - 1] = ss1(i, C, D) - ss1(i - 1, C, D)
    return ord_uh1


@njit(cache=True)
def uh2(C, D):
    """Ordonnees de HU2 par differences successives sur SS2."""
    ord_uh2 = np.zeros(2 * NH)
    for i in range(1, 2 * NH + 1):
        ord_uh2[i - 1] = ss2(i, C, D) - ss2(i - 1, C, D)
    return ord_uh2


# ---------------------------------------------------------------------------
# GR4J (frun_GR4J.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr4j(precip, pe, param, state_start):
    """GR4J : 18 sorties MISC, meme ordre que .FortranOutputs GR."""
    n = precip.shape[0]
    outputs = np.zeros((n, 18))
    state_end = np.zeros(state_start.shape[0])

    st = np.zeros(2)
    st[0] = state_start[0]
    st[1] = state_start[1]
    st_uh1 = np.zeros(NH)
    st_uh2 = np.zeros(2 * NH)
    for i in range(NH):
        st_uh1[i] = state_start[7 + i]
    for i in range(2 * NH):
        st_uh2[i] = state_start[7 + NH + i]

    D = 2.5
    ord_uh1 = uh1(param[3], D)
    ord_uh2 = uh2(param[3], D)

    A = param[0]
    B = _F32_09
    nuh1 = max(1, min(NH - 1, int(param[3] + 1.0)))
    nuh2 = max(1, min(2 * NH - 1, 2 * int(param[3] + 1.0)))

    for k in range(n):
        P1 = precip[k]
        E = pe[k]

        # --- Interception et reservoir de production
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
            PS = 0.0
            PR = 0.0
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

        # --- Percolation
        if st[0] < 0.0:
            st[0] = 0.0
        Sr = st[0] / param[0]
        Sr = Sr * Sr
        Sr = Sr * Sr
        PERC = st[0] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Sr / _STORED_VAL)))
        st[0] = st[0] - PERC
        PR = PR + PERC

        # --- Repartition entre les deux branches de routage
        PRHU1 = PR * B
        PRHU2 = PR * (1.0 - B)

        for i in range(nuh1):
            st_uh1[i] = st_uh1[i + 1] + ord_uh1[i] * PRHU1
        st_uh1[NH - 1] = ord_uh1[NH - 1] * PRHU1

        for i in range(nuh2):
            st_uh2[i] = st_uh2[i + 1] + ord_uh2[i] * PRHU2
        st_uh2[2 * NH - 1] = ord_uh2[2 * NH - 1] * PRHU2

        # --- Echange souterrain potentiel
        Rr = st[1] / param[2]
        EXCH = param[1] * Rr * Rr * Rr * np.sqrt(Rr)

        # --- Reservoir de routage
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

        # --- Branche directe
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
    for i in range(NH):
        state_end[7 + i] = st_uh1[i]
    for i in range(2 * NH):
        state_end[7 + NH + i] = st_uh2[i]
    return outputs, state_end


# ---------------------------------------------------------------------------
# GR5J (frun_GR5J.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr5j(precip, pe, param, state_start):
    """GR5J : 18 sorties MISC."""
    n = precip.shape[0]
    outputs = np.zeros((n, 18))
    state_end = np.zeros(state_start.shape[0])

    st = np.zeros(2)
    st[0] = state_start[0]
    st[1] = state_start[1]
    st_uh2 = np.zeros(2 * NH)
    for i in range(2 * NH):
        st_uh2[i] = state_start[7 + NH + i]

    D = 2.5
    ord_uh2 = uh2(param[3], D)

    A = param[0]
    B = _F32_09
    nuh2 = max(1, min(2 * NH - 1, 2 * int(param[3] + 1.0)))

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
            PS = 0.0
            PR = 0.0
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
        PERC = st[0] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Sr / _STORED_VAL)))
        st[0] = st[0] - PERC
        PR = PR + PERC

        for i in range(nuh2):
            st_uh2[i] = st_uh2[i + 1] + ord_uh2[i] * PR
        st_uh2[2 * NH - 1] = ord_uh2[2 * NH - 1] * PR

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
        outputs[k, 2] = st[0]
        outputs[k, 3] = PN
        outputs[k, 4] = PS
        outputs[k, 5] = AE
        outputs[k, 6] = PERC
        outputs[k, 7] = PR
        outputs[k, 8] = Q9
        outputs[k, 9] = Q1
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
    for i in range(2 * NH):
        state_end[7 + NH + i] = st_uh2[i]
    return outputs, state_end


# ---------------------------------------------------------------------------
# GR6J (frun_GR6J.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr6j(precip, pe, param, state_start):
    """GR6J : 20 sorties MISC (avec reservoir exponentiel)."""
    n = precip.shape[0]
    outputs = np.zeros((n, 20))
    state_end = np.zeros(state_start.shape[0])

    st = np.zeros(3)
    st[0] = state_start[0]
    st[1] = state_start[1]
    st[2] = state_start[2]
    st_uh1 = np.zeros(NH)
    st_uh2 = np.zeros(2 * NH)
    for i in range(NH):
        st_uh1[i] = state_start[7 + i]
    for i in range(2 * NH):
        st_uh2[i] = state_start[7 + NH + i]

    D = 2.5
    ord_uh1 = uh1(param[3], D)
    ord_uh2 = uh2(param[3], D)

    A = param[0]
    B = _F32_09
    C = _F32_04
    nuh1 = max(1, min(NH - 1, int(param[3] + 1.0)))
    nuh2 = max(1, min(2 * NH - 1, 2 * int(param[3] + 1.0)))

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
            PS = 0.0
            PR = 0.0
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
        PERC = st[0] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Sr / _STORED_VAL)))
        st[0] = st[0] - PERC
        PR = PR + PERC

        PRUH1 = PR * B
        PRUH2 = PR * (1.0 - B)

        for i in range(nuh1):
            st_uh1[i] = st_uh1[i + 1] + ord_uh1[i] * PRUH1
        st_uh1[NH - 1] = ord_uh1[NH - 1] * PRUH1

        for i in range(nuh2):
            st_uh2[i] = st_uh2[i + 1] + ord_uh2[i] * PRUH2
        st_uh2[2 * NH - 1] = ord_uh2[2 * NH - 1] * PRUH2

        EXCH = param[1] * (st[1] / param[2] - param[4])

        AEXCH1 = EXCH
        if (st[1] + (1.0 - C) * st_uh1[0] + EXCH) < 0.0:
            AEXCH1 = -st[1] - (1.0 - C) * st_uh1[0]
        st[1] = st[1] + (1.0 - C) * st_uh1[0] + EXCH
        if st[1] < 0.0:
            st[1] = 0.0
        Rr = st[1] / param[2]
        Rr = Rr * Rr
        Rr = Rr * Rr
        QR = st[1] * (1.0 - 1.0 / np.sqrt(np.sqrt(1.0 + Rr)))
        st[1] = st[1] - QR

        # --- Reservoir exponentiel
        st[2] = st[2] + C * st_uh1[0] + EXCH
        AR = st[2] / param[5]
        if AR > 33.0:
            AR = 33.0
        if AR < -33.0:
            AR = -33.0
        if AR > 7.0:
            QRExp = st[2] + param[5] / np.exp(AR)
        elif AR < -7.0:
            QRExp = param[5] * np.exp(AR)
        else:
            QRExp = param[5] * np.log(np.exp(AR) + 1.0)
        st[2] = st[2] - QRExp

        AEXCH2 = EXCH
        if (st_uh2[0] + EXCH) < 0.0:
            AEXCH2 = -st_uh2[0]
        QD = max(0.0, st_uh2[0] + EXCH)

        Q = QR + QD + QRExp
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
        outputs[k, 14] = AEXCH1 + AEXCH2 + EXCH
        outputs[k, 15] = QR
        outputs[k, 16] = QRExp
        outputs[k, 17] = st[2]
        outputs[k, 18] = QD
        outputs[k, 19] = Q

    state_end[0] = st[0]
    state_end[1] = st[1]
    state_end[2] = st[2]
    for i in range(NH):
        state_end[7 + i] = st_uh1[i]
    for i in range(2 * NH):
        state_end[7 + NH + i] = st_uh2[i]
    return outputs, state_end


# ---------------------------------------------------------------------------
# GR2M (frun_GR2M.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr2m(precip, pe, param, state_start):
    """GR2M : 11 sorties MISC."""
    n = precip.shape[0]
    outputs = np.zeros((n, 11))
    state_end = np.zeros(state_start.shape[0])

    st = np.zeros(2)
    st[0] = state_start[0]
    st[1] = state_start[1]

    for k in range(n):
        P = precip[k]
        E = pe[k]

        WS = P / param[0]
        if WS > 13.0:
            WS = 13.0
        expWS = np.exp(2.0 * WS)
        TWS = (expWS - 1.0) / (expWS + 1.0)
        S1 = (st[0] + param[0] * TWS) / (1.0 + st[0] / param[0] * TWS)

        P1 = P + st[0] - S1
        PS = P - P1

        WS = E / param[0]
        if WS > 13.0:
            WS = 13.0
        expWS = np.exp(2.0 * WS)
        TWS = (expWS - 1.0) / (expWS + 1.0)
        S2 = S1 * (1.0 - TWS) / (1.0 + (1.0 - S1 / param[0]) * TWS)
        AE = S1 - S2

        Sr = S2 / param[0]
        Sr = Sr * Sr * Sr + 1.0
        st[0] = S2 / Sr ** _F32_THIRD

        P2 = S2 - st[0]
        P3 = P1 + P2

        R1 = st[1] + P3
        R2 = param[1] * R1
        AEXCH = R2 - R1

        Q = R2 * R2 / (R2 + 60.0)
        st[1] = R2 - Q

        outputs[k, 0] = E
        outputs[k, 1] = P
        outputs[k, 2] = st[0]
        outputs[k, 3] = P1
        outputs[k, 4] = PS
        outputs[k, 5] = AE
        outputs[k, 6] = P2
        outputs[k, 7] = P3
        outputs[k, 8] = st[1]
        outputs[k, 9] = AEXCH
        outputs[k, 10] = Q

    state_end[0] = st[0]
    state_end[1] = st[1]
    return outputs, state_end


# ---------------------------------------------------------------------------
# GR1A (frun_GR1A.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_gr1a(precip, pe, param):
    """GR1A : 3 sorties MISC. Le premier pas de temps n'est pas calcule."""
    n = precip.shape[0]
    outputs = np.full((n, 3), -99e9)
    for k in range(1, n):
        P0 = precip[k - 1]
        P1 = precip[k]
        E1 = pe[k]
        tt = (_F32_07 * P1 + _F32_03 * P0) / param[0] / E1
        Q = P1 * (1.0 - 1.0 / np.sqrt(1.0 + tt * tt))
        outputs[k, 0] = E1
        outputs[k, 1] = P1
        outputs[k, 2] = Q
    return outputs


# ---------------------------------------------------------------------------
# CemaNeige (frun_CEMANEIGE.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_cemaneige(precip, frac_solid, temp, mean_an_solid_precip,
                  param, state_start, is_hyst):
    """CemaNeige, avec ou sans hysteresis lineaire. 11 sorties."""
    n = precip.shape[0]
    outputs = np.zeros((n, 11))
    state_end = np.zeros(4)

    Tmelt = 0.0
    MinSpeed = _F32_01

    G = state_start[0]
    eTG = state_start[1]
    Gratio = 0.0

    CTG = param[0]
    Kf = param[1]

    if is_hyst:
        Gthreshold = state_start[2]
        Glocalmax = state_start[3]
        Gacc = param[2]
        prct = param[3]
        if Gthreshold == 0.0:
            Gthreshold = prct * mean_an_solid_precip
        if Glocalmax == 0.0:
            Glocalmax = Gthreshold
    else:
        Gthreshold = _F32_09 * mean_an_solid_precip
        Glocalmax = _F32_M999
        Gacc = _F32_M999
        prct = _F32_M999

    for k in range(n):
        Pliq = (1.0 - frac_solid[k]) * precip[k]
        Psol = frac_solid[k] * precip[k]

        Ginit = G
        G = G + Psol

        eTG = CTG * eTG + (1.0 - CTG) * temp[k]
        if eTG > 0.0:
            eTG = 0.0

        if eTG == 0.0 and temp[k] > Tmelt:
            PotMelt = Kf * (temp[k] - Tmelt)
            if PotMelt > G:
                PotMelt = G
        else:
            PotMelt = 0.0

        if is_hyst:
            if PotMelt > 0.0:
                if G < Glocalmax and Gratio == 1.0:
                    Glocalmax = G
                Gratio = min(G / Glocalmax, 1.0)
        else:
            if G < Gthreshold:
                Gratio = G / Gthreshold
            else:
                Gratio = 1.0

        Melt = ((1.0 - MinSpeed) * Gratio + MinSpeed) * PotMelt
        G = G - Melt

        if is_hyst:
            dG = G - Ginit
            if dG > 0.0:
                Gratio = min(Gratio + (Psol - Melt) / Gacc, 1.0)
                if Gratio == 1.0:
                    Glocalmax = Gthreshold
            else:
                Gratio = min(G / Glocalmax, 1.0)
        else:
            if G < Gthreshold:
                Gratio = G / Gthreshold
            else:
                Gratio = 1.0

        PliqAndMelt = Pliq + Melt

        outputs[k, 0] = Pliq
        outputs[k, 1] = Psol
        outputs[k, 2] = G
        outputs[k, 3] = eTG
        outputs[k, 4] = Gratio
        outputs[k, 5] = PotMelt
        outputs[k, 6] = Melt
        outputs[k, 7] = PliqAndMelt
        outputs[k, 8] = temp[k]
        outputs[k, 9] = Gthreshold
        outputs[k, 10] = Glocalmax

    state_end[0] = G
    state_end[1] = eTG
    state_end[2] = Gthreshold
    state_end[3] = Glocalmax
    return outputs, state_end


# ---------------------------------------------------------------------------
# Evapotranspiration potentielle d'Oudin (frun_PE.f90)
# ---------------------------------------------------------------------------


@njit(cache=True)
def run_pe_oudin(lat_rad, temp, julian_day):
    """ETP journaliere d'Oudin et al. (2005) [mm/j]."""
    n = temp.shape[0]
    pe = np.zeros(n)
    for k in range(n):
        FI = lat_rad[k]
        DT = temp[k]
        JD = julian_day[k]

        COSFI = np.cos(FI)
        TETA = 0.4093 * np.sin(JD / 58.1 - 1.405)
        COSTETA = np.cos(TETA)
        COSGZ = max(0.001, np.cos(FI - TETA))
        # Le Fortran d'airGR calcule ici GZ et SINGZ, qui ne servent pas
        # dans la suite de la formule d'Oudin : on les omet.

        COSOM = 1.0 - COSGZ / COSFI / COSTETA
        if COSOM < -1.0:
            COSOM = -1.0
        if COSOM > 1.0:
            COSOM = 1.0
        COSOM2 = COSOM * COSOM
        if COSOM2 >= 1.0:
            SINOM = 0.0
        else:
            SINOM = np.sqrt(1.0 - COSOM2)

        OM = np.arccos(COSOM)
        COSPZ = COSGZ + COSFI * COSTETA * (SINOM / OM - 1.0)
        if COSPZ < 0.001:
            COSPZ = 0.001
        ETA = 1.0 + np.cos(JD / 58.1) / 30.0
        GE = 446.0 * OM * COSPZ * ETA

        pe[k] = max(0.0, GE * (DT + 5.0) / 100.0 / 28.5)
    return pe
