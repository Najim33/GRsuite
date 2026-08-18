# References

Every model, module, criterion and algorithm in GRsuite comes from published work.
This page lists the source of each one, so that a result produced with GRsuite can
be traced back to the paper or thesis that defines it.

GRsuite itself introduces no new hydrology: it is an unofficial, independent
Python translation of **airGR** (INRAE, HYCAR), which is the reference
implementation of all of the below. GRsuite is a personal project and is not
endorsed by INRAE. Work published with GRsuite should cite the model reference
*and* airGR — not GRsuite alone.

*Cette page existe en un seul exemplaire : une bibliographie est la même dans les
deux langues. Le README français y renvoie directement.*

---

## 1. The package this one reimplements

> **Coron, L., Thirel, G., Delaigue, O., Perrin, C., Andréassian, V.** (2017).
> The suite of lumped GR hydrological models in an R package.
> *Environmental Modelling & Software*, 94, 166–171.
> doi:[10.1016/j.envsoft.2017.05.002](https://doi.org/10.1016/j.envsoft.2017.05.002)

> **Coron, L., Delaigue, O., Thirel, G., Dorchies, D., Perrin, C., Michel, C.**
> airGR: Suite of GR Hydrological Models for Precipitation-Runoff Modelling.
> R package version 1.7.9. INRAE, HYCAR Research Unit, Antony, France.
> doi:[10.15454/EX11NA](https://doi.org/10.15454/EX11NA) ·
> [CRAN](https://cran.r-project.org/package=airGR) ·
> [website](https://hydrogr.github.io/airGR/)

The GR models themselves were developed at Cemagref / Irstea / INRAE over three
decades, under the scientific direction of **Claude Michel**, and are the work of
the authors listed for each model below.

---

## 2. Rainfall-runoff models

| Model | Time step | Reference |
|---|---|---|
| **GR4J** | daily | Perrin et al. (2003) |
| **GR5J** | daily | Le Moine (2008); see also Pushpalatha et al. (2011) |
| **GR6J** | daily | Pushpalatha et al. (2011) |
| **GR4H** | hourly | Mathevet (2005) |
| **GR5H** | hourly | Ficchì (2017); Ficchì et al. (2019) |
| **GR5H** interception store | hourly | Ficchì et al. (2019) |
| **GR2M** | monthly | Mouelhi et al. (2006a); Mouelhi (2003) |
| **GR1A** | annual | Mouelhi et al. (2006b); Mouelhi (2003) |

**Perrin, C., Michel, C., Andréassian, V.** (2003). Improvement of a parsimonious
model for streamflow simulation. *Journal of Hydrology*, 279(1–4), 275–289.
doi:[10.1016/S0022-1694(03)00225-7](https://doi.org/10.1016/S0022-1694(03)00225-7)

**Le Moine, N.** (2008). *Le bassin versant de surface vu par le souterrain : une
voie d'amélioration des performances et du réalisme des modèles pluie-débit ?*
PhD thesis, Université Pierre et Marie Curie (Paris 6) / Cemagref, Antony.

**Pushpalatha, R., Perrin, C., Le Moine, N., Mathevet, T., Andréassian, V.**
(2011). A downward structural sensitivity analysis of hydrological models to
improve low-flow simulation. *Journal of Hydrology*, 411(1–2), 66–76.
doi:[10.1016/j.jhydrol.2011.09.034](https://doi.org/10.1016/j.jhydrol.2011.09.034)

**Mathevet, T.** (2005). *Quels modèles pluie-débit globaux au pas de temps
horaire ? Développements empiriques et comparaison de modèles sur un large
échantillon de bassins versants.* PhD thesis, ENGREF / Cemagref, Antony.

**Ficchì, A.** (2017). *An adaptive hydrological model for multiple time steps:
diagnostics and improvements based on fluxes consistency.* PhD thesis,
Université Pierre et Marie Curie (Paris 6) / Irstea, Antony.

**Ficchì, A., Perrin, C., Andréassian, V.** (2019). Hydrological modelling at
multiple sub-daily time steps: model improvement via flux-matching.
*Journal of Hydrology*, 575, 1308–1327.
doi:[10.1016/j.jhydrol.2019.05.084](https://doi.org/10.1016/j.jhydrol.2019.05.084)

**Mouelhi, S., Michel, C., Perrin, C., Andréassian, V.** (2006a). Stepwise
development of a two-parameter monthly water balance model.
*Journal of Hydrology*, 318(1–4), 200–214.
doi:[10.1016/j.jhydrol.2005.06.014](https://doi.org/10.1016/j.jhydrol.2005.06.014)

**Mouelhi, S., Michel, C., Perrin, C., Andréassian, V.** (2006b). Linking stream
flow to rainfall at the annual time step: the Manabe bucket model revisited.
*Journal of Hydrology*, 328(1–2), 283–296.
doi:[10.1016/j.jhydrol.2005.12.022](https://doi.org/10.1016/j.jhydrol.2005.12.022)

**Mouelhi, S.** (2003). *Vers une chaîne cohérente de modèles pluie-débit
conceptuels globaux aux pas de temps pluriannuel, annuel, mensuel et journalier.*
PhD thesis, ENGREF / Cemagref, Antony.

---

## 3. Snow: the CemaNeige module

| Component | Reference |
|---|---|
| CemaNeige snow accounting routine | Valéry (2010); Valéry et al. (2014a, 2014b) |
| Linear SWE–SCA hysteresis (`Hyst`) | Riboust et al. (2019) |
| Elevation-band extrapolation, temperature gradients | Valéry (2010) |
| Solid-fraction partition (USACE, Hydrotel) | Valéry (2010); Turcotte et al. (2007) |

**Valéry, A.** (2010). *Modélisation précipitations – débit sous influence
nivale. Élaboration d'un module neige et évaluation sur 380 bassins versants.*
PhD thesis, AgroParisTech / Cemagref, Antony.

**Valéry, A., Andréassian, V., Perrin, C.** (2014a). 'As simple as possible but
not simpler': what is useful in a temperature-based snow-accounting routine?
Part 1 – Comparison of six snow accounting routines on 380 catchments.
*Journal of Hydrology*, 517, 1166–1175.
doi:[10.1016/j.jhydrol.2014.04.059](https://doi.org/10.1016/j.jhydrol.2014.04.059)

**Valéry, A., Andréassian, V., Perrin, C.** (2014b). 'As simple as possible but
not simpler': what is useful in a temperature-based snow-accounting routine?
Part 2 – Sensitivity analysis of the Cemaneige snow accounting routine on 380
catchments. *Journal of Hydrology*, 517, 1176–1187.
doi:[10.1016/j.jhydrol.2014.04.058](https://doi.org/10.1016/j.jhydrol.2014.04.058)

**Riboust, P., Thirel, G., Le Moine, N., Ribstein, P.** (2019). Revisiting a
simple degree-day model for integrating satellite data: implementation of
SWE-SCA hystereses. *Journal of Hydrology and Hydromechanics*, 67(1), 70–81.
doi:[10.2478/johh-2018-0004](https://doi.org/10.2478/johh-2018-0004)

**Turcotte, R., Fortin, L.-G., Fortin, V., Fortin, J.-P., Villeneuve, J.-P.**
(2007). Operational analysis of the spatial distribution and the temporal
evolution of the snowpack water equivalent in southern Québec, Canada.
*Nordic Hydrology*, 38(3), 211–234.
doi:[10.2166/nh.2007.009](https://doi.org/10.2166/nh.2007.009)

The daily temperature-gradient table shipped in
[`src/grsuite/data/gradT_valery2010.csv`](../src/grsuite/data/gradT_valery2010.csv)
is Valéry's (2010), as distributed with airGR.

---

## 4. Potential evapotranspiration

**Oudin, L., Hervieu, F., Michel, C., Perrin, C., Andréassian, V., Anctil, F.,
Loumagne, C.** (2005). Which potential evapotranspiration input for a lumped
rainfall-runoff model? Part 2 – Towards a simple and efficient potential
evapotranspiration model for rainfall-runoff modelling.
*Journal of Hydrology*, 303(1–4), 290–306.
doi:[10.1016/j.jhydrol.2004.08.026](https://doi.org/10.1016/j.jhydrol.2004.08.026)

---

## 5. Error criteria

| Criterion | GRsuite name | Reference |
|---|---|---|
| Nash–Sutcliffe efficiency | `"NSE"` | Nash & Sutcliffe (1970) |
| Kling–Gupta efficiency | `"KGE"` | Gupta et al. (2009) |
| Modified Kling–Gupta efficiency (KGE′) | `"KGE2"` | Kling et al. (2012) |
| Root mean square error | `"RMSE"` | — |
| Box–Cox transformation of discharge | `transfo="boxcox"` | Box & Cox (1964); Santos et al. (2018) |

**Nash, J. E., Sutcliffe, J. V.** (1970). River flow forecasting through
conceptual models part I — a discussion of principles.
*Journal of Hydrology*, 10(3), 282–290.
doi:[10.1016/0022-1694(70)90255-6](https://doi.org/10.1016/0022-1694(70)90255-6)

**Gupta, H. V., Kling, H., Yilmaz, K. K., Martinez, G. F.** (2009). Decomposition
of the mean squared error and NSE performance criteria: implications for
improving hydrological modelling. *Journal of Hydrology*, 377(1–2), 80–91.
doi:[10.1016/j.jhydrol.2009.08.003](https://doi.org/10.1016/j.jhydrol.2009.08.003)

**Kling, H., Fuchs, M., Paulin, M.** (2012). Runoff conditions in the upper
Danube basin under an ensemble of climate change scenarios.
*Journal of Hydrology*, 424–425, 264–277.
doi:[10.1016/j.jhydrol.2012.01.011](https://doi.org/10.1016/j.jhydrol.2012.01.011)

**Box, G. E. P., Cox, D. R.** (1964). An analysis of transformations.
*Journal of the Royal Statistical Society: Series B*, 26(2), 211–243.
doi:[10.1111/j.2517-6161.1964.tb00553.x](https://doi.org/10.1111/j.2517-6161.1964.tb00553.x)

**Santos, L., Thirel, G., Perrin, C.** (2018). Technical note: pitfalls in using
log-transformed flows within the KGE criterion.
*Hydrology and Earth System Sciences*, 22(8), 4583–4591.
doi:[10.5194/hess-22-4583-2018](https://doi.org/10.5194/hess-22-4583-2018)

---

## 6. Calibration and evaluation

| Component | Reference |
|---|---|
| Calibration algorithm (grid screening, then local steepest descent) | Michel (1991) |
| Parameter transformations (bounded search space) | Coron et al. (2017), after the model references above |
| Split-sample test | Klemeš (1986) |

**Michel, C.** (1991). *Hydrologie appliquée aux petits bassins ruraux.*
Hydrology handbook (in French), Cemagref, Antony.

**Klemeš, V.** (1986). Operational testing of hydrological simulation models.
*Hydrological Sciences Journal*, 31(1), 13–24.
doi:[10.1080/02626668609491024](https://doi.org/10.1080/02626668609491024)

---

## 7. Semi-distributed routing

**Lobligeois, F.** (2014). *Mieux connaître la distribution spatiale des pluies
améliore-t-il la modélisation des crues ? Diagnostic sur 181 bassins versants
français.* PhD thesis, AgroParisTech / Irstea, Antony.

**de Lavenne, A., Thirel, G., Andréassian, V., Perrin, C., Ramos, M.-H.** (2016).
Spatial variability of the parameters of a semi-distributed hydrological model.
*Proceedings of the IAHS*, 373, 87–94.
doi:[10.5194/piahs-373-87-2016](https://doi.org/10.5194/piahs-373-87-2016)

---

## 8. Deliberately not implemented

`CreateErrorCrit_GAPX`, the a-priori parameter criterion of

**de Lavenne, A., Andréassian, V., Thirel, G., Ramos, M.-H., Perrin, C.** (2019).
A regularization approach to improve the sequential calibration of a
semidistributed hydrological model. *Water Resources Research*, 55(11),
8821–8839.
doi:[10.1029/2018WR024266](https://doi.org/10.1029/2018WR024266)

---

## 9. Data assimilation

The data assimilation component (`InputsPert`, `run_model_da`,
`Catchment.assimilate`) re-implements airGRdatassim 0.1.4; the EnKF and
particle-filter schemes it realises, and the state-perturbation scheme, are
the work of their authors.

**Piazzi, G., Thirel, G., Perrin, C., Delaigue, O.** (2021). Sequential data
assimilation for streamflow forecasting: assessing the sensitivity to
uncertainties and updated variables of a conceptual hydrological model at
basin scale. *Water Resources Research*, 57, e2020WR028390.
doi:[10.1029/2020WR028390](https://doi.org/10.1029/2020WR028390)

**Salamon, P., Feyen, L.** (2009). Assessing parameter, precipitation, and
predictive uncertainty in a distributed hydrological model using sequential
data assimilation with the particle filter. *Journal of Hydrology*,
376(3–4), 428–442.
doi:[10.1016/j.jhydrol.2009.07.051](https://doi.org/10.1016/j.jhydrol.2009.07.051)

**Piazzi, G., Delaigue, O.** (2025). airGRdatassim: Ensemble-Based Data
Assimilation with GR Hydrological Models. R package version 0.1.4, INRAE,
HYCAR Research Unit, Antony, France.
doi:[10.32614/CRAN.package.airGRdatassim](https://doi.org/10.32614/CRAN.package.airGRdatassim)

---

## 10. Data used for validation

**Delaigue, O., Brigode, P., Andréassian, V., Perrin, C., Etchevers, P.,
Soubeyroux, J.-M., Janet, B., Addor, N.** (2024). CAMELS-FR dataset: a
large-sample hydroclimatic dataset for France to explore hydrological diversity
and support model benchmarking. Recherche Data Gouv.
doi:[10.57745/WH7FJR](https://doi.org/10.57745/WH7FJR) —
meteorological forcing from the SAFRAN/ISBA reanalysis (Météo-France),
streamflow observations from Hydro / Eaufrance.

The three demonstration catchments used by the test suite (L0123001 daily,
L0123002 snow, L0123003 hourly) are the ones distributed with airGR.

---

*Attribution and licensing are in [NOTICE](../NOTICE); how to cite GRsuite itself
is in [CITATION.cff](../CITATION.cff).*
