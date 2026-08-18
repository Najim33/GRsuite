## ============================================================================
## da_instrumented.R - instrumented copies of airGRdatassim 0.1.4 (CRAN)
##
## Provenance: the four functions below (CreateInputsPert, RunModel_DA,
## DA_EnKF, DA_PF) are VERBATIM copies of the R sources of the CRAN package
## airGRdatassim 0.1.4 (GPL-2; authors Gaia Piazzi, Olivier Delaigue;
## (c) INRAE), extracted from
##   https://cran.r-project.org/src/contrib/airGRdatassim_0.1.4.tar.gz
## with ONLY the following modifications, all marked in place:
##
##   1. RNG logging: every rnorm()/runif() call is wrapped in .da_log(), which
##      appends the drawn values to .da_draws$calls (in call order) and
##      returns them unchanged. tools/oracle/export_da_fixtures.R uses this to
##      export the exact random draws so that GRsuite's test suite can replay
##      them (R's Mersenne-Twister + inversion normals cannot be reproduced
##      bit-for-bit by numpy).
##   2. GRSUITE-FIX in DA_EnKF (marked below): CRAN 0.1.4 computes
##      `IndDa <- which(StateEnKF == 1)`. RunModel_DA passes StateEnKF as a
##      character vector (e.g. c("Prod", "Rout"), the documented usage), the
##      comparison is therefore always FALSE, IndDa is empty and the Kalman
##      update never runs: the EnKF is a silent no-op in airGRdatassim 0.1.4
##      (verified empirically: EnsStateA == EnsStateBkg end-to-end). The fix
##      implements the documented behavior (character names select the state
##      variables to update) and keeps the original numeric-indicator branch.
##   3. GRSUITE-ADD in RunModel_DA (marked below): ObsPert is added to the
##      returned list (transposed like QsimEns) so that the EnKF observation
##      perturbations can be validated value by value. CRAN 0.1.4 computes
##      ObsPert but does not return it.
##
## airGRdatassim is free software under GPL-2; these derivative copies carry
## the same license. See https://cran.r-project.org/package=airGRdatassim
## ============================================================================

## ---- RNG draw log ------------------------------------------------------------

## list of list(kind = "rnorm"|"runif", value = numeric) in exact call order
.da_draws <- new.env(parent = emptyenv())
.da_draws$calls <- list()

.da_log <- function(kind, value) {
  .da_draws$calls[[length(.da_draws$calls) + 1L]] <-
    list(kind = kind, value = as.numeric(value))
  value
}

.da_reset_draws <- function() {
  .da_draws$calls <- list()
  invisible(NULL)
}


## =============================================================================
## CreateInputsPert.R - verbatim CRAN 0.1.4 + marked modifications


CreateInputsPert <- function(FUN_MOD, DatesR, Precip = NULL, PotEvap = NULL, TempMean = NULL,
                             ZInputs = NULL, HypsoData = NULL, NLayers = 5, NbMbr = 50, Seed = NULL) {

  # ------ Checks

  FUN_MODList <- c("RunModel_GR4J",
                   "RunModel_GR5J",
                   "RunModel_GR6J")

  FUN_MODSnowList <- c("RunModel_CemaNeigeGR4J",
                       "RunModel_CemaNeigeGR5J",
                       "RunModel_CemaNeigeGR6J")

  FUN_MOD <- match.fun(FUN_MOD)

  if (!any(sapply(c(FUN_MODList, FUN_MODSnowList), function(x) identical(FUN_MOD, match.fun(x))))) {
    stop(sprintf("incorrect 'FUN_MOD' for use in 'CreateInputsPerturb'. Only %s can be used",
                 paste(c(FUN_MODList, FUN_MODSnowList), collapse = ", ")))
  }

  if (!(is.atomic(NbMbr) && is.numeric(NbMbr) && length(NbMbr) == 1 && NbMbr >= 2)) {
    stop("'NbMbr' should be a single vector of a numeric value >= 2")
  }

  Seed0 <- as.numeric(Seed)



  # ------ Settings

  # variables to perturbate
  MeteoNames <- c("Precip", "PotEvap")

  # length of the input vectors
  NbTime <- length(DatesR)

  # model time step [day]
  Dt <- 1

  # temporal decorrelation length for Precip and PotEvap [day]
  Tao  <- c(Precip = 1, PotEvap = 2)
  Alfa <- 1 - (Dt/Tao)

  # fractional error parameter for Precip and PotEvap
  Eps  <- c(Precip = 0.65, PotEvap = 0.65)

  # complementary error function
  Erfc <- function(x) {
    2 * pnorm(q = x * sqrt(2), lower.tail = FALSE)
  }

  # check variables to perturbate
  if (is.null(Precip) & is.null(PotEvap)) {
    stop("'Precip' and 'PotEvap' are both missing. You must provide at least one of them")
  }
  skipVarMeteo <- NULL
  if (is.null(Precip)) {
    skipVarMeteo <- c(skipVarMeteo, "Precip")
    Precip <- rep(0, times = NbTime)
  }
  if (is.null(PotEvap)) {
    skipVarMeteo <- c(skipVarMeteo, "PotEvap")
    PotEvap <- rep(0, times = NbTime)
  }

  # select MeteoNames and Eps in function of the provided variables
  MeteoNames <- MeteoNames[!MeteoNames %in% skipVarMeteo]
  Eps <- Eps[MeteoNames]
  Alfa <- Alfa[MeteoNames]

  # number of variable to perturbate
  NbMeteo <- length(MeteoNames)



  # ------ Generation of the InputsModel object

  if (any(sapply(c(FUN_MODList), function(x) identical(FUN_MOD, match.fun(x)))))  {
    InputsPert <- airGR::CreateInputsModel(FUN_MOD = FUN_MOD, DatesR = DatesR,
                                           Precip = Precip, PotEvap = PotEvap)

  } else if (any(sapply(c(FUN_MODSnowList), function(x) identical(FUN_MOD, match.fun(x))))) {
    InputsPert <- airGR::CreateInputsModel(FUN_MOD = FUN_MOD, DatesR = DatesR,
                                           Precip = Precip, PotEvap = PotEvap,
                                           TempMean = TempMean, ZInputs = ZInputs,
                                           HypsoData = HypsoData, NLayers = NLayers)
  }

  #***************************************************************************************************************
  # Initialisation of the parameters of the first-order autoregressive model

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))

  # time names
  TimeNames <- sprintf("Time_%s", seq_len(NbTime))


  # error time evolution
  S <- array(data = rep(NaN, times = NbMeteo*NbMbr*NbTime),
             dim = c(NbMeteo, NbMbr, NbTime),
             dimnames = list(MeteoNames, MbrNames, TimeNames))
  # uniform random numbers
  U <- array(data = rep(NaN, times = NbMeteo*NbMbr*NbTime),
             dim = c(NbMeteo, NbMbr, NbTime),
             dimnames = list(MeteoNames, MbrNames, TimeNames))

  # multiplicative perturbations of model inputs
  Fi <- array(data = rep(NaN, times = NbMeteo*NbMbr*NbTime),
              dim = c(NbMeteo, NbMbr, NbTime),
              dimnames = list(MeteoNames, MbrNames, TimeNames))



  # ------ Initialisation of meteorological ensembles

  MeteoEns <- sapply(MeteoNames, function(iMeteo) {
    matrix(data = rep(InputsPert[[iMeteo]], each = NbMbr),
           nrow = NbTime, byrow = TRUE,
           dimnames = list(TimeNames, MbrNames))
  }, simplify = "array")

  # same dimensions as Fi
  MeteoEns <- aperm(MeteoEns, perm = 3:1)



  # ------ Perturbation of the time series of model forcings (i.e., precipitation and potential evapotranspiration)

  for (iTime in seq_len(NbTime)) { # olivier, can we do apply here?

    # white noise
    if (!is.null(Seed)) {
      Seed <- Seed0 + iTime
      set.seed(seed = Seed)
      on.exit(set.seed(seed = NULL))
    }

    W <- matrix(data = .da_log("rnorm", rnorm(NbMeteo*NbMbr, mean = 0, sd = 1)), byrow = TRUE,
                nrow = NbMeteo, ncol = NbMbr,
                dimnames = list(MeteoNames, MbrNames))

    if (iTime == 1) {
      # error initialisation
      S0 <- matrix(data = .da_log("runif", runif(NbMeteo*NbMbr, min = 0, max = 1)), byrow = TRUE,
                   nrow = NbMeteo, ncol = NbMbr,
                   dimnames = list(MeteoNames, MbrNames))

      S[, , iTime] <- Alfa * S0 + sqrt(1-Alfa^2) * W

    } else {
      # error time evolution
      S[, , iTime] <- Alfa * S[, , iTime-1] + sqrt(1-Alfa^2) * W
    }

    # generation of uniform random numbers
    U[, , iTime] <- 0.5 * Erfc(S[, , iTime] / sqrt(2))

    # generation of the multiplicative perturbations of model inputs
    Fi[, , iTime] <- (1-Eps) + (2 * U[, , iTime] * Eps)

    MeteoEns[, , iTime] <- MeteoEns[, , iTime]  * Fi[, , iTime]


  } # END FOR time

  # split MeteoEns array into a list of matrix
  MeteoEns <- asplit(MeteoEns, MARGIN = 1)

  # update InputsPert
  for (iName in c(names(MeteoEns), skipVarMeteo)) {
    if (iName %in% MeteoNames) {
      InputsPert[[iName]] <- t(MeteoEns[[iName]])
    } else {
      InputsPert[[iName]] <- NULL
    }
  }
  InputsPert$NbMbr <- NbMbr



  # ------ Class

  class(InputsPert) <- c("InputsPert", class(InputsPert))
  return(InputsPert)
}

## =============================================================================
## RunModel_DA.R - verbatim CRAN 0.1.4 + marked modifications


RunModel_DA <- function(InputsModel, InputsPert = NULL, Qobs = NULL,
                        IndRun,
                        FUN_MOD, Param,
                        DaMethod = c("EnKF", "PF", "none"), NbMbr = NULL,
                        StateEnKF = NULL, StatePert = NULL,
                        Seed = NULL) {

  # ------ Checks

  # FUN_MOD
  TimeUnit <- "daily"

  FUN_MODList <- c("RunModel_GR4J",
                   "RunModel_GR5J",
                   "RunModel_GR6J")

  FUN_MODSnowList <- c("RunModel_CemaNeigeGR4J",
                       "RunModel_CemaNeigeGR5J",
                       "RunModel_CemaNeigeGR6J")

  FUN_MOD <- match.fun(FUN_MOD)
  if (!any(sapply(c(FUN_MODList, FUN_MODSnowList), function(x) identical(FUN_MOD, match.fun(x))))) {
    stop(sprintf("incorrect 'FUN_MOD' for use in 'CreateInputsPerturb'. Only %s can be used",
                 paste(c(FUN_MODList, FUN_MODSnowList), collapse = ", ")))
  }

  if (identical (FUN_MOD, RunModel_GR5J) | identical (FUN_MOD, RunModel_CemaNeigeGR5J)) {
    StateNames <- c("Prod", "Rout", "UH2")
  } else {
    StateNames <- c("Prod", "Rout", "UH1", "UH2")
  }

  # DaMethod
  DaMethod <- match.arg(DaMethod)

  # Seed
  Seed0 <- as.numeric(Seed)

  # StateEnKF & StatePert
  if (DaMethod == "none" && (!is.null(StateEnKF) | !is.null(StatePert))) {
    warning("'StateEnKF' and/or 'StatePert' not taken into account when 'DaMethod' is \"none\"")
  }
  if (DaMethod == "PF" && !is.null(StateEnKF)) {
    warning("'StateEnKF' not taken into account when 'DaMethod' is \"PF\"")
  }
  if (DaMethod == "EnKF" && is.null(StateEnKF)) {
    stop("'StateEnKF' must be defined when 'DaMethod' is \"EnKF\"")
  }
  if (DaMethod != "none") {
    if (!is.null(StateEnKF)) {
      StateEnKF <- match.arg(StateEnKF, choices = StateNames, several.ok = TRUE)
    }
    if (!is.null(StatePert)) {
      StatePert <- match.arg(StatePert, choices = StateNames, several.ok = TRUE)
    }
  }
  if (DaMethod == "EnKF" && any(!StatePert %in% StateEnKF)) {
    stop(sprintf("Perturbation is allowed only for the state variables updated via EnKF (%s). Please check the consistency between 'StatePert' and 'StateEnKF'",
                 sQuote(paste(StateEnKF, collapse = ", "))))
  }
  if (DaMethod != "none" && any(!StatePert %in% StateNames)) {
    warning(StatePert[!StatePert %in% StateNames])
  }

  # InputsModel
  if (!inherits(InputsModel, "InputsModel")) {
    stop("'InputsModel' must of class 'InputsModel'")
  }

  # InputsPert
  if (is.null(InputsPert)) {
    IsMeteo <- FALSE
    if (is.null(NbMbr)) {
      NbMbr <- 50
      message("'InputsPert' and 'NbMbr' not defined: number of ensemble members automatically set to 50")
    }
  } else {
    if (!inherits(InputsPert, "InputsPert")) {
      stop("'InputsPert' must of class 'InputsPert' or NULL")
    } else {
      if (length(InputsPert$DatesR) != length(InputsModel$DatesR)) {
        stop("'InputsPert' elements must have the same length as the 'InputsModel'elements")
      }
      IsMeteo <- TRUE
      ClassInputsModel <- class(InputsModel)
      ClassInputsPert  <- class(InputsPert)[-grep("InputsPert", class(InputsPert))]
      ClassDiffInputsModel <- setdiff(ClassInputsModel, ClassInputsPert)
      ClassDiffInputsPert  <- setdiff(ClassInputsPert, ClassInputsModel)
      if (length(ClassDiffInputsModel) != 0 | length(ClassDiffInputsPert) != 0) {
        msgClassInputs <- "'InputsModel' and 'InputsPert' classes are not consistent:"
        if (length(ClassDiffInputsModel) != 0) {
          msgClassInputs <- sprintf("%s\n\tInputsModel: %s", msgClassInputs, paste(dQuote(ClassDiffInputsModel), collapse = "\t"))
        }
        if (length(ClassDiffInputsPert) != 0) {
          msgClassInputs <- sprintf("%s\n\tInputsPert:  %s", msgClassInputs, paste(dQuote(ClassDiffInputsPert), collapse = "\t"))
        }
        stop(msgClassInputs)
      }
    }

    # NbMbr
    if (is.null(NbMbr)) {
      message(sprintf("'NbMbr' not defined: number of ensemble members automatically set to 'InputsPert$NbMbr' (%i)", InputsPert$NbMbr))
      NbMbr <- InputsPert$NbMbr
    }
    if (!(is.atomic(NbMbr) && is.numeric(NbMbr) && length(NbMbr) == 1 && NbMbr >= 2)) {
      stop("'NbMbr' should be a single vector of a numeric value >= 2")
    } else {
      NbMbr <- as.integer(NbMbr)
      if (IsMeteo) {
        NbMbrMeteo <- ncol(InputsPert[[2]])
        if (NbMbr > NbMbrMeteo) {
          stop(sprintf("cannot take a number of ensemble members (%i) larger than the number available for the perturbed meteorological variables (%i)",
                       NbMbr, NbMbrMeteo))
        }
        if (NbMbr < NbMbrMeteo) {
          warning(sprintf("only %i ensemble members are taken, whereas the number available for the perturbed meteorological variables is equal to %i",
                          NbMbr, NbMbrMeteo))
        }
      }
    }

    # Qobs
    if (is.null(Qobs) || all(is.na(Qobs)) || all(Qobs < 0,  na.rm = TRUE)) {
      DaMethod <- "none"
      warning("'DaMethod' is automatically set to 'none'. All Qobs may be 'NULL', 'NA' or negative")
    } else {
      if (length(Qobs) != length(InputsModel$DatesR)) {
        stop("'Qobs' must have the same length as the 'InputsModel' elements")
      }
      if (any(Qobs < 0,  na.rm = TRUE)) {
        warning("negative value(s) of Qobs are automatically set to 'NA'")
      }
    }

  }


  # ------ Settings

  # data assimilation method used (not open-loop simulation)
  IsDa <- DaMethod != "none"

  NbTime <- length(IndRun)

  NbState <- length(StateNames)

  Qobs[Qobs < 0] <- NaN
  VarThr <- quantile(Qobs, probs = 0.1, na.rm = TRUE)

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))

  # time names
  TimeNames <- sprintf("Time_%s", seq_len(NbTime))

  # InputsModel
  InputsModel <- InputsModel[IndRun]

  # InputsPert
  InputsPert <- InputsPert[IndRun]

  # Qobs
  Qobs <- Qobs[IndRun]



  # ------ Ensemble initializations

  ObsPert <- matrix(data = NA,
                    nrow = NbMbr, ncol = NbTime,
                    dimnames = list(MbrNames,
                                    TimeNames))

  QsimEns <- ObsPert

  IniStatesEns   <- list()
  IniStatesEnsNbTime <- list()

  EnsStateBkg <- array(data = rep(NaN, times = NbState*NbMbr*NbTime),
                       dim = c(NbState, NbMbr, NbTime),
                       dimnames = list(StateNames,
                                       MbrNames,
                                       TimeNames))

  EnsStateA <- EnsStateBkg

  ItAssim <- 0

  # fake RunOptions
  RunOptionsIni <- airGR::CreateRunOptions(FUN_MOD = FUN_MOD,
                                           InputsModel = InputsModel,
                                           IndPeriod_Run = 1L,
                                           warning = FALSE, verbose = FALSE)
  RunOptionsIter <- airGR::CreateRunOptions(FUN_MOD = FUN_MOD,
                                            InputsModel = InputsModel,
                                            IndPeriod_Run = 1L,
                                            IndPeriod_WarmUp = 0L,
                                            IniStates = NULL,
                                            warning = FALSE, verbose = FALSE)



  # ------ Run

  if (IsMeteo) {
    if (is.null(InputsPert$Precip)) {
      InputsPert$Precip <- replicate(n = ncol(InputsPert$PotEvap),
                                     expr = InputsModel$Precip)
      dimnames(InputsPert$Precip) <- dimnames(InputsPert$PotEvap)
    }
    if (is.null(InputsPert$PotEvap)) {
      InputsPert$PotEvap <- replicate(n = ncol(InputsPert$Precip),
                                      expr = InputsModel$PotEvap)
      dimnames(InputsPert$PotEvap) <- dimnames(InputsPert$Precip)
    }
  }

  for (iTime in seq_along(IndRun)) {
    if (!is.null(Seed)) {
      Seed <- Seed0 + iTime
      set.seed(Seed)
      on.exit(set.seed(seed = NULL))
    }
    for (iMbr in seq_len(NbMbr)) {

      if (iTime == 1) { # default (one year by default) warmup

        RunOptionsIni$IndPeriod_Run <- iTime
        OutputsModel <- FUN_MOD(InputsModel = InputsModel,
                                RunOptions = RunOptionsIni, Param = Param)

      } else { # IF iTime > 1

        IniStates <- IniStatesEns[[iMbr]]
        IniStates$Store$Rest <- rep(NA, times = 3)
        IniStates <- unlist(IniStates)
        IniStates[is.na(IniStates)] <- 0
        RunOptionsIter$IniStates <- IniStates
        RunOptionsIter$IniResLevels <- NULL

        # definition of run options
        if (IsMeteo) {
          InputsPertMbr <- InputsPert
          InputsPertMbr$Precip <- InputsPert$Precip[, iMbr]
          InputsPertMbr$PotEvap <- InputsPert$PotEvap[, iMbr]
          RunOptionsIter$IndPeriod_Run <- as.integer(iTime)
          InputsModel <- InputsPertMbr
        } else {
          RunOptionsIter$IndPeriod_Run <- as.integer(iTime)
        } # END IF(IsMeteo)
        OutputsModel <- FUN_MOD(InputsModel = InputsModel,
                                RunOptions = RunOptionsIter, Param = Param)
      } # END IF(t == 1)

      IniStatesEns[[iMbr]] <- OutputsModel$StateEnd
      names(IniStatesEns)[iMbr] <- sprintf("Mbr_%s", iMbr)

      EnsStateBkg["Prod", iMbr, iTime] <- OutputsModel$Prod
      EnsStateBkg["Rout", iMbr, iTime] <- OutputsModel$Rout
      EnsStateBkg["UH2" , iMbr, iTime] <- OutputsModel$StateEnd$UH$UH2[1]
      if ("UH1" %in% StateNames) {
        EnsStateBkg["UH1", iMbr, iTime] <- OutputsModel$StateEnd$UH$UH1[1]
      }

      QsimEns[iMbr, iTime]  <- OutputsModel$Qsim

    } # END FOR particles



    # ------ Assimilation [if an observation is available]

    if (IsDa & is.finite(Qobs[iTime])) {

      ItAssim <- ItAssim + 1

      if (DaMethod == "EnKF") {
        ans <- DA_EnKF(Obs = Qobs[iTime], Qsim = QsimEns[, iTime], EnsState = EnsStateBkg[, , iTime],
                       Param = Param, StateNames = StateNames,
                       StatePert = StatePert,
                       NbMbr = NbMbr,
                       StateEnKF = StateEnKF, VarThr = VarThr)

        for (iMbr in seq_len(NbMbr)) { # olivier, it is possible to write the following 3 loops without loops?
          IniStatesEns[[iMbr]]$Store$Prod <- ans$EnsStateEnkf["Prod", iMbr]
          IniStatesEns[[iMbr]]$Store$Rout <- ans$EnsStateEnkf["Rout", iMbr]
          IniStatesEns[[iMbr]]$UH$UH2[1]  <- ans$EnsStateEnkf["UH2" , iMbr]
          if ("UH1" %in% StateNames) {
            IniStatesEns[[iMbr]]$UH$UH1[1] <- ans$EnsStateEnkf["UH1", iMbr]
          }
        }

        if (iTime < NbTime) {
          IniStatesEnsNbTime[[iTime+1]] <- IniStatesEns
          names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s",iTime+1)
        }

        if (!is.null(StatePert)) {
          for (iMbr in seq_len(NbMbr)) {
            IniStatesEns[[iMbr]]$Store$Prod <- ans$EnsStatePert["Prod", iMbr]
            IniStatesEns[[iMbr]]$Store$Rout <- ans$EnsStatePert["Rout", iMbr]
            IniStatesEns[[iMbr]]$UH$UH2[1]  <- ans$EnsStatePert["UH2" , iMbr]
            if ("UH1" %in% StateNames) {
              IniStatesEns[[iMbr]]$UH$UH1[1] <- ans$EnsStatePert["UH1", iMbr]
            }
          }
        }

        EnsStateA[, , iTime] <- ans$EnsStateEnkf

        if (iTime < NbTime) {
          if (!is.null(StatePert)) {
            EnsStateBkg[, , iTime+1] <- ans$EnsStatePert
          } else {
            EnsStateBkg[, , iTime+1] <- ans$EnsStateEnkf
          }
        }
        ObsPert[, iTime] <- ans$ObsPert

      } else if (DaMethod == "PF") {
        ans <- DA_PF(Obs = Qobs[iTime], Qsim = QsimEns[, iTime], States = IniStatesEns,
                     Param = Param, StateNames = StateNames,
                     NbMbr = NbMbr,
                     StatePert = StatePert, VarThr = VarThr)

        if (!is.null(StatePert)) {
          IniStatesEns <- ans$EnsStatePert
        } else {
          IniStatesEns <- ans$EnsStatePf
        }

        EnsStateA["Prod", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$Store$Prod)
        EnsStateA["Rout", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$Store$Rout)
        EnsStateA["UH2" , , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$UH$UH2[1])
        if ("UH1" %in% StateNames) {
          EnsStateA["UH1", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$UH$UH1[1])
        }

        if (iTime < NbTime) { # olivier?
          IniStatesEnsNbTime[[iTime+1]] <- ans$EnsStatePf
          names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s", iTime+1)
        }
      }

    } else { # IF no assimilation

      if (iTime < NbTime) {
        IniStatesEnsNbTime[[iTime+1]] <- IniStatesEns
        names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s", iTime+1)

        EnsStateBkg[, , iTime+1] <- EnsStateBkg[, , iTime]
      }

      EnsStateA [, , iTime] <- EnsStateBkg[, , iTime]
      if (DaMethod == "EnKF") {
        ObsPert[, iTime] <- rep(Qobs[iTime], times = NbMbr)
      }

    } # END IF assimilation

  } # END FOR time



  # ------ Outputs and class

  res <- list(DatesR = InputsModel$DatesR,
              QsimEns = t(QsimEns),
              EnsStateBkg = aperm(EnsStateBkg),
              EnsStateA = aperm(EnsStateA),
              NbTime = NbTime,
              NbMbr = NbMbr,
              NbState = NbState,
              ## GRSUITE-ADD: also export ObsPert (not part of the CRAN 0.1.4
              ## outputs) so that the EnKF observation perturbations can be
              ## validated value by value; transposed like QsimEns.
              ObsPert = t(ObsPert))
  class(res) <- c("OutputsModelDA", "OutputsModel", DaMethod, TimeUnit)
  return(res)

}

## =============================================================================
## DA_EnKF.R - verbatim CRAN 0.1.4 + marked modifications


DA_EnKF <- function(Obs, Qsim, EnsState,
                   Param, StateNames,
                   NbMbr, StateEnKF = NULL,
                   StatePert = NULL, VarThr) {

  # ------ Settings

  NbState  <- nrow(EnsState)
## GRSUITE-FIX >>> documented behavior for the (documented) character
  ## StateEnKF; CRAN 0.1.4 line was `IndDa <- which(StateEnKF == 1)`, empty
  ## for character input -> silent no-op EnKF (see file header).
  if (is.numeric(StateEnKF)) {
    IndDa    <- which(StateEnKF == 1)
  } else {
    IndDa    <- which(StateNames %in% StateEnKF)
  }
  ## <<< END GRSUITE-FIX
  NbVarDa  <- length(IndDa)
  StateBkg <- EnsState[IndDa, , drop = FALSE]

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))



  # ------ Observation error covariance matrix

  VarObs <- max(VarThr^2, (0.1*Obs)^2)

  Pert    <- .da_log("rnorm", rnorm(NbMbr, mean = 0, sd = sqrt(VarObs)))
  ObsPert <- rep(Obs, times = NbMbr)
  ObsPert <- ObsPert + Pert
  ObsPert[ObsPert < 0] <- 0

  ObsErr <- var(Pert)



  # ------ Innovations

  Innov <- ObsPert - Qsim
  names(Innov) <- MbrNames



  # ------ Kalman Gain

  EnsMeanBkg <- rowMeans(StateBkg)
  EnsMeanQ <- mean(Qsim)

  # evaluation of anomalies
  Anom  <- as.matrix(StateBkg - EnsMeanBkg)
  AnomQ <- t(as.matrix(Qsim - EnsMeanQ))

  # evaluation of Kalman Gain
  BhtMbr  <- matrix(data = NA, nrow = NbVarDa, ncol = NbMbr,
                    dimnames = list(unlist(dimnames(StateBkg)[1], use.names = FALSE),
                                    sprintf("Mbr_%s",seq_len(NbMbr))))
  HbhtMbr <- matrix(data = NA, nrow = 1, ncol = NbMbr,
                    dimnames = list(NULL, MbrNames))

  for (iMbr in seq_len(NbMbr)) {
    BhtMbr[, iMbr] <- Anom[, iMbr] %*% t(AnomQ[iMbr])
    HbhtMbr[iMbr]  <- AnomQ[iMbr] %*% t(AnomQ[iMbr])
  }
  Bht  <- as.matrix((1/(NbMbr-1)) * rowSums(BhtMbr))
  Hbht <- as.matrix((1/(NbMbr-1)) * rowSums(HbhtMbr))

  K <- Bht %*% ((Hbht+ObsErr)^(-1))



  # ------ Analysis

  StateA <- StateBkg + K %*% Innov

  EnsStateEnkf          <- EnsState
  EnsStateEnkf[IndDa, ] <- StateA



  # ------ Constraints on analysis states

  EnsStateEnkf["Prod", EnsStateEnkf["Prod", ] < 0.05*Param[1]] <- 0.05 * Param[1]
  EnsStateEnkf["Rout", EnsStateEnkf["Rout", ] <= 0] <- 1e-3
  EnsStateEnkf["UH2" , EnsStateEnkf["UH2" , ] <  0] <- 1e-3
  if ("UH1" %in% StateNames) {
    EnsStateEnkf["UH1" , EnsStateEnkf["UH1" , ] <  0] <- 1e-3
  }

  EnsStateEnkf["Prod", EnsStateEnkf["Prod", ] > Param[1]] <- Param[1] # if Prod > X1 -> Prod = X1
  EnsStateEnkf["Rout", EnsStateEnkf["Rout", ] > Param[3]] <- Param[3] # if Rout > X3 -> Rout = X3



  # ------ States perturbation

  if (!is.null(StatePert)) {
    IndPert <- as.numeric(StateNames %in% StatePert)
    Sd0 <- apply(EnsStateEnkf, 1, sd)
    SdState <- pmin(3, pmax(1.2, Sd0))
    names(SdState) <- StateNames
    MuState  <- rep(0, times = NbState)
    names(MuState) <- StateNames
    TaoState <- matrix(data = .da_log("rnorm", rnorm(NbState*NbMbr, mean = MuState, sd = SdState)),
                       nrow = NbState, ncol = NbMbr, byrow = TRUE,
                       dimnames = list(StateNames,
                                       MbrNames))

    TaoState[IndPert == 0, ] <- 0

    EnsStatePert <- EnsStateEnkf + TaoState

    # Positive state variables
    EnsStatePert["Prod", EnsStatePert["Prod", ] < 0.05*Param[1]] <- 0.05 * Param[1]
    EnsStatePert["Rout", EnsStatePert["Rout", ] <= 0] <- 1e-3
    EnsStatePert["UH2" , EnsStatePert["UH2" , ] <  0] <- 1e-3
    if ("UH1" %in% StateNames) {
      EnsStatePert["UH1", EnsStatePert["UH1" , ] <  0] <- 1e-3
    }

    EnsStatePert["Prod", EnsStatePert["Prod", ] > Param[1]] <- Param[1] # if Prod > X1 -> Prod = X1
    EnsStatePert["Rout", EnsStatePert["Rout", ] > Param[3]] <- Param[3] # if Rout > X3 -> Rout = X3
  }



  # ------ Outputs

  ans <- list(EnsStateEnkf = EnsStateEnkf, ObsPert = ObsPert)
  if (!is.null(StatePert)) {
    ans$EnsStatePert <- EnsStatePert
  }
  return(ans)
}

## =============================================================================
## DA_PF.R - verbatim CRAN 0.1.4 + marked modifications


DA_PF <- function(Obs, Qsim, States,
                  Param, StateNames,
                  NbMbr,
                  StatePert = NULL, VarThr) {

  # ------ Settings

  NbState <- length(StateNames)

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))

  Weights <- rep(0, times = NbMbr)
  names(Weights) <- MbrNames

  Indices <- rep(NA, times = NbMbr)

  CurrentState     <- list()
  CurrentStatePert <- list()



  # ------ Particles weighting

  VarObs <- max(VarThr^2, (0.1*Obs)^2)

  # evaluation of innovations
  Innov <- Obs - Qsim
  names(Innov) <- MbrNames

  # evaluation of weights
  Weights <- dnorm(Innov, mean = 0, sd = sqrt(VarObs))

  # normalisation of weights
  Weights <- Weights/sum(Weights)

  # if the ensemble is squeezed, NA-valued weights can be generated --> all particles are assigned equal weights
  if (!all(is.finite(Weights))) {
    Weights <- rep(1/NbMbr, times = NbMbr)
  }



  # ------ Particles resampling

  # evaluation of the cumulative density function of weights
  CdfW <- cumsum(Weights)

  A <- CdfW[1]
  B <- min(1, A + (1/(NbMbr+1)))

  Urand0    <- A + ((B - A)*(.da_log("runif", runif(1))))
  StepUrand <- (1-Urand0)/NbMbr

  Urand <- seq(from = Urand0, length.out = NbMbr, by = StepUrand)
  Indices <- findInterval(Urand, vec = CdfW, rightmost.closed = TRUE) + 1 # indices of particles to be resampled

  # if the ensemble is squeezed, almost equivalent and low weight values
  # --> indices of most likely particles can not be identified
  # --> all the particles are assigned equal weights and the resampling is re-evaluated according to the new weights

  if (!all(is.finite(Indices))) {
    Weights <- rep(1/NbMbr, times = NbMbr)
    CdfW    <- cumsum(Weights)

    A <- CdfW[1]
    B <- min(1, A + (1/(NbMbr+1)))

    Urand0    <- A + ((B - A)*(.da_log("runif", runif(1))))
    StepUrand <- (1-Urand0)/NbMbr

    Urand <- seq(from = Urand0, length.out = NbMbr, by = StepUrand)
    Indices <- findInterval(Urand, vec = CdfW, rightmost.closed = TRUE) + 1 # indices of particles to be resampled
  }

  Repeats <- as.data.frame(table(Indices)) # it indicates how many time each selected particle must be replicated



  # ------ State perturbation

  if (!is.null(StatePert)) {
    EnsState <- matrix(data = NA, nrow = NbState, ncol = NbMbr,
                       dimnames = list(StateNames,
                                       MbrNames))

    EnsState["Prod", ] <- sapply(seq_along(States), function(x) States[[x]]$Store$Prod)
    EnsState["Rout", ] <- sapply(seq_along(States), function(x) States[[x]]$Store$Rout)
    EnsState["UH2" , ] <- sapply(seq_along(States), function(x) States[[x]]$UH$UH2[1L])
    if ("UH1" %in% StateNames) {
      EnsState["UH1" , ] <- sapply(seq_along(States), function(x) States[[x]]$UH$UH1[1L])
    }

    MuState <- rep(0, times = NbState)
    names(MuState) <- StateNames

    # evaluation of the variance of state variables
    SdState <- pmin(3,
                     pmax(1.2,
                          apply(EnsState[, as.numeric(as.character(Repeats$Indices)), drop = FALSE], MARGIN = 1, sd, na.rm = TRUE),
                          na.rm = TRUE))
    names(SdState) <- StateNames
  }   # END IF



  # ------ Resampling

  for (iPart in seq_len(nrow(Repeats))) {

    IndexParticle <- as.numeric(as.character(Repeats$Indices[iPart]))
    RepParticle   <- as.numeric(as.character(Repeats$Freq[iPart]))

    # the particle identified by the 'IndexParticle' is replicated 'RepParticle' times
    TempState        <- rep(States[IndexParticle], times = RepParticle)
    names(TempState) <- sprintf("Rep%s_Part%s", seq_len(RepParticle), IndexParticle)

    if (!is.null(StatePert)) { # state perturbation
      IndPert <- as.numeric(StateNames %in% StatePert)
      TempStatePert <- TempState  # if the selected particle is NOT replicated --> it is not perturbed

      if (RepParticle > 1) {      # if the selected particle is replicated --> its replications are perturbed

        StateRep <- matrix(data = NA, nrow = NbState, ncol = RepParticle,
                           dimnames = list(StateNames,
                                           sprintf("Rep_%s", seq_len(RepParticle))))

        StateRep["Prod", ] <- rep(States[[IndexParticle]]$Store$Prod, RepParticle)
        StateRep["Rout", ] <- rep(States[[IndexParticle]]$Store$Rout, RepParticle)
        StateRep["UH2" , ] <- rep(States[[IndexParticle]]$UH$UH2[1L], RepParticle)
        if ("UH1" %in% StateNames) {
          StateRep["UH1", ] <- rep(States[[IndexParticle]]$UH$UH1[1L], RepParticle)
        }

        # noise generation
        NoiseState <- matrix(data = .da_log("rnorm", rnorm(RepParticle*NbState, mean = MuState, sd = SdState)),
                             nrow = NbState, ncol = RepParticle, byrow = FALSE,
                             dimnames = list(StateNames,
                                             sprintf("Rep_%s", seq_len(RepParticle))))

        NoiseState[IndPert == 0, ] <- 0      # null noise for the i-th variable if its uncertainty is not considered
        StateRepPert <- StateRep + NoiseState

        # perturbation constraints
        StateRepPert["Prod", StateRepPert["Prod", ] > Param[1]] <- Param[1]
        StateRepPert["Rout", StateRepPert["Rout", ] > Param[3]] <- Param[3]
        for (iRep in seq_len(RepParticle)) {
          TempStatePert[[iRep]]$Store$Prod <- StateRepPert["Prod", iRep]
          TempStatePert[[iRep]]$Store$Rout <- StateRepPert["Rout", iRep]
          TempStatePert[[iRep]]$UH$UH2[1L] <- StateRepPert["UH2" , iRep]
          if ("UH1" %in% StateNames) {
            TempStatePert[[iRep]]$UH$UH1[1L] <- StateRepPert["UH1", iRep]
          }
        }
      } # END IF RepParticle

      CurrentStatePert <- c(CurrentStatePert, TempStatePert)

    } # END IF
    CurrentState <- c(CurrentState, TempState)

    rm(TempState)
    if (exists("TempStatePert")) {
      rm(TempStatePert)
    }

  } # END FOR repeats

  EnsStatePf <- CurrentState

  if (!is.null(StatePert)) {
    EnsStatePert <- CurrentStatePert
  }



  # ------ Outputs

  ans <- list(EnsStatePf = EnsStatePf)
  if (!is.null(StatePert)) {
    ans$EnsStatePert <- EnsStatePert
  }
  return(ans)
}
