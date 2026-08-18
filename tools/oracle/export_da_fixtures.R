# Generate the airGRdatassim reference files used by the GRsuite test suite.
#
# This script extends the chain of trust described in docs/VALIDATION.md to
# data assimilation: the numbers below are produced by R, not by GRsuite.
# The assimilation itself is run by the instrumented copies of airGRdatassim
# 0.1.4 in tools/oracle/da_instrumented.R - read its header: the copies are
# verbatim CRAN except (1) RNG draw logging and (2) the marked GRSUITE-FIX to
# DA_EnKF that repairs the documented-but-broken character StateEnKF (the
# EnKF is a silent no-op in official airGRdatassim 0.1.4).
#
# Because R's random draws cannot be reproduced bit-for-bit by numpy, the
# exact draws are exported alongside the outputs (files *_draws.csv.gz);
# GRsuite's tests replay them through the same call sequence.
#
# The script also cross-checks the instrumented copies against the official
# airGRdatassim package where both are supposed to agree (CreateInputsPert
# and the PF path) and documents the official EnKF no-op.
#
# Usage:        Rscript tools/oracle/export_da_fixtures.R
# Requirement:  install.packages(c("airGR", "airGRdatassim"))

suppressPackageStartupMessages({
  library(airGR)
  library(airGRdatassim)
})

source("tools/oracle/da_instrumented.R")

OUT <- "tests/data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# gzipped output, 17 significant digits (exact round trip for a double)
wz <- function(df, name) {
  path <- file.path(OUT, paste0(name, ".csv.gz"))
  con <- gzfile(path, "wt")
  on.exit(close(con))
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) vapply(x, function(v)
    if (is.na(v)) "NA" else sprintf("%.17g", v), character(1)))
  write.csv(df, con, row.names = FALSE, quote = FALSE)
  invisible(path)
}

jstr <- function(x) paste(x, collapse = "|")   # join vector fields for meta

# RNG draws logged by da_instrumented.R -> long data.frame (call, kind, value)
draws_df <- function() {
  calls <- .da_draws$calls
  if (length(calls) == 0) {
    return(data.frame(call = integer(0), kind = character(0),
                      value = numeric(0)))
  }
  data.frame(
    call = rep(seq_along(calls), vapply(calls, function(x) length(x$value), 1L)),
    kind = rep(vapply(calls, function(x) x$kind, ""), 
               vapply(calls, function(x) length(x$value), 1L)),
    value = unlist(lapply(calls, function(x) x$value)))
}

# ensemble matrix (NbTime x NbMbr) -> wide data.frame with a Date column
ens_df <- function(dates, mat) {
  df <- as.data.frame(mat)
  names(df) <- sprintf("Mbr_%i", seq_len(ncol(mat)))
  data.frame(Date = format(dates, "%Y-%m-%d"), df, check.names = FALSE)
}

# state cube (NbTime x NbMbr x NbState) -> long data.frame
# NB: t() avant as.vector pour aplatir en ligne-major, coherent avec les
# colonnes Date/Mbr (Date en each, Mbr en times)
state_df <- function(dates, cube, state_names) {
  nb_time <- dim(cube)[1]; nb_mbr <- dim(cube)[2]
  df <- data.frame(Date = rep(format(dates, "%Y-%m-%d"), each = nb_mbr),
                   Mbr = rep(seq_len(nb_mbr), times = nb_time))
  for (s in state_names) df[[s]] <- as.vector(t(cube[, , s]))
  df
}

meta_df <- function(...) {
  args <- list(...)
  args <- lapply(args, function(x) if (is.null(x)) "" else jstr(x))
  as.data.frame(args, stringsAsFactors = FALSE)
}

cat("== Daily basin (L0123001) ==\n")
data(L0123001)
BI <- BasinObs
PARAM <- list(GR4J = c(257.238, 1.012, 88.235, 2.208),
              GR5J = c(245.918, 1.027, 90.017, 2.198, 0.318),
              GR6J = c(250, 0.8, 80, 2.1, 0.2, 30))
FUN <- list(GR4J = RunModel_GR4J, GR5J = RunModel_GR5J, GR6J = RunModel_GR6J)

iRun <- seq(which(format(BI$DatesR, "%Y-%m-%d") == "1994-01-01"),
            which(format(BI$DatesR, "%Y-%m-%d") == "1994-12-31"))
IM <- CreateInputsModel(RunModel_GR4J, BI$DatesR, BI$P, PotEvap = BI$E)

# ---------------------------------------------------------------------------
# CreateInputsPert on a two-year window
# ---------------------------------------------------------------------------
cat("== CreateInputsPert ==\n")
iWin <- seq(which(format(BI$DatesR, "%Y-%m-%d") == "1994-01-01"),
            which(format(BI$DatesR, "%Y-%m-%d") == "1995-12-31"))

.da_reset_draws()
IP <- CreateInputsPert(FUN_MOD = RunModel_GR4J, DatesR = BI$DatesR[iWin],
                       Precip = BI$P[iWin], PotEvap = BI$E[iWin],
                       NbMbr = 8, Seed = 42)
IPo <- airGRdatassim::CreateInputsPert(FUN_MOD = RunModel_GR4J,
                       DatesR = BI$DatesR[iWin], Precip = BI$P[iWin],
                       PotEvap = BI$E[iWin], NbMbr = 8, Seed = 42)
stopifnot(isTRUE(all.equal(IP$Precip, IPo$Precip)),
          isTRUE(all.equal(IP$PotEvap, IPo$PotEvap)))
wz(meta_df(basin = "daily", model = "GR4J", nb_mbr = 8, seed = 42,
           variables = "Precip|PotEvap",
           ind_start = iWin[1], ind_end = tail(iWin, 1)), "da_pert_GR4J_meta")
wz(ens_df(BI$DatesR[iWin], IP$Precip), "da_pert_GR4J_precip")
wz(ens_df(BI$DatesR[iWin], IP$PotEvap), "da_pert_GR4J_potevap")
wz(draws_df(), "da_pert_GR4J_draws")
cat("   perturb GR4J OK (official cross-check passed)\n")

# ---------------------------------------------------------------------------
# RunModel_DA without meteorological perturbation
# ---------------------------------------------------------------------------
run_da <- function(cfg, model, method, nb_mbr, state_enkf, state_pert, seed,
                   inputs_model, fun_mod, qobs, ind_run) {
  .da_reset_draws()
  om <- suppressWarnings(RunModel_DA(
    InputsModel = inputs_model, Qobs = qobs, IndRun = ind_run,
    FUN_MOD = fun_mod, Param = PARAM[[model]], DaMethod = method,
    NbMbr = nb_mbr, StateEnKF = state_enkf, StatePert = state_pert,
    Seed = seed))
  dates <- inputs_model$DatesR[ind_run]
  sn <- if (model == "GR5J") c("Prod", "Rout", "UH2") else
    c("Prod", "Rout", "UH1", "UH2")
  wz(meta_df(basin = "daily", model = model, da_method = method,
             nb_mbr = nb_mbr, seed = seed, state_enkf = state_enkf,
             state_pert = state_pert, param = PARAM[[model]],
             ind_start = ind_run[1], ind_end = tail(ind_run, 1)),
     paste0("da_", cfg, "_meta"))
  wz(ens_df(dates, om$QsimEns), paste0("da_", cfg, "_qsimens"))
  wz(state_df(dates, om$EnsStateBkg, sn), paste0("da_", cfg, "_bkg"))
  wz(state_df(dates, om$EnsStateA, sn), paste0("da_", cfg, "_ana"))
  if (method == "EnKF") {
    wz(ens_df(dates, om$ObsPert), paste0("da_", cfg, "_obspert"))
  }
  wz(draws_df(), paste0("da_", cfg, "_draws"))
  cat(sprintf("   %s OK\n", cfg))
  invisible(om)
}

OMe4 <- run_da("enkf_GR4J", "GR4J", "EnKF", 8,
               c("Prod", "Rout", "UH1", "UH2"), c("Prod", "Rout"), 1,
               IM, RunModel_GR4J, BI$Qmm, iRun)
run_da("enkf_GR5J", "GR5J", "EnKF", 8,
       c("Prod", "Rout", "UH2"), NULL, 2,
       IM, RunModel_GR5J, BI$Qmm, iRun)
run_da("enkf_GR6J", "GR6J", "EnKF", 8,
       c("Prod", "Rout", "UH1", "UH2"), c("UH1", "UH2"), 3,
       IM, RunModel_GR6J, BI$Qmm, iRun)
OMp4 <- run_da("pf_GR4J", "GR4J", "PF", 8,
               NULL, c("Prod", "Rout"), 5,
               IM, RunModel_GR4J, BI$Qmm, iRun)
run_da("none_GR4J", "GR4J", "none", 8,
       NULL, NULL, 7,
       IM, RunModel_GR4J, BI$Qmm, iRun)

# official package cross-checks ------------------------------------------------
# PF path: the instrumented copy must give the very same numbers
OMp4o <- suppressWarnings(airGRdatassim::RunModel_DA(
  InputsModel = IM, Qobs = BI$Qmm, IndRun = iRun, FUN_MOD = RunModel_GR4J,
  Param = PARAM$GR4J, DaMethod = "PF", NbMbr = 8,
  StatePert = c("Prod", "Rout"), Seed = 5))
stopifnot(isTRUE(all.equal(OMp4$QsimEns, OMp4o$QsimEns)),
          isTRUE(all.equal(OMp4$EnsStateA, OMp4o$EnsStateA)))
cat("   PF official cross-check passed\n")

# EnKF path: the official package is a documented no-op (GRSUITE-FIX applies)
OMe4o <- suppressWarnings(airGRdatassim::RunModel_DA(
  InputsModel = IM, Qobs = BI$Qmm, IndRun = iRun, FUN_MOD = RunModel_GR4J,
  Param = PARAM$GR4J, DaMethod = "EnKF", NbMbr = 8,
  StateEnKF = c("Prod", "Rout", "UH1", "UH2"),
  StatePert = c("Prod", "Rout"), Seed = 1))
cat(sprintf(paste0(
  "   official airGRdatassim 0.1.4 EnKF no-op check: max|EnsStateA-EnsStateBkg|",
  " = %.3g (instrumented with fix: %.3g)\n"),
  max(abs(OMe4o$EnsStateA - OMe4o$EnsStateBkg), na.rm = TRUE),
  max(abs(OMe4$EnsStateA - OMe4$EnsStateBkg), na.rm = TRUE)))

# ---------------------------------------------------------------------------
# RunModel_DA with meteorological perturbation (InputsPert), 3-year window
# ---------------------------------------------------------------------------
cat("== RunModel_DA with InputsPert ==\n")
iWin3 <- seq(which(format(BI$DatesR, "%Y-%m-%d") == "1993-01-01"),
             which(format(BI$DatesR, "%Y-%m-%d") == "1995-12-31"))
IM3 <- CreateInputsModel(RunModel_GR4J, BI$DatesR[iWin3], BI$P[iWin3],
                         PotEvap = BI$E[iWin3])
iRun3 <- seq(which(format(BI$DatesR[iWin3], "%Y-%m-%d") == "1994-01-01"),
             which(format(BI$DatesR[iWin3], "%Y-%m-%d") == "1994-12-31"))

.da_reset_draws()
IP3 <- CreateInputsPert(FUN_MOD = RunModel_GR4J, DatesR = BI$DatesR[iWin3],
                        Precip = BI$P[iWin3], PotEvap = BI$E[iWin3],
                        NbMbr = 8, Seed = 44)
wz(ens_df(BI$DatesR[iWin3], IP3$Precip), "da_enkfmet_GR4J_precip")
wz(ens_df(BI$DatesR[iWin3], IP3$PotEvap), "da_enkfmet_GR4J_potevap")
wz(draws_df(), "da_enkfmet_GR4J_pdraws")

.da_reset_draws()
OMm <- suppressWarnings(RunModel_DA(
  InputsModel = IM3, InputsPert = IP3, Qobs = BI$Qmm[iWin3], IndRun = iRun3,
  FUN_MOD = RunModel_GR4J, Param = PARAM$GR4J, DaMethod = "EnKF", NbMbr = 8,
  StateEnKF = c("Prod", "Rout", "UH1", "UH2"),
  StatePert = c("Prod", "Rout"), Seed = 8))
dates3 <- IM3$DatesR[iRun3]
wz(meta_df(basin = "daily", model = "GR4J", da_method = "EnKF",
           nb_mbr = 8, seed = 8,
           state_enkf = "Prod|Rout|UH1|UH2", state_pert = "Prod|Rout",
           param = PARAM$GR4J, pert_seed = 44,
           win_start = iWin3[1], win_end = tail(iWin3, 1),
           ind_start = iRun3[1], ind_end = tail(iRun3, 1)),
   "da_enkfmet_GR4J_meta")
wz(ens_df(dates3, OMm$QsimEns), "da_enkfmet_GR4J_qsimens")
wz(state_df(dates3, OMm$EnsStateBkg, c("Prod", "Rout", "UH1", "UH2")),
   "da_enkfmet_GR4J_bkg")
wz(state_df(dates3, OMm$EnsStateA, c("Prod", "Rout", "UH1", "UH2")),
   "da_enkfmet_GR4J_ana")
wz(ens_df(dates3, OMm$ObsPert), "da_enkfmet_GR4J_obspert")
wz(draws_df(), "da_enkfmet_GR4J_draws")
cat("   enkfmet_GR4J OK\n")

# ---------------------------------------------------------------------------
# Snow basin (L0123002), CemaNeige configurations
# ---------------------------------------------------------------------------
cat("== Snow basin (L0123002) ==\n")
data(L0123002)
BS <- BasinObs
hyp <- BasinInfo$HypsoData
ZIn <- median(hyp)
iSCN <- seq(which(format(BS$DatesR, "%Y-%m-%d") == "1994-01-01"),
            which(format(BS$DatesR, "%Y-%m-%d") == "1994-12-31"))
PARAMCN <- list(CemaNeigeGR4J = c(408.774, 2.646, 131.264, 1.174, 0.962, 2.249),
                CemaNeigeGR6J = c(250, 0.8, 80, 2.1, 0.2, 30, 0.962, 2.249))

IMs <- suppressWarnings(CreateInputsModel(
  RunModel_CemaNeigeGR4J, BS$DatesR, BS$P, PotEvap = BS$E, TempMean = BS$T,
  ZInputs = ZIn, HypsoData = hyp, NLayers = 5, verbose = FALSE))

# CreateInputsPert on the snow model (two-year window)
iWinS <- seq(which(format(BS$DatesR, "%Y-%m-%d") == "1994-01-01"),
             which(format(BS$DatesR, "%Y-%m-%d") == "1995-12-31"))
.da_reset_draws()
IPs <- CreateInputsPert(FUN_MOD = RunModel_CemaNeigeGR4J,
                        DatesR = BS$DatesR[iWinS], Precip = BS$P[iWinS],
                        PotEvap = BS$E[iWinS], TempMean = BS$T[iWinS],
                        ZInputs = ZIn, HypsoData = hyp, NLayers = 5,
                        NbMbr = 8, Seed = 43)
IPso <- airGRdatassim::CreateInputsPert(FUN_MOD = RunModel_CemaNeigeGR4J,
                        DatesR = BS$DatesR[iWinS], Precip = BS$P[iWinS],
                        PotEvap = BS$E[iWinS], TempMean = BS$T[iWinS],
                        ZInputs = ZIn, HypsoData = hyp, NLayers = 5,
                        NbMbr = 8, Seed = 43)
stopifnot(isTRUE(all.equal(IPs$Precip, IPso$Precip)),
          isTRUE(all.equal(IPs$PotEvap, IPso$PotEvap)))
wz(meta_df(basin = "snow", model = "CemaNeigeGR4J", nb_mbr = 8, seed = 43,
           variables = "Precip|PotEvap",
           ind_start = iWinS[1], ind_end = tail(iWinS, 1)),
   "da_pert_CemaNeigeGR4J_meta")
wz(ens_df(BS$DatesR[iWinS], IPs$Precip), "da_pert_CemaNeigeGR4J_precip")
wz(ens_df(BS$DatesR[iWinS], IPs$PotEvap), "da_pert_CemaNeigeGR4J_potevap")
wz(draws_df(), "da_pert_CemaNeigeGR4J_draws")
cat("   perturb CemaNeigeGR4J OK (official cross-check passed)\n")

run_da_snow <- function(cfg, model, fun_mod, method, nb_mbr, state_enkf,
                        state_pert, seed) {
  .da_reset_draws()
  om <- suppressWarnings(RunModel_DA(
    InputsModel = IMs, Qobs = BS$Qmm, IndRun = iSCN, FUN_MOD = fun_mod,
    Param = PARAMCN[[model]], DaMethod = method, NbMbr = nb_mbr,
    StateEnKF = state_enkf, StatePert = state_pert, Seed = seed))
  dates <- BS$DatesR[iSCN]
  sn <- c("Prod", "Rout", "UH1", "UH2")
  wz(meta_df(basin = "snow", model = model, da_method = method,
             nb_mbr = nb_mbr, seed = seed, state_enkf = state_enkf,
             state_pert = state_pert, param = PARAMCN[[model]],
             ind_start = iSCN[1], ind_end = tail(iSCN, 1)),
     paste0("da_", cfg, "_meta"))
  wz(ens_df(dates, om$QsimEns), paste0("da_", cfg, "_qsimens"))
  wz(state_df(dates, om$EnsStateBkg, sn), paste0("da_", cfg, "_bkg"))
  wz(state_df(dates, om$EnsStateA, sn), paste0("da_", cfg, "_ana"))
  if (method == "EnKF") {
    wz(ens_df(dates, om$ObsPert), paste0("da_", cfg, "_obspert"))
  }
  wz(draws_df(), paste0("da_", cfg, "_draws"))
  cat(sprintf("   %s OK\n", cfg))
  invisible(om)
}

run_da_snow("enkf_CemaNeigeGR4J", "CemaNeigeGR4J", RunModel_CemaNeigeGR4J,
            "EnKF", 8, c("Prod", "Rout", "UH1", "UH2"), c("Prod", "Rout"), 4)
run_da_snow("pf_CemaNeigeGR6J", "CemaNeigeGR6J", RunModel_CemaNeigeGR6J,
            "PF", 8, NULL, c("Prod", "Rout", "UH1", "UH2"), 6)

sz <- sum(file.info(list.files(OUT, pattern = "^da_", full.names = TRUE))$size)
cat(sprintf("\nWrote %i DA reference files to %s (%.1f MB)\n",
            length(list.files(OUT, pattern = "^da_")), OUT, sz / 1048576))
cat("Produced by ", R.version.string, " with airGR ",
    as.character(packageVersion("airGR")), " and airGRdatassim ",
    as.character(packageVersion("airGRdatassim")), "\n", sep = "")
