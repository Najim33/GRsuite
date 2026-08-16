# Generate the airGR reference files used by the GRsuite test suite.
#
# This script is the first link in the chain of trust described in
# docs/VALIDATION.md section 1.1: it is airGR, not GRsuite, that produces every
# number under tests/data/. Re-run it and compare the decompressed content of
# the files with the committed ones - they must be identical.
#
# Output is gzipped and deliberately compact: it has to fit in the repository
# while covering the whole suite.
#
# Usage:        Rscript tools/oracle/export_fixtures.R
# Requirement:  install.packages("airGR")

suppressPackageStartupMessages(library(airGR))

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

drop_meta <- function(om) {
  as.data.frame(om[setdiff(names(om),
    c("DatesR", "StateEnd", "RunOptions", "CemaNeigeLayers"))])
}

cat("== Daily models ==\n")
data(L0123001)
BI <- BasinObs
iRun <- seq(which(format(BI$DatesR, "%Y-%m-%d") == "1994-01-01"),
            which(format(BI$DatesR, "%Y-%m-%d") == "1999-12-31"))
wz(data.frame(Date = format(BI$DatesR, "%Y-%m-%d"), P = BI$P, E = BI$E,
              T = BI$T, Qmm = BI$Qmm), "basin_daily")
wz(data.frame(ind_start = iRun[1], ind_end = tail(iRun, 1)), "idx_daily")

IM <- CreateInputsModel(RunModel_GR4J, BI$DatesR, BI$P, PotEvap = BI$E)
daily <- list(
  GR4J = list(f = RunModel_GR4J, p = c(257.238, 1.012, 88.235, 2.208)),
  GR5J = list(f = RunModel_GR5J, p = c(245.918, 1.027, 90.017, 2.198, 0.318)),
  GR6J = list(f = RunModel_GR6J, p = c(250, 0.8, 80, 2.1, 0.2, 30)))
for (m in names(daily)) {
  RO <- suppressWarnings(CreateRunOptions(daily[[m]]$f, IM, IndPeriod_Run = iRun,
                                          warnings = FALSE))
  OM <- RunModel(IM, RO, daily[[m]]$p, FUN_MOD = daily[[m]]$f)
  wz(drop_meta(OM), paste0("sim_", m))
  ROs <- suppressWarnings(CreateRunOptions(daily[[m]]$f, IM, IndPeriod_Run = iRun,
                                           Outputs_Sim = "StateEnd", warnings = FALSE))
  OMs <- RunModel(IM, ROs, daily[[m]]$p, FUN_MOD = daily[[m]]$f)
  wz(data.frame(state = unlist(OMs$StateEnd)), paste0("state_", m))
  cat("  ", m, "OK\n")
}

cat("== GR2M and GR1A ==\n")
BM <- SeriesAggreg(BI[, c("DatesR", "P", "E", "Qmm")], Format = "%Y%m",
                   ConvertFun = c("sum", "sum", "sum"))
wz(data.frame(Date = format(BM$DatesR, "%Y-%m-%d"), P = BM$P, E = BM$E, Qmm = BM$Qmm),
   "basin_monthly")
IMm <- CreateInputsModel(RunModel_GR2M, BM$DatesR, BM$P, PotEvap = BM$E)
iM <- seq(which(format(BM$DatesR, "%Y-%m") == "1994-01"),
          which(format(BM$DatesR, "%Y-%m") == "1999-12"))
ROm <- suppressWarnings(CreateRunOptions(RunModel_GR2M, IMm, IndPeriod_Run = iM,
                                         warnings = FALSE))
wz(drop_meta(RunModel(IMm, ROm, c(265.072, 1.007), FUN_MOD = RunModel_GR2M)), "sim_GR2M")
wz(data.frame(ind_start = iM[1], ind_end = tail(iM, 1)), "idx_monthly")

BY <- SeriesAggreg(BI[, c("DatesR", "P", "E", "Qmm")], Format = "%Y",
                   ConvertFun = c("sum", "sum", "sum"))
wz(data.frame(Date = format(BY$DatesR, "%Y-%m-%d"), P = BY$P, E = BY$E, Qmm = BY$Qmm),
   "basin_yearly")
IMy <- CreateInputsModel(RunModel_GR1A, BY$DatesR, BY$P, PotEvap = BY$E)
iY <- seq(which(format(BY$DatesR, "%Y") == "1994"),
          which(format(BY$DatesR, "%Y") == "1999"))
ROy <- suppressWarnings(CreateRunOptions(RunModel_GR1A, IMy, IndPeriod_Run = iY,
                                         warnings = FALSE))
wz(drop_meta(RunModel(IMy, ROy, 0.840, FUN_MOD = RunModel_GR1A)), "sim_GR1A")
wz(data.frame(ind_start = iY[1], ind_end = tail(iY, 1)), "idx_yearly")

# extra aggregations (means, hydrological year)
wz(data.frame(Date = format(SeriesAggreg(BI[, c("DatesR", "P", "E", "Qmm")],
      Format = "%Y%m", ConvertFun = c("mean", "mean", "mean"))$DatesR, "%Y-%m-%d"),
   SeriesAggreg(BI[, c("DatesR", "P", "E", "Qmm")], Format = "%Y%m",
      ConvertFun = c("mean", "mean", "mean"))[, -1]), "aggreg_monthly_mean")
BYs <- SeriesAggreg(BI[, c("DatesR", "P", "E", "Qmm")], Format = "%Y",
                    ConvertFun = c("sum", "sum", "sum"), YearFirstMonth = 9)
wz(data.frame(Date = format(BYs$DatesR, "%Y-%m-%d"), P = BYs$P, E = BYs$E, Qmm = BYs$Qmm),
   "aggreg_yearly_sept")

cat("== Oudin potential evapotranspiration ==\n")
JD <- as.POSIXlt(BI$DatesR)$yday + 1
wz(data.frame(JD = JD, Temp = BI$T,
              PE = PE_Oudin(JD = JD, Temp = BI$T, Lat = 0.8, LatUnit = "rad")),
   "pe_oudin")

cat("== CemaNeige (5 elevation bands) ==\n")
data(L0123002)
BS <- BasinObs
hyp <- BasinInfo$HypsoData
ZIn <- median(hyp)
iCN <- seq(which(format(BS$DatesR, "%Y-%m-%d") == "1994-01-01"),
           which(format(BS$DatesR, "%Y-%m-%d") == "1999-12-31"))
wz(data.frame(Date = format(BS$DatesR, "%Y-%m-%d"), P = BS$P, E = BS$E,
              T = BS$T, Qmm = BS$Qmm), "basin_snow")
wz(data.frame(hypso = hyp), "hypso_snow")
wz(data.frame(ind_start = iCN[1], ind_end = tail(iCN, 1), zinputs = ZIn), "idx_snow")

cn <- list(
  list(n = "CemaNeigeGR4J", f = RunModel_CemaNeigeGR4J, h = FALSE,
       p = c(408.774, 2.646, 131.264, 1.174, 0.962, 2.249), layers = TRUE),
  list(n = "CemaNeigeGR5J", f = RunModel_CemaNeigeGR5J, h = FALSE,
       p = c(245.918, 1.027, 90.017, 2.198, 0.318, 0.962, 2.249), layers = FALSE),
  list(n = "CemaNeigeGR6J", f = RunModel_CemaNeigeGR6J, h = FALSE,
       p = c(250, 0.8, 80, 2.1, 0.2, 30, 0.962, 2.249), layers = FALSE),
  list(n = "CemaNeigeGR4J_Hyst", f = RunModel_CemaNeigeGR4J, h = TRUE,
       p = c(408.774, 2.646, 131.264, 1.174, 0.962, 2.249, 80, 0.4), layers = TRUE))

for (m in cn) {
  Inp <- suppressWarnings(CreateInputsModel(m$f, BS$DatesR, BS$P, PotEvap = BS$E,
    TempMean = BS$T, ZInputs = ZIn, HypsoData = hyp, NLayers = 5, verbose = FALSE))
  RO <- suppressWarnings(CreateRunOptions(m$f, Inp, IndPeriod_Run = iCN,
                                          IsHyst = m$h, warnings = FALSE))
  OM <- RunModel(Inp, RO, m$p, FUN_MOD = m$f)
  wz(drop_meta(OM), paste0("sim_", m$n))
  wz(data.frame(masp = RO$MeanAnSolidPrecip[1]), paste0("masp_", m$n))
  if (m$layers) {
    for (i in 1:5) {
      wz(as.data.frame(OM$CemaNeigeLayers[[i]]),
         sprintf("sim_%s_layer%02i", m$n, i))
    }
  }
  cat("  ", m$n, "OK\n")
}
# input layers, to check the elevation extrapolation on its own
InpA <- suppressWarnings(CreateInputsModel(RunModel_CemaNeigeGR4J, BS$DatesR, BS$P,
  PotEvap = BS$E, TempMean = BS$T, ZInputs = ZIn, HypsoData = hyp, NLayers = 5,
  verbose = FALSE))
for (i in 1:5) {
  wz(data.frame(P = InpA$LayerPrecip[[i]], T = InpA$LayerTemp[[i]],
                FS = InpA$LayerFracSolidPrecip[[i]]),
     sprintf("inputs_layer%02i", i))
}

cat("== Hourly models ==\n")
data(L0123003)
BH <- BasinObs
iH <- seq(which(format(BH$DatesR, "%Y-%m-%d %H") == "2005-01-01 00"),
          which(format(BH$DatesR, "%Y-%m-%d %H") == "2005-12-31 23"))
# keep only the useful window (warm-up + simulation) to save space
keepH <- seq(max(1, iH[1] - 24 * 400), tail(iH, 1))
wz(data.frame(Date = format(BH$DatesR[keepH], "%Y-%m-%d %H:%M:%S"),
              P = BH$P[keepH], E = BH$E[keepH], Qmm = BH$Qmm[keepH]), "basin_hourly")
wz(data.frame(ind_start = iH[1] - keepH[1] + 1, ind_end = tail(iH, 1) - keepH[1] + 1),
   "idx_hourly")

IMh <- CreateInputsModel(RunModel_GR4H, BH$DatesR[keepH], BH$P[keepH],
                         PotEvap = BH$E[keepH])
iHl <- seq(iH[1] - keepH[1] + 1, tail(iH, 1) - keepH[1] + 1)
pH <- c(756.930, -0.773, 138.638, 5.247)
ROh <- suppressWarnings(CreateRunOptions(RunModel_GR4H, IMh, IndPeriod_Run = iHl,
                                         warnings = FALSE))
wz(drop_meta(RunModel(IMh, ROh, pH, FUN_MOD = RunModel_GR4H)), "sim_GR4H")

IMh5 <- CreateInputsModel(RunModel_GR5H, BH$DatesR[keepH], BH$P[keepH],
                          PotEvap = BH$E[keepH])
pH5 <- c(756.930, -0.773, 138.638, 5.247, 0.400)
ROh5 <- suppressWarnings(CreateRunOptions(RunModel_GR5H, IMh5, IndPeriod_Run = iHl,
                                          warnings = FALSE))
wz(drop_meta(RunModel(IMh5, ROh5, pH5, FUN_MOD = RunModel_GR5H)), "sim_GR5H")

IMAX <- 0.7
ROh5i <- suppressWarnings(CreateRunOptions(RunModel_GR5H, IMh5, IndPeriod_Run = iHl,
                                           Imax = IMAX, warnings = FALSE))
wz(drop_meta(RunModel(IMh5, ROh5i, pH5, FUN_MOD = RunModel_GR5H)), "sim_GR5H_interception")
wz(data.frame(imax = IMAX), "imax_value")
cat("   GR4H, GR5H, GR5H+interception OK\n")

cat("== Error criteria and transformations ==\n")
RO4 <- suppressWarnings(CreateRunOptions(RunModel_GR4J, IM, IndPeriod_Run = iRun,
                                         warnings = FALSE))
OM4 <- RunModel(IM, RO4, c(257.238, 1.012, 88.235, 2.208), FUN_MOD = RunModel_GR4J)
crits <- list(NSE = ErrorCrit_NSE, KGE = ErrorCrit_KGE, KGE2 = ErrorCrit_KGE2,
              RMSE = ErrorCrit_RMSE)
transfos <- c("", "sqrt", "log", "inv", "sort", "boxcox", "^0.5", "^2")
rows <- NULL
for (cr in names(crits)) for (tr in transfos) {
  eps <- if (tr %in% c("log", "inv")) 0.01 else NULL
  IC <- suppressWarnings(CreateInputsCrit(crits[[cr]], IM, RO4, Obs = BI$Qmm[iRun],
    transfo = tr, epsilon = eps, warnings = FALSE))
  OC <- suppressWarnings(crits[[cr]](IC, OM4, verbose = FALSE, warnings = FALSE))
  rows <- rbind(rows, data.frame(crit = cr, transfo = tr,
    epsilon = ifelse(is.null(eps), NA, eps),
    value = sprintf("%.17g", OC$CritValue), stringsAsFactors = FALSE))
}
ICc <- suppressWarnings(CreateInputsCrit(list(ErrorCrit_NSE, ErrorCrit_NSE), IM, RO4,
  Obs = list(BI$Qmm[iRun], BI$Qmm[iRun]), VarObs = list("Q", "Q"),
  transfo = list("", "log"), epsilon = list(NULL, 0.01),
  Weights = list(0.6, 0.4), warnings = FALSE))
rows <- rbind(rows, data.frame(crit = "Composite_NSE", transfo = "|log", epsilon = 0.01,
  value = sprintf("%.17g", suppressWarnings(ErrorCrit(ICc, OM4, verbose = FALSE))$CritValue),
  stringsAsFactors = FALSE))
con <- gzfile(file.path(OUT, "error_crits.csv.gz"), "wt")
write.csv(rows, con, row.names = FALSE); close(con)
wz(data.frame(Date = format(BI$DatesR, "%Y-%m-%d"), Qmm = BI$Qmm), "qobs_daily")

tp <- list(list("GR4J", TransfoParam_GR4J, 4), list("GR5J", TransfoParam_GR5J, 5),
           list("GR6J", TransfoParam_GR6J, 6), list("GR4H", TransfoParam_GR4H, 4),
           list("GR5H", TransfoParam_GR5H, 5), list("GR2M", TransfoParam_GR2M, 2),
           list("GR1A", TransfoParam_GR1A, 1),
           list("CemaNeige", TransfoParam_CemaNeige, 2),
           list("CemaNeigeHyst", TransfoParam_CemaNeigeHyst, 4))
set.seed(42)
tr <- NULL
for (d in tp) for (k in 1:5) {
  pt <- runif(d[[3]], -9.99, 9.99)
  pr <- d[[2]](matrix(pt, nrow = 1), "TR")
  tr <- rbind(tr, data.frame(model = d[[1]], k = k,
    paramT = paste(sprintf("%.17g", pt), collapse = "|"),
    paramR = paste(sprintf("%.17g", as.vector(pr)), collapse = "|"),
    paramT2 = paste(sprintf("%.17g", as.vector(d[[2]](pr, "RT"))), collapse = "|"),
    stringsAsFactors = FALSE))
}
con <- gzfile(file.path(OUT, "transfo_param.csv.gz"), "wt")
write.csv(tr, con, row.names = FALSE); close(con)

cat("== Calibration, Michel algorithm ==\n")
res <- NULL
for (m in names(daily)) {
  RO <- suppressWarnings(CreateRunOptions(daily[[m]]$f, IM, IndPeriod_Run = iRun,
                                          warnings = FALSE))
  for (cr in c("NSE", "KGE")) for (tf in c("", "log")) {
    if (cr == "KGE" && tf != "") next
    IC <- suppressWarnings(CreateInputsCrit(crits[[cr]], IM, RO, Obs = BI$Qmm[iRun],
      transfo = tf, epsilon = if (tf == "log") 0.01 else NULL, warnings = FALSE))
    CO <- suppressWarnings(CreateCalibOptions(daily[[m]]$f))
    OC <- suppressWarnings(Calibration_Michel(IM, RO, IC, CO,
      FUN_MOD = daily[[m]]$f, verbose = FALSE))
    res <- rbind(res, data.frame(model = m, crit = cr, transfo = tf,
      crit_final = sprintf("%.17g", OC$CritFinal), n_iter = OC$NIter,
      param = paste(sprintf("%.17g", OC$ParamFinalR), collapse = "|"),
      stringsAsFactors = FALSE))
    cat(sprintf("   %s %s %-4s crit=%8.5f\n", m, cr, ifelse(tf == "", "-", tf),
                OC$CritFinal))
  }
}
con <- gzfile(file.path(OUT, "calib_reference.csv.gz"), "wt")
write.csv(res, con, row.names = FALSE); close(con)

cat("== Semi-distributed routing ==\n")
areas <- c(180, 240, 360); lens <- c(30, 55)
Qup <- cbind(BI$Qmm * areas[1] * 1e3 * 0.9, BI$Qmm * areas[2] * 1e3 * 1.1)
ISD <- suppressWarnings(CreateInputsModel(RunModel_GR4J, BI$DatesR, BI$P,
  PotEvap = BI$E, Qupstream = Qup, LengthHydro = lens, BasinAreas = areas,
  QupstrUnit = "m3", verbose = FALSE))
ROsd <- suppressWarnings(CreateRunOptions(RunModel_GR4J, ISD, IndPeriod_Run = iRun,
                                          warnings = FALSE))
OMsd <- suppressWarnings(RunModel(ISD, ROsd, c(1.2, 257.238, 1.012, 88.235, 2.208),
                                  FUN_MOD = RunModel_GR4J))
wz(data.frame(Qsim = OMsd$Qsim, Qsim_m3 = OMsd$Qsim_m3, QsimDown = OMsd$QsimDown),
   "sim_SD_lag")
wz(data.frame(Qup1 = Qup[, 1], Qup2 = Qup[, 2]), "sd_qupstream")
wz(data.frame(speed = 1.2, length1 = lens[1], length2 = lens[2],
              area1 = areas[1], area2 = areas[2], area3 = areas[3]), "sd_config")

sz <- sum(file.info(list.files(OUT, full.names = TRUE))$size)
cat(sprintf("\nWrote %i reference files to %s (%.1f MB)\n",
            length(list.files(OUT)), OUT, sz / 1048576))
cat("Produced by ", R.version.string, " with airGR ",
    as.character(packageVersion("airGR")), "\n", sep = "")
