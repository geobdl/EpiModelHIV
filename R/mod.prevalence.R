
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  ## Variables

  # Attributes

  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepElig <- dat$attr$prepElig
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt


  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA
    dat$epi$num.B <- rNA
    dat$epi$num.W <- rNA
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    dat$epi$i.num.B <- rNA
    dat$epi$i.num.W <- rNA
    dat$epi$i.prev <- rNA
    dat$epi$i.prev.B <- rNA
    dat$epi$i.prev.W <- rNA
    dat$epi$incid <- rNA
    dat$epi$ir100 <- rNA

    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA

    dat$epi$prev.rgc <- rNA
    dat$epi$prev.ugc <- rNA
    dat$epi$prev.gc <- rNA
    dat$epi$prev.gc.sympt <- rNA
    dat$epi$prev.gc.dual <- rNA

    dat$epi$prev.rct <- rNA
    dat$epi$prev.uct <- rNA
    dat$epi$prev.ct <- rNA
    dat$epi$prev.ct.sympt <- rNA
    dat$epi$prev.ct.dual <- rNA

    dat$epi$prev.rgcct <- rNA
    dat$epi$prev.ugcct <- rNA

    dat$epi$incid.rgc <- rNA
    dat$epi$incid.ugc <- rNA
    dat$epi$incid.gc <- rNA
    dat$epi$incid.rct <- rNA
    dat$epi$incid.uct <- rNA
    dat$epi$incid.ct <- rNA

    dat$epi$ir100.rgc <- rNA
    dat$epi$ir100.ugc <- rNA
    dat$epi$ir100.gc <- rNA
    dat$epi$ir100.rct <- rNA
    dat$epi$ir100.uct <- rNA
    dat$epi$ir100.ct <- rNA

    dat$epi$ir100.sti <- rNA
    dat$epi$incid.gcct.prep <- rNA

    dat$epi$recov.rgc <- rNA
    dat$epi$recov.ugc <- rNA
    dat$epi$recov.rct <- rNA
    dat$epi$recov.uct <- rNA

    dat$epi$trans.main <- rNA
    dat$epi$trans.casl <- rNA
    dat$epi$trans.inst <- rNA

    dat$epi$txGC <- rNA
    dat$epi$txCT <- rNA
  }

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]
  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, na.rm = TRUE)) * 5200

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) & status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
    sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] / sum(prepStat == 1, na.rm = TRUE)
  }

  dat$epi$prev.rgc[at] <- sum(rGC == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ugc[at] <- sum(uGC == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.gc.sympt[at] <- sum((rGC.sympt == 1 | uGC.sympt == 1)) / dat$epi$num[at]
  dat$epi$prev.gc.dual[at] <- sum((rGC == 1 & uGC == 1), na.rm = TRUE) / dat$epi$num[at]

  dat$epi$prev.rct[at] <- sum(rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.uct[at] <- sum(uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct.sympt[at] <- sum((rCT.sympt == 1 | uCT.sympt == 1)) / dat$epi$num[at]
  dat$epi$prev.ct.dual[at] <- sum((rCT == 1 & uCT == 1), na.rm = TRUE) / dat$epi$num[at]

  dat$epi$prev.rgcct[at] <- sum(rGC == 1 | rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ugcct[at] <- sum(uGC == 1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at]

  dat$epi$ir100.rgc[at] <- (dat$epi$incid.rgc[at] / sum(rGC == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.ugc[at] <- (dat$epi$incid.ugc[at] / sum(uGC == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.gc[at] <- (dat$epi$incid.gc[at] /
                             (sum(rGC == 0, na.rm = TRUE) +
                                sum(uGC == 0, na.rm = TRUE))) * 5200

  dat$epi$ir100.rct[at] <- (dat$epi$incid.rct[at] / sum(rCT == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.uct[at] <- (dat$epi$incid.uct[at] / sum(uCT == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.ct[at] <- (dat$epi$incid.ct[at] /
                             (sum(rCT == 0, na.rm = TRUE) +
                                sum(uCT == 0, na.rm = TRUE))) * 5200

  dat$epi$prev.sti[at] <- sum(rGC == 1 | uGC == 1 |
                                rCT ==1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$ir100.sti[at] <- ((dat$epi$incid.ct[at] + dat$epi$incid.gc[at]) /
                              (sum(rGC == 0, na.rm = TRUE) +
                                 sum(uGC == 0, na.rm = TRUE) +
                                 sum(rCT == 0, na.rm = TRUE) +
                                 sum(uCT == 0, na.rm = TRUE))) * 5200

  dat$epi$ir100.sti.prep[at] <- (dat$epi$incid.gcct.prep[at] /
                                  (sum(rGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(rCT == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(uCT == 0 & prepStat == 1, na.rm = TRUE))) * 5200

  return(dat)
}


#' @export
#' @rdname prevalence_msm
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
