
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at >= dat$param$riskh.start) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < dat$param$prep.start) {
    return(dat)
  }

  # Set attributes ----------------------------------------------------------

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepElig.la <- dat$attr$prepElig.la

  prepStat <- dat$attr$prepStat
  prepStat.la <- dat$attr$prepStat.la

  prepClass <- dat$attr$prepClass
  prepClass.la <- dat$attr$prepClass.la

  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen

  prepTimeLastInj <- dat$attr$prepTimeLastInj

  # Parameters
  prep.replace.mod <- dat$param$prep.replace.mod

  prep.coverage <- dat$param$prep.coverage
  prep.coverage.la <- dat$param$prep.coverage.la

  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method

  prep.adhr.dist <- dat$param$prep.adhr.dist
  prep.adhr.dist.la <- dat$param$prep.adhr.dist.la

  prep.discont.rate <- dat$param$prep.discont.rate

  prep.hadr.int <- dat$param$prep.hadr.int
  prep.ladr.int <- dat$param$prep.ladr.int


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 &
                        status == 0 &
                        prepStat == 0 &
                        prepStat.la == 0 &
                        lnt == at)

  if (prep.replace.mod == "all") {
    idsEligStart.la <- which(active == 1 &
                             status == 0 &
                             prepStat == 0 &
                             prepStat.la == 0 &
                             lnt == at)
  } else if (prep.replace.mod == "curr.oral") {
    idsEligStart.la <- which(active == 1 &
                             status == 0 &
                             prepStat == 1 &
                             prepStat.la == 0 &
                             lnt == at)
  } else if (prep.replace.mod == "curr.oral.ladhr") {
    idsEligStart.la <- which(active == 1 &
                             status == 0 &
                             prepStat == 1 &
                             prepStat.la == 0 &
                             prepClass == 1 &
                             lnt == at)
  }


  # Core eligiblity
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind | ind4 >= twind)

  idsEligStart <- intersect(idsIndic, idsEligStart)
  idsEligStart.la <- intersect(idsIndic, idsEligStart.la)

  prepElig[idsEligStart] <- 1
  prepElig.la[idsEligStart.la] <- 1


  ## Stoppage ------------------------------------------------------------------

  # No indications
  idsNoIndic <- which((ind1 < twind | is.na(ind1)) &
                      (ind2 < twind | is.na(ind2)) &
                      (ind3 < twind | is.na(ind3)) &
                      (ind4 < twind | is.na(ind4)))

  prepElig[idsNoIndic] <- 0
  prepElig.la[idsNoIndic] <- 0

  # Risk reassessment rule
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at &
                             (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random (memoryless) discontinuation
  ## TO ADD: different discontinuation rates by formulation
  idsEligStpRand <- which(active == 1 & (prepStat == 1 | prepStat.la == 1))
  vecStpRand <- rbinom(length(idsEligStpRand), 1, prep.discont.rate)
  idsStpRand <- idsEligStpRand[which(vecStpRand == 1)]

  # Diagnosis
  idsStpDx <- which(active == 1 & (prepStat == 1 | prepStat.la == 1) & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & (prepStat == 1 | prepStat.la == 1))

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)
  idsStp.oral <- intersect(idsStp, which(prepStat == 1))
  idsStp.la <- intersect(idsStp, which(prepStat.la == 1))

  prepStat[idsStp.oral] <- 0
  prepElig[idsStp.oral] <- 0

  prepStat.la[idsStp.la] <- 0
  prepElig.la[idsStp.la] <- 0

  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  # Oral
  prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  nEligSt <- length(idsEligStart)
  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    idsStart <- ssample(idsEligStart, nStart)
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP adherence class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 1:3, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist)
  }


  # Injectable
  prepCov.la <- sum(prepStat.la == 1, na.rm = TRUE)/sum(prepElig.la == 1, na.rm = TRUE)
  prepCov.la <- ifelse(is.nan(prepCov.la), 0, prepCov.la)

  nEligSt.la <- length(idsEligStart.la)
  nStart.la <- max(0, min(nEligSt.la, round((prep.coverage.la - prepCov.la) *
                                        sum(prepElig.la == 1, na.rm = TRUE))))
  idsStart.la <- NULL
  if (nStart.la > 0) {
    idsStart.la <- ssample(idsEligStart.la, nStart.la)
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStartTime[idsStart.la] <- at
    prepLastRisk[idsStart.la] <- at

    # PrEP adherence class
    needPC <- which(is.na(prepClass.la[idsStart.la]))
    prepClass.la[idsStart.la[needPC]] <- sample(x = 1:2, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist.la)
  }


  # Injection Process -------------------------------------------------------

  # Started Today
  start.today <- which(prepStat.la == 1 & prepStartTime == at)
  prepTimeLastInj[start.today] <- at

  last.inj <- at - prepTimeLastInj

  # High Adherence
  idsLA.hadr <- which(prepStat.la == 1 & prepClass.la == 2)
  get.injection.hadr <- intersect(idsLA.hadr, which(last.inj >= prep.hadr.int))

  # Low Adherence
  idsLA.ladr <- which(prepStat.la == 1 & prepClass.la == 1)
  get.injection.ladr <- intersect(idsLA.ladr, which(last.inj >= prep.ladr.int))

  get.injection <- union(get.injection.hadr, get.injection.ladr)

  prepTimeLastInj[get.injection] <- at



  # Drug Level --------------------------------------------------------------

  # attributes
  prepLA.dlevel <- dat$attr$prepLA.dlevel
  prepLA.dlevel.int <- dat$attr$prepLA.dlevel.int

  # parameters
  intcept <- dat$param$prepla.dlevel.int # 4.5
  intcept.err <- dat$param$prepla.dlevel.int.err # 2.5/3
  slope <- dat$param$prepla.dlevel.slope # 25

  # set dlevel.int for newly injected
  ## TODO: constrain as positive
  prepLA.dlevel.int[start.today] <- rnorm(length(start.today), intcept, intcept.err)

  # update dlevel for all active users
  prepLA.dlevel <- prepLA.dlevel.int * 10^(-(1/slope)*last.inj)

  ## TODO: make sure these get set to NA when LA PrEP stops


  ## Output --------------------------------------------------------------------

  # Attributes

  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  dat$attr$prepStat.la <- prepStat.la
  dat$attr$prepClass.la <- prepClass.la

  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  dat$attr$prepTimeLastInj <- prepTimeLastInj

  # Summary stats

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected persons
#'              for purpose of PrEP targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Attributes
  n <- length(dat$attr$active)
  uid <- dat$attr$uid
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.uai.mono)) {
    dat$attr$prep.ind.uai.mono <- rep(NA, n)
    dat$attr$prep.ind.uai.nmain <- rep(NA, n)
    dat$attr$prep.ind.ai.sd <- rep(NA, n)
    dat$attr$prep.ind.sti <- rep(NA, n)
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.uai.mono[part.not.tested.6mo] <- at

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$prep.ind.uai.nmain[uai.nmain] <- at

  ## Condition 3a: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- el2.cond3$p2[discl == TRUE]
  dat$attr$prep.ind.ai.sd[ai.sd] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

  return(dat)
}
