
## MSM Het Functions

#' @export
param_msmhet <- function(nwstats,
                      race.method = 1,
                      last.neg.test.B.int = 301,
                      last.neg.test.W.int = 315,
                      mean.test.B.int = 301,
                      mean.test.W.int = 315,
                      testing.pattern = "memoryless",
                      test.window.int = 21,

                      tt.traj.B.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.prob = c(0.052, 0.000, 0.331, 0.617),

                      tx.init.B.prob = 0.092,
                      tx.init.W.prob = 0.127,
                      tx.halt.B.prob = 0.0102,
                      tx.halt.W.prob = 0.0071,
                      tx.reinit.B.prob = 0.00066,
                      tx.reinit.W.prob = 0.00291,

                      max.time.off.tx.full.int = 520 * 7,
                      max.time.on.tx.part.int = 52 * 15 * 7,
                      max.time.off.tx.part.int = 520 * 7,
                      vl.acute.rise.int = 45,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 45,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520 * 7,
                      vl.aids.int = 52 * 2 * 7,
                      vl.fatal = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      full.supp.down.slope = 0.25,
                      full.supp.up.slope = 0.25,
                      part.supp.down.slope = 0.25,
                      part.supp.up.slope = 0.25,

                      b.B.rate = 1e-3 / 7,
                      b.W.rate = 1e-3 / 7,
                      birth.age = 18,
                      b.method = "fixed",

                      URAI.prob = 0.0082 * 1.09,
                      UIAI.prob = 0.0031 * 1.09,
                      acute.rr = 6,
                      circ.rr = 0.4,
                      condom.rr = 0.295,

                      disc.outset.main.B.prob = 0.685,
                      disc.outset.main.W.prob = 0.889,
                      disc.at.diag.main.B.prob = 1,
                      disc.at.diag.main.W.prob = 1,
                      disc.post.diag.main.B.prob = 0,
                      disc.post.diag.main.W.prob = 0,
                      disc.outset.pers.B.prob = 0.527,
                      disc.outset.pers.W.prob = 0.828,
                      disc.at.diag.pers.B.prob = 1,
                      disc.at.diag.pers.W.prob = 1,
                      disc.post.diag.pers.B.prob = 0,
                      disc.post.diag.pers.W.prob = 0,
                      disc.inst.B.prob = 0.445,
                      disc.inst.W.prob = 0.691,

                      circ.B.prob = 0.874,
                      circ.W.prob = 0.918,

                      ccr5.B.prob = c(0, 0.034),
                      ccr5.W.prob = c(0.021, 0.176),
                      ccr5.heteroz.rr = 0.3,

                      num.inst.ai.classes = 1,
                      base.ai.main.BB.rate = 0.17,
                      base.ai.main.BW.rate = 0.26,
                      base.ai.main.WW.rate = 0.23,
                      base.ai.pers.BB.rate = 0.11,
                      base.ai.pers.BW.rate = 0.16,
                      base.ai.pers.WW.rate = 0.14,
                      ai.scale = 1,

                      cond.main.BB.prob = 0.38,
                      cond.main.BW.prob = 0.10,
                      cond.main.WW.prob = 0.15,
                      cond.pers.always.prob = 0.216,
                      cond.pers.BB.prob = 0.26,
                      cond.pers.BW.prob = 0.26,
                      cond.pers.WW.prob = 0.26,
                      cond.inst.always.prob = 0.326,
                      cond.inst.BB.prob = 0.27,
                      cond.inst.BW.prob = 0.27,
                      cond.inst.WW.prob = 0.27,
                      cond.always.prob.corr = 0.5,
                      cond.rr.BB = 1,
                      cond.rr.BW = 1,
                      cond.rr.WW = 1,
                      cond.diag.main.beta = -0.67,
                      cond.discl.main.beta = -0.85,
                      cond.diag.pers.beta = -0.67,
                      cond.discl.pers.beta = -0.85,
                      cond.diag.inst.beta = -0.67,
                      cond.discl.inst.beta = -0.85,

                      vv.iev.BB.prob = 0.42,
                      vv.iev.BW.prob = 0.56,
                      vv.iev.WW.prob = 0.49,

                      prep.start = Inf,
                      prep.elig.model = "base",
                      prep.class.prob = c(0.211, 0.07, 0.1, 0.619),
                      prep.class.hr = c(1, 0.69, 0.19, 0.05),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
         call. = FALSE)
  }

  if (race.method == 1) {
    p$last.neg.test.B.int = (last.neg.test.B.int + last.neg.test.W.int)/2
    p$last.neg.test.W.int = (last.neg.test.B.int + last.neg.test.W.int)/2
    p$mean.test.B.int = (mean.test.W.int + mean.test.B.int)/2
    p$mean.test.W.int = (mean.test.W.int + mean.test.B.int)/2
    p$tt.traj.B.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
    p$tt.traj.W.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
    p$tx.init.B.prob = (tx.init.B.prob + tx.init.W.prob)/2
    p$tx.init.W.prob = (tx.init.B.prob + tx.init.W.prob)/2
    p$tx.halt.B.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
    p$tx.halt.W.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
    p$tx.reinit.B.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
    p$tx.reinit.W.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
    p$disc.outset.main.B.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
    p$disc.outset.main.W.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
    p$disc.outset.pers.B.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
    p$disc.outset.pers.W.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
    p$disc.inst.B.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
    p$disc.inst.W.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
    p$circ.B.prob = (circ.B.prob + circ.W.prob)/2
    p$circ.W.prob = (circ.B.prob + circ.W.prob)/2
    p$ccr5.B.prob = (ccr5.B.prob + ccr5.W.prob)/2
    p$ccr5.W.prob = (ccr5.B.prob + ccr5.W.prob)/2
    p$base.ai.main.BB.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.main.BW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.main.WW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.pers.BB.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$base.ai.pers.BW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$base.ai.pers.WW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$cond.main.BB.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.main.BW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.main.WW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.pers.BB.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.pers.BW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.pers.WW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.inst.BB.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$cond.inst.BW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$cond.inst.WW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$vv.iev.BB.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
    p$vv.iev.BW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
    p$vv.iev.WW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
  }

  p$time.unit <- nwstats$time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)

  p$role.B.prob <- nwstats$role.B.prob
  p$role.W.prob <- nwstats$role.W.prob

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)


  p$riskh.start <- max(1, prep.start - prep.risk.int - 1)

  p$method <- nwstats$method
  p$modes <- 1

  p$asmr.B <- nwstats$asmr.B
  p$asmr.W <- nwstats$asmr.W

  p$nwstats <- NULL

  class(p) <- "param.net"
  return(p)
}


#' @export
init_msmhet <- function(nwstats,
                     prev.B.male = 0.15,
                     prev.W.male = 0.15,
                     prev.B.feml = 0.01,
                     prev.W.feml = 0.01,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  p$num.B <- nwstats$num.B + nwstats$num.B.het
  p$num.W <- nwstats$num.W + nwstats$num.W.het

  p$ages <- nwstats$ages

  p$init.prev.age.slope.B.male <- 0.05 / 12
  p$init.prev.age.slope.W.male <- 0.05 / 12
  p$init.prev.age.slope.B.feml <- 0.05 / 12
  p$init.prev.age.slope.W.feml <- 0.05 / 12

  p$nwstats <- NULL

  class(p) <- "init.net"
  return(p)
}


#' @export
control_msmhet <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_msm,
                        aging.FUN = aging_msm,
                        deaths.FUN = deaths_msm,
                        births.FUN = births_msm,
                        test.FUN = test_msm,
                        tx.FUN = tx_msm,
                        prep.FUN = NULL,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_msm,
                        aiclass.FUN = NULL,
                        roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_msm,
                        disclose.FUN = disclose_msm,
                        acts.FUN = acts_msm,
                        condoms.FUN = condoms_msm,
                        riskhist.FUN = NULL,
                        position.FUN = position_msm,
                        trans.FUN = trans_msm,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose_msm,
                        save.nwstats = FALSE,
                        verbose = TRUE,
                        verbose.int = 1,
                        ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  p$save.other = c("attr", "temp", "el", "p")

  p$save.network = FALSE

  class(p) <- "control.net"
  return(p)
}


#' @export
initialize_msmhet <- function(x, param, init, control, s) {

  # Master data list
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  ## Network simulation ##
  nw <- list()
  for (i in 1:3) {
    nw[[i]] <- simulate(x[[i]]$fit)
    nw[[i]] <- remove_bad_roles_msm(nw[[i]])
  }

  ## ergm_prep here
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    attributes(dat$el[[i]])$vnames <- NULL
    p <- tergmLite::stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                                x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- tergmLite::ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p


  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }


  ## Nodal attributes ##

  # Male
  dat$attr$male <- get.vertex.attribute(nw[[1]], "male")

  # Race
  dat$attr$race <- get.vertex.attribute(nw[[1]], "race")
  num.B <- sum(dat$attr$race == "B")
  num.W <- sum(dat$attr$race == "W")
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")

  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$sqrt.age <- get.vertex.attribute(nw[[1]], "sqrt.age")
  dat$attr$age <- dat$attr$sqrt.age^2

  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[3]], "riskg")

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(num, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers <- uai.always[, 1]
  dat$attr$cond.always.inst <- uai.always[, 2]

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  idsMale.B <- which(dat$attr$male == 1 & dat$attr$race == "B")
  idsMale.W <- which(dat$attr$male == 1 & dat$attr$race == "W")

  circ <- rep(NA, num)
  circ[idsMale.B] <- sample(apportion_lr(length(idsMale.B), 0:1, 1 - param$circ.B.prob))
  circ[idsMale.W] <- sample(apportion_lr(length(idsMale.W), 0:1, 1 - param$circ.W.prob))
  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)

  # MSM class
  msm.class <- get.vertex.attribute(nw[[1]], "msm.class")
  dat$attr$msm.class <- msm.class

  # Role class
  role.class <- get.vertex.attribute(nw[[1]], "role.class")
  dat$attr$role.class <- role.class

  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[which(role.class == "I")]  <- 1
  ins.quot[which(role.class == "R")]  <- 0
  ins.quot[which(role.class == "V")]  <- runif(sum(role.class == "V", na.rm = TRUE))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msmhet(dat)

  # CCR5
  dat <- init_ccr5_msmhet(dat)


  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- matrix(NA, nrow = 0, ncol = 3)
  colnames(dat$temp$discl.list) <- c("pos", "neg", "discl.time")

  if (control$save.nwstats == TRUE) {
    dat$stats <- list()
    dat$stats$nwstats <- list()

  }

  dat <- prevalence_msmhet(dat, at = 1)

  class(dat) <- "dat"
  return(dat)
}


#' @export
init_status_msmhet <- function(dat) {

  race <- dat$attr$race
  male <- dat$attr$male

  num.B.male <- sum(race == "B" & male == 1)
  num.W.male <- sum(race == "W" & male == 1)
  num.B.feml <- sum(race == "B" & male == 0)
  num.W.feml <- sum(race == "W" & male == 0)

  ids.B.male <- which(race == "B" & male == 1)
  ids.W.male <- which(race == "W" & male == 1)
  ids.B.feml <- which(race == "B" & male == 0)
  ids.W.feml <- which(race == "W" & male == 0)

  num.B <- sum(race == "B")
  num.W <- sum(race == "W")
  ids.B <- union(ids.B.male, ids.B.feml)
  ids.W <- union(ids.W.male, ids.W.feml)

  num <- length(race)
  age <- dat$attr$age

  # Infection Status
  nInfB.male <- round(dat$init$prev.B.male * num.B.male)
  nInfW.male <- round(dat$init$prev.W.male * num.W.male)
  nInfB.feml <- round(dat$init$prev.B.feml * num.B.feml)
  nInfW.feml <- round(dat$init$prev.W.feml * num.W.feml)

  # Infection probability
  probInfB.male <- dat$init$prev.B.male
  probInfW.male <- dat$init$prev.W.male
  probInfB.feml <- dat$init$prev.B.feml
  probInfW.feml <- dat$init$prev.W.feml

  # Infection status
  status <- rep(0, num)
  while (sum(status[ids.B.male]) != nInfB.male) {
    status[ids.B.male] <- rbinom(num.B.male, 1, probInfB.male)
  }
  while (sum(status[ids.W.male]) != nInfW.male) {
    status[ids.W.male] <- rbinom(num.W.male, 1, probInfW.male)
  }
  while (sum(status[ids.B.feml]) != nInfB.feml) {
    status[ids.B.feml] <- rbinom(num.B.feml, 1, probInfB.feml)
  }
  while (sum(status[ids.W.feml]) != nInfW.feml) {
    status[ids.W.feml] <- rbinom(num.W.feml, 1, probInfW.feml)
  }
  dat$attr$status <- status

  # Treatment trajectory
  tt.traj <- rep(NA, num)

  tt.traj[ids.B] <- sample(apportion_lr(num.B, 1:4, dat$param$tt.traj.B.prob))
  tt.traj[ids.W] <- sample(apportion_lr(num.W, 1:4, dat$param$tt.traj.W.prob))
  dat$attr$tt.traj <- tt.traj


  ## Infection-related attributes
  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  last.neg.test <- rep(NA, num)
  tx.status <- rep(NA, num)
  tx.init.time <- rep(NA, num)
  cum.time.on.tx <- rep(NA, num)
  cum.time.off.tx <- rep(NA, num)

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                  min(dat$init$ages), 0))

  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vlds <- (vlf - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int


  ### Non-treater type: tester and non-tester
  selected <- which(status == 1 & tt.traj %in% c(1, 2))
  max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  tx.status[selected] <- 0
  cum.time.on.tx[selected] <- 0
  cum.time.off.tx[selected] <- time.since.inf

  stage[selected[time.since.inf <= vlar.int]] <- 1
  stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- 2
  stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- 3
  stage[selected[time.since.inf > vldo.int]] <- 4

  stage.time[selected][stage[selected] == 1] <- time.since.inf[stage[selected] == 1]
  stage.time[selected][stage[selected] == 2] <- time.since.inf[stage[selected] == 2] -
    vlar.int
  stage.time[selected][stage[selected] == 3] <- time.since.inf[stage[selected] == 3] -
    vl.acute.int
  stage.time[selected][stage[selected] == 4] <- time.since.inf[stage[selected] == 4] -
    vldo.int

  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
    (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

  selected <- which(status == 1 & tt.traj == 1)
  diag.status[selected] <- 0

  selected <- which(status == 1 & tt.traj == 2)

  # Time to next test
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                               dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                            dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  twind.int <- dat$param$test.window.int
  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1


  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

  tx.init.time.B <- twind.int + dat$param$last.neg.test.B.int + 1 / dat$param$tx.init.B.prob
  tx.init.time.W <- twind.int + dat$param$last.neg.test.W.int + 1 / dat$param$tx.init.W.prob

  # Stage for Blacks
  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
    (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  numsteps.B <- (dat$param$max.time.off.tx.full.int - tx.init.time.B) /
    (1 - prop.time.on.tx.B)
  offon.B <- rbind(offon.B,
                   cbind(tx.init.time.B + (1 - prop.time.on.tx.B) * 1:numsteps.B,
                         prop.time.on.tx.B * 1:numsteps.B))
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  # Stage for Whites
  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
    (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)
  numsteps.W <- (dat$param$max.time.off.tx.full.int - tx.init.time.W) /
    (1 - prop.time.on.tx.W)
  offon.W <- rbind(offon.W,
                   cbind(tx.init.time.W + (1 - prop.time.on.tx.W) * 1:numsteps.W,
                         prop.time.on.tx.W * 1:numsteps.W))
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # Vl for Blacks
  selected <- which(status == 1 & tt.traj == 4 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B) * (vlsp) +
    (time.since.inf > exp.onset.aids.B) *
    (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 4 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W) * (vlsp) +
    (time.since.inf > exp.onset.aids.W) *
    (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # Diagnosis
  selected <- which(status == 1 & tt.traj == 4)
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                               dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                            dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]
  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  ### Part adherent type

  # Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
    (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  while (offon.B[nrow(offon.B), 1] / dat$param$max.time.off.tx.part.int +
         offon.B[nrow(offon.B), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B <- rbind(offon.B,
                     offon.B[nrow(offon.B), ] + c(1 - prop.time.on.tx.B,
                                                  prop.time.on.tx.B))
  }
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
    (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)

  while (offon.W[nrow(offon.W), 1] / dat$param$max.time.off.tx.part.int +
         offon.W[nrow(offon.W), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W <- rbind(offon.W,
                     offon.W[nrow(offon.W), ] + c(1 - prop.time.on.tx.W,
                                                  prop.time.on.tx.W))
  }
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # VL for Blacks
  selected <- which(status == 1 & tt.traj == 3 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B) * (vlsp) +
    (time.since.inf > exp.onset.aids.B) *
    (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 3 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W) * (vlsp) +
    (time.since.inf > exp.onset.aids.W) *
    (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == 3)
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                               dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                            dat$param$mean.test.W.int * (race[selected] == "W")))
  }


  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(2, 3, 4))

  if (dat$param$testing.pattern == "interval") {
    tslt <- ceiling(runif(length(selected),
                          min = 0,
                          max = dat$param$mean.test.B.int * (race[selected] == "B") +
                            dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    tslt <- rgeom(length(selected),
                  1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                         dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  last.neg.test[selected] <- -tslt


  ## Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$tx.status <- tx.status
  dat$attr$tx.init.time <- tx.init.time
  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx

  return(dat)

}


#' @export
init_ccr5_msmhet <- function(dat) {

  race <- dat$attr$race
  num.B <- sum(race == "B")
  num.W <- sum(race == "W")
  num <- num.B + num.W
  status <- dat$attr$status

  nInfB <- sum(race == "B" & status == 1)
  nInfW <- sum(race == "W" & status == 1)

  ##  CCR5 genotype
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  ccr5 <- rep("WW", num)

  # homozygotes for deletion
  num.ccr5.DD.B <- dat$param$ccr5.B.prob[1] * num.B
  # heterozygotes
  num.ccr5.DW.B <- dat$param$ccr5.B.prob[2] * num.B
  # homozygotes for deletion
  num.ccr5.WW.B <- num.B - num.ccr5.DD.B - num.ccr5.DW.B
  # DD's can't be infected
  num.uninf.ccr5.DD.B <- round(num.ccr5.DD.B)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B <- round(num.ccr5.DW.B * nInfB * ccr5.heteroz.rr /
                               (num.ccr5.WW.B + num.ccr5.DW.B * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B <- round(num.ccr5.DW.B - num.inf.ccr5.DW.B)
  inf.B <- which(status == 1 & race == "B")
  inf.ccr5.DW.B <- sample(inf.B, num.inf.ccr5.DW.B, replace = FALSE)
  ccr5[inf.ccr5.DW.B] <- "DW"
  uninf.B <- which(status == 0 & race == "B")
  uninf.ccr5.DWDD.B <- sample(uninf.B, num.uninf.ccr5.DW.B + num.uninf.ccr5.DD.B)
  uninf.ccr5.DW.B <- sample(uninf.ccr5.DWDD.B, num.uninf.ccr5.DW.B)
  uninf.ccr5.DD.B <- setdiff(uninf.ccr5.DWDD.B, uninf.ccr5.DW.B)
  ccr5[uninf.ccr5.DW.B] <- "DW"
  ccr5[uninf.ccr5.DD.B] <- "DD"

  num.ccr5.DD.W <- dat$param$ccr5.W.prob[1] * num.W
  num.ccr5.DW.W <- dat$param$ccr5.W.prob[2] * num.W
  num.ccr5.WW.W <- num.W - num.ccr5.DD.W - num.ccr5.DW.W
  num.uninf.ccr5.DD.W <- round(num.ccr5.DD.W)
  num.inf.ccr5.DW.W <- round(num.ccr5.DW.W * nInfW * ccr5.heteroz.rr /
                               (num.ccr5.WW.W + num.ccr5.DW.W * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.W <- round(num.ccr5.DW.W - num.inf.ccr5.DW.W)
  inf.W <- which(status == 1 & race == "W")
  inf.ccr5.DW.W <- sample(inf.W, num.inf.ccr5.DW.W)
  ccr5[inf.ccr5.DW.W] <- "DW"
  uninf.W <- which(status == 0 & race == "W")
  uninf.ccr5.DWDD.W <- sample(uninf.W, num.uninf.ccr5.DW.W + num.uninf.ccr5.DD.W)
  uninf.ccr5.DW.W <- sample(uninf.ccr5.DWDD.W, num.uninf.ccr5.DW.W)
  uninf.ccr5.DD.W <- setdiff(uninf.ccr5.DWDD.W, uninf.ccr5.DW.W)
  ccr5[uninf.ccr5.DW.W] <- "DW"
  ccr5[uninf.ccr5.DD.W] <- "DD"

  dat$attr$ccr5 <- ccr5

  return(dat)
}


#' @export
reinit_msmhet <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}


#' @export
deaths_msmhet <- function(dat, at) {

  ## General deaths
  age <- floor(dat$attr$age)
  race <- dat$attr$race

  alive.B <- which(race == "B")
  age.B <- age[alive.B]
  death.B.prob <- dat$param$asmr.B[age.B]
  deaths.B <- alive.B[rbinom(length(death.B.prob), 1, death.B.prob) == 1]

  alive.W <- which(race == "W")
  age.W <- age[alive.W]
  death.W.prob <- dat$param$asmr.W[age.W]
  deaths.W <- alive.W[rbinom(length(death.W.prob), 1, death.W.prob) == 1]

  dth.gen <- c(deaths.B, deaths.W)


  ## Disease deaths
  dth.dis <- which(dat$attr$stage == 4 &
                     dat$attr$vl >= dat$param$vl.fatal)

  dth.all <- NULL
  dth.all <- unique(c(dth.gen, dth.dis))

  if (length(dth.all) > 0) {
    dat$attr$active[dth.all] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], dth.all)
    }
    dat$attr <- deleteAttr(dat$attr, dth.all)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in death mod")
    }
  }


  ## Summary Output
  dat$epi$dth.gen[at] <- length(dth.gen)
  dat$epi$dth.dis[at] <- length(dth.dis)

  return(dat)
}


#' @export
births_msmhet <- function(dat, at){

  ## Variables

  # Parameters
  b.B.rate <- dat$param$b.B.rate
  b.W.rate <- dat$param$b.W.rate
  b.method <- dat$param$b.method


  ## Process
  if (b.method == "fixed") {
    numB <- dat$epi$num.B[1]
    numW <- dat$epi$num.W[1]
  }
  if (b.method == "varying") {
    numB <- dat$epi$num.B[at - 1]
    numW <- dat$epi$num.W[at - 1]
  }

  nBirths.B <- rpois(1, b.B.rate * numB)
  nBirths.W <- rpois(1, b.W.rate * numW)
  nBirths <- nBirths.B + nBirths.W


  ## Update Attr
  if (nBirths > 0) {
    dat <- setBirthAttr_msmhet(dat, at, nBirths.B, nBirths.W)
  }
  stopifnot(length(unique(vapply(dat$attr, length, FUN.VALUE = 1L))) == 1)

  # Update Networks
  if (nBirths > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::add_vertices(dat$el[[i]], nBirths)
    }
  }


  ## Output
  dat$epi$nBirths[at] <- nBirths

  return(dat)
}

#' @export
setBirthAttr_msmhet <- function(dat, at, nBirths.B, nBirths.W) {

  nBirths <- nBirths.B + nBirths.W

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths)
  dat$temp$max.uid <- dat$temp$max.uid + nBirths

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)

  male <- sample(0:1, nBirths, TRUE)
  dat$attr$male[newIds] <- male

  race <- sample(rep(c("B", "W"), c(nBirths.B, nBirths.W)))
  newB <- which(race == "B")
  newW <- which(race == "W")
  dat$attr$race[newIds] <- race

  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])

  # Disease status and related
  dat$attr$status[newIds] <- rep(0, nBirths)

  dat$attr$tt.traj[newIds[newB]] <- sample(1:4,
                                           nBirths.B,
                                           replace = TRUE,
                                           prob = dat$param$tt.traj.B.prob)
  dat$attr$tt.traj[newIds[newW]] <- sample(1:4,
                                           nBirths.W, replace = TRUE,
                                           prob = dat$param$tt.traj.W.prob)

  # Circumcision
  newB.male <- which(male == 1 & race == "B")
  newW.male <- which(male == 1 & race == "W")

  dat$attr$circ[newIds[newB.male]] <- rbinom(length(newB.male), 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[newW.male]] <- rbinom(length(newW.male), 1, dat$param$circ.W.prob)

  # msm.class
  msm.class <- rep(NA, nBirths)
  msm.class[male == 0] <- 4
  msm.class[male == 1] <- sample(1:3, sum(male == 1), TRUE, c(0.1, 0.1, 0.8))
  dat$attr$msm.class[newIds] <- msm.class

  # Role
  table(dat$attr$msm.class, dat$attr$role.class)
  newB.msm <- which(msm.class %in% 1:2 & race == "B")
  newW.msm <- which(msm.class %in% 1:2 & race == "W")

  dat$attr$role.class[newIds[newB.msm]] <- sample(c("I", "R", "V"),
                                              length(newB.msm), replace = TRUE,
                                              prob = dat$param$role.B.prob)
  dat$attr$role.class[newIds[newW.msm]] <- sample(c("I", "R", "V"),
                                              length(newW.msm), replace = TRUE,
                                              prob = dat$param$role.W.prob)

  ins.quot <- rep(NA, nBirths)
  ins.quot[which(dat$attr$role.class[newIds] == "I")]  <- 1
  ins.quot[which(dat$attr$role.class[newIds] == "R")]  <- 0
  ins.quot[which(dat$attr$role.class[newIds] == "V")]  <- runif(sum(dat$attr$role.class[newIds] == "V", na.rm = TRUE))
  dat$attr$ins.quot[newIds] <- ins.quot

  # CCR5
  ccr5.B.prob <- dat$param$ccr5.B.prob
  ccr5.W.prob <- dat$param$ccr5.W.prob
  dat$attr$ccr5[newIds[newB]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.prob),
                                                 ccr5.B.prob[2], ccr5.B.prob[1]))
  dat$attr$ccr5[newIds[newW]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W, replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.prob),
                                                 ccr5.W.prob[2], ccr5.W.prob[1]))


  # One-off risk group
  table(dat$attr$riskg, dat$attr$msm.class)

  dat$attr$riskg[newIds[msm.class %in% 1:2]] <- sample(1:5, sum(msm.class %in% 1:2), TRUE)

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(nBirths, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers[newIds] <- uai.always[, 1]
  dat$attr$cond.always.inst[newIds] <- uai.always[, 2]

  # PrEP
  dat$attr$prepStat[newIds] <- 0

  return(dat)
}


#' @export
test_msmhet <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.int <- dat$param$mean.test.B.int
  mean.test.W.int <- dat$param$mean.test.W.int
  twind.int <- dat$param$test.window.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  ## Process

  if (testing.pattern == "memoryless") {
    elig.B <- which(race == "B" &
                      tt.traj != 1 &
                      (diag.status == 0 | is.na(diag.status)) &
                      prepStat == 0)
    rates.B <- rep(1/mean.test.B.int, length(elig.B))
    tst.B <- elig.B[rbinom(length(elig.B), 1, rates.B) == 1]

    elig.W <- which(race == "W" &
                      tt.traj != 1 &
                      (diag.status == 0 | is.na(diag.status)) &
                      prepStat == 0)
    rates.W <- rep(1/mean.test.W.int, length(elig.W))
    tst.W <- elig.W[rbinom(length(elig.W), 1, rates.W) == 1]
    tst.nprep <- c(tst.B, tst.W)
  }

  if (testing.pattern == "interval") {
    tst.B <- which(race == "B" &
                     tt.traj != 1 &
                     (diag.status == 0 | is.na(diag.status)) &
                     tsincelntst >= 2*(mean.test.B.int) &
                     prepStat == 0)

    tst.W <- which(race == "W" &
                     tt.traj != 1 &
                     (diag.status == 0 | is.na(diag.status)) &
                     tsincelntst >= 2*(mean.test.W.int) &
                     prepStat == 0)
    tst.nprep <- c(tst.B, tst.W)
  }

  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                      prepStat == 1 &
                      tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.prep)

  tst.pos <- tst.all[status[tst.all] == 1 & inf.time[tst.all] <= at - twind.int]
  tst.neg <- setdiff(tst.all, tst.pos)

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  return(dat)
}


#' @export
tx_msmhet <- function(dat, at) {

  ## Variables

  # Attributes
  race <- dat$attr$race
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage

  # Parameters
  tx.init.B.prob <- dat$param$tx.init.B.prob
  tx.init.W.prob <- dat$param$tx.init.W.prob
  tx.halt.B.prob <- dat$param$tx.halt.B.prob
  tx.halt.W.prob <- dat$param$tx.halt.W.prob
  tx.reinit.B.prob <- dat$param$tx.reinit.B.prob
  tx.reinit.W.prob <- dat$param$tx.reinit.W.prob


  ## Initiation
  tx.init.elig.B <- which(race == "B" & status == 1 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.B <- tx.init.elig.B[rbinom(length(tx.init.elig.B), 1,
                                     tx.init.B.prob) == 1]

  tx.init.elig.W <- which(race == "W" & status == 1 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.W <- tx.init.elig.W[rbinom(length(tx.init.elig.W), 1,
                                     tx.init.W.prob) == 1]

  tx.init <- c(tx.init.B, tx.init.W)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  tx.halt.elig.B <- which(race == "B" & tx.status == 1)
  tx.halt.B <- tx.halt.elig.B[rbinom(length(tx.halt.elig.B), 1,
                                     tx.halt.B.prob) == 1]

  tx.halt.elig.W <- which(race == "W" & tx.status == 1)
  tx.halt.W <- tx.halt.elig.W[rbinom(length(tx.halt.elig.W),
                                     1, tx.halt.W.prob) == 1]
  tx.halt <- c(tx.halt.B, tx.halt.W)
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  tx.reinit.elig.B <- which(race == "B" & tx.status == 0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B <- tx.reinit.elig.B[rbinom(length(tx.reinit.elig.B),
                                         1, tx.reinit.B.prob) == 1]

  tx.reinit.elig.W <- which(race == "W" & tx.status == 0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W <- tx.reinit.elig.W[rbinom(length(tx.reinit.elig.W),
                                         1, tx.reinit.W.prob) == 1]

  tx.reinit <- c(tx.reinit.B, tx.reinit.W)
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  dat$attr$cum.time.on.tx <- dat$attr$cum.time.on.tx +
    ((dat$attr$tx.status == 1) %in% TRUE)
  dat$attr$cum.time.off.tx <- dat$attr$cum.time.off.tx +
    ((dat$attr$tx.status == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)

  return(dat)
}


#' @export
progress_msmhet <- function(dat, at) {

  ## Variables

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  time.since.inf <- at - dat$attr$inf.time
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  cum.time.off.tx <- dat$attr$cum.time.off.tx
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  tt.traj <- dat$attr$tt.traj
  tx.status <- dat$attr$tx.status

  # Parameters
  vl.acute.rise.int <- dat$param$vl.acute.rise.int
  vl.acute.fall.int <- dat$param$vl.acute.fall.int
  vl.aids.onset <- dat$param$vl.aids.onset
  max.time.off.tx.part <- dat$param$max.time.off.tx.part
  max.time.on.tx.part <- dat$param$max.time.on.tx.part

  max.time.off.tx.full <- dat$param$max.time.off.tx.full


  ## Process

  # Increment day
  stage.time[active == 1] <- stage.time[active == 1] + 1

  # Change stage to Acute Falling
  toAF <- which(active == 1 & time.since.inf == (vl.acute.rise.int + 1))
  stage[toAF] <- 2
  stage.time[toAF] <- 1

  # Change stage to Chronic
  toC <- which(active == 1 & time.since.inf == (vl.acute.rise.int +
                                                  vl.acute.fall.int + 1))
  stage[toC] <- 3
  stage.time[toC] <- 1

  # Change stage to AIDS
  aids.tx.naive <- which(active == 1 & status == 1 & cum.time.on.tx == 0 &
                           (time.since.inf >= vl.aids.onset) & stage != 4)

  part.tx.score <- (cum.time.off.tx / max.time.off.tx.part) +
    (cum.time.on.tx / max.time.on.tx.part)

  aids.part.escape <- which(active == 1 & cum.time.on.tx > 0 & tt.traj == 3 &
                              stage == 3 & part.tx.score >= 1 & stage != 4)

  aids.off.tx.full.escape <- which(active == 1 & tx.status == 0 & tt.traj == 4 &
                                     cum.time.on.tx > 0 &
                                     cum.time.off.tx >= max.time.off.tx.full &
                                     stage != 4)

  isAIDS <- c(aids.tx.naive, aids.part.escape, aids.off.tx.full.escape)
  stage[isAIDS] <- 4
  stage.time[isAIDS] <- 1


  ## Output
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time

  return(dat)
}


#' @export
vl_msmhet <- function(dat, at) {

  ## Variables

  # Attributes
  inf.time.bp <- at - dat$attr$inf.time
  # cum.time.off.tx <- dat$attr$cum.time.off.tx
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  status <- dat$attr$status
  tt.traj <- dat$attr$tt.traj
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  tx.status <- dat$attr$tx.status

  # Parameters
  vlard <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlafd <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo <- dat$param$vl.aids.onset
  vldd <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vl.full.supp <- dat$param$vl.full.supp
  full.supp.down.slope <- dat$param$full.supp.down.slope
  vl.part.supp <- dat$param$vl.part.supp
  part.supp.down.slope <- dat$param$part.supp.down.slope
  full.supp.up.slope <- dat$param$full.supp.up.slope
  part.supp.up.slope <- dat$param$part.supp.up.slope
  # max.time.off.tx.part <- dat$param$max.time.off.tx.part
  # max.time.on.tx.part <- dat$param$max.time.on.tx.part

  # Calculations
  vlds <- (vlf - vlsp) / vldd
  # part.tx.score <-  (cum.time.off.tx / max.time.off.tx.part) +
  #                   (cum.time.on.tx / max.time.on.tx.part)


  ## Process

  # 1. tx-naive men
  target <- which(status == 1 & cum.time.on.tx == 0)
  inf.time.bp.tn <- inf.time.bp[target]
  new.vl <- (inf.time.bp.tn <= vlard) * (vlap * inf.time.bp.tn / vlard) +
    (inf.time.bp.tn > vlard) * (inf.time.bp.tn <= vlard + vlafd) *
    ((vlsp - vlap) * (inf.time.bp.tn - vlard) / vlafd + vlap) +
    (inf.time.bp.tn > vlard + vlafd) * (inf.time.bp.tn <= vldo) * (vlsp) +
    (inf.time.bp.tn > vldo) * (vlsp + (inf.time.bp.tn - vldo) * vlds)
  vl[target] <- new.vl

  # 2. men on tx, tt.traj=full, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 4 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - full.supp.down.slope, vl.full.supp)
  vl[target] <- new.vl

  # 3. men on tx, tt.traj=part, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 3 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - part.supp.down.slope, vl.part.supp)
  vl[target] <- new.vl

  # 4. men off tx, not naive, tt.traj=full, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                    cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + full.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 5. men off tx, not naive, tt.traj=part, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                    cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + part.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 6. men on tx, tt.traj=full, escaped
  # Doesn't exist.

  # 7. men on tx, tt.traj=part, escaped
  target <- which(tx.status == 1 &
                    tt.traj == 3 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 8. men off tx, tt.traj=full, and escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                    cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 9. men off tx, tt.traj=part, and escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                    cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl


  ## Output
  dat$attr$vl <- vl

  return(dat)
}


#' @export
simnet_msmhet <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.m$coef.form,
                                             coef.diss = nwparam.m$coef.diss$coef.adj,
                                             save.changes = TRUE)

  dat$temp$new.edges <- NULL
  if (at == 2) {
    new.edges.m <- matrix(dat$el[[1]], ncol = 2)
  } else {
    new.edges.m <- attributes(dat$el[[1]])$changes
    new.edges.m <- new.edges.m[new.edges.m[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.m], ncol = 2)


  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

  dat <- tergmLite::updateModelTermInputs(dat, network = 2)

  dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
                                             el = dat$el[[2]],
                                             coef.form = nwparam.p$coef.form,
                                             coef.diss = nwparam.p$coef.diss$coef.adj,
                                             save.changes = TRUE)

  if (at == 2) {
    new.edges.p <- matrix(dat$el[[2]], ncol = 2)
  } else {
    new.edges.p <- attributes(dat$el[[2]])$changes
    new.edges.p <- new.edges.p[new.edges.p[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- rbind(dat$temp$new.edges,
                              matrix(dat$attr$uid[new.edges.p], ncol = 2))


  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)

  dat <- tergmLite::updateModelTermInputs(dat, network = 3)

  dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
                                          el = dat$el[[3]],
                                          coef = nwparam.i$coef.form)

  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_resim_nwstats(dat, at)
  }

  return(dat)
}


#' @export
disclose_msmhet <- function(dat, at){

  for (type in c("main", "pers", "inst")) {

    # Variables --------------------------------------------------------------

    # Attributes
    status <- dat$attr$status
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    diag.time <- dat$attr$diag.time
    race <- dat$attr$race

    # Parameters and network
    if (type == "main") {
      disc.outset.B.prob <- dat$param$disc.outset.main.B.prob
      disc.at.diag.B.prob <- dat$param$disc.at.diag.main.B.prob
      disc.post.diag.B.prob <- dat$param$disc.post.diag.main.B.prob
      disc.outset.W.prob <- dat$param$disc.outset.main.W.prob
      disc.at.diag.W.prob <- dat$param$disc.at.diag.main.W.prob
      disc.post.diag.W.prob <- dat$param$disc.post.diag.main.W.prob
      el <- dat$el[[1]]
    }

    if (type == "pers") {
      disc.outset.B.prob <- dat$param$disc.outset.pers.B.prob
      disc.at.diag.B.prob <- dat$param$disc.at.diag.pers.B.prob
      disc.post.diag.B.prob <- dat$param$disc.post.diag.pers.B.prob
      disc.outset.W.prob <- dat$param$disc.outset.pers.W.prob
      disc.at.diag.W.prob <- dat$param$disc.at.diag.pers.W.prob
      disc.post.diag.W.prob <- dat$param$disc.post.diag.pers.W.prob
      el <- dat$el[[2]]
    }

    if (type == "inst") {
      disc.inst.B.prob <- dat$param$disc.inst.B.prob
      disc.inst.W.prob <- dat$param$disc.inst.W.prob
      el <- dat$el[[3]]
    }


    # Processes --------------------------------------------------------------

    # Check for discordant rels
    posneg <- el[which(status[el[, 1]] - status[el[, 2]] == 1), , drop = FALSE]
    negpos <- el[which(status[el[, 2]] - status[el[, 1]] == 1), , drop = FALSE]
    disc.el <- rbind(posneg, negpos[, 2:1])

    # Check for not already disclosed
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    discord.cdl <- uid[disc.el[, 1]] * 1e7 + uid[disc.el[, 2]]
    notdiscl <- !(discord.cdl %in% disclose.cdl)


    # data frame of non-disclosed pairs
    nd <- disc.el[notdiscl, , drop = FALSE]

    # Check for positive diagnosis
    notdiscl.dx <- which(diag.status[nd[, 1]] == 1)

    # data frame of non-disclosed pairs where infected is dx'ed
    nd.dx <- nd[notdiscl.dx, , drop = FALSE]

    # If there are any eligible pairs
    if (nrow(nd.dx) > 0) {

      # Split by race of pos node
      pos.race <- race[nd.dx[, 1]]

      if (type %in% c("main", "pers")) {

        # Check that rel is new
        # new.edges matrix is expressed in uid, so need to transform nd.dx
        new.edges <- dat$temp$new.edges
        new.rel <- ((uid[nd.dx[, 1]] * 1e7 + uid[nd.dx[, 2]]) %in%
                      (new.edges[, 1] * 1e7 + new.edges[, 2])) |
          ((uid[nd.dx[, 2]] * 1e7 + uid[nd.dx[, 1]]) %in%
             (new.edges[, 1] * 1e7 + new.edges[, 2]))

        # Check if diag is new
        new.dx <- diag.time[nd.dx[, 1]] == at

        # Assign disclosure probs
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[pos.race == "B" & new.rel == TRUE] <- disc.outset.B.prob
        dl.prob[pos.race == "B" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.prob
        dl.prob[pos.race == "B" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.prob

        dl.prob[pos.race == "W" & new.rel == TRUE] <- disc.outset.W.prob
        dl.prob[pos.race == "W" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.prob
        dl.prob[pos.race == "W" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.prob
      }

      if (type == "inst") {
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[pos.race == "B"] <- disc.inst.B.prob
        dl.prob[pos.race == "W"] <- disc.inst.W.prob
      }

      # Determine disclosers
      discl <- which(rbinom(length(dl.prob), 1, dl.prob) == 1)

      # Write output
      if (length(discl) > 0) {
        discl.mat <- cbind(pos = uid[nd.dx[discl, 1]],
                           neg = uid[nd.dx[discl, 2]],
                           discl.time = at)
        dat$temp$discl.list <- rbind(dat$temp$discl.list, discl.mat)
      }
    }
  }

  if (at > 2) {
    discl.list <- dat$temp$discl.list
    master.el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])
    m <- which(match(discl.list[, 1] * 1e7 + discl.list[, 2],
                     uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]) |
                 match(discl.list[, 2] * 1e7 + discl.list[, 1],
                       uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]))
    dat$temp$discl.list <- discl.list[m, ]
  }

  return(dat)
}

#' @export
acts_msmhet <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    status <- dat$attr$status
    race <- dat$attr$race

    # Parameters
    ai.scale <- dat$param$ai.scale
    if (type == "main") {
      base.ai.BB.rate <- dat$param$base.ai.main.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.main.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.main.WW.rate
      fixed <- FALSE
      ptype <- 1
      el <- dat$el[[1]]
    }
    if (type == "pers") {
      base.ai.BB.rate <- dat$param$base.ai.pers.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.pers.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.pers.WW.rate
      fixed <- FALSE
      ptype <- 2
      el <- dat$el[[2]]
    }
    if (type == "inst") {
      base.ai.BB.rate <- 1
      base.ai.BW.rate <- 1
      base.ai.WW.rate <- 1
      fixed <- ifelse(ai.scale != 1, FALSE, TRUE)
      ptype <- 3
      el <- dat$el[[3]]
    }

    ## Processes ##

    # Construct edgelist

    st1 <- status[el[, 1]]
    st2 <- status[el[, 2]]
    disc <- abs(st1 - st2) == 1
    el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1]
    el <- cbind(el, status[el[, 1]], status[el[, 2]])
    colnames(el) <- c("p1", "p2", "st1", "st2")

    if (nrow(el) > 0) {

      # Base AI rates
      ai.rate <- rep(NA, nrow(el))
      race.p1 <- race[el[, 1]]
      race.p2 <- race[el[, 2]]
      num.B <- (race.p1 == "B") + (race.p2 == "B")
      ai.rate <- (num.B == 2) * base.ai.BB.rate +
        (num.B == 1) * base.ai.BW.rate +
        (num.B == 0) * base.ai.WW.rate
      ai.rate <- ai.rate * ai.scale

      # Final act number
      if (fixed == FALSE) {
        ai <- rpois(length(ai.rate), ai.rate)
      } else {
        ai <- round(ai.rate)
      }

      # Full edge list
      el <- cbind(el, ptype, ai)
      colnames(el)[5:6] <- c("ptype", "ai")

      if (type == "main") {
        dat$temp$el <- el
      } else {
        dat$temp$el <- rbind(dat$temp$el, el)
      }
    }

  } # loop over type end

  # Remove inactive edges from el
  dat$temp$el <- dat$temp$el[-which(dat$temp$el[, "ai"] == 0), ]

  return(dat)
}


#' @export
condoms_msmhet <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race

    # Parameters
    cond.rr.BB <- dat$param$cond.rr.BB
    cond.rr.BW <- dat$param$cond.rr.BW
    cond.rr.WW <- dat$param$cond.rr.WW

    if (type == "main") {
      cond.BB.prob <- dat$param$cond.main.BB.prob
      cond.BW.prob <- dat$param$cond.main.BW.prob
      cond.WW.prob <- dat$param$cond.main.WW.prob
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      cond.always <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.BB.prob <- dat$param$cond.pers.BB.prob
      cond.BW.prob <- dat$param$cond.pers.BW.prob
      cond.WW.prob <- dat$param$cond.pers.WW.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "inst") {
      cond.BB.prob <- dat$param$cond.inst.BB.prob
      cond.BW.prob <- dat$param$cond.inst.BW.prob
      cond.WW.prob <- dat$param$cond.inst.WW.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
      ptype <- 3
    }

    el <- dat$temp$el
    elt <- el[el[, "ptype"] == ptype, ]

    ## Process ##

    # Base condom probs
    race.p1 <- race[elt[, 1]]
    race.p2 <- race[elt[, 2]]
    num.B <- (race.p1 == "B") + (race.p2 == "B")
    cond.prob <- (num.B == 2) * (cond.BB.prob * cond.rr.BB) +
      (num.B == 1) * (cond.BW.prob * cond.rr.BW) +
      (num.B == 0) * (cond.WW.prob * cond.rr.WW)

    # Transform to UAI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt[, 1]]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta

    # Disclosure modifier
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1)
    delt <- elt[isDiscord, ]
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))

    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1

    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always[elt[, 1]]
      ca2 <- cond.always[elt[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
    }

    ai.vec <- elt[, "ai"]
    pos <- rep(elt[, "p1"], ai.vec)
    neg <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(pos), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(pos, neg, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(pos, neg, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  return(dat)
}


#' @export
position_msmhet <- function(dat, at) {

  ## Variables
  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  status <- dat$attr$status
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dat$temp$al <- NULL

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  vv.iev.BB.prob <- dat$param$vv.iev.BB.prob
  vv.iev.BW.prob <- dat$param$vv.iev.BW.prob
  vv.iev.WW.prob <- dat$param$vv.iev.WW.prob


  ## Process
  pos.role.class <- role.class[dal[, 1]]
  neg.role.class <- role.class[dal[, 2]]

  ins <- rep(NA, length(pos.role.class))
  ins[which(pos.role.class == "I")] <- 1  # "P"
  ins[which(pos.role.class == "R")] <- 0  # "N"
  ins[which(neg.role.class == "I")] <- 0  # "N"
  ins[which(neg.role.class == "R")] <- 1  # "P"

  vv <- which(pos.role.class == "V" & neg.role.class == "V")
  vv.race.combo <- paste0(race[dal[, 1]][vv], race[dal[, 2]][vv])
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.BB.prob +
    (vv.race.combo == "BW") * vv.iev.BW.prob +
    (vv.race.combo == "WW") * vv.iev.WW.prob

  iev <- rbinom(length(vv), 1, vv.iev.prob)
  ins[vv[iev == 1]] <- 2 # "B"
  vv.remaining <- vv[iev == 0]

  inspos.prob <- ins.quot[dal[, 1][vv.remaining]] /
    (ins.quot[dal[, 1][vv.remaining]] + ins.quot[dal[, 2][vv.remaining]])
  inspos <- rbinom(length(vv.remaining), 1, inspos.prob)
  ins[vv.remaining[inspos == 1]] <- 1  # "P"
  ins[vv.remaining[inspos == 0]] <- 0  # "N"


  ## Output
  dat$temp$dal <- cbind(dal, ins)

  return(dat)
}


#' @export
trans_msmhet <- function(dat, at){

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr

  # Data
  dal <- dat$temp$dal
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]

  if (nrow(dal) == 0) {
    return(dat)
  }

  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- c("i", "r")
  colnames(disc.rp)[1:2] <- c("i", "r")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1-ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr)

  # CCR5
  ip.tlo[ip.ccr5 == "DD"] <- ip.tlo[ip.ccr5 == "DD"] + -Inf
  ip.tlo[ip.ccr5 == "DW"] <- ip.tlo[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo[temp.ids] <- ip.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c(1, 2))
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  # Retransformation to probability
  ip.tprob <- plogis(ip.tlo)
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)


  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1-rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr)

  # CCR5
  rp.tlo[rp.ccr5 == "DD"] <- rp.tlo[rp.ccr5 == "DD"] + -Inf
  rp.tlo[rp.ccr5 == "DW"] <- rp.tlo[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo[temp.ids] <- rp.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c(1, 2))
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  # Retransformation to probability
  rp.tprob <- plogis(rp.tlo)
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)

  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------

  # Update attributes

  infected <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])

    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- 1
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
  }

  # Summary Output
  dat$epi$incid[at] <- length(infected)


  return(dat)
}


#' @export
prevalence_msmhet <- function(dat, at) {

  race <- dat$attr$race
  status <- dat$attr$status
  male <- dat$attr$male

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA
    dat$epi$num.B <- rNA
    dat$epi$num.W <- rNA
    dat$epi$num.B.male <- rNA
    dat$epi$num.W.male <- rNA
    dat$epi$num.B.feml <- rNA
    dat$epi$num.W.feml <- rNA
    dat$epi$i.num <- rNA
    dat$epi$i.num.B <- rNA
    dat$epi$i.num.W <- rNA
    dat$epi$i.num.B.male <- rNA
    dat$epi$i.num.W.male <- rNA
    dat$epi$i.num.B.feml <- rNA
    dat$epi$i.num.W.feml <- rNA

    dat$epi$nBirths <- rNA
    dat$epi$dth.gen <- rNA
    dat$epi$dth.dis <- rNA

    dat$epi$incid <- rNA
  }

  dat$epi$num[at] <- length(status)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$num.B.male <- sum(race == "B" & male == 1, na.rm = TRUE)
  dat$epi$num.W.male <- sum(race == "W" & male == 1, na.rm = TRUE)
  dat$epi$num.B.feml <- sum(race == "B" & male == 0, na.rm = TRUE)
  dat$epi$num.W.feml <- sum(race == "W" & male == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.num.B.male <- sum(race == "B" & male == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.num.W.male <- sum(race == "W" & male == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.num.B.feml <- sum(race == "B" & male == 0 & status == 1, na.rm = TRUE)
  dat$epi$i.num.W.feml <- sum(race == "W" & male == 0 & status == 1, na.rm = TRUE)

  return(dat)
}
