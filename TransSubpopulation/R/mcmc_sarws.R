mcmc.sarws.biass <- function(dc, iteration, seed = NULL)
{
  require(data.table)
  require(Boom)	# dirichlet

  # set up mcmc objects
  mc <- list()
  mc$n <- iteration
  mc$pars <- list()
  mc$pars$LAMBDA <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = 1)		#prior for proportions
  mc$pars$S <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n)	#sampling probabilities (product)
  mc$pars$S_DTL <- data.table()														#sampling probabilities (detailed terms that make up product)
  mc$pars$Z <- matrix(NA_integer_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n) #augmented data
  mc$pars$NU <- NA_real_															#prior for N
  mc$pars$N <- matrix(NA_integer_, ncol = 1, nrow = mc$n)							#total number of counts on augmented data
  mc$pars$PI <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n)	#proportions
  mc$it.info <- data.table(
    IT = seq.int(1, mc$n),
    PAR_ID = rep(NA_integer_, mc$n),
    BLOCK = rep(NA_character_, mc$n),
    MHRATIO = rep(NA_real_, mc$n),
    ACCEPT = rep(NA_integer_, mc$n),
    LOG_LKL = rep(NA_real_, mc$n),
    LOG_PRIOR = rep(NA_real_, mc$n)
  )

  # define helper functions
  lddirichlet_vector <- function(x, nu) {
    ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
    stopifnot(is.finite(ans))
    ans
  }

  # initialise MCMC
  setkey(dc, COUNT_ID)
  mc$verbose <- 0L
  mc$curr.it <- 1L
  mc$seed <- seed
  set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
  set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
  set.seed(mc$seed)
  #	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
  #	(https://projecteuclid.org/euclid.ba/1422556416)
  mc$pars$LAMBDA[1,]	<- 0.8 / nrow(dc)
  # 	sampling: set to draw from prior
  tmp	<- dc[, list(TR_SEQ_P = rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),REC_SEQ_P = rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)),
            by = c('COUNT_ID','TR_P_SEQ_ALPHA','TR_P_SEQ_BETA','REC_P_SEQ_ALPHA','REC_P_SEQ_BETA')]
  tmp[, TR_SEQ_LOGD := dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log = TRUE)]
  tmp[, REC_SEQ_LOGD := dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log = TRUE)]
  mc$pars$S_DTL	<- copy(tmp)
  tmp	<- mc$pars$S_DTL[, list(S = TR_SEQ_P * REC_SEQ_P), by = 'COUNT_ID']
  setkey(tmp, COUNT_ID)
  mc$pars$S[1,] <- tmp$S
  #	augmented data: proposal draw under sampling probability
  #  mc$pars$Z[1,] <- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, mc$pars$S[1,])]
  mc$pars$Z[1,] <- dc$TR_OBS + rnbinom(length(dc$TR_OBS), dc$TR_OBS, mc$pars$S[1,])
  #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
  mc$pars$NU <- sum(dc$TR_OBS) / mean(dc$S)
  #	total count: that s just the sum of Z
  mc$pars$N[1,] <- sum(mc$pars$Z[1,])
  #	proportions: draw from full conditional
  mc$pars$PI[1,] <- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])
  #	store log likelihood
  tmp	<- sum(dbinom(
      dc$TR_OBS,
      size = mc$pars$Z[1,],
      prob = mc$pars$S[1,],
      log = TRUE)) +
    dmultinom(
      mc$pars$Z[1,],
      size = mc$pars$N[1,],
      prob = mc$pars$PI[1,],
      log = TRUE)
  set(mc$it.info, 1L, 'LOG_LKL', tmp)
  # 	store log prior
  tmp	<- dpois(mc$pars$N[1,], lambda = mc$pars$NU, log = TRUE) +
    lddirichlet_vector(mc$pars$PI[1,], nu = mc$pars$LAMBDA[1,]) +
    mc$pars$S_DTL[, sum(TR_SEQ_LOGD) + sum(REC_SEQ_LOGD)]
  set(mc$it.info, 1L, 'LOG_PRIOR', tmp)

  #
  # run mcmc
  options(warn = 0)
  mc$sweep <- ncol(mc$pars$Z) + 1L
  for (i in 1L:(mc$n - 1L))
  {
    mc$curr.it <- i
    # determine source-recipient combination that will be updated in this iteration
    update.count <- (i - 1L) %% mc$sweep + 1L
    # update S, Z, N for the source-recipient combination 'update.count'
    if (update.count < mc$sweep) {
      #	propose
      S.prop <- mc$pars$S[mc$curr.it,]
      S_DTL.prop <- copy(mc$pars$S_DTL)
      tmp <- dc[COUNT_ID == update.count, list(
        TR_SEQ_P = rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
        REC_SEQ_P = rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
      ),
      by = c(
        'COUNT_ID',
        'TR_P_SEQ_ALPHA',
        'TR_P_SEQ_BETA',
        'REC_P_SEQ_ALPHA',
        'REC_P_SEQ_BETA'
      )]
      tmp[, TR_SEQ_LOGD := dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log = TRUE)]
      tmp[, REC_SEQ_LOGD := dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log = TRUE)]
      S_DTL.prop <-
        rbind(subset(S_DTL.prop, COUNT_ID != update.count), tmp)
      setkey(S_DTL.prop, COUNT_ID)
      S.prop[update.count] <- tmp[, TR_SEQ_P * REC_SEQ_P]
      Z.prop <- mc$pars$Z[mc$curr.it,]
      Z.prop[update.count] <-
        dc$TR_OBS[update.count] + rnbinom(1, dc$TR_OBS[update.count], S.prop[update.count])
      #Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S.prop)]
      N.prop <- sum(Z.prop)
      #	calculate MH ratio
      log.prop.ratio <-
        sum(dnbinom(
          mc$pars$Z[mc$curr.it, update.count] - dc$TR_OBS[update.count],
          size = dc$TR_OBS[update.count],
          prob = mc$pars$S[mc$curr.it, update.count],
          log = TRUE
        )) -
        sum(dnbinom(
          Z.prop[update.count] - dc$TR_OBS[update.count],
          size = dc$TR_OBS[update.count],
          prob = S.prop[update.count],
          log = TRUE
        ))
      log.fc <-
        sum(dbinom(
          dc$TR_OBS[update.count],
          size = mc$pars$Z[mc$curr.it, update.count],
          prob = mc$pars$S[mc$curr.it, update.count],
          log = TRUE
        )) +
        dmultinom(mc$pars$Z[mc$curr.it,], prob = mc$pars$PI[mc$curr.it,], log =
                    TRUE) +
        dpois(mc$pars$N[mc$curr.it,], lambda = mc$pars$NU, log = TRUE)
      log.fc.prop <- sum(dbinom(
        dc$TR_OBS[update.count],
        size = Z.prop[update.count],
        prob = S.prop[update.count],
        log = TRUE
      )) +
        dmultinom(Z.prop, prob = mc$pars$PI[mc$curr.it,], log = TRUE) +
        dpois(N.prop, lambda = mc$pars$NU, log = TRUE)
      log.mh.ratio			<- log.fc.prop - log.fc + log.prop.ratio
      mh.ratio				<- min(1, exp(log.mh.ratio))
      #	update
      mc$curr.it				<- mc$curr.it + 1L
      set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
      set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
      set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
      set(mc$it.info,
          mc$curr.it,
          'ACCEPT',
          as.integer(runif(1) < mh.ratio))
      mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it - 1L,]
      if (mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
      {
        print(paste0('it ', mc$curr.it, ' ACCEPT S-Z-N block ', update.count))
      }
      if (mc$it.info[mc$curr.it, ACCEPT])
      {
        mc$pars$Z[mc$curr.it,]	<- Z.prop
        mc$pars$N[mc$curr.it,]	<- N.prop
        mc$pars$S[mc$curr.it,]	<- S.prop
        mc$pars$S_DTL			<- copy(S_DTL.prop)
      }
      if (mc$it.info[mc$curr.it, !ACCEPT])
      {
        mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it - 1L,]
        mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it - 1L,]
        mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it - 1L,]
      }
    }
    # update PI
    if (update.count == mc$sweep)
    {
      #	propose
      PI.prop					<-
        rdirichlet(1L, nu = mc$pars$Z[mc$curr.it,] + mc$pars$LAMBDA[1,])
      #	this is the full conditional of PI given S, N, Z
      #	always accept
      #	update
      mc$curr.it				<- mc$curr.it + 1L
      set(mc$it.info, mc$curr.it, 'BLOCK', 'PI')
      set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
      set(mc$it.info, mc$curr.it, 'MHRATIO', 1L)
      set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)
      mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it - 1L,]
      mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it - 1L,]
      mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it - 1L,]
      mc$pars$PI[mc$curr.it,]	<- PI.prop
    }
    # store log likelihood
    tmp	<-
      sum(dbinom(
        dc$TR_OBS,
        size = mc$pars$Z[mc$curr.it,],
        prob = mc$pars$S[mc$curr.it,],
        log = TRUE
      )) +
      dmultinom(
        mc$pars$Z[mc$curr.it,],
        size = mc$pars$N[mc$curr.it,],
        prob = mc$pars$PI[mc$curr.it,],
        log = TRUE
      )
    set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
    # store log prior
    tmp	<-
      dpois(mc$pars$N[mc$curr.it,], lambda = mc$pars$NU, log = TRUE) +
      lddirichlet_vector(mc$pars$PI[1,], nu = mc$pars$LAMBDA[1,]) +
      mc$pars$S_DTL[, sum(TR_SEQ_LOGD) + sum(REC_SEQ_LOGD)]
    set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)
  }
  save(dc, mc, iteration, file = 'core_inference_SNPIZ_mcmcEachCount.rda')
}

mcmc.sarws.biassap<-function(dc, iteration, seed = NULL)
  {
    require(data.table)
    require(Boom)	# dirichlet

    # set up mcmc objects
    mc <- list()
    mc$n <- iteration
    mc$pars <- list()
    mc$pars$LAMBDA <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = 1)		#prior for proportions
    mc$pars$S <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n)	#sampling probabilities (product)
    mc$pars$S_DTL <- data.table()														#sampling probabilities (detailed terms that make up product)
    mc$pars$Z <- matrix(NA_integer_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n) #augmented data
    mc$pars$NU <- NA_real_															#prior for N
    mc$pars$N <- matrix(NA_integer_, ncol = 1, nrow = mc$n)							#total number of counts on augmented data
    mc$pars$PI <- matrix(NA_real_, ncol = length(unique(dc$COUNT_ID)), nrow = mc$n)	#proportions
    mc$it.info <- data.table(
      IT = seq.int(1, mc$n),
      PAR_ID = rep(NA_integer_, mc$n),
      BLOCK = rep(NA_character_, mc$n),
      MHRATIO = rep(NA_real_, mc$n),
      ACCEPT = rep(NA_integer_, mc$n),
      LOG_LKL = rep(NA_real_, mc$n),
      LOG_PRIOR = rep(NA_real_, mc$n)
    )


    # define helper functions
    lddirichlet_vector	<- function(x, nu) {
      ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
      stopifnot(is.finite(ans))
      ans
    }

    # initialise MCMC
    setkey(dc, COUNT_ID)
    mc$verbose			<- 0L
    mc$curr.it			<- 1L
    mc$seed				<- seed
    set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
    set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
    set.seed(mc$seed)
    #	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
    #	(https://projecteuclid.org/euclid.ba/1422556416)
    mc$pars$LAMBDA[1, ]	<- 0.8 / nrow(dc)
    # 	sampling: set to draw from prior
    tmp					<-
      dc[, list(
        TR_PART_P = rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
        TR_SEQ_P = rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
        REC_PART_P = rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
        REC_SEQ_P = rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
      ),
      by = c(
        'COUNT_ID',
        'TR_P_PART_ALPHA',
        'TR_P_PART_BETA',
        'REC_P_PART_ALPHA',
        'REC_P_PART_BETA',
        'TR_P_SEQ_ALPHA',
        'TR_P_SEQ_BETA',
        'REC_P_SEQ_ALPHA',
        'REC_P_SEQ_BETA'
      )]
    tmp[, TR_PART_LOGD := dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log =
                                  TRUE)]
    tmp[, TR_SEQ_LOGD := dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log =
                                 TRUE)]
    tmp[, REC_PART_LOGD := dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log =
                                   TRUE)]
    tmp[, REC_SEQ_LOGD := dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log =
                                  TRUE)]
    mc$pars$S_DTL	<- copy(tmp)
    tmp				<-
      mc$pars$S_DTL[, list(S = TR_PART_P * TR_SEQ_P * REC_PART_P * REC_SEQ_P), by =
                      'COUNT_ID']
    setkey(tmp, COUNT_ID)
    mc$pars$S[1, ]		<- tmp$S
    #	augmented data: proposal draw under sampling probability
    #  mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, mc$pars$S[1,])]
    mc$pars$Z[1, ]		<-
      dc$TR_OBS + rnbinom(nrow(dc), dc$TR_OBS, mc$pars$S[1, ])
    #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
    mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
    #	total count: that s just the sum of Z
    mc$pars$N[1, ]		<- sum(mc$pars$Z[1, ])
    #	proportions: draw from full conditional
    mc$pars$PI[1, ]		<-
      rdirichlet(1, mc$pars$Z[1, ] + mc$pars$LAMBDA[1, ])
    #	store log likelihood
    tmp	<-
      sum(dbinom(
        dc$TR_OBS,
        size = mc$pars$Z[1, ],
        prob = mc$pars$S[1, ],
        log = TRUE
      )) +
      dmultinom(
        mc$pars$Z[1, ],
        size = mc$pars$N[1, ],
        prob = mc$pars$PI[1, ],
        log = TRUE
      )
    set(mc$it.info, 1L, 'LOG_LKL', tmp)
    # 	store log prior
    tmp	<- dpois(mc$pars$N[1, ], lambda = mc$pars$NU, log = TRUE) +
      lddirichlet_vector(mc$pars$PI[1, ], nu = mc$pars$LAMBDA[1, ]) +
      mc$pars$S_DTL[, sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]
    set(mc$it.info, 1L, 'LOG_PRIOR', tmp)

    #
    # run mcmc
    options(warn = 0)
    mc$sweep			<- ncol(mc$pars$Z) + 1L
    for (i in 1L:(mc$n - 1L))
    {
      mc$curr.it		<- i
      # determine source-recipient combination that will be updated in this iteration
      update.count	<- (i - 1L) %% mc$sweep + 1L
      # update S, Z, N for the source-recipient combination 'update.count'
      if (update.count < mc$sweep)
      {
        #	propose
        S.prop					<- mc$pars$S[mc$curr.it, ]
        S_DTL.prop				<- copy(mc$pars$S_DTL)
        tmp						<-
          dc[COUNT_ID == update.count, list(
            TR_PART_P = rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
            TR_SEQ_P = rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
            REC_PART_P = rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
            REC_SEQ_P = rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
          ),
          by = c(
            'COUNT_ID',
            'TR_P_PART_ALPHA',
            'TR_P_PART_BETA',
            'REC_P_PART_ALPHA',
            'REC_P_PART_BETA',
            'TR_P_SEQ_ALPHA',
            'TR_P_SEQ_BETA',
            'REC_P_SEQ_ALPHA',
            'REC_P_SEQ_BETA'
          )]
        tmp[, TR_PART_LOGD := dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log =
                                      TRUE)]
        tmp[, TR_SEQ_LOGD := dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log =
                                     TRUE)]
        tmp[, REC_PART_LOGD := dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log =
                                       TRUE)]
        tmp[, REC_SEQ_LOGD := dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log =
                                      TRUE)]
        S_DTL.prop				<-
          rbind(subset(S_DTL.prop, COUNT_ID != update.count), tmp)
        setkey(S_DTL.prop, COUNT_ID)
        S.prop[update.count]	<-
          tmp[, TR_PART_P * TR_SEQ_P * REC_PART_P * REC_SEQ_P]
        Z.prop					<- mc$pars$Z[mc$curr.it, ]
        #      Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S.prop)]
        Z.prop[update.count]	<-
          dc$TR_OBS[update.count] + rnbinom(1, dc$TR_OBS[update.count], S.prop[update.count])

        N.prop					<- sum(Z.prop)
        #	calculate MH ratio
        log.prop.ratio			<-
          sum(dnbinom(
            mc$pars$Z[mc$curr.it, update.count] - dc$TR_OBS[update.count],
            size = dc$TR_OBS[update.count],
            prob = mc$pars$S[mc$curr.it, update.count],
            log = TRUE
          )) -
          sum(dnbinom(
            Z.prop[update.count] - dc$TR_OBS[update.count],
            size = dc$TR_OBS[update.count],
            prob = S.prop[update.count],
            log = TRUE
          ))
        log.fc					<-
          sum(dbinom(
            dc$TR_OBS[update.count],
            size = mc$pars$Z[mc$curr.it, update.count],
            prob = mc$pars$S[mc$curr.it, update.count],
            log = TRUE
          )) +
          dmultinom(mc$pars$Z[mc$curr.it, ], prob = mc$pars$PI[mc$curr.it, ], log =
                      TRUE) +
          dpois(mc$pars$N[mc$curr.it, ], lambda = mc$pars$NU, log = TRUE)
        log.fc.prop				<-
          sum(dbinom(
            dc$TR_OBS[update.count],
            size = Z.prop[update.count],
            prob = S.prop[update.count],
            log = TRUE
          )) +
          dmultinom(Z.prop, prob = mc$pars$PI[mc$curr.it, ], log = TRUE) +
          dpois(N.prop, lambda = mc$pars$NU, log = TRUE)
        log.mh.ratio			<- log.fc.prop - log.fc + log.prop.ratio
        mh.ratio				<- min(1, exp(log.mh.ratio))
        #	update
        mc$curr.it				<- mc$curr.it + 1L
        set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
        set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
        set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
        set(mc$it.info,
            mc$curr.it,
            'ACCEPT',
            as.integer(runif(1) < mh.ratio))
        mc$pars$PI[mc$curr.it, ]	<- mc$pars$PI[mc$curr.it - 1L, ]
        if (mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
        {
          print(paste0('it ', mc$curr.it, ' ACCEPT S-Z-N block ', update.count))
        }
        if (mc$it.info[mc$curr.it, ACCEPT])
        {
          mc$pars$Z[mc$curr.it, ]	<- Z.prop
          mc$pars$N[mc$curr.it, ]	<- N.prop
          mc$pars$S[mc$curr.it, ]	<- S.prop
          mc$pars$S_DTL			<- copy(S_DTL.prop)
        }
        if (mc$it.info[mc$curr.it,!ACCEPT])
        {
          mc$pars$Z[mc$curr.it, ]	<- mc$pars$Z[mc$curr.it - 1L, ]
          mc$pars$N[mc$curr.it, ]	<- mc$pars$N[mc$curr.it - 1L, ]
          mc$pars$S[mc$curr.it, ]	<- mc$pars$S[mc$curr.it - 1L, ]
        }
      }
      # update PI
      if (update.count == mc$sweep)
      {
        #	propose
        PI.prop					<-
          rdirichlet(1L, nu = mc$pars$Z[mc$curr.it, ] + mc$pars$LAMBDA[1, ])
        #	this is the full conditional of PI given S, N, Z
        #	always accept
        #	update
        mc$curr.it				<- mc$curr.it + 1L
        set(mc$it.info, mc$curr.it, 'BLOCK', 'PI')
        set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
        set(mc$it.info, mc$curr.it, 'MHRATIO', 1L)
        set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)
        mc$pars$S[mc$curr.it, ]	<- mc$pars$S[mc$curr.it - 1L, ]
        mc$pars$Z[mc$curr.it, ]	<- mc$pars$Z[mc$curr.it - 1L, ]
        mc$pars$N[mc$curr.it, ]	<- mc$pars$N[mc$curr.it - 1L, ]
        mc$pars$PI[mc$curr.it, ]	<- PI.prop
      }
      # store log likelihood
      tmp	<-
        sum(dbinom(
          dc$TR_OBS,
          size = mc$pars$Z[mc$curr.it, ],
          prob = mc$pars$S[mc$curr.it, ],
          log = TRUE
        )) +
        dmultinom(
          mc$pars$Z[mc$curr.it, ],
          size = mc$pars$N[mc$curr.it, ],
          prob = mc$pars$PI[mc$curr.it, ],
          log = TRUE
        )
      set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
      # store log prior
      tmp	<-
        dpois(mc$pars$N[mc$curr.it, ], lambda = mc$pars$NU, log = TRUE) +
        lddirichlet_vector(mc$pars$PI[1, ], nu = mc$pars$LAMBDA[1, ]) +
        mc$pars$S_DTL[, sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]
      set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)
    }
    save(dc, mc, iteration, file = 'core_inference_SNPIZ_mcmcEachCount.rda')
}




mcmc.sarws.diagnostics <- function(mc, iteration, burnin) {
  require(coda)
  require(data.table)
  diagnostic <- list()

  diagnostic$acc.rate <- matrix(NA_real_, nrow = 1, ncol = mc$sweep)
  for (i in 1L:mc$sweep) {
    diagnostic$acc.rate[1, i] <- mean(mc$it.info[IT > burnin & PAR_ID == i, ACCEPT])
  }

  #sort(diagnostic$acc.rate)

  it.pi <- mc$it.info[BLOCK == 'PI', IT]
  it.pi <- it.pi[it.pi > burnin]

  S <- mc$pars$S[it.pi, ]
  Z <- mc$pars$Z[it.pi, ]
  N <- mc$pars$N[it.pi, ]
  PI <- mc$pars$PI[it.pi, ]

  # effetive sampling size
  eff.size.PI <- matrix(NA_real_, nrow = 1, ncol = ncol(mc$pars$Z))
  for (i in 1:ncol(mc$pars$Z)) {
    eff.size.PI[i] <- effectiveSize(PI[, i])
  }

  sort.effsize.PI <- sort(eff.size.PI, index.return = TRUE)
  if (length(sort.effsize.PI$x) < 6) {
    worst.effsize.PI.index <- sort.effsize.PI$ix
    worst.effsize.PI <- sort.effsize.PI$x
    print(
      paste0(
        "The worst effective sample sizes for PI ",
        worst.effsize.PI.index,
        " are ",
        worst.effsize.PI,
        " out of ",
        length(PI[, 1]),
        " samples."
      )
    )
  } else{
    worst.effsize.PI.index <- sort.effsize.PI$ix[1:6]
    worst.effsize.PI <- sort.effsize.PI$x[1:6]
    print(
      paste0(
        "The worst effective sample sizes for PI ",
        worst.effsize.PI.index,
        " are ",
        worst.effsize.PI,
        " out of ",
        length(PI[, 1]),
        " samples."
      )
    )
  }

  diagnostic$PI$eff.size.index <- worst.effsize.PI.index
  diagnostic$PI$eff.size <- worst.effsize.PI

  # traceplot
  png("traceplot_PI.png", width = 1400, height = 800)
  par(ps = 12,
      cex = 1,
      cex.main = 2.5)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.PI.index) {
    traceplot(
      mcmc(PI[, i]),
      main = paste('Traceplot of Pi', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  # autocorrelation plot:
  png("autocorrelation_PI.png",
      width = 1400,
      height = 800)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.PI.index) {
    acf(mcmc(PI[, i]), main = '')
    mtext(
      sprintf(paste('Autocorrelation plot of PI', i)),
      side = 3,
      cex = 2,
      line = 1,
      font = 2
    )
    #acf(mcmc(PI[,i]),main=paste('Autocorrelation plot of Pi',i),ylab='',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # density plot:
  png("density_PI.png", width = 1400, height = 800)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.PI.index) {
    densplot(
      mcmc(PI[, i]),
      main = paste('Density plot of Pi', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  #########################################S###################################################
  # effetive sampling size
  eff.size.S <- matrix(NA_real_, nrow = 1, ncol = ncol(mc$pars$S))
  for (i in 1:ncol(mc$pars$S)) {
    eff.size.S[i] <- effectiveSize(S[, i])
  }

  sort.effsize.S <- sort(eff.size.S, index.return = TRUE)
  if (length(sort.effsize.S$x) < 6) {
    worst.effsize.S.index <- sort.effsize.S$ix
    worst.effsize.S <- sort.effsize.S$x
    print(
      paste0(
        "The worst effective sample sizes for S ",
        worst.effsize.S.index,
        " are ",
        worst.effsize.S,
        " out of ",
        length(S[, 1]),
        " samples."
      )
    )
  } else{
    worst.effsize.S.index <- sort.effsize.S$ix[1:6]
    worst.effsize.S <- sort.effsize.S$x[1:6]
    print(
      paste0(
        "The worst effective sample sizes for S are ",
        worst.effsize.S.index,
        " are ",
        worst.effsize.S,
        " out of ",
        length(S[, 1]),
        " samples."
      )
    )
  }

  diagnostic$S$eff.size.index <- worst.effsize.S.index
  diagnostic$S$eff.size <- worst.effsize.S

  # traceplot
  png("traceplot_S.png", width = 1400, height = 800)
  par(ps = 12,
      cex = 1,
      cex.main = 2.5)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.S.index) {
    traceplot(
      mcmc(S[, i]),
      main = paste('Traceplot of S', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  # autocorrelation plot:
  png("autocorrelation_S.png",
      width = 1400,
      height = 800)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.S.index) {
    acf(mcmc(S[, i]), main = '')
    mtext(
      sprintf(paste('Autocorrelation plot of S', i)),
      side = 3,
      cex = 2,
      line = 1,
      font = 2
    )
    #acf(mcmc(S[,i]),main=paste('Autocorrelation plot of S',i),ylab='',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  png("density_S.png", width = 1400, height = 800)
  # density plot:
  par(mfrow = c(2, 3))
  for (i in worst.effsize.S.index) {
    densplot(
      mcmc(S[, i]),
      main = paste('Density plot of S', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  #########################################Z###################################################
  # effetive sampling size
  eff.size.Z <- matrix(NA_real_, nrow = 1, ncol = ncol(mc$pars$Z))
  for (i in 1:ncol(mc$pars$Z)) {
    eff.size.Z[i] <- effectiveSize(Z[, i])
  }

  sort.effsize.Z <- sort(eff.size.Z, index.return = TRUE)
  if (length(sort.effsize.Z$x) < 6) {
    worst.effsize.Z.index <- sort.effsize.Z$ix
    worst.effsize.Z <- sort.effsize.Z$x
    print(
      paste0(
        "The worst effective sample sizes for Z are ",
        worst.effsize.Z.index,
        " are ",
        worst.effsize.Z,
        " out of ",
        length(Z[, 1]),
        " samples."
      )
    )
  } else{
    worst.effsize.Z.index <- sort.effsize.Z$ix[1:6]
    worst.effsize.Z <- sort.effsize.Z$x[1:6]
    print(
      paste0(
        "The worst effective sample sizes for Z are ",
        worst.effsize.Z.index,
        " are ",
        worst.effsize.Z,
        " out of ",
        length(Z[, 1]),
        " samples."
      )
    )
  }

  diagnostic$Z$eff.size.index <- worst.effsize.Z.index
  diagnostic$Z$eff.size <- worst.effsize.Z

  # traceplot
  png("traceplot_Z.png", width = 1400, height = 800)
  par(ps = 12,
      cex = 1,
      cex.main = 2.5)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.Z.index) {
    traceplot(
      mcmc(Z[, i]),
      main = paste('Traceplot of Z', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  # autocorrelation plot:
  png("autocorrelation_Z.png",
      width = 1400,
      height = 800)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.Z.index) {
    acf(mcmc(Z[, i]), main = '')
    mtext(
      sprintf(paste('Autocorrelation plot of Z', i)),
      side = 3,
      cex = 2,
      line = 1,
      font = 2
    )
    #acf(mcmc(Z[,i]),main=paste('Autocorrelation plot of Z',i),ylab='',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # density plot:
  png("density_Z.png", width = 1400, height = 800)
  par(mfrow = c(2, 3))
  for (i in worst.effsize.Z.index) {
    densplot(
      mcmc(Z[, i]),
      main = paste('Density plot of Z', i),
      cex.lab = 2.5,
      cex.axis = 1.5,
      cex.main = 2.5,
      cex.sub = 2.5
    )
  }
  dev.off()

  #########################################N###################################################
  # effetive sampling size
  eff.size.N <- effectiveSize(N)
  diagnostic$N$eff.size <- eff.size.N

  print(paste0("The worst effective sample size is ", eff.size.N))


  # traceplot
  png("traceplot_N.png", width = 600, height = 500)
  par(ps = 12,
      cex = 1,
      cex.main = 2.5)
  par(mfrow = c(1, 1))
  traceplot(
    mcmc(N),
    main = 'Traceplot of N',
    cex.lab = 2.5,
    cex.axis = 1.5,
    cex.main = 2.5,
    cex.sub = 2.5
  )
  dev.off()

  # autocorrelation plot:
  png("autocorrelation_N.png",
      width = 600,
      height = 500)
  par(mfrow = c(1, 1))
  acf(mcmc(N), main = '')
  mtext(
    sprintf('Autocorrelation plot of N'),
    side = 3,
    cex = 2,
    line = 1,
    font = 2
  )
  dev.off()

  # density plot:
  png("density_N.png", width = 600, height = 500)
  par(mfrow = c(1, 1))
  densplot(
    mcmc(N),
    main = 'Density plot of N',
    cex.lab = 2.5,
    cex.axis = 1.5,
    cex.main = 2.5,
    cex.sub = 2.5
  )
  dev.off()
  save(diagnostic, file = "core_inference_diagnostic.rda")
}
