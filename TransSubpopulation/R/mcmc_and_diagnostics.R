mcmc.core.inference<-function(s.dtl.prior,tr.obs,seed){
  library(data.table)
  library(gtools)
  ptm<-Sys.time()

  set.seed(seed)
  TR_OBS<-tr.obs$OBS
  pair<-cbind(tr.obs$TR_CATEGORY_ID,tr.obs$REC_CATEGORY_ID)

  update.dim.list<-list()
  for (i in 1:ncol(s.dtl.prior$SAM_P)){
    tmp<-which(pair[,1]==i|pair[,2]==i)
    update.dim.list[[i]]<-tmp
  }

  # set up mcmc objects
  mc<-list()
  nsweep<-nrow(s.dtl.prior$SAM_P) # number of samples after thining
  mc$n<-(nrow(s.dtl.prior$SAM_P)-1)*(ncol(s.dtl.prior$SAM_P)+1)+1 # the number of iterations

  mc$pars<-list()
  mc$pars$LAMBDA<-matrix(NA_real_, ncol=length(TR_OBS), nrow=1)		#prior for proportions

  mc$pars$SAM_P<-matrix(NA_real_,ncol = ncol(s.dtl.prior$SAM_P),nrow=nsweep)
  mc$pars$SAM_P_LOGD<-matrix(NA_real_,ncol = ncol(s.dtl.prior$SAM_P),nrow=nsweep)
  mc$pars$S<-matrix(NA_integer_, ncol=length(TR_OBS), nrow=nsweep)
  mc$pars$Z<-matrix(NA_integer_, ncol=length(TR_OBS), nrow=nsweep) #augmented data
  mc$pars$NU<-NA_real_															#prior for N
  mc$pars$N<-matrix(NA_integer_, ncol=1, nrow=nsweep)							#total number of counts on augmented data
  mc$pars$PI<-matrix(NA_real_, ncol=length(TR_OBS), nrow=nsweep)	#proportions

  mc$it.info<-data.table(	IT= seq.int(1,nsweep),
                          LOG_LKL=rep(NA_real_, nsweep),
                          LOG_PRIOR=rep(NA_real_, nsweep))

  mc$z.acc.count<-matrix(NA_integer_, ncol=1, nrow=length(TR_OBS))
  mc$xi.acc.count<-matrix(NA_integer_, ncol=1, nrow=ncol(s.dtl.prior$SAM_P))
  mc$z.total.count<-matrix(NA_integer_, ncol=1, nrow=length(TR_OBS))
  mc$xi.total.count<-matrix(NA_integer_, ncol=1, nrow=ncol(s.dtl.prior$SAM_P))

  # define helper functions
  lddirichlet_vector	<- function(x, nu){
    ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
    stopifnot(is.finite(ans))
    ans
  }


  # initialise MCMC
  #  mc$verbose			<- 0L
  mc$curr.it			<- 1L
  #mc$seed				<- seed
  #set.seed( mc$seed )

  #	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
  #	(https://projecteuclid.org/euclid.ba/1422556416)
  mc$pars$LAMBDA[1,]	<- 0.8/length(TR_OBS)

  mc$pars$SAM_P[1,] <- s.dtl.prior$SAM_P[1,]
  mc$pars$SAM_P_LOGD[1,] <- s.dtl.prior$SAM_P_LOGD[1,]
  mc$pars$S[1,]	<- mc$pars$SAM_P[1,][pair[,1]]*mc$pars$SAM_P[1,][pair[,2]]
  #	augmented data: proposal draw under sampling probability
  mc$pars$Z[1,]		<- TR_OBS + rnbinom(length(TR_OBS),TR_OBS,mc$pars$S[1,])
  #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
  mc$pars$NU			<- sum(TR_OBS) / mean(mc$pars$S[1,])
  #	total count: that s just the sum of Z
  mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
  #	proportions: draw from full conditional
  mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])


  #	store log likelihood
  tmp	<- sum( dbinom(TR_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
    dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
  set(mc$it.info, 1L, 'LOG_LKL', tmp)
  # 	store log prior
  tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
    lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
    sum(mc$pars$SAM_P_LOGD[1,])
  set(mc$it.info, 1L, 'LOG_PRIOR', tmp)

  mc$z.acc.count<-rep(0,length(TR_OBS))
  mc$xi.acc.count<-rep(0,ncol(s.dtl.prior$SAM_P))
  mc$z.total.count<-rep(0,length(TR_OBS))
  mc$xi.total.count<-rep(0,ncol(s.dtl.prior$SAM_P))

  # parameter value at the current step
  SAM_P_CURR<-mc$pars$SAM_P[1,]
  SAM_P_LOGD_CURR<-mc$pars$SAM_P_LOGD[1,]
  S_CURR<-mc$pars$S[1,]
  PI_CURR<-mc$pars$PI[1,]
  N_CURR<-mc$pars$N[1,]
  Z_CURR<-mc$pars$Z[1,]

  # run mcmc
  options(warn=0)
  mc$sweep<-ncol(s.dtl.prior$SAM_P) + 1L
  for(i in 1L:(mc$n-1L)){
    mc$curr.it		<- i
    # determine source-recipient combination that will be updated in this iteration
    update.count	<- (i-1L) %% mc$sweep + 1L
    update.round <- (i-1L) %/% mc$sweep + 1L
    # update S, Z, N for the source-recipient combination 'update.count'
    if(update.count<mc$sweep)
    {
      update.dim<-update.dim.list[[update.count]]
      # update S
      SAM_P.prop<-SAM_P_CURR
      SAM_P.prop[update.count]<-s.dtl.prior$SAM_P[update.round+1L,update.count]
      SAM_P_LOGD.prop<-SAM_P_LOGD_CURR
      SAM_P_LOGD.prop[update.count]<-s.dtl.prior$SAM_P_LOGD[update.round+1L,update.count]

      S.prop<-S_CURR
      S.prop[update.dim]<-SAM_P.prop[pair[update.dim,1]]*SAM_P.prop[pair[update.dim,2]]

      #	calculate MH ratio
      log.fc					<- sum(dbinom(TR_OBS[update.dim], size=Z_CURR[update.dim], prob=S_CURR[update.dim], log=TRUE))
      log.fc.prop				<- sum(dbinom(TR_OBS[update.dim], size=Z_CURR[update.dim], prob=S.prop[update.dim], log=TRUE))
      log.mh.ratio			<- log.fc.prop - log.fc

      mh.ratio				<- min(1,exp(log.mh.ratio))

      mc$xi.total.count[update.count]<-mc$xi.total.count[update.count]+1

      #	update
      if(runif(1) < mh.ratio)
      {
        mc$xi.acc.count[update.count]<-mc$xi.acc.count[update.count]+1
        SAM_P_CURR<-SAM_P.prop
        SAM_P_LOGD_CURR<-SAM_P_LOGD.prop
        S_CURR<- S.prop
      }


      # update z for update.dim
      Z.prop					<- Z_CURR

      for (j in 1:length(update.dim)){
        Z.prop[update.dim[j]]<-TR_OBS[update.dim[j]]+rnbinom(1L, TR_OBS[update.dim[j]], S_CURR[update.dim[j]])
        N.prop<-sum(Z.prop)
        log.prop.ratio <- sum(dnbinom(Z_CURR[update.dim[j]]-TR_OBS[update.dim[j]], size=TR_OBS[update.dim[j]], prob= S_CURR[update.dim[j]], log=TRUE)) -
          sum(dnbinom(Z.prop[update.dim[j]]-TR_OBS[update.dim[j]], size=TR_OBS[update.dim[j]], prob=S_CURR[update.dim[j]], log=TRUE))
        log.fc	<- sum(dbinom(TR_OBS[update.dim[j]], size=Z_CURR[update.dim[j]], prob=S_CURR[update.dim[j]], log=TRUE)) +
          dmultinom(Z_CURR, prob=PI_CURR, log=TRUE) +
          dpois(N_CURR, lambda=mc$pars$NU, log=TRUE)
        log.fc.prop				<- sum(dbinom(TR_OBS[update.dim[j]], size=Z.prop[update.dim[j]], prob=S_CURR[update.dim[j]], log=TRUE)) +
          dmultinom(Z.prop, prob=PI_CURR, log=TRUE) +
          dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
        log.mh.ratio			<- log.fc.prop - log.fc + log.prop.ratio
        mh.ratio				<- min(1,exp(log.mh.ratio))
        #	update

        mc$z.total.count[update.dim[j]]<-mc$z.total.count[update.dim[j]]+1

        if(runif(1) < mh.ratio)
        {
          Z_CURR	<- Z.prop
          N_CURR	<- N.prop
          mc$z.acc.count[update.dim[j]]<-mc$z.acc.count[update.dim[j]]+1
        }
      }
    }
    # update PI
    if(update.count==mc$sweep)
    {
      #	propose
      PI.prop					<- rdirichlet(1L, Z_CURR + mc$pars$LAMBDA[1,])
      #	this is the full conditional of PI given S, N, Z
      #	always accept
      #	update
      PI_CURR	<- PI.prop

      # store values
      mc$pars$SAM_P[update.round+1L,]<-SAM_P_CURR
      mc$pars$SAM_P_LOGD[update.round+1L,]<-SAM_P_LOGD_CURR
      mc$pars$Z[update.round+1L,]<-Z_CURR
      mc$pars$S[update.round+1L,]<-S_CURR
      mc$pars$N[update.round+1L,]<-N_CURR
      mc$pars$PI[update.round+1L,]<-PI_CURR

      tmp	<- sum(dbinom(TR_OBS, size=Z_CURR, prob=S_CURR, log=TRUE) ) +
        dmultinom(Z_CURR, size=N_CURR, prob=PI_CURR, log=TRUE)
      set(mc$it.info, update.round+1L, 'LOG_LKL', tmp)
      # 	store log prior
      tmp	<- dpois(N_CURR, lambda=mc$pars$NU, log=TRUE) +
        lddirichlet_vector(PI_CURR, nu=mc$pars$LAMBDA[1,]) +
        sum(SAM_P_LOGD_CURR)
      set(mc$it.info, update.round+1L, 'LOG_PRIOR', tmp)
    }
  }
  mc$method<-s.dtl.prior$method
  time<-Sys.time()-ptm
  save(time,mc,file=paste0('core_inference_SNPIZ_mcmcEachCount_',mc$method,'.rda'))
}

mcmc.core.inference.diagnostics<-function(mc,burnin){
  library(coda)
  library(data.table)
  diagnostic<-list()

  # take indexes
  it<-1:nrow(mc$pars$S)
  it.rm.burnin<-it[it>burnin]

  S<-mc$pars$S[it.rm.burnin,]
  Z<-mc$pars$Z[it.rm.burnin,]
  N<-mc$pars$N[it.rm.burnin,]
  PI<-mc$pars$PI[it.rm.burnin,]

  diagnostic$xi.acc.rate<-mc$xi.acc.count/mc$xi.total.count
  diagnostic$xi.acc.rate.sorted<-sort(diagnostic$xi.acc.rate,index.return=TRUE)
  diagnostic$z.acc.rate<-mc$z.acc.count/mc$z.total.count
  diagnostic$z.acc.rate.sorted<-sort(diagnostic$z.acc.rate,index.return=TRUE)


  # effetive sampling size for PI
  eff.size.PI<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$Z))
  for(i in 1:ncol(mc$pars$Z)){
    eff.size.PI[i]<-effectiveSize(PI[,i])
  }

  sort.effsize.PI<-sort(eff.size.PI,index.return=TRUE)
  if (length(sort.effsize.PI$x)<6){
    worst.effsize.PI.index<-sort.effsize.PI$ix
    worst.effsize.PI<-sort.effsize.PI$x
    print(paste0("The worst effective sample sizes for PI ", worst.effsize.PI.index, " are ", worst.effsize.PI," out of ", length(PI[,1])," samples."))
  }else{
    worst.effsize.PI.index<-sort.effsize.PI$ix[1:6]
    worst.effsize.PI<-sort.effsize.PI$x[1:6]
    print(paste0("The worst effective sample sizes for PI ", worst.effsize.PI.index, " are ", worst.effsize.PI," out of ", length(PI[,1])," samples."))
  }

  diagnostic$PI$eff.size.index<-worst.effsize.PI.index
  diagnostic$PI$eff.size<-worst.effsize.PI

  # traceplot for PI
  png(paste0('traceplot_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    coda::traceplot(mcmc(PI[,i]),main=paste('Traceplot of Pi',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for PI
  png(paste0('autocorrelation_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    acf(mcmc(PI[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of PI',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for PI
  png(paste0('density_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    densplot(mcmc(PI[,i]),main=paste('Density plot of Pi',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # effetive sampling size for S
  eff.size.S<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$S))
  for(i in 1:ncol(mc$pars$S)){
    eff.size.S[i]<-effectiveSize(S[,i])
  }

  sort.effsize.S<-sort(eff.size.S,index.return=TRUE)
  if (length(sort.effsize.S$x)<6){
    worst.effsize.S.index<-sort.effsize.S$ix
    worst.effsize.S<-sort.effsize.S$x
    print(paste0("The worst effective sample sizes for S ", worst.effsize.S.index, " are ", worst.effsize.S," out of ", length(S[,1])," samples."))
  }else{
    worst.effsize.S.index<-sort.effsize.S$ix[1:6]
    worst.effsize.S<-sort.effsize.S$x[1:6]
    print(paste0("The worst effective sample sizes for S are ", worst.effsize.S.index, " are ", worst.effsize.S," out of ", length(S[,1])," samples."))
  }

  diagnostic$S$eff.size.index<-worst.effsize.S.index
  diagnostic$S$eff.size<-worst.effsize.S

  # traceplot for S
  png(paste0('traceplot_S_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    coda::traceplot(mcmc(S[,i]),main=paste('Traceplot of S',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for S
  png(paste0('autocorrelation_S_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    acf(mcmc(S[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of S',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for S
  png(paste0('density_S_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    densplot(mcmc(S[,i]),main=paste('Density plot of S',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # effetive sampling size for Z
  eff.size.Z<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$Z))
  for(i in 1:ncol(mc$pars$Z)){
    eff.size.Z[i]<-effectiveSize(Z[,i])
  }

  sort.effsize.Z<-sort(eff.size.Z,index.return=TRUE)
  if (length(sort.effsize.Z$x)<6){
    worst.effsize.Z.index<-sort.effsize.Z$ix
    worst.effsize.Z<-sort.effsize.Z$x
    print(paste0("The worst effective sample sizes for Z are ", worst.effsize.Z.index, " are ", worst.effsize.Z," out of ", length(Z[,1])," samples."))
  }else{
    worst.effsize.Z.index<-sort.effsize.Z$ix[1:6]
    worst.effsize.Z<-sort.effsize.Z$x[1:6]
    print(paste0("The worst effective sample sizes for Z are ", worst.effsize.Z.index, " are ", worst.effsize.Z," out of ", length(Z[,1])," samples."))
  }

  diagnostic$Z$eff.size.index<-worst.effsize.Z.index
  diagnostic$Z$eff.size<-worst.effsize.Z


  # traceplot for Z
  png(paste0('traceplot_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    coda::traceplot(mcmc(Z[,i]),main=paste('Traceplot of Z',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for Z
  png(paste0('autocorrelation_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    acf(mcmc(Z[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of Z',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for Z:
  png(paste0('density_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    densplot(mcmc(Z[,i]),main=paste('Density plot of Z',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  #########################################N###################################################
  # effetive sampling size for N
  eff.size.N<-effectiveSize(N)
  diagnostic$N$eff.size<-eff.size.N

  print(paste0("The worst effective sample size is ",eff.size.N))


  # traceplot for N
  png(paste0('traceplot_N_',mc$method,'.png'), width = 600, height = 500)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(1,1))
  coda::traceplot(mcmc(N),main='Traceplot of N',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  dev.off()

  # autocorrelation plot for N
  png(paste0('autocorrelation_N_',mc$method,'.png'), width = 600, height = 500)
  par(mfrow=c(1,1))
  acf(mcmc(N),main='')
  mtext(sprintf('Autocorrelation plot of N'), side=3,cex=2,line = 1,font=2)
  dev.off()

  # density plot for N
  png(paste0('density_N_',mc$method,'.png'), width = 600, height = 500)
  par(mfrow=c(1,1))
  densplot(mcmc(N),main='Density plot of N',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  dev.off()

  save(diagnostic,file=paste0('core_inference_diagnostic_',mc$method,'.rda'))
}
