GP1D.Gaussian <- function()
{
	home <- '~/git/phyloscanner/phyloflows/inst/misc'
	dir.mb.tutorial <- file.path(home,'MB_gp_tutorial')
	#
	#	this is Mike Betancourts GP tutorial 1-3
	#	with some modifications
	#	https://github.com/betanalpha/knitr_case_studies/tree/master/gaussian_processes
	#
	library(rstan)
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	source(file.path(dir.mb.tutorial,"gp_utility.R"))
	source(file.path(dir.mb.tutorial,"stan_utility.R"))
	
	alpha_true <- 3
	rho_true <- 5.5
	sigma_true <- 2	
	N_total <- 501
	x_total <- 20 * (0:(N_total - 1)) / (N_total - 1) - 10	
	simu_data <- list(	alpha=alpha_true, 
						rho=rho_true, 
						sigma=sigma_true,
						N=N_total, x=x_total)
		
	#	simulate data set	
	simu_fit <- stan(	file=file.path(dir.mb.tutorial,"simu_gauss.stan"), 
						data=simu_data, iter=1,
						chains=1, seed=494838, algorithm="Fixed_param")
	f_total <- extract(simu_fit)$f[1,]
	y_total <- extract(simu_fit)$y[1,]	
	true_realization <- data.frame(f_total, x_total)
	names(true_realization) <- c("f_total", "x_total")	
	observed_idx <- c(50*(0:10)+1)
	N = length(observed_idx)
	x <- x_total[observed_idx]
	y <- y_total[observed_idx]	
	#	plot simulated data
	plot(x_total, f_total, type="l", lwd=2, xlab="x", ylab="y",
			xlim=c(-10, 10), ylim=c(-10, 10))
	points(x_total, y_total, col="white", pch=16, cex=0.6)
	points(x_total, y_total, col=c_mid_teal, pch=16, cex=0.4)
	points(x, y, col="white", pch=16, cex=1.2)
	points(x, y, col="black", pch=16, cex=0.8)
	#	save simulated data and truth
	N_predict <- N_total
	x_predict <- x_total
	y_predict <- y_total	
	sample_idx <- observed_idx
	stan_rdump(c("N", "x", "y",
					"N_predict", "x_predict", "y_predict",
					"sample_idx"), file=file.path(dir.mb.tutorial,"gp.data.R"))
	data <- read_rdump(file.path(dir.mb.tutorial,"gp.data.R"))	
	stan_rdump(c("f_total", "x_total", "sigma_true"), file=file.path(dir.mb.tutorial,"gp.truth.R"))
	#	simulate 1e3 times and plot data generating process quantiles
	f_data <- list(sigma=sigma_true, N=N_total, f=f_total)
	dgp_fit <- stan(file=file.path(dir.mb.tutorial,'simu_gauss_dgp.stan'), 
					data=f_data, iter=1000, warmup=0,
					chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")	
	plot_gp_pred_quantiles(dgp_fit, data, true_realization,
			"True Data Generating Process Quantiles")
	
	
	#
	#	simulate from GP conditional on Gaussian observations
	pred_data <- list(	alpha=alpha_true, rho=rho_true, 
						sigma=sigma_true, N=N, x=x, y=y,
						N_predict=N_predict, x_predict=x_predict
						)
	pred_fit <- stan(file=file.path(dir.mb.tutorial,'predict_gauss.stan'), 
						data=pred_data, iter=1000, warmup=0,
						chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")
	plot_gp_quantiles(pred_fit, data, true_realization,"Posterior Quantiles")
	
	#
	#	estimate GP and hyperparameters with informative pior on length scale
	fit <- stan(file=file.path(dir.mb.tutorial,'gp3.stan'), data=data, seed=5838298)
	#	about 20 seconds for each of 4 chains
	check_all_diagnostics(fit)
	params <- extract(fit)
	#	plot 1D posteriors and true param
	par(mfrow=c(1, 3))	
	alpha_breaks=10 * (0:50) / 50 - 5
	hist(params$alpha, main="", xlab="alpha", col=c_dark, border=c_dark_highlight, yaxt='n')
	abline(v=3, col=c_light, lty=1, lwd=3)	
	beta_breaks=10 * (0:50) / 50 - 5
	hist(params$rho, main="", xlab="rho", col=c_dark, border=c_dark_highlight, yaxt='n')
	abline(v=5.5, col=c_light, lty=1, lwd=3)	
	sigma_breaks=5 * (0:50) / 50
	hist(params$sigma, main="", xlab="sigma", col=c_dark, border=c_dark_highlight, yaxt='n')
	abline(v=2, col=c_light, lty=1, lwd=3)
	#	plot divergent transitions in 2D posterior
	partition <- partition_div(fit)
	div_params <- partition[[1]]
	nondiv_params <- partition[[2]]	
	par(mfrow=c(1, 3))	
	par(mar = c(4, 4, 0.5, 0.5))
	plot(nondiv_params$rho, nondiv_params$alpha, log="xy", col=c_dark_trans, pch=16, cex=0.8, xlab="rho", ylab="alpha")
	points(div_params$rho, div_params$alpha, col='green', pch=16, cex=0.8)	
	par(mar = c(4, 4, 0.5, 0.5))
	plot(nondiv_params$rho, nondiv_params$sigma, col=c_dark_trans, pch=16, cex=0.8, xlab="rho", ylab="sigma")
	points(div_params$rho, div_params$sigma, col='green', pch=16, cex=0.8)	
	par(mar = c(4, 4, 0.5, 0.5))
	plot(nondiv_params$alpha, nondiv_params$sigma, col=c_dark_trans, pch=16, cex=0.8, xlab="alpha", ylab="sigma")
	points(div_params$alpha, div_params$sigma, col='green', pch=16, cex=0.8)
	#
	par(mfrow=c(1, 1))	
	plot_gp_quantiles(fit, data, true_realization, "Posterior Quantiles")	
}

GP2D.Gaussian <- function()
{
	home <- '~/git/phyloscanner/phyloflows/inst/misc'
	
	library(rstan)
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	
	alpha.true <- c(3)
	rho.true <- c(5.5)
	sigma.true <- c(2)
	gp.dim <- 2
	
	x.all.inputs <- expand.grid(x1= seq(15,50,1), x2=seq(15,50,1))
	n.all <- nrow( x.all.inputs )
	#x.all.inputs <- lapply(seq_len(ncol(x.all.inputs)), function(d) x.all.inputs[,d])	 			
	simu.pars <- list(	N=n.all, D=gp.dim, 
						alpha=alpha.true, rho=rho.true, sigma=sigma.true,
						x=x.all.inputs)
	
	#	simulate data set	
	simu.fit <- stan(	file=file.path(home,"191022_simu_gauss_v1.stan"), 
			data=simu.pars, iter=1,
			chains=1, seed=494838, algorithm="Fixed_param")
}

GP1D.Poisson <- function()
{
	home <- '~/git/phyloscanner/phyloflows/inst/misc'
	dir.mb.tutorial <- file.path(home,'MB_gp_tutorial')
	#
	#	this is Mike Betancourts GP tutorial 1-3
	#	with some modifications
	#	https://github.com/betanalpha/knitr_case_studies/tree/master/gaussian_processes
	#
	library(rstan)
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	source(file.path(dir.mb.tutorial,"gp_utility.R"))
	source(file.path(dir.mb.tutorial,"stan_utility.R"))
	
	alpha_true <- 3
	rho_true <- 5.5
	sigma_true <- 2	
	N_total <- 501
	x_total <- 20 * (0:(N_total - 1)) / (N_total - 1) - 10	
	simu_data <- list(	alpha=alpha_true, 
						rho=rho_true, 
						sigma=sigma_true,
						N=N_total, x=x_total)
	
	#	simulate data set		
	simu_fit <- stan(	file=file.path(dir.mb.tutorial,'simu_poisson.stan'), 
						data=simu_data, iter=1,
						chains=1, seed=494838, algorithm="Fixed_param")
	f_total <- extract(simu_fit)$f[1,]
	y_total <- extract(simu_fit)$y[1,]	
	true_realization <- data.frame(exp(f_total), x_total)
	names(true_realization) <- c("f_total", "x_total")	
	sample_idx <- c(50*(0:10)+1)
	N = length(sample_idx)
	x <- x_total[sample_idx]
	y <- y_total[sample_idx]	
	data <- list("N"=N, "x"=x, "y"=y, "N_predict"=N_predict, "x_predict"=x_total, "y_predict"=y_total)
	#	plot simulated data
	plot(x_total, exp(f_total), type="l", lwd=2, xlab="x", ylab="y", xlim=c(-10, 10), ylim=c(0, 10))
	points(x_total, y_total, col="white", pch=16, cex=0.6)
	points(x_total, y_total, col=c_mid_teal, pch=16, cex=0.4)
	points(x, y, col="white", pch=16, cex=1.2)
	points(x, y, col="black", pch=16, cex=0.8)
	#	save simulated data and truth
	N_predict <- N_total
	x_predict <- x_total
	y_predict <- y_total	
	stan_rdump(c("N", "x", "y", "N_predict", "x_predict", "y_predict", "sample_idx"), file=file.path(dir.mb.tutorial,"gppois.data.R"))
	data <- read_rdump(file.path(dir.mb.tutorial,"gppois.data.R"))	
	stan_rdump(c("f_total", "x_total", "sigma_true"), file=file.path(dir.mb.tutorial,"gppois.truth.R"))
	#	simulate 1e3 times and plot data generating process quantiles
	pred_data <- list(	alpha=alpha_true, rho=rho_true,
						N_predict=N_predict, x_predict=x_predict,
						N_observed=N, y_observed=y, observed_idx=observed_idx)
	pred_fit <- stan(file=file.path(dir.mb.tutorial,'predict_poisson.stan'), data=pred_data, seed=5838298, refresh=1000)
	#	this takes 139 seconds
	plot_gp_quantiles(pred_fit, data, true_realization, "Posterior Quantiles")	
	plot_gp_pred_quantiles(pred_fit, data, true_realization, "Posterior Predictive Quantiles")
}

GP2D.Poisson<- function()
{
	home <- '/Users/apple/Desktop/phyloscanner.utility/gp/'
	home <- '/Users/xx4515/Desktop/gp/'
	home <- '~/git/phyloscanner/phyloflows/inst/misc'
	dir.mb.tutorial <- file.path(home,'MB_gp_tutorial')
	out.file <-file.path(home,'figure')
	
	library(rstan)
	library(data.table)
	library(ggplot2)
	library(metR)
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	
	alpha_true <- c(3)
	rho_true <- c(5.5)
	gp_dim <- 2
	
	x_total <- expand.grid(x1= seq(15,50,1), x2=seq(15,50,1))
	N_total <- nrow(x_total)
	#x_total <- lapply(seq_len(ncol(x_total)), function(d) x_total[,d])     
	simu_data <- list( alpha=alpha_true, 
			rho=rho_true, 
			N=N_total, x=x_total,
			D=gp_dim)
	
	# simulate data set  
	simu_fit <- stan( file=file.path(dir.mb.tutorial,'simu_poisson_2d.stan'), 
			data=simu_data, iter=1,
			chains=1, seed=494838, algorithm="Fixed_param")
	
	f_total <- extract(simu_fit)$f[1,]
	y_total <- extract(simu_fit)$y[1,] 
	true_realization <- data.frame(exp(f_total), x_total)
	names(true_realization) <- c("f_total", "x_total") 
	sample_idx <- seq(from=1,to=N_total,by=50)
	N = length(sample_idx)
	x <- x_total[sample_idx,]
	y <- y_total[sample_idx] 
	
	# plot simulated data
	df <- data.table(x=x_total[,1],y=x_total[,2],
			z=f_total)
	range(df$z)
	ggplot(df, aes(x, y, z = z))+
			geom_contour(bins = 10)+
			geom_raster(aes(fill = z)) +
			geom_contour(colour = "white")+
			# geom_text_contour(aes(z = z), stroke = 0.2)+ 
			scale_fill_gradientn(colours = terrain.colors(10),limits=c(-7,7))+
			scale_x_continuous(expand = c(0,0), limits = c(15,50))+
			scale_y_continuous(expand = c(0,0), limits = c(15,50))+
			theme_classic()
	ggsave(file=paste0(out.file,'/pos2d/simulated_f.pdf'),width = 7.5, height = 6)
	
	df <- data.table(x=x_total[,1],y=x_total[,2],
			z=y_total)
	df2 <- data.table(x=x[,1],y=x[,2],z=y)
	range(df$z)
	range(df2$z)
	ggplot(df, aes(x, y, z = z))+
			geom_contour()+
			geom_raster(aes(fill = z)) +
			geom_contour(colour = "white")+
			# geom_text_contour(aes(z = z), stroke = 0.2)+ 
			scale_fill_gradientn(colours = terrain.colors(10),limits=c(0,40),
					name='all the data')+
			scale_x_continuous(expand = c(0,0), limits = c(15,50))+
			scale_y_continuous(expand = c(0,0), limits = c(15,50))+
			theme_classic()+
			geom_point(df2,mapping = aes(x=x,y=y,colour=z),size=3)+
			scale_colour_gradient(low = "white", high = "black",
					limits=c(0,40),name='observed data')
	ggsave(file=paste0(out.file,'/pos2d/simulated_y.pdf'),width = 7.5, height = 6)
	
	# save simulated data and truth
	N_predict <- N_total
	x_predict <- x_total
	y_predict <- y_total 
	stan_rdump(c("N", "x", "y", "N_predict", "x_predict", "y_predict", "sample_idx"), 
			file=file.path(dir.mb.tutorial,"gppois.data.2d.R"))
	data <- read_rdump(file.path(dir.mb.tutorial,"gppois.data.2d.R")) 
	stan_rdump(c("f_total", "x_total"), file=file.path(dir.mb.tutorial,"gppois.truth.2d.R"))
	
	# simulate from GP conditional on observations
	pred_data <- list( alpha=alpha_true, rho=rho_true,
			N_predict=N_predict, D=gp_dim, x_predict=x_predict,
			N_observed=N, y_observed=y, observed_idx=sample_idx)
	pred_fit <- stan(file=file.path(dir.mb.tutorial,'predict_poisson_2d.stan'), 
			data=pred_data, seed=5838298, refresh=1000,chains = 1)
	
	save(pred_fit,file=file.path(dir.mb.tutorial,'fit_predict_poisson_2d.rda'))
	params <- extract(pred_fit)
	cred <- sapply(1:N_total,
			function(n) quantile(params$y_predict[,n], probs=c(0.025,0.5,0.975)))
	df <- data.table(x=x_total[,1],y=x_total[,2],
			z=cred[1,])
	range(df$z)
	ggplot(df, aes(x, y, z = z))+
			geom_contour(bins = 10)+
			geom_raster(aes(fill = z)) +
			geom_contour(colour = "white")+
			# geom_text_contour(aes(z = z), stroke = 0.2)+ 
			scale_fill_gradientn(colours = terrain.colors(10),limits=c(0,5))+
			scale_x_continuous(expand = c(0,0), limits = c(15,50))+
			scale_y_continuous(expand = c(0,0), limits = c(15,50))+
			theme_classic()
	ggsave(file=paste0(out.file,'/pos2d/predicted_ylb.pdf'),width = 7.5, height = 6)
	
	df <- data.table(x=x_total[,1],y=x_total[,2],
			z=cred[2,])
	range(df$z)
	ggplot(df, aes(x, y, z = z))+
			geom_contour(bins = 10)+
			geom_raster(aes(fill = z)) +
			geom_contour(colour = "white")+
			# geom_text_contour(aes(z = z), stroke = 0.2)+ 
			scale_fill_gradientn(colours = terrain.colors(10),limits=c(0,20))+
			scale_x_continuous(expand = c(0,0), limits = c(15,50))+
			scale_y_continuous(expand = c(0,0), limits = c(15,50))+
			theme_classic()
	ggsave(file=paste0(out.file,'/pos2d/predicted_ym.pdf'),width = 7.5, height = 6)
	
	df <- data.table(x=x_total[,1],y=x_total[,2],
			z=cred[3,])
	range(df$z)
	ggplot(df, aes(x, y, z = z))+
			geom_contour(bins = 10)+
			geom_raster(aes(fill = z)) +
			geom_contour(colour = "white")+
			# geom_text_contour(aes(z = z), stroke = 0.2)+ 
			scale_fill_gradientn(colours = terrain.colors(10),limits=c(0,460))+
			scale_x_continuous(expand = c(0,0), limits = c(15,50))+
			scale_y_continuous(expand = c(0,0), limits = c(15,50))+
			theme_classic()
	ggsave(file=paste0(out.file,'/pos2d/predicted_yub.pdf'),width = 7.5, height = 6)
	
	cred <- data.table(t(cred))
	cred[,id:=seq_len(nrow(cred))]
	cred[,y:=y_predict]
	cred[,x1:=x_predict[,1]]
	cred[,x2:=x_predict[,2]]
	setnames(cred,c('2.5%','50%','97.5%'),c('lb','m','ub'))
	cred[,observed:=0]
	cred$observed[sample_idx] <- 1
	
	ggplot(cred, aes(x=x1, y=m,color=factor(observed))) + 
			geom_point(cred,mapping=aes(x=x1,y=y,color=factor(observed)))+
			facet_grid(x2~.)+
			geom_errorbar(aes(ymin=lb, ymax=ub, color=factor(observed)), width=.2,
					position=position_dodge(0.05))+
			theme_classic()+
			theme(legend.position = 'none')+
			scale_color_manual(values=c('black','red'))
	ggsave(file=paste0(out.file,'/pos2d/prediction.pdf'),
			width = 50, height=50,limitsize = FALSE)
	
	# estimate hyperparameter
	pred_fit <- stan(file=file.path(dir.mb.tutorial,'predict_poisson2_2d.stan'), 
			data=pred_data, seed=5838298, refresh=1000,chains = 1)
}


Rakai.age.no.adjustment<- function(infile.inference=NULL,
                                   opt=NULL
)
{
  
  require(data.table)	
  require(phyloflows)
  require(rstan)
  
  #
  #	input args
  #
  if(is.null(opt))
  {
    opt									<- list()
    opt$adjust.sequencing.bias			<- 1
    opt$adjust.participation.bias		<- 1
    opt$migration.def.code				<- '24'
    opt$set.missing.migloc.to.inland	<- 0
    opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
  }
  
  if(is.null(infile.inference))
  {
    infile.inference	<- "RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
  }

  cat('\ninfile.inference=',infile.inference)
  cat('\nopt=',unlist(opt))			
  load(infile.inference)
  
  #
  #	prepare data on observed transmission flows
  #
  #	subset to variables needed, using RTR3	
  rtr	<- copy(rtr3)
  if(opt$migration.def.code=='06')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
    setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='12')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
    setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='24')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
    setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='36')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
    setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='48')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
    setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
  }
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
                              'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
                              'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
                              'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
  # inmigrant status
  rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
  rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
  set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
  
  #	set unknown origin to either fishing or inland
  tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
  if(opt$set.missing.migloc.to.inland)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
  }		
  if(opt$set.missing.migloc.to.fishing)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
  }
  
  # add age 
  rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
  rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
  
  # impute age
  tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
  set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
  tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
  set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
  
  # fixup from latest surveillance data
  set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
  set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
  set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
  set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
  
  #	stratify age
  range(c(rtr$TR_AGE_AT_MID,rtr$REC_AGE_AT_MID)) #16.867 50.458
  rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(0,17:50,60), labels=paste0(16:50,'-',17:51), right=FALSE))]
  rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(0,17:50,60), labels=paste0(16:50,'-',17:51), right=FALSE))]
  stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
  stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
  
  # define TR_COMM_TYPE_F, REC_COMM_TYPE_F (i: inland; f: fishing) 
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # define TR_COMM_TYPE_F_MIG (i: inland; f: fishing; e: external) 
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_SEX,':',TR_AGE_AT_MID_C)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_SEX,':',TR_AGE_AT_MID_C)]
  
  # make all combinations of variables
  dac <- expand.grid( SEX= c('M','F'),
                      AGE_AT_MID_C= paste0(16:50,'-',17:51))
  
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(SEX, ':', AGE_AT_MID_C)]  					
  dac <- as.data.table(expand.grid(TR_CATEGORY= dac$CATEGORY, REC_CATEGORY= dac$CATEGORY))
  # ignore Male-Male and Female-Female combinations
  dac <- subset(dac, !(grepl('F',TR_CATEGORY)&grepl('F',REC_CATEGORY)) &
                  !(grepl('M',TR_CATEGORY)&grepl('M',REC_CATEGORY))  
  )
  # add transmission categories
  dac[, REC_TRM_CATEGORY:= REC_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_CATEGORY] 
  dac[, REC_SAMPLING_CATEGORY:=REC_TRM_CATEGORY]
  dac[, TR_SAMPLING_CATEGORY:=TR_TRM_CATEGORY]
  
  #	calculate observed number of transmissions
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  dac[, DUMMY:= 1]
  dobs <- merge(dac, dobs, by=c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'), all=TRUE)
  stopifnot( dobs[, !any(is.na(DUMMY))] )
  set(dobs, NULL, 'DUMMY', NULL)
  set(dobs, dobs[, which(is.na(TRM_OBS))], 'TRM_OBS', 0L)
  
  # take f-m transmissions
  dobs <- dobs[grepl('F',TR_TRM_CATEGORY),]
  #	make PAIR_ID
  setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY,TR_SAMPLING_CATEGORY,REC_SAMPLING_CATEGORY)	
  dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
  dobs[,TR_SMOOTH_CATEGORY:= as.numeric(substr(TR_TRM_CATEGORY,3,4))+0.5]
  dobs[,REC_SMOOTH_CATEGORY:= as.numeric(substr(REC_TRM_CATEGORY,3,4))+0.5]
  
  #	estimate GP and hyperparameters with informative pior on length scale
  range(dist(dobs[TRM_OBS!=0,TR_SMOOTH_CATEGORY]))
  range(dist(dobs[TRM_OBS!=0,REC_SMOOTH_CATEGORY]))
  
  # [0,30]
  fit <- stan(file=file.path(outfile.base,'gp_prior_tune_30.stan'), 
              iter=1, warmup=0, chains=1,
              seed=5838298, algorithm="Fixed_param")
  # a = 1.78207
  # b = 3.11324
  
  # set up stan 
  M <- 20
  D <- 2 
  indices <- matrix(NA, M^D, D)
  mm=0;
  for (m1 in 1:M){
    for (m2 in 1:M){
      mm = mm+1
      indices[mm,] = c(m1, m2)
    }
  }
  
  standata_bf12 <- list( M= M, M_nD= M^D, 
                         L= c(5/2*max(dobs$TR_SMOOTH_CATEGORY),5/2*max(dobs$REC_SMOOTH_CATEGORY)), 
                         N = nrow(dobs),
                         x = cbind(dobs$TR_SMOOTH_CATEGORY,dobs$REC_SMOOTH_CATEGORY),
                         D = D,
                         y = dobs$TRM_OBS,
                         indices= indices)
  
  
  fit <- stan(file = file.path('gpa_ard_nonzero.stan'),
              data = standata_bf12,
              iter = 10000,  warmup = 2000, chains=1, thin=1, seed = 1234,
              algorithm = "NUTS", verbose = FALSE,
              control = list(adapt_delta = 0.999))
  
}


Rakai.sampling.glm <- function(){
  require(data.table)
  require(rethinking)
  
  indir <- '/Users/xx4515/Desktop/data'
  infile.data <- file.path(indir,"190327_sampling_by_gender_age.rda")
  outfile.base			<- '/Users/xx4515/Desktop/gp2/rakai'
  
  #	set up variables for STAN
  load(infile.data)
  des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, 
                              MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
  
  # remove individuals in the communities where no sequences were obtained successfully
  tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
  des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')

  # age group
  des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(16:51), labels = paste0(16:50,'-',17:51))]
  
  des <- des[!is.na(AGE_AT_MID_C),]
  
  #	binarize age, sex
  des[, AGE1:= as.integer(AGE_AT_MID_C=="16-17")]
  des[, AGE2:= as.integer(AGE_AT_MID_C=="17-18")]
  des[, AGE3:= as.integer(AGE_AT_MID_C=="18-19")]
  des[, AGE4:= as.integer(AGE_AT_MID_C=="19-20")]
  des[, AGE5:= as.integer(AGE_AT_MID_C=="20-21")]
  des[, AGE6:= as.integer(AGE_AT_MID_C=="21-22")]
  des[, AGE7:= as.integer(AGE_AT_MID_C=="22-23")]
  des[, AGE8:= as.integer(AGE_AT_MID_C=="23-24")]
  des[, AGE9:= as.integer(AGE_AT_MID_C=="24-25")]
  des[, AGE10:= as.integer(AGE_AT_MID_C=="25-26")]
  des[, AGE11:= as.integer(AGE_AT_MID_C=="26-27")]
  des[, AGE12:= as.integer(AGE_AT_MID_C=="27-28")]
  des[, AGE13:= as.integer(AGE_AT_MID_C=="28-29")]
  des[, AGE14:= as.integer(AGE_AT_MID_C=="29-30")]
  des[, AGE15:= as.integer(AGE_AT_MID_C=="30-31")]
  des[, AGE16:= as.integer(AGE_AT_MID_C=="31-32")]
  des[, AGE17:= as.integer(AGE_AT_MID_C=="32-33")]
  des[, AGE18:= as.integer(AGE_AT_MID_C=="33-34")]
  des[, AGE19:= as.integer(AGE_AT_MID_C=="34-35")]
  des[, AGE20:= as.integer(AGE_AT_MID_C=="35-36")]
  des[, AGE21:= as.integer(AGE_AT_MID_C=="36-37")]
  des[, AGE22:= as.integer(AGE_AT_MID_C=="37-38")]
  des[, AGE23:= as.integer(AGE_AT_MID_C=="38-39")]
  des[, AGE24:= as.integer(AGE_AT_MID_C=="39-40")]
  des[, AGE25:= as.integer(AGE_AT_MID_C=="40-41")]
  des[, AGE26:= as.integer(AGE_AT_MID_C=="41-42")]
  des[, AGE27:= as.integer(AGE_AT_MID_C=="42-43")]
  des[, AGE28:= as.integer(AGE_AT_MID_C=="43-44")]
  des[, AGE29:= as.integer(AGE_AT_MID_C=="44-45")]
  des[, AGE30:= as.integer(AGE_AT_MID_C=="45-46")]
  des[, AGE31:= as.integer(AGE_AT_MID_C=="46-47")]
  des[, AGE32:= as.integer(AGE_AT_MID_C=="47-48")]
  des[, AGE33:= as.integer(AGE_AT_MID_C=="48-49")]
  des[, AGE34:= as.integer(AGE_AT_MID_C=="49-50")]
  
  
  des[, MALE:= as.integer(SEX=='M')]	
  
  #	aggregate participation
  dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), 
             by=c('AGE_AT_MID_C','SEX','MALE',
                  'AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','AGE7','AGE8','AGE9','AGE10',
                  'AGE11','AGE12','AGE13','AGE14','AGE15','AGE16','AGE17','AGE18','AGE19','AGE20',
                  'AGE21','AGE22','AGE23','AGE24','AGE25','AGE26','AGE27','AGE28','AGE29','AGE30',
                  'AGE31','AGE32','AGE33','AGE34')]
  dp[, CATEGORY:= paste0(SEX,':',AGE_AT_MID_C)]
  
  opt.exclude.onART.from.denominator<-1
  #	aggregate sequencing
  des		<- subset(des, HIV_1517==1)
  if(opt.exclude.onART.from.denominator)
  {
    stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
    des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
  }
  
  ds<-des[,list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))),
          by=c('AGE_AT_MID_C','SEX','MALE',
               'AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','AGE7','AGE8','AGE9','AGE10',
               'AGE11','AGE12','AGE13','AGE14','AGE15','AGE16','AGE17','AGE18','AGE19','AGE20',
               'AGE21','AGE22','AGE23','AGE24','AGE25','AGE26','AGE27','AGE28','AGE29','AGE30',
               'AGE31','AGE32','AGE33','AGE34')]
  ds[, CATEGORY:= paste0(SEX,':',AGE_AT_MID_C)]
  
  
  library(rethinking)
  mp1 	<- map2stan(
    alist(
      SUC ~ dbinom(TRIAL, p_part),
      logit(p_part) <- a +  male*MALE + 
        age1*AGE1 + age2*AGE2 + age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 + age7 * AGE7 + age8 * AGE8 + age9 * AGE9 +  age10 * AGE10 +
        age11*AGE11 + age12*AGE12 + age13 * AGE13 + age14 * AGE14 + age15 * AGE15 + age16 * AGE16 + age17 * AGE17 + age18 * AGE18 + age19 * AGE19 +  age20 * AGE20 +
        age21*AGE21 + age22*AGE22 + age23 * AGE23 + age24 * AGE24 + age25 * AGE25 + age26 * AGE26 + age27 * AGE27 + age28 * AGE28 + age29 * AGE29 +  age30 * AGE30 +
        age31*AGE31 + age32*AGE32 + age33 * AGE33 + age34 * AGE34,
      a ~ dnorm(0, 100),
      c(male, age1, age2, age3, age4, age5, age6, age7, age8, age9, age10,
        age11, age12, age13, age14, age15, age16, age17, age18, age19, age20,
        age21, age22, age23, age24, age25, age26, age27, age28, age29, age30,
        age31, age32, age33, age34
      ) ~ dnorm(0,10)
    ),
    data=as.data.frame(dp), 
    start=list(a=0, male=0,
               age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, age7=0, age8=0, age9=0, age10=0,
               age11=0, age12=0, age13=0, age14=0, age15=0, age16=0, age17=0, age18=0, age19=0, age20=0,
               age21=0, age22=0, age23=0, age24=0, age25=0, age26=0, age27=0, age28=0, age29=0, age30=0,
               age31=0, age32=0, age33=0, age34=0),
    warmup=5e2, iter=5e3, chains=1, cores=4
  )	
  
  plot( precis(mp1, depth=2, prob=0.95) )
  #	posterior check to see if this makes sense
  sims 	<- sim(mp1, n=1e3)
  tmp		<- apply(sims, 2, median)
  dpp 	<- apply(sims, 2, PI, prob=0.95)
  dpp		<- rbind(tmp, dpp)
  rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
  dpp		<- as.data.table(t(dpp))
  dpp		<- cbind(dp, dpp)		
  dpp[, c(mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95), sum(SUC<predicted_obs_l95 | SUC>predicted_obs_u95))]
  #[1] 0.2285714 16.0000000
  
  
  dpp[, OFF:= as.integer(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
  dpp[OFF!=0]
  table(dpp[OFF!=0,paste0(AGE_AT_MID_C,':',SEX)])
  
  
  mp2 	<- map2stan(
    alist(
      SUC ~ dbinom(TRIAL, p_part),
      logit(p_part) <- a +  male*MALE + 
        age1*AGE1 + age2*AGE2 + age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 + age7 * AGE7 + age8 * AGE8 + age9 * AGE9 +  age10 * AGE10 +
        age11*AGE11 + age12*AGE12 + age13 * AGE13 + age14 * AGE14 + age15 * AGE15 + age16 * AGE16 + age17 * AGE17 + age18 * AGE18 + age19 * AGE19 +  age20 * AGE20 +
        age21*AGE21 + age22*AGE22 + age23 * AGE23 + age24 * AGE24 + age25 * AGE25 + age26 * AGE26 + age27 * AGE27 + age28 * AGE28 + age29 * AGE29 +  age30 * AGE30 +
        age31*AGE31 + age32*AGE32 + age33 * AGE33 + age34 * AGE34 +
        age1_female * AGE1 * (1-MALE) + age1_male * AGE1 * MALE +
        age2_female * AGE2 * (1-MALE) + age2_male * AGE2 * MALE +
        age3_male * AGE3 * MALE +
        age4_female * AGE4 * (1-MALE) + age4_male * AGE4 * MALE +
        age13_female * AGE13 * (1-MALE) + 
        age15_female * AGE15 * (1-MALE) + 
        age17_female * AGE17 * (1-MALE) + 
        age18_female * AGE18 * (1-MALE) + 
        age19_female * AGE19 * (1-MALE) + age19_male * AGE19 * MALE +
        age20_female * AGE20 * (1-MALE) +
        age22_female * AGE22 * (1-MALE) +
        age28_female * AGE28 * (1-MALE),
      a ~ dnorm(0, 100),
      c(male, age1, age2, age3, age4, age5, age6, age7, age8, age9, age10,
        age11, age12, age13, age14, age15, age16, age17, age18, age19, age20,
        age21, age22, age23, age24, age25, age26, age27, age28, age29, age30,
        age31, age32, age33, age34,
        age1_female, age1_male, age2_female, age2_male, 
        age3_male, age4_female, age4_male, age13_female, 
        age15_female, age17_female, age18_female,
        age19_female, age19_male, age20_female, age22_female,
        age28_female
      ) ~ dnorm(0,10)
    ),
    data=as.data.frame(dp), 
    start=list(a=0, male=0,
               age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, age7=0, age8=0, age9=0, age10=0,
               age11=0, age12=0, age13=0, age14=0, age15=0, age16=0, age17=0, age18=0, age19=0, age20=0,
               age21=0, age22=0, age23=0, age24=0, age25=0, age26=0, age27=0, age28=0, age29=0, age30=0,
               age31=0, age32=0, age33=0, age34=0, 
               age1_female=0, age1_male=0, age2_female=0, age2_male=0, 
               age3_male=0, age4_female=0, age4_male=0, age13_female=0, 
               age15_female=0, age17_female=0, age18_female=0,
               age19_female=0, age19_male=0, age20_female=0, age22_female=0,
               age28_female=0),
    warmup=5e2, iter=5e3, chains=1, cores=4
  )	
  
  plot( precis(mp2, depth=2, prob=0.95) )
  #	posterior check to see if this makes sense
  sims 	<- sim(mp2, n=1e3)
  tmp		<- apply(sims, 2, median)
  dpp 	<- apply(sims, 2, PI, prob=0.95)
  dpp		<- rbind(tmp, dpp)
  rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
  dpp		<- as.data.table(t(dpp))
  dpp		<- cbind(dp, dpp)		
  dpp[, c(mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95), sum(SUC<predicted_obs_l95 | SUC>predicted_obs_u95))]
  #0 0
  
  
  dpp[, OFF:= as.integer(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
  dpp[OFF!=0]
  table(dpp[OFF!=0,paste0(AGE_AT_MID_C,':',SEX)])
  
  stancode(mp2)
  
  save(mp1,mp2,file=file.path(outfile.base,'partm.rda'))
  require(data.table)
  require(rstan)
  require(extraDistr)
  require(bayesplot)
  
  #	run STAN 
  tmp			<- as.list(subset(dp,select=c('TRIAL','SUC', 
                                      'MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','AGE7','AGE8','AGE9','AGE10',
                                      'AGE11','AGE12','AGE13','AGE14','AGE15','AGE16','AGE17','AGE18','AGE19','AGE20',
                                      'AGE21','AGE22','AGE23','AGE24','AGE25','AGE26','AGE27','AGE28','AGE29','AGE30',
                                      'AGE31','AGE32','AGE33','AGE34')))
  tmp$N		<- nrow(dp)
  
  fit.par 	<- stan(	file = file.path(outfile.base,'part.stan'), 
                    data = tmp, 
                    iter = 10e3,
                    warmup = 5e2,
                    cores = 1,
                    chains = 1,
                    control=list(max_treedepth=15,adapt_delta=0.9),
                    init = list(list(a=0, male=0,age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, age7=0, age8=0, age9=0, age10=0,
                                     age11=0, age12=0, age13=0, age14=0, age15=0, age16=0, age17=0, age18=0, age19=0, age20=0,
                                     age21=0, age22=0, age23=0, age24=0, age25=0, age26=0, age27=0, age28=0, age29=0, age30=0,
                                     age31=0, age32=0, age33=0, age34=0,
                                     age1_female=0, age1_male=0, age2_female=0, age2_male=0, 
                                     age3_male=0, age4_female=0, age4_male=0, age13_female=0, 
                                     age15_female=0, age17_female=0, age18_female=0,
                                     age19_female=0, age19_male=0, age20_female=0, age22_female=0,
                                     age28_female=0)))
  
  # assess convergence
  
  fit.pars	<- c('a','male', 'age1', 'age2', 'age3', 'age4', 'age5', 'age6', 'age7', 'age8', 'age9', 'age10',
                'age11', 'age12', 'age13', 'age14', 'age15', 'age16', 'age17', 'age18', 'age19', 'age20',
                'age21', 'age22', 'age23', 'age24', 'age25', 'age26', 'age27', 'age28', 'age29', 'age30',
                'age31', 'age32', 'age33', 'age34', 
                'age1_female', 'age1_male', 'age2_female', 'age2_male', 
                'age3_male', 'age4_female', 'age4_male', 'age13_female', 
                'age15_female', 'age17_female', 'age18_female',
                'age19_female', 'age19_male', 'age20_female', 'age22_female',
                'age28_female')
  
  any(rhat(fit.par, pars=fit.pars)>1.02)
  any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
  po              <- as.matrix(fit.par) 
  po              <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]	
  p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
    geom_vline(xintercept = 0)
  pdf(file=file.path(outfile.base,'part_marginal_posterior.pdf'), w=7, h=20)
  p
  dev.off()
  p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
  pdf(file=file.path(outfile.base,'part_marginal_trace.pdf'), w=7, h=100)
  p
  dev.off()
  
  # summary csv file
  write.csv(  summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
              file=file.path(outfile.base,'part_summary.csv')
  )
  
  # extract samples for unique strata levels
  nprior		<- 1e3
  dps <- dp
  
  fit.e		<- extract(fit.par)
  set.seed(42)
  tmp			<- sample(length(fit.e$a), nprior)
  dps			<- dps[,	
               {
                 z<- with(fit.e, a +  male*MALE + 
                            age1*AGE1 + age2*AGE2 + age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 + age7 * AGE7 + age8 * AGE8 + age9 * AGE9 +  age10 * AGE10 +
                            age11*AGE11 + age12*AGE12 + age13 * AGE13 + age14 * AGE14 + age15 * AGE15 + age16 * AGE16 + age17 * AGE17 + age18 * AGE18 + age19 * AGE19 +  age20 * AGE20 +
                            age21*AGE21 + age22*AGE22 + age23 * AGE23 + age24 * AGE24 + age25 * AGE25 + age26 * AGE26 + age27 * AGE27 + age28 * AGE28 + age29 * AGE29 +  age30 * AGE30 +
                            age31*AGE31 + age32*AGE32 + age33 * AGE33 + age34 * AGE34 +
                            age1_female * AGE1 * (1-MALE) + age1_male * AGE1 * MALE +
                            age2_female * AGE2 * (1-MALE) + age2_male * AGE2 * MALE +
                            age3_male * AGE3 * MALE +
                            age4_female * AGE4 * (1-MALE) + age4_male * AGE4 * MALE +
                            age13_female * AGE13 * (1-MALE) + 
                            age15_female * AGE15 * (1-MALE) + 
                            age17_female * AGE17 * (1-MALE) + 
                            age18_female * AGE18 * (1-MALE) + 
                            age19_female * AGE19 * (1-MALE) + age19_male * AGE19 * MALE +
                            age20_female * AGE20 * (1-MALE) +
                            age22_female * AGE22 * (1-MALE) +
                            age28_female * AGE28 * (1-MALE))
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dps[, P:= exp(ETA)/(1+exp(ETA))]
  
  require(bde)
  tmp	<- dps[, {
    bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  }, by='CATEGORY']
  
  dps	<- merge(dps, tmp, by=c('CATEGORY','SAMPLE'))
  set(dps, NULL, c('ETA'), NULL)		
  
  save(dp, dps, fit.par, file=file.path(outfile.base,'part_samples2.rda'))
  
  
  ms1 	<- map2stan(
    alist(
      SUC ~ dbinom(TRIAL, p_seq),
      logit(p_seq) <- a + male*MALE + 
        age1*AGE1 + age2*AGE2 + age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 + age7 * AGE7 + age8 * AGE8 + age9 * AGE9 +  age10 * AGE10 +
        age11*AGE11 + age12*AGE12 + age13 * AGE13 + age14 * AGE14 + age15 * AGE15 + age16 * AGE16 + age17 * AGE17 + age18 * AGE18 + age19 * AGE19 +  age20 * AGE20 +
        age21*AGE21 + age22*AGE22 + age23 * AGE23 + age24 * AGE24 + age25 * AGE25 + age26 * AGE26 + age27 * AGE27 + age28 * AGE28 + age29 * AGE29 +  age30 * AGE30 +
        age31*AGE31 + age32*AGE32 + age33 * AGE33 + age34 * AGE34,
      a ~ dnorm(0, 100),
      c(male, age1, age2, age3, age4, age5, age6, age7, age8, age9, age10,
        age11, age12, age13, age14, age15, age16, age17, age18, age19, age20,
        age21, age22, age23, age24, age25, age26, age27, age28, age29, age30,
        age31, age32, age33, age34) ~ dnorm(0,10)
    ),
    data=as.data.frame(ds), 
    start=list(a=0, male=0,age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, age7=0, age8=0, age9=0, age10=0,
               age11=0, age12=0, age13=0, age14=0, age15=0, age16=0, age17=0, age18=0, age19=0, age20=0,
               age21=0, age22=0, age23=0, age24=0, age25=0, age26=0, age27=0, age28=0, age29=0, age30=0,
               age31=0, age32=0, age33=0, age34=0),
    warmup=5e2, iter=5e3, chains=1, cores=4
  )	
  
  plot( precis(ms1, depth=2, prob=0.95) )
  #	posterior check to see if this makes sense
  sims 	<- sim(ms1, n=1e3)
  tmp		<- apply(sims, 2, median)
  dpp 	<- apply(sims, 2, PI, prob=0.95)
  dpp		<- rbind(tmp, dpp)
  rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
  dpp		<- as.data.table(t(dpp))
  dpp		<- cbind(ds, dpp)		
  dpp[, c(mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95), sum(SUC<predicted_obs_l95 | SUC>predicted_obs_u95))]
  
  dpp[, OFF:= as.integer(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
  
  #  0 0
  stancode(ms1)
  save(ms1,file=file.path(outfile.base,'seqm.rda'))
  #	run STAN 
  tmp			<- as.list(subset(ds,select=c('AGE_AT_MID_C', 'SEX', 'TRIAL', 'SUC',
                                      'MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','AGE7','AGE8','AGE9','AGE10',
                                      'AGE11','AGE12','AGE13','AGE14','AGE15','AGE16','AGE17','AGE18','AGE19','AGE20',
                                      'AGE21','AGE22','AGE23','AGE24','AGE25','AGE26','AGE27','AGE28','AGE29','AGE30',
                                      'AGE31','AGE32','AGE33','AGE34')))
  tmp$N		<- nrow(ds)
  fit.seq 	<- stan(file = file.path(outfile.base,'seq.stan'), 
                   data = tmp, 
                   iter = 30e3,
                   warmup = 5e2,
                   cores = 1,
                   chains = 1,
                   init = list(list(a=0, male=0,age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, age7=0, age8=0, age9=0, age10=0,
                                    age11=0, age12=0, age13=0, age14=0, age15=0, age16=0, age17=0, age18=0, age19=0, age20=0,
                                    age21=0, age22=0, age23=0, age24=0, age25=0, age26=0, age27=0, age28=0, age29=0, age30=0,
                                    age31=0, age32=0, age33=0, age34=0)))
  
  # assess convergence
  fit.pars	<- c('a', 'male', 'age1', 'age2', 'age3', 'age4', 'age5', 'age6', 'age7', 'age8', 'age9', 'age10',
                'age11', 'age12', 'age13', 'age14', 'age15', 'age16', 'age17', 'age18', 'age19', 'age20',
                'age21', 'age22', 'age23', 'age24', 'age25', 'age26', 'age27', 'age28', 'age29', 'age30',
                'age31', 'age32', 'age33', 'age34')
  any(rhat(fit.seq, pars=fit.pars)>1.02)
  any(neff_ratio(fit.seq, pars=fit.pars) * 2.95e4 < 500)
  po              <- as.matrix(fit.seq) 
  po              <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
  p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
    geom_vline(xintercept = 0)
  pdf(file=file.path(outfile.base,'seq_marginal_posterior.pdf'), w=7, h=12)
  p
  dev.off()
  p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
  pdf(file=file.path(outfile.base,'seq_marginal_traces.pdf'), w=7, h=100)
  p
  dev.off()
  
  # summary csv file
  write.csv(  summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
              file=file.path(outfile.base,'seq_summary.csv')
  )	
  
  # extract samples for unique strata levels
  nprior		<- 1e3
  dss   <- data.table(AGE_AT_MID_C=rep(paste0(16:50,'-',17:51),2),SEX=c(rep('M',35),rep('F',35)))
  dss[, AGE1:= as.integer(AGE_AT_MID_C=="16-17")]
  dss[, AGE2:= as.integer(AGE_AT_MID_C=="17-18")]
  dss[, AGE3:= as.integer(AGE_AT_MID_C=="18-19")]
  dss[, AGE4:= as.integer(AGE_AT_MID_C=="19-20")]
  dss[, AGE5:= as.integer(AGE_AT_MID_C=="20-21")]
  dss[, AGE6:= as.integer(AGE_AT_MID_C=="21-22")]
  dss[, AGE7:= as.integer(AGE_AT_MID_C=="22-23")]
  dss[, AGE8:= as.integer(AGE_AT_MID_C=="23-24")]
  dss[, AGE9:= as.integer(AGE_AT_MID_C=="24-25")]
  dss[, AGE10:= as.integer(AGE_AT_MID_C=="25-26")]
  dss[, AGE11:= as.integer(AGE_AT_MID_C=="26-27")]
  dss[, AGE12:= as.integer(AGE_AT_MID_C=="27-28")]
  dss[, AGE13:= as.integer(AGE_AT_MID_C=="28-29")]
  dss[, AGE14:= as.integer(AGE_AT_MID_C=="29-30")]
  dss[, AGE15:= as.integer(AGE_AT_MID_C=="30-31")]
  dss[, AGE16:= as.integer(AGE_AT_MID_C=="31-32")]
  dss[, AGE17:= as.integer(AGE_AT_MID_C=="32-33")]
  dss[, AGE18:= as.integer(AGE_AT_MID_C=="33-34")]
  dss[, AGE19:= as.integer(AGE_AT_MID_C=="34-35")]
  dss[, AGE20:= as.integer(AGE_AT_MID_C=="35-36")]
  dss[, AGE21:= as.integer(AGE_AT_MID_C=="36-37")]
  dss[, AGE22:= as.integer(AGE_AT_MID_C=="37-38")]
  dss[, AGE23:= as.integer(AGE_AT_MID_C=="38-39")]
  dss[, AGE24:= as.integer(AGE_AT_MID_C=="39-40")]
  dss[, AGE25:= as.integer(AGE_AT_MID_C=="40-41")]
  dss[, AGE26:= as.integer(AGE_AT_MID_C=="41-42")]
  dss[, AGE27:= as.integer(AGE_AT_MID_C=="42-43")]
  dss[, AGE28:= as.integer(AGE_AT_MID_C=="43-44")]
  dss[, AGE29:= as.integer(AGE_AT_MID_C=="44-45")]
  dss[, AGE30:= as.integer(AGE_AT_MID_C=="45-46")]
  dss[, AGE31:= as.integer(AGE_AT_MID_C=="46-47")]
  dss[, AGE32:= as.integer(AGE_AT_MID_C=="47-48")]
  dss[, AGE33:= as.integer(AGE_AT_MID_C=="48-49")]
  dss[, AGE34:= as.integer(AGE_AT_MID_C=="49-50")]
  dss[, MALE := as.integer(SEX=='M')]
  dss[, CATEGORY:=paste0(SEX,':',AGE_AT_MID_C)]
  
  
  fit.e		<- extract(fit.seq)
  set.seed(42)
  tmp			<- sample(length(fit.e$a), nprior)
  dss			<- dss[,	
               {
                 z<- with(fit.e, a + male*MALE + 
                            age1*AGE1 + age2*AGE2 + age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 + age7 * AGE7 + age8 * AGE8 + age9 * AGE9 +  age10 * AGE10 +
                            age11*AGE11 + age12*AGE12 + age13 * AGE13 + age14 * AGE14 + age15 * AGE15 + age16 * AGE16 + age17 * AGE17 + age18 * AGE18 + age19 * AGE19 +  age20 * AGE20 +
                            age21*AGE21 + age22*AGE22 + age23 * AGE23 + age24 * AGE24 + age25 * AGE25 + age26 * AGE26 + age27 * AGE27 + age28 * AGE28 + age29 * AGE29 +  age30 * AGE30 +
                            age31*AGE31 + age32*AGE32 + age33 * AGE33 + age34 * AGE34)
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dss[, P:= exp(ETA)/(1+exp(ETA))]
  
  require(bde)
  tmp	<- dss[, {
    bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  }, by='CATEGORY']
  
  dss	<- merge(dss, tmp, by=c('CATEGORY','SAMPLE'))
  set(dss, NULL, c('ETA'), NULL)		
  
  save(ds, dss, fit.seq, file=file.path(outfile.base,'seq_samples.rda'))
  
  set.seed(42)
  outfile.base <- '/Users/xx4515/Desktop/gp2/rakai'
  load(file.path(outfile.base,'part_samples2.rda'))
  load(file.path(outfile.base,'seq_samples.rda'))
  dprior	<- merge(dps, dss, by=c('CATEGORY','SAMPLE'))
  if(opt$adjust.sequencing.bias==0)
  {
    set(dprior, NULL, 'P.y', 1)
    set(dprior, NULL, 'LP.y', 0)
  }
  if(opt$adjust.participation.bias==0)
  {
    set(dprior, NULL, 'P.x', 1)
    set(dprior, NULL, 'LP.x', 0)
  }
  dprior[, P:= P.x*P.y]		# multiply participation and sequencing probabilities 
  dprior[, LP:= LP.x+LP.y]	# add log posterior densities
  set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
  setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
  
  ggplot(dprior, aes(x=P)) +
    geom_histogram()+
    facet_grid(SAMPLING_CATEGORY~.)+
    theme_classic()
  ggsave(file=paste0(outfile.base,'/sampling_fractions.pdf'),width = 6, height = 1.2*70,
         limitsize = FALSE)
  
  dprior.fit <- dprior[,{tmp <- fitdistr(P, "beta", start=list(shape1=1,shape2=1),lower=c(0,0))
  list(SHAPE1=tmp$estimate[1],SHAPE2=tmp$estimate[2])
  },by='SAMPLING_CATEGORY']
  save(dprior.fit,file=paste0(outfile.base,'samples_fit.rda'))
  
  dprior.fit.samples <- dprior.fit[,{
    list(P = rbeta(10000,SHAPE1,SHAPE2))
  },by='SAMPLING_CATEGORY']
  
  dprior.fit.plot <- subset(dprior,select = c('SAMPLING_CATEGORY','P'))
  dprior.fit.plot[,TYPE:='prior']
  dprior.fit.samples[,TYPE:='fit']
  dprior.fit.plot <- rbind(dprior.fit.plot, dprior.fit.samples)
  ggplot(dprior.fit.plot, aes(x=P,color=TYPE)) +
    geom_density()+
    theme_classic()+
    facet_grid(SAMPLING_CATEGORY~.)+
    theme(legend.position="bottom")
  ggsave(file=paste0(outfile.base,'/sampling_fractions_fit.pdf'),width = 6, height = 1.2*70,
         limitsize = FALSE)
  
  # check
  dprior.fit.plot.quantile <- dprior.fit.plot[,list(
    PLB = quantile(P,0.025),
    PUB = quantile(P,0.975),
    PM = quantile(P,0.5)
  ),by=c('SAMPLING_CATEGORY','TYPE')]
  
  ggplot(dprior.fit.plot.quantile, aes(x=SAMPLING_CATEGORY, y=PM, 
                                       group=TYPE, color=TYPE)) + 
    geom_point(position=position_dodge(0.9))+
    geom_errorbar(aes(ymin=PLB, ymax=PUB), width=.2,
                  position=position_dodge(0.9))+
    theme_classic()+
    theme(legend.position="bottom")+ 
    coord_flip() +
    labs(x='age and gender groups \n', y='\n sampling probabilities',color='type')
  ggsave(file=paste0(outfile.base,'/sampling_fractions_fit_ci.pdf'),width = 4, height = 10)
  
  dprior.fit.plot.quantile <- dcast(dprior.fit.plot.quantile,
                                    SAMPLING_CATEGORY~TYPE, value.var=c('PM','PLB','PUB'))
  dprior.fit.plot.quantile.diff <- dprior.fit.plot.quantile[,list(O = max(PUB_fit,PUB_prior)-min(PLB_fit,PLB_prior)-
                                                                    min(PUB_fit,PUB_prior)+max(PLB_fit,PLB_prior)),
                                                            by='SAMPLING_CATEGORY']
  # SAMPLING_CATEGORY            O
  # 1:           F:16-17 0.0257351133
  # 2:           F:17-18 0.0288658728
  # 3:           F:18-19 0.0018241417
  # 4:           F:19-20 0.0069969013
  # 5:           F:20-21 0.0029778113
  # 6:           F:21-22 0.0019753799
  # 7:           F:22-23 0.0049651749
  # 8:           F:23-24 0.0039252678
  # 9:           F:24-25 0.0025458977
  # 10:           F:25-26 0.0029932017
  # 11:           F:26-27 0.0031336749
  # 12:           F:27-28 0.0039138612
  # 13:           F:28-29 0.0027128447
  # 14:           F:29-30 0.0030755388
  # 15:           F:30-31 0.0018986706
  # 16:           F:31-32 0.0056885963
  # 17:           F:32-33 0.0030613273
  # 18:           F:33-34 0.0025019226
  # 19:           F:34-35 0.0033655579
  # 20:           F:35-36 0.0100678373
  # 21:           F:36-37 0.0079761620
  # 22:           F:37-38 0.0026454563
  # 23:           F:38-39 0.0038065852
  # 24:           F:39-40 0.0076068994
  # 25:           F:40-41 0.0075368182
  # 26:           F:41-42 0.0068427402
  # 27:           F:42-43 0.0069261303
  # 28:           F:43-44 0.0006046767
  # 29:           F:44-45 0.0031918049
  # 30:           F:45-46 0.0128447123
  # 31:           F:46-47 0.0056958667
  # 32:           F:47-48 0.0039531433
  # 33:           F:48-49 0.0087503231
  # 34:           F:49-50 0.0089547945
  # 35:           F:50-51 0.0119370597
  # 36:           M:16-17 0.0288795777
  # 37:           M:17-18 0.0234119617
  # 38:           M:18-19 0.0069168122
  # 39:           M:19-20 0.0056760381
  # 40:           M:20-21 0.0036523345
  # 41:           M:21-22 0.0047503859
  # 42:           M:22-23 0.0028848129
  # 43:           M:23-24 0.0032446476
  # 44:           M:24-25 0.0006390216
  # 45:           M:25-26 0.0012226868
  # 46:           M:26-27 0.0022026948
  # 47:           M:27-28 0.0006101664
  # 48:           M:28-29 0.0017168258
  # 49:           M:29-30 0.0037392963
  # 50:           M:30-31 0.0024530950
  # 51:           M:31-32 0.0050381700
  # 52:           M:32-33 0.0009600297
  # 53:           M:33-34 0.0019649197
  # 54:           M:34-35 0.0033922763
  # 55:           M:35-36 0.0014809766
  # 56:           M:36-37 0.0089759653
  # 57:           M:37-38 0.0020885551
  # 58:           M:38-39 0.0050795669
  # 59:           M:39-40 0.0017969722
  # 60:           M:40-41 0.0092960226
  # 61:           M:41-42 0.0046218961
  # 62:           M:42-43 0.0112116878
  # 63:           M:43-44 0.0083687619
  # 64:           M:44-45 0.0126878636
  # 65:           M:45-46 0.0063629972
  # 66:           M:46-47 0.0065972156
  # 67:           M:47-48 0.0143950763
  # 68:           M:48-49 0.0089652537
  # 69:           M:49-50 0.0155531349
  # 70:           M:50-51 0.0163946721
  
  # beta distribution fits data well
}

Rakai.age.adjustment<- function(  infile.inference=NULL,
                                  opt=NULL
)
{
  require(data.table)	
  require(phyloflows)
  require(rstan)
  
  #
  #	input args
  #
  if(is.null(opt))
  {
    opt									<- list()
    opt$adjust.sequencing.bias			<- 1
    opt$adjust.participation.bias		<- 1
    opt$migration.def.code				<- '24'
    opt$set.missing.migloc.to.inland	<- 0
    opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
  }
  if(is.null(infile.inference))
  {
    infile.inference	<- "/Users/xx4515/Desktop/data/RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
  }

  cat('\ninfile.inference=',infile.inference)
  cat('\ninfile.participation.prior.samples=',infile.participation.prior.samples)
  cat('\ninfile.sequencing.prior.samples=',infile.sequencing.prior.samples)
  cat('\nopt=',unlist(opt))			
  indir					<- dirname(infile.inference)	
  outfile.base			<- '/Users/xx4515/Desktop/gp2/rakai/gpa_nonzero_adj/'
  load(infile.inference)
  
  #
  #	prepare data on observed transmission flows
  #
  #	subset to variables needed, using RTR3	
  rtr	<- copy(rtr3)
  if(opt$migration.def.code=='06')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
    setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='12')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
    setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='24')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
    setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='36')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
    setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='48')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
    setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
  }
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
                              'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
                              'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
                              'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
  # inmigrant status
  rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
  rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
  set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
  
  #	set unknown origin to either fishing or inland
  tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
  if(opt$set.missing.migloc.to.inland)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
  }		
  if(opt$set.missing.migloc.to.fishing)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
  }
  
  
  # add age 
  rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
  rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
  
  # impute age
  tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
  set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
  tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
  set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
  
  # fixup from latest surveillance data
  set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
  set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
  set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
  set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
  
  #	stratify age
  range(c(rtr$TR_AGE_AT_MID,rtr$REC_AGE_AT_MID)) #16.867 50.458
  rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(0,17:50,60), labels=paste0(16:50,'-',17:51), right=FALSE))]
  rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(0,17:50,60), labels=paste0(16:50,'-',17:51), right=FALSE))]
  stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
  stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
  
  # define TR_COMM_TYPE_F, REC_COMM_TYPE_F (i: inland; f: fishing) 
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # define TR_COMM_TYPE_F_MIG (i: inland; f: fishing; e: external) 
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_SEX,':',TR_AGE_AT_MID_C)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_SEX,':',TR_AGE_AT_MID_C)]
  
  # make all combinations of variables
  dac <- expand.grid( SEX= c('M','F'),
                      AGE_AT_MID_C= paste0(16:50,'-',17:51))
  
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(SEX, ':', AGE_AT_MID_C)]  					
  dac <- as.data.table(expand.grid(TR_CATEGORY= dac$CATEGORY, REC_CATEGORY= dac$CATEGORY))
  # ignore Male-Male and Female-Female combinations
  dac <- subset(dac, !(grepl('F',TR_CATEGORY)&grepl('F',REC_CATEGORY)) &
                  !(grepl('M',TR_CATEGORY)&grepl('M',REC_CATEGORY))  
  )
  # add transmission categories
  dac[, REC_TRM_CATEGORY:= REC_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_CATEGORY] 
  dac[, REC_SAMPLING_CATEGORY:=REC_TRM_CATEGORY]
  dac[, TR_SAMPLING_CATEGORY:=TR_TRM_CATEGORY]
  
  #	calculate observed number of transmissions
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  dac[, DUMMY:= 1]
  dobs <- merge(dac, dobs, by=c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'), all=TRUE)
  stopifnot( dobs[, !any(is.na(DUMMY))] )
  set(dobs, NULL, 'DUMMY', NULL)
  set(dobs, dobs[, which(is.na(TRM_OBS))], 'TRM_OBS', 0L)
  
  # take f-m transmissions
  dobs <- dobs[grepl('F',TR_TRM_CATEGORY),]
  #	make PAIR_ID
  setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY,TR_SAMPLING_CATEGORY,REC_SAMPLING_CATEGORY)	
  dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
  dobs[,TR_SMOOTH_CATEGORY:= as.numeric(substr(TR_TRM_CATEGORY,3,4))+0.5]
  dobs[,REC_SMOOTH_CATEGORY:= as.numeric(substr(REC_TRM_CATEGORY,3,4))+0.5]
  
  load("~/Desktop/gp2/rakai/samples_fit.rda")
  dprior.fit
  shape <- cbind(dprior.fit$SHAPE1,dprior.fit$SHAPE2)
  dprior.id <- subset(dprior.fit,select = 'SAMPLING_CATEGORY')
  setkey(dprior.id, SAMPLING_CATEGORY)
  dprior.id[,ID:= seq_len(nrow(dprior.id))]
  
  tmp <- subset(dobs, select = c('TR_SAMPLING_CATEGORY',
                                 'REC_SAMPLING_CATEGORY',
                                 'TRM_CAT_PAIR_ID'))
  
  setnames(dprior.id,colnames(dprior.id),paste0('TR_',colnames(dprior.id)))
  tmp <- merge(tmp,dprior.id,by='TR_SAMPLING_CATEGORY',all.x = TRUE)
  setnames(dprior.id,colnames(dprior.id),gsub('TR_','REC_',colnames(dprior.id)))
  tmp <- merge(tmp,dprior.id,by='REC_SAMPLING_CATEGORY',all.x = TRUE)
  setkey(tmp, TRM_CAT_PAIR_ID)
  xi_id <- cbind(tmp$TR_ID, tmp$REC_ID)
  
  #	estimate GP and hyperparameters with informative pior on length scale
  range(dist(dobs$TR_SMOOTH_CATEGORY))
  range(dist(dobs$REC_SMOOTH_CATEGORY))
  
  range(dist(dobs[TRM_OBS!=0,TR_SMOOTH_CATEGORY]))
  range(dist(dobs[TRM_OBS!=0,REC_SMOOTH_CATEGORY]))
  
  # [0,30]
  fit <- stan(file=file.path(outfile.base,'gp_prior_tune_ard_nonzero.stan'), 
              iter=1, warmup=0, chains=1,
              seed=5838298, algorithm="Fixed_param")
  # a = 1.78207
  # b = 3.11324
  
  # # [0,15]
  # fit <- stan(file=file.path(outfile.base,'gp_prior_tune_ard_div.stan'), 
  #             iter=1, warmup=0, chains=1,
  #             seed=5838298, algorithm="Fixed_param")
  # a = 2.38497
  # b = 3.6696
  
  M <- 20
  D <- 2 
  indices <- matrix(NA, M^D, D)
  mm=0;
  for (m1 in 1:M){
    for (m2 in 1:M){
      mm = mm+1
      indices[mm,] = c(m1, m2)
    }
  }
  
  standata_bf12 <- list( M= M, M_nD= M^D, 
                         L= c(5/2*max(dobs$TR_SMOOTH_CATEGORY),5/2*max(dobs$REC_SMOOTH_CATEGORY)), 
                         N = nrow(dobs),
                         x = cbind(dobs$TR_SMOOTH_CATEGORY,dobs$REC_SMOOTH_CATEGORY),
                         D = D,
                         y = dobs$TRM_OBS,
                         indices= indices,
                         N_xi = nrow(dprior.fit),
                         shape = shape,
                         xi_id = xi_id)
  
  fit <- stan(file = file.path(outfile.base,'gpa_ard_nonzero_adj.stan'),
              data = standata_bf12,
              iter = 10000,  warmup = 2000, chains=1, thin=1,
              algorithm = "NUTS", verbose = FALSE,
              control = list(adapt_delta = 0.99))
}