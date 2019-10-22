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