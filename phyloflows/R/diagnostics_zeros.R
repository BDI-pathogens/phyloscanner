#' @title MCMC diagnostics for the source attribution algorithm
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom bayesplot mcmc_trace mcmc_acf_bar mcmc_hist
#' @importFrom coda mcmc effectiveSize
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param mcmc.file Full file name to MCMC output from function \code{source.attribution.mcmc}
#' @param mc MCMC object from function \code{source.attribution.mcmc}; this can be specified alternatively to `mcmc.file`.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"burnin.p"}{Proportion of MCMC iterations that are removed as burn-in period.}
#'  \item{"regex_pars"}{Regular expression to select the parameter names for which diagnostics are computed. The default is all parameters, which is '*'.}
#'  \item{"credibility.interval"}{Width of the marginal posterior credibility intervals.}
#'  \item{"pdf.plot.n.worst.case.parameters"}{Integer which specifies the number of parameters with smallest effective sample size that are inspected in detail. If set to 0, worst case analyses are not performed.}
#'  \item{"pdf.plot.all.parameters"}{Flag which specifies if traces shall be plotted for all parameters. If set to TRUE, very large pdfs may be created.}
#'  \item{"pdf.height.per.par"}{Most plots show diagnostics with parameters listed on the y-axis. This value controls the plot height in inches for each free parameter.}
#'  \item{"outfile.base"}{Start of the full file name for all output files.}
#' }
#' @return Nothing. Diagnostic plots and csv files are written to disk.
#' @seealso \link{\code{source.attribution.mcmc}}
#' @examples
#' require(data.table)
#' require(phyloflows)
#' # load input data set
#' data(twoGroupFlows1, package="phyloflows")
#' dobs <- twoGroupFlows1$dobs
#' dprior <- twoGroupFlows1$dprior
#' # load MCMC output obtained with
#' # control <- list(seed=42, mcmc.n=5e4, verbose=0)
#' # mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
#' data(twoGroupFlows1_mcmc, package="phyloflows")
#'
#' #
#' # calculate standard diagnostics:
#' # trace plots, histograms, autocorrelation plots,
#' # remove burn-in, thin,
#' # calculate effective sample size,
#' # get mean, 50% quantile, and credibility intervals for marginal posterior densities
#' #
#' \dontrun{
#' outfile.base <- 'mcmc_diagnostics_'
#' control <- list( burnin.p=0.05,                # percentage of chain to remove as burn-in
#'                  regex_pars='*',                # '*' for all parameters, 'PI' for the transmission flow proportions
#'                  credibility.interval=0.95,     # width of credibility intervals
#'                  pdf.plot.all.parameters=FALSE, # plot traces for all parameters selected with \code{regex_pars}
#'                  pdf.plot.n.worst.case.parameters=10, # plot traces, histograms, autocorellation plots for the 10 parameters with lowest effective sample size. Set this to 0 if you don t want to plot these.
#'                  pdf.height.per.par=1.2,        # change this to adjust the height of pdf figures
#'                  outfile.base=outfile.base      # file name base for all output
#'                  )
#' phyloflows:::source.attribution.mcmc.diagnostics(mc=mc, control=control)
#' }
#'
source.attribution.mcmc.diagnostics    <- function(mcmc.file, mc=NULL, control=list(burnin.p=0.2, regex_pars='*', credibility.interval=0.95, pdf.plot.n.worst.case.parameters=10, pdf.plot.all.parameters=FALSE, pdf.height.per.par=1.2, outfile.base=gsub('\\.rda','',mcmc.file))){
    #
    #library(coda); library(data.table); library(bayesplot); library(ggplot2)
    if(!is.null(mc))
    {
        cat('\nUsing MCMC output specified as input...')
    }
    if(is.null(mc))
    {
        cat('\nLoading MCMC output from file...')
		cat('\nLoading ', mcmc.file[1])
        load(mcmc.file[1])
        XI <- mc$pars$XI
        S <- mc$pars$S
        LOG_LAMBDA <- mc$pars$LOG_LAMBDA
        it.info <- mc$it.info
        PI <-  t(apply(exp(mc$pars$LOG_LAMBDA), 1, function(rw) rw/sum(rw)))
        
        if (length(mcmc.file)>1)
		{
            for (k in 2:length(mcmc.file))
			{
				cat('\nLoading ', mcmc.file[k])
                load(mcmc.file[k])
                XI <- rbind(XI, mc$pars$XI)
                S <- rbind(S, mc$pars$S)
                LOG_LAMBDA <- rbind(LOG_LAMBDA, mc$pars$LOG_LAMBDA)
                it.info <- rbind(it.info, mc$it.info)
                PI <- rbind(PI, t(apply(exp(mc$pars$LOG_LAMBDA), 1, function(rw) rw/sum(rw)))    )
            }			
        }
        
        mc$pars$XI <- XI
        mc$pars$S <- S
        mc$pars$LOG_LAMBDA <- LOG_LAMBDA
        mc$it.info <- it.info
        mc$pars$PI <- PI
		cat('\nTotal number of MCMC iterations loaded from file ', nrow(mc$pars$LOG_LAMBDA))
		XI <- S <- LOG_LAMBDA <- PI <- NULL
		gc()
    }
    
    
    burnin.n <- floor(control$burnin.p*nrow(mc$pars$LOG_LAMBDA))
    
    cat('\nCollecting parameters...')
    pars <- matrix(NA,nrow=nrow(mc$pars$XI),ncol=0)
    if(grepl(control$regex_pars,'XI'))
    {
        tmp <- mc$pars$XI
        colnames(tmp)    <- paste0('XI-',1:ncol(tmp))
        pars <- cbind(pars, tmp)
    }
    if(grepl(control$regex_pars,'LOG_LAMBDA'))
    {
        tmp <- mc$pars$LOG_LAMBDA
        colnames(tmp) <- paste0('LOG_LAMBDA-',1:ncol(tmp))
        pars <- cbind(pars, tmp)
		
    }
	if(grepl(control$regex_pars,'PI') & ('PI' %in% names(mc$pars)))
	{
<<<<<<< HEAD
		tmp <- exp(mc$pars$LOG_LAMBDA)
		tmp <- t(apply(tmp, 1, function(rw) rw/sum(rw)))		
=======

		tmp <- mc$pars$PI
>>>>>>> phyloflow_v120
		colnames(tmp) <- paste0('PI-',1:ncol(tmp))
		pars <- cbind(pars, tmp)		
	}
	
    #    traces for parameters
    if(control$pdf.plot.all.parameters)
    {
        cat('\nPlotting traces for all parameters...')
        p <- mcmc_trace(pars, pars=colnames(pars), facet_args = list(ncol = 1), n_warmup=burnin.n)
        pdf(file=paste0(control$outfile.base,'_marginaltraces.pdf'), w=7, h=control$pdf.height.per.par*ncol(pars))
        print(p)
        dev.off()
    }
    
    #    traces for log likelihood and log posterior
    if(1)
    {
        pars2 <- as.matrix(subset(mc$it.info, BLOCK=='LOG_LAMBDA', select=c(LOG_LKL, LOG_PRIOR)))
        pars2[,2] <- pars2[,1]+pars2[,2]
        colnames(pars2) <- c('log likelihood','log posterior')
        cat('\nPlotting traces for log likelihood and log posterior...')
        p <- mcmc_trace(pars2, pars=colnames(pars2), facet_args = list(ncol = 1), n_warmup=burnin.n)
        pdf(file=paste0(control$outfile.base,'_loglklpotrace.pdf'), w=10, h=control$pdf.height.per.par*ncol(pars2)*2)
        print(p)
        dev.off()
		cat('\nPlotting histograms for log likelihood and log posterior...')
        p <- mcmc_hist(pars2, pars=colnames(pars2), facet_args = list(ncol=4))
        pdf(file=paste0(control$outfile.base,'_loglklpohist.pdf'), w=10, h=control$pdf.height.per.par*ncol(pars2))
		suppressMessages(print(p))
        dev.off()
    }
    
    #    acceptance rate per MCMC update ID
    cat('\nCalculating acceptance rates...')
    da <- subset(mc$it.info, BLOCK=='XI' & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
    setnames(da, 'PAR_ID', 'UPDATE_ID')
    # tmp <- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
    # da <- merge(da, tmp, by='UPDATE_ID')
    # ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) +
    # 	geom_point() +
    # 	theme_bw() +
    # 	scale_y_continuous(label=scales::percent) +
    # 	labs(    x='\nNumber of transmission pair categories updated per sampling category',
    # 	y='Acceptance rate\n')
    # ggsave(file=paste0(control$outfile.base,'_acceptance_per_updateID.pdf'), w=6, h=6)
    cat('\nAverage acceptance rate= ',subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, round(mean(ACCEPT), d=3)])
    cat('\nUpdate IDs with lowest acceptance rates')
    print( da[order(ACC_RATE)[1:min(10, nrow(da))],] )
    
    # remove burn-in
    cat('\nRemoving burnin in set to ', 100*control$burnin.p,'% of chain, corresponding to the first iterations=',burnin.n)
    tmp    <- seq.int(burnin.n+1L,nrow(mc$pars$LOG_LAMBDA))
    pars    <- pars[tmp,,drop=FALSE]
    
    # effective sampling sizes
    cat('\nCalculating effective sample size for all parameters...')
    tmp    <- mcmc(pars)
    ans    <- data.table(ID= seq_len(ncol(pars)), VAR= colnames(pars), NEFF=as.numeric(effectiveSize(tmp)))
    set(ans, NULL, 'ID', ans[, factor(ID, labels=VAR)])
    if(control$pdf.plot.all.parameters)
    {
        ggplot(ans, aes(x=NEFF, y=ID)) +
        geom_point() +
        theme_bw() +
        labs(x='\neffective sample size', y='')
        ggsave(file=paste0(control$outfile.base,'_neff.pdf'), w=6, h=control$pdf.height.per.par*ncol(pars)*0.15, limitsize=FALSE)
    }
    
    # summarise mean, sd, quantiles
    cat('\nCalculating posterior summaries for all parameters...')
    tmp    <- apply(pars, 2, function(x) quantile(x, p=c((1-control$credibility.interval)/2, 0.5, control$credibility.interval+(1-control$credibility.interval)/2)))
    tmp    <- data.table(   VAR= colnames(pars),
    						MEAN= apply(pars, 2, mean),
    						SD= apply(pars, 2, sd),
    						MEDIAN=tmp[2,],
    						CI_L=tmp[1,],
    						CI_U=tmp[3,])
    ans    <- merge(tmp, ans, by='VAR')
    cat('\nSummary of parameters with lowest effective samples\n')
    print( ans[order(NEFF)[1:min(10, nrow(ans))],] )
    
    # write to file
    cat('\nWriting summary file to',paste0(control$outfile.base,'_summary.csv'))
    setkey(ans, ID)
    write.csv(ans, file=paste0(control$outfile.base,'_summary.csv'))
    
    # plots for worst case parameters
    if(control$pdf.plot.n.worst.case.parameters>0)
    {
		control$pdf.plot.n.worst.case.parameters <- min(control$pdf.plot.n.worst.case.parameters, nrow(ans))
		control$pdf.plot.n.worst.case.parameters <- max(control$pdf.plot.n.worst.case.parameters, 2)
        worst.pars <- pars[, ans[order(NEFF)[1:control$pdf.plot.n.worst.case.parameters], VAR]]
        #    traces
        cat('\nPlotting traces for worst parameters...')
        p <- mcmc_trace(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=1))
        pdf(file=paste0(control$outfile.base,'_worst_traces.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
        print(p)
        dev.off()
        
        #    histograms
        cat('\nPlotting marginal posterior densities for worst parameters...')
        p <- mcmc_hist(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=4))
        pdf(file=paste0(control$outfile.base,'_worst_marginalposteriors.pdf'), w=10, h=control$pdf.height.per.par*ncol(worst.pars)/4)
		suppressMessages(print(p))
        dev.off()
        
        #    autocorrelations
        cat('\nPlotting autocorrelations for worst parameters...')
        p <- mcmc_acf_bar(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol = 1))
        pdf(file=paste0(control$outfile.base,'_worst_acf.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
        print(p)
        dev.off()
    }
}

