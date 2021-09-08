#' @title Source attribution while adjusting for sampling bias
#' @export
#' @import data.table rstan
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param dobs Data.table of observed number of transmission events for each transmission pair category.
#' @param dprior Data.table of parameters from the prior distribution of sampling probabilities.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"seed"}{Random number seed, for reproducibility.}
#'  \item{"mcmc.n"}{Guide on the number of MCMC iterations. The actual number of iterations will be slightly larger, and a multiple of the number of iterations needed to complete one MCMC sweep through all free parameters.}
#'  \item{"method"}{Choose flow prior among 'gamma'(default) and 'gp'.}
#'  \item{"outfile"}{Full path name to which the MCMC output is saved to in RDA format, in an object called \code{fit}.}
#' }
#' @return If `outfile` is not specified, phyloflows MCMC output is returned, which is a list object. If `outfile` is specified, phyloflows MCMC output is written to an RDA file.
#' @seealso \link{\code{source.attribution.mcmc.diagnostics}}, \link{\code{source.attribution.mcmc.aggregateToTarget}}, \link{\code{source.attribution.mcmc.getKeyQuantities}}
#' @examples
#' require(data.table)
#' require(phyloflows)
#' # load input data set
#' data(twoGroupFlows1, package="phyloflows")
#' dobs <- twoGroupFlows1$dobs
#' dprior <- twoGroupFlows1$dprior
#'
#' #
#' # run MCMC and return output to R
#' # this is usually fine for simple source attribution models
#' # that don t take long to run
#' #
#' \dontrun{
#' control <- list(seed=42,    # seed for reproducibility
#'                 mcmc.n=5e4, # guide on the number of MCMC iterations, the actual number of MCMC iterations may be a little bit larger
#'                 method='gamma'   # switch on for output to each MCMC iteration
#'                 )
#' mc <- phyloflows:::source.attribution.mcmc.stan(dobs, dprior, control)
#' }
#'
#' #
#' # run MCMC and save output to file
#' # this is recommended when the source attribution is more complex
#' # so you can deploy the MCMC on a high performance cluster
#' #
#' \dontrun{
#' mcmc.file <- 'mcmcm_output.rda'
#' control <- list(seed=42,    # seed for reproducibility
#'                 mcmc.n=5e4, # guide on the number of MCMC iterations, the actual number of MCMC iterations may be a little bit larger
#'                  method='gamma',   # switch on for output to each MCMC iteration
#'                 outfile= mcmc.file # store MCMC output as rda file here
#'                 )
#' mc <- phyloflows:::source.attribution.mcmc.stan(dobs, dprior, control)
#' }
#'
source.attribution.mcmc.stan  <- function(dobs, dprior.fit, control=list(seed=42, mcmc.n=1e3, method='gamma', outfile='StanMCMC.rda'))
{
  pkg.dir <- system.file(package = "phyloflows" )
  setkey(dprior.fit,  SAMPLING_CATEGORY)
  dprior.fit[,ID:=seq_len(nrow(dprior.fit))]
  tmp = subset(dprior.fit,select = c('SAMPLING_CATEGORY','ID'))
  setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
  dobs <- merge(dobs,tmp,by='TR_SAMPLING_CATEGORY')
  setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
  dobs <- merge(dobs,tmp,by='REC_SAMPLING_CATEGORY')
  setnames(tmp, colnames(tmp), gsub('REC_','',colnames(tmp)))
  if(is.null(control$method))control$method='gamma'
  if(control$method=='gamma'){
    data.gamma <- list(N=nrow(dobs),
                       Y=dobs$TRM_OBS,
                       N_xi = nrow(dprior.fit),
                       shape = cbind(dprior.fit$ALPHA,dprior.fit$BETA),
                       xi_id = cbind(dobs$TR_ID,dobs$REC_ID),
                       alpha = 0.8/nrow(dobs))
    
    fit<- stan(file = file.path(paste0(pkg.dir,'/misc'),'gamma.stan'),
                      data = data.gamma,
                      iter = control$mcmc.n,  warmup = min(500,floor(control$mcmc.n/2)), chains=1, thin=1, seed = control$seed,
                      algorithm = "NUTS", verbose = FALSE,
                      control = list(adapt_delta = 0.8, max_treedepth=10))
  }else if(control$method=='gp'){
    M <- 30
    D <- 2
    indices <- matrix(NA, M^D, D)
    mm=0;
    for (m1 in 1:M){
      for (m2 in 1:M){
        mm = mm+1
        indices[mm,] = c(m1, m2)
      }
    }
    
    data.gp <- list( M= M, M_nD= M^D,
                     L= c(3/2*max(dobs$TR_SMOOTH_CATEGORY),3/2*max(dobs$REC_SMOOTH_CATEGORY)),
                     N = nrow(dobs),
                     x = cbind(dobs$TR_SMOOTH_CATEGORY,dobs$REC_SMOOTH_CATEGORY),
                     D = D,
                     y = dobs$TRM_OBS,
                     indices= indices,
                     N_xi = nrow(dprior.fit),
                     shape = cbind(dprior.fit$ALPHA,dprior.fit$BETA),
                     xi_id = cbind(dobs$TR_ID,dobs$REC_ID))
    
    fit <- stan(file =file.path(paste0(pkg.dir,'/misc'),'gp.stan'),
                   data = data.gp,
                   iter = control$mcmc.n,  warmup = min(500,floor(control$mcmc.n/2)), chains=1, thin=1, seed = control$seed,
                   algorithm = "NUTS", verbose = FALSE,
                   control = list(adapt_delta = 0.8, max_treedepth=10))
    
  }
  if(!'outfile'%in%names(control))
    return(fit)
  else
    save(object = fit, file = paste0(control$method,'_',control$outfile))

}
