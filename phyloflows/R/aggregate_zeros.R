#' @title Aggregate MCMC output to target parameters
#' @export
#' @import data.table
#' @param mcmc.file Full file name to MCMC output from function \code{source.attribution.mcmc}
#' @param mc MCMC object from function \code{source.attribution.mcmc}; this can be specified alternatively to `mcmc.file`.
#' @param daggregateTo Data.table that maps the categories of transmission pairs used in the MCMC to lower-dimensional categories that are of primary interest.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"burnin.p"}{Proportion of MCMC iterations that are removed as burn-in period.}
#'  \item{"thin"}{Thin MCMC output to every nth MCMC iteration.}
#'  \item{"regex_pars"}{Regular expression to select the parameters that are to be aggregated.}
#'  \item{"outfile"}{Full file name of the output csv file.}
#' }
#' @return If `outfile` is not specified, aggregated Monte Carlo samples are returned, as a data.table. If `outfile` is specified and ends in `csv`, the output is written to a csv file. If `outfile` is specified and ends in `rda`, the output is written to a rda file.
#' @seealso \link{\code{source.attribution.mcmc}}, \link{\code{source.attribution.mcmc.getKeyQuantities}}
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
#' # aggregate MCMC output and return output to R
#' #
#' daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
#' daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]#>    TRM_CAT_PAIR_ID TR_TRM_CATEGORY REC_TRM_CATEGORY TR_TARGETCAT
#' daggregateTo[, REC_TARGETCAT:= 'Any']
#' set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
#' control <- list( burnin.p=0.05, thin=NA_integer_, regex_pars='PI')
#' pars <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mc=mc, daggregateTo=daggregateTo, control=control)#>
#'
#' #
#' # if the source attribution model is complex
#' # the MCMC output tends to be large, and it is
#' # better to load the MCMC output internally, so
#' # pass a file name to the function
#' #
#' \dontrun{
#' mcmc.file <- 'mcmc_output.rda'  # MCMC output
#' daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
#' daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]
#' daggregateTo[, REC_TARGETCAT:= 'Any']
#' set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
#' control <- list( burnin.p=0.05, thin=NA_integer_, regex_pars='PI')
#' pars <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mcmc.file=mcmc.file, daggregateTo=daggregateTo, control=control)
#' }
#'
#' #
#' # you can also write the aggregated output to file
#' #
#' \dontrun{
#' mcmc.file <- 'mcmc_output.rda'  # MCMC output
#' agg.mcmc.file <-  'agg_mcmc_output.rda' # store results here; this can also be a csv file
#' daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
#' daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]
#' daggregateTo[, REC_TARGETCAT:= 'Any']
#' set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
#' control <- list( burnin.p=0.05, thin=NA_integer_, regex_pars='PI', outfile=agg.mcmc.file)
#' pars <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mcmc.file=mcmc.file, daggregateTo=daggregateTo, control=control)
#' }
source.attribution.mcmc.aggregateToTarget    <- function(mcmc.file=NULL, mc=NULL, daggregateTo=NULL, control=list(burnin.p=NA_real_, thin=NA_integer_, regex_pars='*', outfile=gsub('\\.rda','_aggregated.csv',mcmc.file))){
    #    basic checks
    if(!'data.table'%in%class(daggregateTo))
        stop('daggregateTo is not a data.table')
    if(!all(c('TRM_CAT_PAIR_ID','TR_TARGETCAT','REC_TARGETCAT')%in%colnames(daggregateTo)))
        stop('daggregateTo does not contain one of the required columns TRM_CAT_PAIR_ID, TR_TARGETCAT, REC_TARGETCAT')
    
    if(!is.null(mc))
    {
        cat('\nUsing MCMC output specified as input...')
    }
    if(is.null(mc))
    {
        #    load MCMC output
        cat('\nLoading MCMC output from file...')
        cat('\nLoading ', mcmc.file[1])
        load(mcmc.file[1])
        LOG_LAMBDA <- mc$pars$LOG_LAMBDA
        
        if (length(mcmc.file)>1)
        {
            for (k in 2:length(mcmc.file))
            {
                cat('\nLoading ', mcmc.file[k])
                load(mcmc.file[k])
                LOG_LAMBDA <- rbind(LOG_LAMBDA, mc$pars$LOG_LAMBDA)
            }
        }
        mc$pars$LOG_LAMBDA <- LOG_LAMBDA
        cat('\nTotal number of MCMC iterations loaded from file ', nrow(mc$pars$LOG_LAMBDA))
        LOG_LAMBDA <- NULL
        gc()
    }
        
    #    define internal control variables
    burnin.p    <- control$burnin.p
    if(is.na(burnin.p))
        burnin.p<- 0
    burnin.n    <- floor(burnin.p*nrow(mc$pars$LOG_LAMBDA))
    thin        <- control$thin
    if(is.na(thin))
        thin    <- 1
    
    #    collect parameters
    cat('\nCollecting parameters...')
    pars        <- matrix(NA,nrow=nrow(mc$pars$LOG_LAMBDA),ncol=0)
    if(grepl(control$regex_pars,'LOG_LAMBDA|LAMBDA|PI'))
    {
        tmp    <- mc$pars$LOG_LAMBDA
        colnames(tmp)    <- paste0('LOG_LAMBDA-',1:ncol(tmp))
        pars    <- cbind(pars, tmp)
    }
        
    # remove burn-in
    if(burnin.n>0)
    {
        cat('\nRemoving burnin in set to ', 100*burnin.p,'% of chain, total iterations=',burnin.n)
        tmp        <- seq.int(burnin.n,nrow(mc$pars$LOG_LAMBDA))
        pars    <- pars[tmp,,drop=FALSE]
    }
    
    # thin
    if(thin>1)
    {
        cat('\nThinning to every', thin,'th iteration')
        tmp        <- seq.int(1,nrow(pars),thin)
        pars    <- pars[tmp,,drop=FALSE]
    }
    
    cat('\nMaking aggregated MCMC output...')
    # make data.table in long format
    pars    <- as.data.table(pars)
    pars[, SAMPLE:= seq_len(nrow(pars))]
    pars    <- melt(pars, id.vars='SAMPLE')
    pars[, VARIABLE:= pars[, gsub('([A-Z]+)-([0-9]+)','\\1',variable)]]
    pars[, TRM_CAT_PAIR_ID:= pars[, as.integer(gsub('([A-Z_]+)-([0-9]+)','\\2',variable))]]
    
    # aggregate MCMC samples
    if(!all(sort(unique(pars$TRM_CAT_PAIR_ID))==sort(unique(daggregateTo$TRM_CAT_PAIR_ID))))
        warning('The transmission count categories in the MCMC output do not match the transmission count categories in the aggregateTo data table.')
    pars    <- merge(pars, daggregateTo, by='TRM_CAT_PAIR_ID')
    
    # return aggregated lambda values
    tmp1    <- pars[,{
            tmp <- min((max(value) - min(value))/2, 700)
            list(VALUE= exp(-tmp)*sum(exp(tmp+value)))
        },
        by=c('VARIABLE','TR_TARGETCAT','REC_TARGETCAT','SAMPLE')]
    # return sum of all lambda values
    tmp2    <- pars[, {
            tmp <- min((max(value) - min(value))/2, 700)
            list(SUM=exp(-tmp)*sum(exp(tmp+value)))
        },
        by=c('VARIABLE','SAMPLE')]
    
    pars <- tmp1[0,]
    if(grepl(control$regex_pars,'LAMBDA'))
    {
        pars <- rbind(pars, tmp1)
        pars[, VARIABLE:= 'LAMBDA']
    }
    if(grepl(control$regex_pars,'LOG_LAMBDA'))
    {
        tmp <- copy(tmp1)
        set(tmp, NULL, 'VALUE', log(tmp$VALUE))
        pars <- rbind(pars, tmp)
    }
    if(grepl(control$regex_pars,'PI'))
    {
        tmp <- merge(tmp1, tmp2, by=c('VARIABLE','SAMPLE'))
        tmp[, VALUE:=VALUE/SUM]
        tmp[, VARIABLE:='PI']
        set(tmp, NULL, 'SUM', NULL)
        pars <- rbind(pars, tmp)
    }
       
    # save or return
    if(!'outfile'%in%names(control))
        return(pars)
    if(grepl('csv$',control$outfile))
    {
        cat('\nWriting csv file to',control$outfile)
        write.csv(pars, row.names=FALSE, file=control$outfile)
    }
    if(grepl('rda$',control$outfile))
    {
        cat('\nSaving rda file to',control$outfile)
        save(pars, file=control$outfile)
    }
}
