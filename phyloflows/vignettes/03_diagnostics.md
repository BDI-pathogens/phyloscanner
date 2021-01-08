This vignette describes how to run a number of diagnostics on
**phyloflows** MCMC output, obtained with the function
`phyloflows:::source.attribution.mcmc`. Please work through the vignette
*phyloflows: Estimating transmission flows under heterogeneous sampling
– a first example* before you go ahead here.

Getting started
---------------

We continue our “First\_Example”. The following code chunk contains all
code needed, up to running **phyloflows** MCMC routine. The only change
is that the number of iterations is now 50, 000. The MCMC should take
about 5 minutes to run.

    require(data.table)
    require(phyloflows)

    data(twoGroupFlows1, package="phyloflows")
    dobs <- twoGroupFlows1$dobs
    dprior <- twoGroupFlows1$dprior
    tmp= copy(dprior)
    tmp[,WHO:='REC_SAMPLING_CATEGORY']
    dprior[,WHO:='TR_SAMPLING_CATEGORY']
    dprior <- rbind(dprior,tmp)
    control <- list(seed=42, mcmc.n=5e4, verbose=0)
    mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)

MCMC: diagnostics
-----------------

**phyloflow** comes with a function to calculate standard MCMC
diagnostics. You can

1.  Make trace plots for all model parameters;
2.  or make trace plots for the model parameters with smallest effective
    sample size. This may be useful to avoid generating very large pdf
    files that you won t be able to open anyway.
3.  Make trace plots for values of the log likelihood and log posterior
    density.
4.  Calculate acceptance rates.
5.  Remove a burn-in period.
6.  Calculate effective sample sizes.
7.  Calculate summary statistics (mean, median, quantiles) of the
    marginal posterior densities.
8.  Plot marginal posterior densities for the model parameters with
    smallest effective sample size.
9.  Make autocorrelation plots for the model parameters with smallest
    effective sample size.

The syntax is as follows. Look up the help page for the diagnostics
function for a full explanation of the control arguments.

    outfile.base <- file.path(getwd(),'twoGroupFlows1_mcmc') 
    control <- list( burnin.p=0.05, 
                     regex_pars='*', 
                     credibility.interval=0.95, 
                     pdf.plot.all.parameters=TRUE, 
                     pdf.plot.n.worst.case.parameters=1, 
                     pdf.height.per.par=1.2, 
                     outfile.base=outfile.base)
    phyloflows:::source.attribution.mcmc.diagnostics(mc=mc, control=control)
    #> 
    #> Using MCMC output specified as input...
    #> Collecting parameters...
    #> Plotting traces for all parameters...
    #> 
    #> Plotting traces for log likelihood and log posterior...
    #> 
    #> Plotting histograms for log likelihood and log posterior...
    #> 
    #> Calculating acceptance rates...
    #> Average acceptance rate=  0.945
    #> Update IDs with lowest acceptance rates   UPDATE_ID ACC_RATE
    #> 1:         2  0.87776
    #> 2:         4  0.88032
    #> 3:         1  0.89840
    #> 4:         3  0.90496
    #> 
    #> Removing burnin in set to  5 % of chain, corresponding to the first iterations= 312
    #> Calculating effective sample size for all parameters...
    #> 
    #> Calculating posterior summaries for all parameters...
    #> Summary of parameters with lowest effective samples
    #>             VAR      MEAN         SD    MEDIAN      CI_L      CI_U           ID     NEFF
    #> 1:         XI-2 0.4500180 0.01008328 0.4498561 0.4294136 0.4705052         XI-2 3506.238
    #> 2:         XI-4 0.4500289 0.01006284 0.4500080 0.4287890 0.4704544         XI-4 3761.276
    #> 3:         XI-3 0.6001421 0.01061498 0.5999786 0.5791859 0.6208616         XI-3 3850.136
    #> 4: LOG_LAMBDA-1 5.9519630 0.08868220 5.9518514 5.7763276 6.1202251 LOG_LAMBDA-1 4180.379
    #> 5:         XI-1 0.5999227 0.01061797 0.6000620 0.5796611 0.6196614         XI-1 4549.380
    #> 6: LOG_LAMBDA-4 6.4498427 0.09417585 6.4520634 6.2655500 6.6334985 LOG_LAMBDA-4 5386.169
    #> 7: LOG_LAMBDA-2 3.9901164 0.26506317 4.0026036 3.4396043 4.4735366 LOG_LAMBDA-2 5604.403
    #> 8: LOG_LAMBDA-3 4.2871386 0.22762554 4.2940786 3.8182778 4.7077017 LOG_LAMBDA-3 5939.000
    #> 
    #> Writing summary file to /Users/xx4515/phyloscanner/phyloflows/vignettes/twoGroupFlows1_mcmc_summary.csv
    #> Plotting traces for worst parameters...
    #> 
    #> Plotting marginal posterior densities for worst parameters...
    #> 
    #> Plotting autocorrelations for worst parameters...
    #> quartz_off_screen 
    #>                 2

That’s it for now. Use your usual R wizadry to process the output
further, and have a look at the other vignettes.
