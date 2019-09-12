This vignette describes how to run a number of diagnostics on
**phyloflows** MCMC output, obtained with the function
`phyloflows:::source.attribution.mcmc`. Please work through the vignette
*phyloflows: Estimating transmission flows under heterogeneous sampling
- a first example* before you go ahead here.

Getting started
---------------

We continue our “First\_Example”. The following code chunk contains all
code needed, up to running **phyloflows** MCMC routine. The only change
is that the number of iterations is now 50, 000. The MCMC should take
about 2 minutes to run.

    require(data.table)
    require(phyloflows)

    data(twoGroupFlows1, package="phyloflows")
    dobs <- twoGroupFlows1$dobs
    dprior <- twoGroupFlows1$dprior
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

    outfile.base <- file.path(getwd(),'twoGroupFlows1_mcmc_') 
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
    #> Plotting acceptance rates...
    #> 
    #> Average acceptance rate=  0.897
    #> Update IDs with lowest acceptance rates   UPDATE_ID ACC_RATE N_TRM_CAT_PAIRS
    #> 1:         2  0.77528               4
    #> 2:         1  0.81440               4
    #> 
    #> Removing burnin in set to  5 % of chain, corresponding to the first iterations= 625
    #> Calculating effective sample size for all parameters...
    #> 
    #> Calculating posterior summaries for all parameters...
    #> Summary of parameters with lowest effective samples
    #>              VAR       MEAN         SD     MEDIAN       CI_L       CI_U
    #>  1:         XI-2 0.44995864 0.01006440 0.44980794 0.42895635 0.47010807
    #>  2:         XI-1 0.59983962 0.01059647 0.59992896 0.57969886 0.61966686
    #>  3: LOG_LAMBDA-4 6.45209566 0.09947787 6.45316086 6.25486850 6.64483695
    #>  4:         PI-4 0.55143756 0.03162361 0.55149555 0.49018444 0.61277277
    #>  5:         PI-1 0.33546857 0.02824658 0.33486974 0.28195780 0.39197102
    #>  6: LOG_LAMBDA-1 5.95319116 0.09225504 5.95354021 5.77111902 6.13208754
    #>  7: LOG_LAMBDA-3 4.28444033 0.22661016 4.29196911 3.82144791 4.70727508
    #>  8:         PI-2 0.04859381 0.01213510 0.04762408 0.02758604 0.07486730
    #>  9:         PI-3 0.06450006 0.01390632 0.06347870 0.04018642 0.09412482
    #> 10: LOG_LAMBDA-2 3.99308267 0.26076190 4.00238510 3.45173220 4.47650217
    #>               ID      NEFF
    #>  1:         XI-2  4769.454
    #>  2:         XI-1  5468.961
    #>  3: LOG_LAMBDA-4  6824.919
    #>  4:         PI-4  7521.251
    #>  5:         PI-1  7589.843
    #>  6: LOG_LAMBDA-1  8859.387
    #>  7: LOG_LAMBDA-3 11315.509
    #>  8:         PI-2 11876.000
    #>  9:         PI-3 11876.000
    #> 10: LOG_LAMBDA-2 11876.000
    #> 
    #> Writing summary file to /Users/Oliver/git/phyloscanner/phyloflows/vignettes/twoGroupFlows1_mcmc__summary.csv
    #> Plotting traces for worst parameters...
    #> 
    #> Plotting marginal posterior densities for worst parameters...
    #> 
    #> Plotting autocorrelations for worst parameters...
    #> quartz_off_screen 
    #>                 2

That’s it for now. Use your usual R wizadry to process the output
further, and have a look at the other vignettes.
