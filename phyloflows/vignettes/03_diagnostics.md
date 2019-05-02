03 - Performing MCMC diagnostic checks
================
2019-04-30

Here we show you how to run a number of diagnostics on **phyloflows**
MCMC output. Please read the sections “Our Job” and “Our Solution” on
the main page before you go ahead here.

## Getting started

We continue our Very\_First\_Example from the vignette with the same
name here. The following code chunk contains all code needed, up to
running **phyloflows** MCMC routine. The only change is that the number
of iterations is now \(50,000\). The MCMC should take about 2 minutes to
run.

``` r
require(data.table)
require(phyloflows)

data(twoGroupFlows1, package="phyloflows")
dobs <- twoGroupFlows1$dobs
dprior <- twoGroupFlows1$dprior
mcmc.file <- file.path(getwd(),'twoGroupFlows1_mcmc.RData')
control <- list(seed=42, mcmc.n=5e4, verbose=0)
mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
```

## MCMC: diagnostics

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

``` r
outfile.base <- file.path(getwd(),'twoGroupFlows1_mcmc_') 
control <- list( burnin.p=0.05, 
                 regex_pars='*', 
                 credibility.interval=0.95, 
                 pdf.plot.all.parameters=FALSE, 
                 pdf.plot.n.worst.case.parameters=10, 
                 pdf.height.per.par=1.2, 
                 outfile.base=outfile.base)
phyloflows:::source.attribution.mcmc.diagnostics(mc=mc, control=control)
```

That’s it for now. Use your usual R wizadry to process the output
further, and have a look at the other vignettes.
