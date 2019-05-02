04 - Aggregating MCMC output to a new set of variables
================
2019-04-30

Here we show you how to aggregate estimated transmission flows to a
higher level. Please read the sections “Our Job” and “Our Solution” on
the main page before you go ahead here.

## Getting started

We continue our Very\_First\_Example from the vignette with the same
name here. The following code chunk contains all code needed, up to
running **phyloflows** MCMC routine. The only change is that the number
of iterations is now
 $50,000$
. The MCMC should take about 2 minutes to run.

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

## Aggregating flows

Why would it be useful to aggregated the estimated transmission flows?
Let us suppose that group “1” are the individuals aged 15-24 and group
“2” are the individuals aged 25 or older in a population. The
estimated flow vector

$$
\pi=(\pi_{11},\pi_{12},\pi_{21},\pi_{22})
$$


describes the transsmission flow within and between the two age
categories. But what is the overall contribution of transmissions from
individuals aged 15-24, and the overall contribution of transmissions
from individuals aged 25+? We want to estimate

$$
\eta= (\eta_{1},\eta_{2})
$$


where
 $\eta_{1}=\pi_{11}+\pi_{12}$
and\(\eta_{2}=\pi_{21}+\pi_{22}\). There are many similar scenorios
like that, and **phyloflow** has a little function to help you with that
task. The syntax is as follows.

``` r
daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]
daggregateTo[, REC_TARGETCAT:= 'Any']
set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
control <- list( burnin.p=0.05,
                 thin=NA_integer_,
                 regex_pars='PI')
mca <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mc=mc, daggregateTo=daggregateTo, control=control)
#>
#> Using MCMC output specified as input...
#> Collecting parameters...
#> Removing burnin in set to  5 % of chain, total iterations= 500
#> Making aggregated MCMC output...
mca
#>        VARIABLE TR_TARGETCAT REC_TARGETCAT SAMPLE     VALUE
#>     1:       PI            1           Any      1 0.3976094
#>     2:       PI            1           Any      2 0.4063326
#>     3:       PI            1           Any      3 0.4326280
#>     4:       PI            1           Any      4 0.4145340
#>     5:       PI            1           Any      5 0.4041081
#>    ---
#> 19000:       PI            2           Any   9498 0.6192508
#> 19001:       PI            2           Any   9499 0.6133096
#> 19002:       PI            2           Any   9500 0.5770948
#> 19003:       PI            2           Any   9501 0.6073157
#> 19004:       PI            2           Any   9502 0.6018683
```

The output is a data.table that contains the aggregated transmission
flows, and other aggregated variables depending on the value of
`control[['regex_pars]]`. In our case, we removed a burnin-period of 5%
of the MCMC chain, and did not thin the remaining iterations, yielding
about 9500 MCMC samples of the aggregated flows.

That’s it for now. Use your usual R wizadry to process the output
further, and have a look at the other vignettes.
