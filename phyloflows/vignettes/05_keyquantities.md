05 - Calculating sources, onward transmissions and flow ratios
================
2019-05-02

Here we explain how a number of key summary statistics of inferred
transmission flows can be easily calculated with **phyloflows**
`source.attribution.mcmc.getKeyQuantities` function. Please read the
sections “Our Job” and “Our Solution” on the main page before you go
ahead here.

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

## Sources of transmission for each group

OK, so now we have samples from posterior distribution of transmission
flows within and between the two population groups,

$$
\pi=(\pi_{11},\pi_{12},\pi_{21},\pi_{22}).
$$


One important summary statistic are the sources of transmissions into
each recipient group, defined by

$$
\eta=(\eta_{11},\eta_{21},\eta_{12},\eta_{22})
$$


where
 $\eta_{ij}=\pi_{ij} /\sum_{s}\pi_{sj}$
.

## Onward transmissions from each group

Another important summary statistic are the proportions of transmissions
that originate from each group, defined by

$$
\nu=(\nu_{11},\nu_{21},\nu_{12},\nu_{22})
$$


where
 $\nu_{ij}=\pi_{ij} /\sum_{s}\pi_{is}$
.

## Transmission flow ratios

Yes another important summary statistic are ratios of transmission
flows, defined by

$$
\rho_{ij}=\pi_{ij} /\pi_{ji}.
$$


## Calculating key quantities

**phyloflow** has a little function to help you with calculating the
above summary statistics. The basic syntax is as follows, although it is
also possible to specify a file name to MCMC output or aggregated MCMC
output, and it is also possible to input aggregated MCMC output.

``` r
control <- list(  quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                  flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
)
ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(mc=mc, dobs=dobs, control=control)
#>
#> Computing flows...
#> Computing WAIFM...
#> Computing sources...
#> Computing flow ratios...
ans
#>     TR_TARGETCAT REC_TARGETCAT         CL         CU         IL         IU
#>  1:            1             1 0.28965836 0.39241184 0.32184217 0.35724829
#>  2:            1             2 0.02862946 0.07513569 0.04062637 0.05678786
#>  3:            2             1 0.04015472 0.09450119 0.05480740 0.07355452
#>  4:            2             2 0.48788546 0.60210041 0.52679696 0.56656661
#>  5:            1             1 0.81047715 0.92417491 0.85529754 0.89445227
#>  6:            1             2 0.07582509 0.18952285 0.10554773 0.14470246
#>  7:            2             1 0.06581727 0.15511930 0.09001016 0.12010877
#>  8:            2             2 0.84488070 0.93418273 0.87989123 0.90998984
#>  9:            1             1 0.77458839 0.89742492 0.81987801 0.86278707
#> 10:            2             1 0.10257508 0.22541161 0.13721293 0.18012199
#> 11:            1             2 0.04815878 0.12517414 0.06839438 0.09534414
#> 12:            2             2 0.87482586 0.95184122 0.90465586 0.93160562
#> 13:           NA            NA 0.40033567 1.42314246 0.60269712 0.95289301
#> 14:           NA            NA 0.70267034 2.49790385 1.04943576 1.65920820
#>              M                  LABEL              LABEL2       STAT
#>  1: 0.33920017   33.9%\n[29% - 39.2%]   33.9% (29%-39.2%)      flows
#>  2: 0.04812999    4.8%\n[2.9% - 7.5%]    4.8% (2.9%-7.5%)      flows
#>  3: 0.06370609      6.4%\n[4% - 9.5%]      6.4% (4%-9.5%)      flows
#>  4: 0.54688995 54.7%\n[48.8% - 60.2%] 54.7% (48.8%-60.2%)      flows
#>  5: 0.87573515   87.6%\n[81% - 92.4%]   87.6% (81%-92.4%)      waifm
#>  6: 0.12426485    12.4%\n[7.6% - 19%]    12.4% (7.6%-19%)      waifm
#>  7: 0.10430204  10.4%\n[6.6% - 15.5%]  10.4% (6.6%-15.5%)      waifm
#>  8: 0.89569796 89.6%\n[84.5% - 93.4%] 89.6% (84.5%-93.4%)      waifm
#>  9: 0.84218401 84.2%\n[77.5% - 89.7%] 84.2% (77.5%-89.7%)    sources
#> 10: 0.15781599 15.8%\n[10.3% - 22.5%] 15.8% (10.3%-22.5%)    sources
#> 11: 0.08100214   8.1%\n[4.8% - 12.5%]   8.1% (4.8%-12.5%)    sources
#> 12: 0.91899786 91.9%\n[87.5% - 95.2%] 91.9% (87.5%-95.2%)    sources
#> 13: 0.75541129     0.76\n[0.4 - 1.42]     0.76 (0.4-1.42) flow_ratio
#> 14: 1.32378218      1.32\n[0.7 - 2.5]      1.32 (0.7-2.5) flow_ratio
#>     FLOWRATIO_CAT
#>  1:          <NA>
#>  2:          <NA>
#>  3:          <NA>
#>  4:          <NA>
#>  5:          <NA>
#>  6:          <NA>
#>  7:          <NA>
#>  8:          <NA>
#>  9:          <NA>
#> 10:          <NA>
#> 11:          <NA>
#> 12:          <NA>
#> 13:           1/2
#> 14:           2/1
```
