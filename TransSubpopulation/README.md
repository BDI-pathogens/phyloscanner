# TransSubpopulation

*TransSubpopulation* analyses transmissions between subpopulations while considering variation in sampling intensity. Inputs for this analysis are transmission networks inferred from [phyloscanner](https://github.com/BDI-pathogens/phyloscanner) and demographic data. The software package comprises

1. R functions for implementing Markov Chain Monte Carlo to estimate probabilities of transmission between subpopulations.
2. R functions for analysing outputs from MCMC including diagnostics, aggregation to the transmission flows of interest and key statistics on the aggregated flows.


### Installing

You can install TransSubpopulation in R

```{r}
install_github("BDI-pathogens/phyloscanner/TransSubpopulation", dependencies=FALSE, build_vignettes=FALSE)
require(TransSubpopulation)
```

This step takes less than XX minutes. If you have issues with installation/running of TransSubpopulation, please report it here and we will get back to you.

## Tutorials for transmission flow analysis.

Transmission networks resolved by phyloscanner could be analysed through TransSubpopulation:

1. generate samples from the prior distribution of the sampling probabilities under SARWS/GLM assumptions. 
   * if you have the number of eligible and sampled individuals, please use SARWS assumption.
   * if you have additional individual covariates data, please use GLM assumption.
2. implement MCMC algorithm through *source.attribution.mcmc*.
3. apply diagnostic tools for the MCMC algorithm through *source.attribution.mcmc.diagnostics*. 
4. aggreagate to the transmission flows of interests through *source.attribution.mcmc.aggregateToTarget*.
5. extract key statistics on the aggregated flows through *source.attribution.aggmcmc.getKeyQuantities*.

We provide two tutorials on simulated data to to demonstrate analysis of reconstructed networks. 

**Source_attribution_under_SARWS_assumption**
[This tutorial](Source_attribution_under_SARWS_assumption.html) uses the generated transmission networks to understand transmissions between subpopulations under the assumption that individuals are sampled at random within subpopulation. 

**Source_attribution_under_GLM_assumption**
[This tutorial](Source_attribution_under_GLM_assumption.html) uses the generated transmission networks to understand transmissions between subpopulations under the assumption that the probability of an individual being sampled depends on his/her characteristics.

