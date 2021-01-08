07 - Age analysis
================
Xiaoyue Xi and Oliver Ratmann
2019-12-10

This vignette provides an extension of the general method of
**phyloflow**. The aim is to understand the transmission flows between
one-year increment age groups. The differences between this example with
the general pipeline of **phyloflow** is the correlation between flows.
To tackle this problem, we impose a Gaussian process (GP) prior on
transmission flows. Also, a more flexible platform, Stan, is introduced
to solve the problem.

# Dataset

**We start with simulating transmission counts** between seven age
groups called “15-19”,“20-24”,“25-29”,“30-34”,“35-39”,“40-44”,“45-49”.
Note that in practice, it would be good to use this method to
investigate transmission dynamics between one-year increment age group,
as the squared exponential kernel is for the continuous input space. We
first set up hyperparameters for GP, the baseline and sampling
fractions.

``` r
library(rstan)
library(data.table)
library(ggplot2)
library(viridis)
library(phyloflows)
set.seed(42)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
alpha_true <- c(2.5)
rho_true <- c(12,9)
mu_true <- -1
gp_dim <- 2
xi <- c(0.35, 0.45, 0.5, 0.55, 0.5, 0.55, 0.4)
```

**Next, we calculate the sampling probabilities of transmission flows
within and between the two groups.**

``` r
dobs <- data.table(expand.grid(TR_TRM_CATEGORY = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
                               REC_TRM_CATEGORY = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")))
```

``` r
ds <- data.table(CATEGORY = c("15-19","20-24","25-29",
                              "30-34","35-39","40-44","45-49"),
                 P = xi,ID = 1:7)
```

``` r
setnames(ds,colnames(ds),paste0('TR_TRM_',colnames(ds)))
dobs <- merge(dobs,ds,by='TR_TRM_CATEGORY')
setnames(ds,colnames(ds),gsub('TR_','REC_',colnames(ds)))
dobs <- merge(dobs,ds,by='REC_TRM_CATEGORY')
setnames(ds,colnames(ds),gsub('REC_TRM_','',colnames(ds)))
dobs[,P:= TR_TRM_P * REC_TRM_P]
dobs[,TR_SMOOTH_CATEGORY:=as.numeric(substr(TR_TRM_CATEGORY,1,2))+2]
dobs[,REC_SMOOTH_CATEGORY:=as.numeric(substr(REC_TRM_CATEGORY,1,2))+2]
dobs[, TR_SAMPLING_CATEGORY:= TR_TRM_CATEGORY]
dobs[, REC_SAMPLING_CATEGORY:= REC_TRM_CATEGORY]
```

**Then we simulate true transmission flows and observed transmission
counts.**

``` r
simu_pars <- list(  N=nrow(dobs), D=gp_dim, x=cbind(dobs$TR_SMOOTH_CATEGORY,dobs$REC_SMOOTH_CATEGORY),
                   alpha=alpha_true, rho=rho_true,
                   mu=mu_true,xi=dobs$P)
#   simulate data set   
simu_fit <- stan(   file="simu_poiss.stan", 
                  data=simu_pars, iter=1,
                  chains=1, seed=424838, algorithm="Fixed_param")
dobs$TRM_OBS <- extract(simu_fit)$y[1,] 
```

The data can be loaded through

``` r
data(sevenGroupFlows1)
```

# Input data: observed transmission flows

Input data of the similar format of **phyloflow** are expected.

``` r
dobs <- subset(dobs, select = c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY',
                                'REC_SAMPLING_CATEGORY', 'TR_SMOOTH_CATEGORY','REC_SMOOTH_CATEGORY',
                                'TRM_OBS'))
head(dobs)
```

    ##    TR_TRM_CATEGORY REC_TRM_CATEGORY TR_SAMPLING_CATEGORY REC_SAMPLING_CATEGORY
    ## 1:           15-19            15-19                15-19                 15-19
    ## 2:           20-24            15-19                20-24                 15-19
    ## 3:           25-29            15-19                25-29                 15-19
    ## 4:           30-34            15-19                30-34                 15-19
    ## 5:           35-39            15-19                35-39                 15-19
    ## 6:           40-44            15-19                40-44                 15-19
    ##    TR_SMOOTH_CATEGORY REC_SMOOTH_CATEGORY TRM_OBS
    ## 1:                 17                  17       0
    ## 2:                 22                  17       1
    ## 3:                 27                  17       0
    ## 4:                 32                  17       1
    ## 5:                 37                  17       0
    ## 6:                 42                  17       0

**`dobs` specifies observed counts of transmissions from a transmitter
age group to a recipient age group.** It must contain the following
columns:

  - *TR\_TRM\_CATEGORY* name of transmitter group.
  - *REC\_TRM\_CATEGORY* name of recipient group.
  - *TR\_SMOOTH\_CATEGORY* midpoint of transmitter age group.
  - *REC\_SMOOTH\_CATEGORY* midpoint of recipient age group.
  - *TRM\_CAT\_PAIR\_ID* identifier of transmitter-recipient pair
  - *TRM\_OBS* observed transmission counts

Let us look at the data. The first row shows zero counts of transmission
flows from age group “15-19” to age group “15-19”. We visualise our
simulated input data in the heatmap :

``` r
ggplot(dobs, aes(TR_SMOOTH_CATEGORY, REC_SMOOTH_CATEGORY))+
  geom_tile(aes(fill = TRM_OBS)) +
  scale_fill_viridis() +
  scale_x_continuous(expand = c(0,0), limits = c(14.5,49.5))+
  scale_y_continuous(expand = c(0,0), limits = c(14.5,49.5))+
  theme_classic()+
  labs(x='transmitter age group \n', y='\n recipient age group',fill='transmission \n counts')
```

![](07_age_analysis_files/figure-gfm/data-1.png)<!-- -->

## Input data: sampling information

**`dobs` also requires information about how each group was sampled. **
This is stored in the following columns:

  - *TR\_SAMPLING\_CATEGORY* sampling strata of transmitter group
  - *REC\_SAMPLING\_CATEGORY* sampling strata of recipient group

**`dprior.fit` specifies the distribution of sampling probability in
each sampling group.** The sampling probability is either characterised
by beta-binomial model or GLM model (See 01 - Simulating data). and
given in the form of samples. In the tutorial, we opted for a
probabilistic programming language Stan, and it cannot take samples as
input. However, it is possible to fit statistical distributions to
samples and input distribution parameters. The information is stored in
data *dprior.fit* in the following columns:

  - *SAMPLING\_CATEGORY* name of sampling strata
  - *ALPHA, BETA* shape parameters of the distribution of sampling
    probability.

Under beta-binomial model, the posterior distribution of the sampling
probability is known analytically, i.e. beta. Let us look at the
sampling information:

``` r
ds$TRIAL <- c(4000, 3700, 3300, 2500, 1700, 1000, 500)
ds[,SUC := round(TRIAL * P)]
dprior.fit <- copy(ds)
dprior.fit[,ALPHA := SUC+1]
dprior.fit[,BETA := TRIAL-SUC+1]
dprior.fit <- subset(dprior.fit, select = c('CATEGORY','ALPHA','BETA'))
colnames(dprior.fit) <- c('SAMPLING_CATEGORY','ALPHA','BETA')
head(dprior.fit)
```

    ##    SAMPLING_CATEGORY ALPHA BETA
    ## 1:             15-19  1401 2601
    ## 2:             20-24  1666 2036
    ## 3:             25-29  1651 1651
    ## 4:             30-34  1376 1126
    ## 5:             35-39   851  851
    ## 6:             40-44   551  451

# Method

Again we use a Bayesian approach to estimate the proportion of
transmissions between the two population groups. Please see the problem
setting in **phyloflows: Estimating transmission flows under
heterogeneous sampling – a first example**. Unlike the simple example, a
probabilistic programming language Stan provides a flexible platform to
impose different kinds of priors. Recall the posterior distribution of
the parameters \((\lambda, s)\) is given by \[
\begin{aligned}
p(\lambda, s | n) & \propto p(n | \lambda, s) p(\lambda, s) \\
              & = \prod_{i=1,\cdots,7;j=1\cdots,7} Poisson(n_{ij};\lambda_{ij}*s_i*s_j) p(\lambda_{ij}) p(s_i) p(s_j).
\end{aligned}
\] Then, we calculate the main quantity of interest, \(\pi\), via \[
\pi_{ij}= \lambda_{ij} / \sum_{k=1,\cdots,7;  l=1,\cdots,7} \lambda_{kl}.
\] for \(i=1,\cdots,7\) and \(j=1,\cdots,7\).

**For the prior distributions**, we specify for \(p(\lambda_{ij})\),
\(i=1,2; j=1,2\) uninformative prior distributions. We use a Gamma
distribution with parameters \(\alpha_i=0.8/7^2\) and \(\beta=0.8/Z\)
with
\(Z= \sum_{ij | n_{ij}>0} n_{ij}/(s_i*s_j) + \sum_{ij | n_{ij}>0} (1-s_i*s_j)/(s_i*s_j)\).
This choice implies for \(\pi\) a Dirichlet prior distribution with
parameters \(\alpha_i\), which is considered to be an objective choice.
For \(p(s_i)\), we use a strongly informative prior distribution, based
on the available data as illustrated above.

An alternative prior distribution is Gaussian process prior, which
penalises large changes between neighbouring age groups, and the kernel
\(k\) is given by

\[
k((a,b),(a',b'))= \sigma^2 \exp(-0.5 ((\frac{a-a'}{\ell_1})^2 + (\frac{b-b'}{\ell_2})^2 ))
\] That is, \(\lambda_{i,j} = \mu + f\) where
\(f\sim \mathcal{GP}(0,k)\).

# MCMC

## MCMC syntax

We use a Markov Chain Monte Carlo algorithm to sample from the posterior
The syntax for running the algorithm is as follows.

``` r
# specify a list of control variables:
#   seed    random number seed
#   mcmc.n  number of MCMC iterations
#   method "gamma" or "gp"
#   outfile output file name if you like to have the results 
#           written to an *.rda* file
control <- list(seed=42, mcmc.n=500, method='gamma')
# run MCMC
ans.gamma <- source.attribution.mcmc.stan(dobs, dprior.fit, control)
```

    ## Error in new_CppObject_xp(fields$.module, fields$.pointer, ...) : 
    ##   Exception: variable does not exist; processing stage=data initialization; variable name=xi_id; base type=int  (in 'model54b1778e3d67_gamma' at line 7)

    ## failed to create the sampler; sampling not done

``` r
save(ans.gamma, file = '~/phyloscanner/phyloflows/data/sevenGroupFlows1_gamma_mcmc.RData')

control <- list(seed=42, mcmc.n=500, method='gp')
# run MCMC
ans.gp <- source.attribution.mcmc.stan(dobs, dprior.fit, control)
```

    ## Error in new_CppObject_xp(fields$.module, fields$.pointer, ...) : 
    ##   Exception: variable does not exist; processing stage=data initialization; variable name=xi_id; base type=int  (in 'model296e35706ff8_gp' at line 42)

    ## failed to create the sampler; sampling not done

``` r
save(ans.gp, file = '~/phyloscanner/phyloflows/data/sevenGroupFlows1_gp_mcmc.RData')
```
