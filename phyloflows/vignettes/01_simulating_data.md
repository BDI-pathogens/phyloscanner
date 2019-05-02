01 - Simulating data
================
2019-04-30

Here is some code to set up simulated data sets. The main aim of these
simulations are to test the performance of phyloflows MCMC algorithm.
Please read the sections “Our Job” and “Our Solution” on the main page
before you go ahead here.

## Data set 1 (simple, SARWS)

We start with simulating transmission flows between two population
groups called “1” and “2”. We will assume a very simple sampling
process, sampling at random within population strata, which we
abbreviate to SARWS.

**We first set up the sampling process within the two population
groups.** Let us assume population group 1 consists of \(2000\)
individuals and group 2 of \(2500\) individuals, and that the sampling
rates are \(0.6\) for group 1 and \(0.45\) for group 2. Here, we will
suppose that sampling is at random within each of the population groups
with these two sampling probabilities:

``` r
library(data.table)
set.seed(42)
ds <- data.table(CATEGORY=c(1,2),TRIAL=c(2000,2500),P_SEQ_EMP=c(0.6,0.45))
ds[,SUC:=TRIAL*P_SEQ_EMP]
ds
```

    ##    CATEGORY TRIAL P_SEQ_EMP  SUC
    ## 1:        1  2000      0.60 1200
    ## 2:        2  2500      0.45 1125

**Next, we calculate the sampling probabilities of transmission flows
within and between the two
groups.**

``` r
dobs <- data.table(TR_TRM_CATEGORY=c(1,1,2,2),REC_TRM_CATEGORY=c(1,2,1,2))
tmp <- subset(ds,select=c('CATEGORY','P_SEQ_EMP'))
setnames(tmp,colnames(tmp),paste0('TR_TRM_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='TR_TRM_CATEGORY')
setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='REC_TRM_CATEGORY')
dobs[, S:= TR_TRM_P_SEQ_EMP * REC_TRM_P_SEQ_EMP]
set(dobs, NULL, c('TR_TRM_P_SEQ_EMP','REC_TRM_P_SEQ_EMP'), NULL)
setkey(dobs,TR_TRM_CATEGORY,REC_TRM_CATEGORY)
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY      S
    ## 1:                1               1 0.3600
    ## 2:                2               1 0.2700
    ## 3:                1               2 0.2700
    ## 4:                2               2 0.2025

**Next, we simulate true transmission flows.** Let us assume 36% and 54%
transmissions are within group 1 and 2 respectively, and 4% are from
group 1 to group 2 and 6% are from group 2 to group 1: \[
\pi=(0.36,0.04,0.06,0.54).
\] We further assume the total number of observed transmissions is
\(N=300\). We will simulate the actual transmission count \(Z\) from a
Poisson distribution. Then we will generate transmission flows between
groups by \[
z \sim \mbox{Multinomial} (Z,\pi).
\] Finally we will generate observed transmissions flows by subsampling
the actual transmission flows by \[
n_{ab} \sim \mbox{Binomial} (z_{ab},\xi_{ab}), \forall a,b,
\] where \(\xi_{ab}\) is the probability of sampling a transmission
event from \(a\) to \(b\).

``` r
TRUE_PI <- c(0.36,0.04,0.06,0.54)
N <- rpois(1,300/mean(dobs$S))
z <- rmultinom(1,size=N,prob=TRUE_PI)
n <- matrix(NA_integer_,ncol=1,nrow=length(TRUE_PI))
for (i in 1:length(TRUE_PI)){
    n[i] <- rbinom(1,size=z[i],dobs$S[i])
}
dobs[, TRM_OBS:= n]
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY      S TRM_OBS
    ## 1:                1               1 0.3600     139
    ## 2:                2               1 0.2700      20
    ## 3:                1               2 0.2700      15
    ## 4:                2               2 0.2025     129

**Next, we will bring the transmission count data into the form needed
for phyloflows MCMC algorithm.** We add an ID to each observation,
called `TRM_CAT_PAIR_ID`. We also define the sampling groups. In this
example, they correspond directly to the transmission groups.

``` r
dobs[, TR_SAMPLING_CATEGORY:=TR_TRM_CATEGORY]
dobs[, REC_SAMPLING_CATEGORY:=REC_TRM_CATEGORY]
setkey(dobs,TR_TRM_CATEGORY,REC_TRM_CATEGORY)
dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
dobs[, S:=NULL]
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY TRM_OBS TR_SAMPLING_CATEGORY
    ## 1:                1               1     139                    1
    ## 2:                2               1      20                    1
    ## 3:                1               2      15                    2
    ## 4:                2               2     129                    2
    ##    REC_SAMPLING_CATEGORY TRM_CAT_PAIR_ID
    ## 1:                     1               1
    ## 2:                     2               2
    ## 3:                     1               3
    ## 4:                     2               4

**We still need to define the prior distribution on the unknown sampling
probabilities, and generate samples from it**. At the very top of this
page, defined the number of infected and sampled individuals in
data.frame `ds`. Let us denote these by \(X_a^i\) and \(X_a^s\) for our
two population groups \(a\). Usually this type of information is
available to us in real-world data analyses, and so we work from these
numbers here also. Under the Binomial sampling model that we assume
throughout, \(X_a^s\sim Binom(X_a^i, \xi_a)\). If we suppose a flat
prior on \(\xi_a\), we obtain the posterior distribution of the sampling
probabilities conditional on the number of total and sampled
individuals, \[
p(\xi_a|X_a^i,X_s^i)= Beta(\xi_a;X_a^s+1,X_a^i-X_a^s+1).
\] This density typically contains a lot of information on the sampling
process. To get this information into the form needed for **phyloflows**
MCMC algorithm, we need to derive samples from that distribution. We
also need to calculate their log density.

``` r
nprior<-1000
dprior<-ds[,list(P=rbeta(nprior,SUC+1,TRIAL-SUC+1),
         SAMPLE=1:nprior),
   by='CATEGORY']
dprior<-merge(dprior,ds,by='CATEGORY')
dprior[,LP:=dbeta(P,SUC+1,TRIAL-SUC+1,log=TRUE)]
set(dprior,NULL,c('TRIAL','SUC'),NULL)
setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
head(dprior)
```

    ##    SAMPLING_CATEGORY         P SAMPLE P_SEQ_EMP       LP
    ## 1:                 1 0.5824160      1       0.6 2.318750
    ## 2:                 1 0.6184042      2       0.6 2.168504
    ## 3:                 1 0.6033518      3       0.6 3.548540
    ## 4:                 1 0.6015475      4       0.6 3.585452
    ## 5:                 1 0.5918721      5       0.6 3.321375
    ## 6:                 1 0.6034198      6       0.6 3.546614

This data set can be loaded through

``` r
data(twoGroupFlows1)
```

To test our MCMC algorithm, we repeated generated 100 such data set, and
they can be loaded through

``` r
data(twoGroupFlows100)
```

## Data set 2 (more complex, GLM)

We then simulate transmission flows between two locations called “1” and
“2”. We will assume sampling depends on genders and is thus not at
random within each location. Due to the dependence of sampling on
genders and locations, sampling rates could be predicted through genders
and locations, which we abbreviate to GLM.

**We first set up the sampling process within the four population
groups.** We consider four population groups, “1M”, “1F”, “2M” and “2F”.
Let us assume there are \(10,000\) individuals in total, \(40\%\)
individuals are in location 1, \(60\%\) individuals are in location 2,
\(60\%\) individuals are women within each location, and that the
sampling rates are \(0.6\) for men in group 1, \(0.8\) for women in
group 1, \(0.3\) for men in group 2 and \(0.7\) for women in group 2.
Here, we will suppose that sampling is at random within each of the
population groups with these four sampling probabilities:

``` r
library(data.table)
set.seed(42)
ds<-data.table(SEX=c(rep(c('M','F'),2)),
               COMM_NUM_B=c(rep(1,2),rep(2,2)),
               P_SEQ=c(0.6,0.8,0.3,0.7)) # set up seq rates

ninf<-1000
ninf1<-round(ninf*0.4)
ninf2<-ninf-ninf1
ninf1f<-round(ninf1*0.7)
ninf1m<-ninf1-ninf1f
ninf2f<-round(ninf2*0.6)
ninf2m<-ninf2-ninf2f
ds[,TRIAL:=c(ninf1m,ninf1f,ninf2m,ninf2f)]
ds[,SUC:=rbinom(nrow(ds),TRIAL,P_SEQ)]

ds[,CATEGORY:=paste0(COMM_NUM_B,':',SEX)]
ds[,COMM_NUM_B:=as.integer(COMM_NUM_B)]
ds[,MALE:=as.integer(SEX=='M')]
```

**Next, we calculate the sampling probabilities of transmission flows
within and between the two groups.**

``` r
dobs<-data.table(TR_COMM_NUM_B=c(rep(1,4),rep(2,4)),
                 REC_COMM_NUM_B=c(rep(1,2),rep(2,2),rep(1,2),rep(2,2)),
                 TR_SEX=c(rep(c('F','M'),4)),
                 REC_SEX=c(rep(c('M','F'),4)))
dobs[,TR_TRM_CATEGORY:=paste0(TR_COMM_NUM_B,":",TR_SEX)]
dobs[,REC_TRM_CATEGORY:=paste0(REC_COMM_NUM_B,":",REC_SEX)]
dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]

tmp <- subset(ds,select=c('CATEGORY','P_SEQ'))
setnames(tmp,colnames(tmp),paste0('TR_TRM_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='TR_TRM_CATEGORY')
setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='REC_TRM_CATEGORY')
dobs[, S:= TR_TRM_P_SEQ * REC_TRM_P_SEQ]
set(dobs, NULL, c('TR_TRM_P_SEQ','REC_TRM_P_SEQ'), NULL)
setkey(dobs,TRM_CAT_PAIR_ID)
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY TR_COMM_NUM_B REC_COMM_NUM_B TR_SEX
    ## 1:              1:M             1:F             1              1      F
    ## 2:              1:F             1:M             1              1      M
    ## 3:              2:M             1:F             1              2      F
    ## 4:              2:F             1:M             1              2      M
    ## 5:              1:M             2:F             2              1      F
    ## 6:              1:F             2:M             2              1      M
    ## 7:              2:M             2:F             2              2      F
    ## 8:              2:F             2:M             2              2      M
    ##    REC_SEX TRM_CAT_PAIR_ID    S
    ## 1:       M               1 0.48
    ## 2:       F               2 0.48
    ## 3:       M               3 0.24
    ## 4:       F               4 0.42
    ## 5:       M               5 0.42
    ## 6:       F               6 0.24
    ## 7:       M               7 0.21
    ## 8:       F               8 0.21

**Next, we simulate true transmission flows.** Let us assume 30% and 45%
transmissions are within location 1 and 2 respectively, and 10% are from
location 1 to location 2 and 15% are from location 2 to location 1.
\(60\%\) of each transmission flow between locations originate from men:

\[
\pi=(0.12,0.18,0.04,0.06,0.06,0.09,0.18,0.27).
\] We further assume the total number of observed transmissions is
\(N=300\). We will simulate the actual transmission count \(Z\) from a
Poisson distribution. Then we will generate transmission flows between
groups by \[
z \sim \mbox{Multinomial} (Z,\pi).
\] Finally we will generate observed transmissions flows by subsampling
the actual transmission flows by \[
n_{ab} \sim \mbox{Binomial} (z_{ab},\xi_{ab}), \forall a,b,
\] where \(\xi_{ab}\) is the probability of sampling a transmission
event from \(a\) to \(b\).

``` r
TRUE_PI <- c(0.12,0.18,0.04,0.06,0.06,0.09,0.18,0.27)
N <- rpois(1,300/mean(dobs$S))
z <- rmultinom(1,size=N,prob=TRUE_PI)
n <- matrix(NA_integer_,ncol=1,nrow=length(TRUE_PI))
for (i in 1:length(TRUE_PI)){
    n[i] <- rbinom(1,size=z[i],dobs$S[i])
}
setkey(dobs,TR_COMM_NUM_B,REC_COMM_NUM_B,TR_SEX,REC_SEX)
dobs[, TRM_OBS:= n]
dobs.zero<-which(dobs$TRM_OBS!=0)
dobs<-dobs[TRM_OBS!=0,]
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY TR_COMM_NUM_B REC_COMM_NUM_B TR_SEX
    ## 1:              1:M             1:F             1              1      F
    ## 2:              1:F             1:M             1              1      M
    ## 3:              2:M             1:F             1              2      F
    ## 4:              2:F             1:M             1              2      M
    ## 5:              1:M             2:F             2              1      F
    ## 6:              1:F             2:M             2              1      M
    ## 7:              2:M             2:F             2              2      F
    ## 8:              2:F             2:M             2              2      M
    ##    REC_SEX TRM_CAT_PAIR_ID    S TRM_OBS
    ## 1:       M               1 0.48      55
    ## 2:       F               2 0.48      77
    ## 3:       M               3 0.24       7
    ## 4:       F               4 0.42      22
    ## 5:       M               5 0.42      17
    ## 6:       F               6 0.24      28
    ## 7:       M               7 0.21      35
    ## 8:       F               8 0.21      56

**Next, we will bring the transmission count data into the form needed
for phyloflows MCMC algorithm.** We add an ID to each observation,
called `TRM_CAT_PAIR_ID`. We also define the sampling groups. In this
example, they correspond directly to the transmission groups.

``` r
dobs[, TR_SAMPLING_CATEGORY:=TR_TRM_CATEGORY]
dobs[, REC_SAMPLING_CATEGORY:=REC_TRM_CATEGORY]
setkey(dobs,TR_TRM_CATEGORY,REC_TRM_CATEGORY)
dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
dobs[, S:=NULL]
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY TR_COMM_NUM_B REC_COMM_NUM_B TR_SEX
    ## 1:              1:M             1:F             1              1      F
    ## 2:              2:M             1:F             1              2      F
    ## 3:              1:F             1:M             1              1      M
    ## 4:              2:F             1:M             1              2      M
    ## 5:              1:M             2:F             2              1      F
    ## 6:              2:M             2:F             2              2      F
    ## 7:              1:F             2:M             2              1      M
    ## 8:              2:F             2:M             2              2      M
    ##    REC_SEX TRM_CAT_PAIR_ID TRM_OBS TR_SAMPLING_CATEGORY
    ## 1:       M               1      55                  1:F
    ## 2:       M               2       7                  1:F
    ## 3:       F               3      77                  1:M
    ## 4:       F               4      22                  1:M
    ## 5:       M               5      17                  2:F
    ## 6:       M               6      35                  2:F
    ## 7:       F               7      28                  2:M
    ## 8:       F               8      56                  2:M
    ##    REC_SAMPLING_CATEGORY
    ## 1:                   1:M
    ## 2:                   2:M
    ## 3:                   1:F
    ## 4:                   2:F
    ## 5:                   1:M
    ## 6:                   2:M
    ## 7:                   1:F
    ## 8:                   2:F

**We still need to define the prior distribution on the unknown sampling
probabilities, and generate samples from it**. Like Data Set 1, we
defined the number of infected and sampled individuals in data.frame
`ds`. Let us denote these by \(X_a^i\) and \(X_a^s\) for our four
population groups \(a\). We have an additional individual covariate,
gender denoted as \(x_a\), for population group \(a\) and we will work
from these data. Under the Binomial sampling model that we assume
throughout, \(X_a^s\sim Binom(X_a^i, \xi_a)\). Because of dependence of
sampling rates on covariates and locations \(b\), \(\xi_a\) could be
linked with a regression with gender,
\(\mbox{logit}(\xi_{ab})=\beta_0+\gamma_b+\beta_1 x_{a}\). If we suppose
broad priors on \(\beta_0\) and \(\beta_1\), we obtain the posterior
distribution of the sampling probabilities conditional on gender, the
number of total and sampled individuals from rstan.

This density typically contains a lot of information on the sampling
process. To get this information into the form needed for **phyloflows**
MCMC algorithm, we need to derive samples from that distribution. We
also need to calculate their log density.

``` r
library(rstan)
```

    ## Loading required package: ggplot2

    ## Need help getting started? Try the cookbook for R:
    ## http://www.cookbook-r.com/Graphs/

    ## Loading required package: StanHeaders

    ## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## options(mc.cores = parallel::detectCores()).
    ## To avoid recompilation of unchanged Stan programs, we recommend calling
    ## rstan_options(auto_write = TRUE)

``` r
library(bayesplot)
```

    ## This is bayesplot version 1.6.0

    ## - Online documentation and vignettes at mc-stan.org/bayesplot

    ## - bayesplot theme set to bayesplot::theme_default()

    ##    * Does _not_ affect other ggplot2 plots

    ##    * See ?bayesplot_theme_set for details on theme setting

``` r
infile.sequencing.stan.model<-'glm_sex2comm2.stan'
#   run STAN 
tmp         <- as.list(subset(ds,select=c('COMM_NUM_B','TRIAL','SUC','MALE')))
tmp$N       <- nrow(ds)
tmp$N_COMM  <- length(unique(ds$COMM_NUM_B))
fit.par     <- stan(    file = infile.sequencing.stan.model, 
                  data = tmp, 
                  iter = 10e3,
                  warmup = 5e2,
                  cores = 1,
                  chains = 1,
                  init = list(list(a=0, comm=rep(0,2), sig_comm=1, male=0)))
```

    ## 
    ## SAMPLING FOR MODEL 'glm_sex2comm2' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 1.7e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 10000 [  0%]  (Warmup)
    ## Chain 1: Iteration:  501 / 10000 [  5%]  (Sampling)
    ## Chain 1: Iteration: 1500 / 10000 [ 15%]  (Sampling)
    ## Chain 1: Iteration: 2500 / 10000 [ 25%]  (Sampling)
    ## Chain 1: Iteration: 3500 / 10000 [ 35%]  (Sampling)
    ## Chain 1: Iteration: 4500 / 10000 [ 45%]  (Sampling)
    ## Chain 1: Iteration: 5500 / 10000 [ 55%]  (Sampling)
    ## Chain 1: Iteration: 6500 / 10000 [ 65%]  (Sampling)
    ## Chain 1: Iteration: 7500 / 10000 [ 75%]  (Sampling)
    ## Chain 1: Iteration: 8500 / 10000 [ 85%]  (Sampling)
    ## Chain 1: Iteration: 9500 / 10000 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 10000 / 10000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 0.081532 seconds (Warm-up)
    ## Chain 1:                1.88541 seconds (Sampling)
    ## Chain 1:                1.96694 seconds (Total)
    ## Chain 1:

``` r
# check
fit.pars    <- c('a','comm','sig_comm','male')
any(rhat(fit.par, pars=fit.pars)>1.02)
```

    ## [1] FALSE

``` r
any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
```

    ## [1] FALSE

``` r
# extract samples from the prior distribution
fit.e       <- extract(fit.par)
tmp         <- sample(length(fit.e$a), nprior)
dprior          <- ds[, 
                   {
                     z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE)
                     list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
                   },   
                   by=c('CATEGORY')]
dprior[, P:= exp(ETA)/(1+exp(ETA))]

# log density
require(bde)
```

    ## Loading required package: bde

    ## Loading required package: shiny

    ## 
    ## Attaching package: 'bde'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     density, quantile

    ## The following object is masked from 'package:methods':
    ## 
    ##     getSubclasses

``` r
tmp <- dprior[, {
  bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1,
              options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
  list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
}, by='CATEGORY']
dprior  <- merge(dprior, tmp, by=c('CATEGORY','SAMPLE'))
set(dprior, NULL, c('ETA'), NULL)   
setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')

head(dprior)
```

    ##    SAMPLING_CATEGORY SAMPLE         P       LP
    ## 1:               1:F      1 0.8279894 2.780059
    ## 2:               1:F      2 0.8316006 2.703062
    ## 3:               1:F      3 0.8056642 2.688324
    ## 4:               1:F      4 0.7994241 2.499298
    ## 5:               1:F      5 0.8252191 2.820963
    ## 6:               1:F      6 0.8211598 2.853287

This data set can be loaded through

``` r
data(fourGroupFlows1)
```
