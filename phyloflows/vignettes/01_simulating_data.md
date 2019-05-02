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
groups.** Let us assume population group 1 consists of
![2000](https://latex.codecogs.com/png.latex?2000 "2000") individuals
and group 2 of ![2500](https://latex.codecogs.com/png.latex?2500 "2500")
individuals, and that the sampling rates are
![0.6](https://latex.codecogs.com/png.latex?0.6 "0.6") for group 1 and
![0.45](https://latex.codecogs.com/png.latex?0.45 "0.45") for group 2.
Here, we will suppose that sampling is at random within each of the
population groups with these two sampling probabilities:

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
within and between the two groups.**

``` r
dobs <- data.table(TR_TRM_CATEGORY=c(1,1,2,2),REC_TRM_CATEGORY=c(1,2,1,2))
tmp <- subset(ds,select=c('CATEGORY','P_SEQ_EMP'))
setnames(tmp,colnames(tmp),paste0('TR_TRM_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='TR_TRM_CATEGORY')
setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
dobs <- merge(dobs,tmp,by='REC_TRM_CATEGORY')
dobs[, S:= TR_TRM_P_SEQ_EMP * REC_TRM_P_SEQ_EMP]
set(dobs, NULL, c('TR_TRM_P_SEQ_EMP','REC_TRM_P_SEQ_EMP'), NULL)
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY      S
    ## 1:                1               1 0.3600
    ## 2:                1               2 0.2700
    ## 3:                2               1 0.2700
    ## 4:                2               2 0.2025

**Next, we simulate true transmission flows.** Let us assume 36% and 54%
transmissions are within group 1 and 2 respectively, and 4% are from
group 1 to group 2 and 6% are from group 2 to group 1:   
![
\\pi=(0.36,0.04,0.06,0.54).
](https://latex.codecogs.com/png.latex?%0A%5Cpi%3D%280.36%2C0.04%2C0.06%2C0.54%29.%0A
"
\\pi=(0.36,0.04,0.06,0.54).
")  
We further assume the total number of observed transmissions is
![N=300](https://latex.codecogs.com/png.latex?N%3D300 "N=300"). We will
simulate the actual transmission count
![Z](https://latex.codecogs.com/png.latex?Z "Z") from a Poisson
distribution. Then we will generate transmission flows between groups by
  
![
z \\sim \\mbox{Multinomial} (Z,\\pi).
](https://latex.codecogs.com/png.latex?%0Az%20%5Csim%20%5Cmbox%7BMultinomial%7D%20%28Z%2C%5Cpi%29.%0A
"
z \\sim \\mbox{Multinomial} (Z,\\pi).
")  
Finally we will generate observed transmissions flows by subsampling the
actual transmission flows by   
![
n\_{ab} \\sim \\mbox{Binomial} (z\_{ab},\\xi\_{ab}), \\forall a,b,
](https://latex.codecogs.com/png.latex?%0An_%7Bab%7D%20%5Csim%20%5Cmbox%7BBinomial%7D%20%28z_%7Bab%7D%2C%5Cxi_%7Bab%7D%29%2C%20%5Cforall%20a%2Cb%2C%0A
"
n_{ab} \\sim \\mbox{Binomial} (z_{ab},\\xi_{ab}), \\forall a,b,
")  
where ![\\xi\_{ab}](https://latex.codecogs.com/png.latex?%5Cxi_%7Bab%7D
"\\xi_{ab}") is the probability of sampling a transmission event from
![a](https://latex.codecogs.com/png.latex?a "a") to
![b](https://latex.codecogs.com/png.latex?b "b").

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
    ## 2:                1               2 0.2700      20
    ## 3:                2               1 0.2700      15
    ## 4:                2               2 0.2025     129

**Next, we will bring the transmission count data into the form needed
for phyloflows MCMC algorithm.** We add an ID to each observation,
called `TRM_CAT_PAIR_ID`. We also define the sampling groups. In this
example, they correspond directly to the transmission groups.

``` r
dobs[, TR_SAMPLING_CATEGORY:=TR_TRM_CATEGORY]
dobs[, REC_SAMPLING_CATEGORY:=REC_TRM_CATEGORY]
dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
dobs[, S:=NULL]
setkey(dobs,TR_TRM_CATEGORY,REC_TRM_CATEGORY)
dobs
```

    ##    REC_TRM_CATEGORY TR_TRM_CATEGORY TRM_OBS TR_SAMPLING_CATEGORY
    ## 1:                1               1     139                    1
    ## 2:                2               1      15                    1
    ## 3:                1               2      20                    2
    ## 4:                2               2     129                    2
    ##    REC_SAMPLING_CATEGORY TRM_CAT_PAIR_ID
    ## 1:                     1               1
    ## 2:                     2               3
    ## 3:                     1               2
    ## 4:                     2               4

**We still need to define the prior distribution on the unknown sampling
probabilities, and generate samples from it**. At the very top of this
page, defined the number of infected and sampled individuals in
data.frame `ds`. Let us denote these by
![X\_a^i](https://latex.codecogs.com/png.latex?X_a%5Ei "X_a^i") and
![X\_a^s](https://latex.codecogs.com/png.latex?X_a%5Es "X_a^s") for our
two population groups ![a](https://latex.codecogs.com/png.latex?a "a").
Usually this type of information is available to us in real-world data
analyses, and so we work from these numbers here also. Under the
Binomial sampling model that we assume throughout, ![X\_a^s\\sim
Binom(X\_a^i,
\\xi\_a)](https://latex.codecogs.com/png.latex?X_a%5Es%5Csim%20Binom%28X_a%5Ei%2C%20%5Cxi_a%29
"X_a^s\\sim Binom(X_a^i, \\xi_a)"). If we suppose a flat prior on
![\\xi\_a](https://latex.codecogs.com/png.latex?%5Cxi_a "\\xi_a"), we
obtain the posterior distribution of the sampling probabilities
conditional on the number of total and sampled individuals,   
![
p(\\xi\_a|X\_a^i,X\_s^i)= Beta(\\xi\_a;X\_a^s+1,X\_a^i-X\_a^s+1).
](https://latex.codecogs.com/png.latex?%0Ap%28%5Cxi_a%7CX_a%5Ei%2CX_s%5Ei%29%3D%20Beta%28%5Cxi_a%3BX_a%5Es%2B1%2CX_a%5Ei-X_a%5Es%2B1%29.%0A
"
p(\\xi_a|X_a^i,X_s^i)= Beta(\\xi_a;X_a^s+1,X_a^i-X_a^s+1).
")  
This density typically contains a lot of information on the sampling
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

xxx
