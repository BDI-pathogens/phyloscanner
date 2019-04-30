**phyloflows**
================
Xiaoyue Xi, Oliver Ratmann
2019-04-30

**phyloflows** provides an MCMC algorithm to infer transmission flows
within and between population groups from observed flow counts.

**phyloflows** was primarily written to estimate transmission flows from
counts of phylogenetically reconstructed likely transmission pairs, in
whom the direction of transmission could be inferred with high support.
These input data can be obtained with
[phyloscanner](https://github.com/BDI-pathogens/phyloscanner), though
different types of input data can also be used, so long as they
represent observed counts of flow events between two discrete groups.
The primary feature of **phyloflows** is that the underlying statistical
model allows adjusting for sampling heterogeneity across population
groups while estimating transmission flows, a common problem in
molecular epidemiologic analyses.

Installation
------------

Install the latest development version from **GitHub**:

``` r
install_github("BDI-pathogens/phyloscanner/phyloflows", dependencies=TRUE, build_vignettes=FALSE)
require(phyloflows)
```

If you have issues with installation/running of phyloflows, please raise
an Issue on the GitHub page and we will get back to you.

Getting started
---------------

If you are just getting started with **phyloflows**, please take a look
at the sections below, and the tutorial vignettes at the bottom of this
page.

Our job
-------

Suppose the true number of transmission flows between two population
groups “1” and “2” is

![
(z\_{11}=50, z\_{12}=10, z\_{21}=20, z\_{22}=20),
](https://latex.codecogs.com/png.latex?%0A%28z_%7B11%7D%3D50%2C%20z_%7B12%7D%3D10%2C%20z_%7B21%7D%3D20%2C%20z_%7B22%7D%3D20%29%2C%0A "
(z_{11}=50, z_{12}=10, z_{21}=20, z_{22}=20),
")

so that the true proportion of transmission flows within and between
groups is

![
(\\pi\_{11}=0.5, \\pi\_{12}=0.1, \\pi\_{21}=0.2, \\pi\_{22}=0.2).
](https://latex.codecogs.com/png.latex?%0A%28%5Cpi_%7B11%7D%3D0.5%2C%20%5Cpi_%7B12%7D%3D0.1%2C%20%5Cpi_%7B21%7D%3D0.2%2C%20%5Cpi_%7B22%7D%3D0.2%29.%0A "
(\pi_{11}=0.5, \pi_{12}=0.1, \pi_{21}=0.2, \pi_{22}=0.2).
")

**The primary aim** is to estimate the vector of proportions
![\\pi=(\\pi\_{11},\\pi\_{12},\\pi\_{21},\\pi\_{22})](https://latex.codecogs.com/png.latex?%5Cpi%3D%28%5Cpi_%7B11%7D%2C%5Cpi_%7B12%7D%2C%5Cpi_%7B21%7D%2C%5Cpi_%7B22%7D%29 "\pi=(\pi_{11},\pi_{12},\pi_{21},\pi_{22})").
This looks simple. **However in real life one main complication is**
that only a proportion of the true transmission flows are observed (or
sampled), and that the sampling probabilities differ. In our experience
this situation not only occurs frequently in observational epidemiologic
studies but also controlled situations (trials). Suppose that population
“1” is sampled with probability
![s\_1=0.6](https://latex.codecogs.com/png.latex?s_1%3D0.6 "s_1=0.6")
and population “2” with probability
![s\_2=1](https://latex.codecogs.com/png.latex?s_2%3D1 "s_2=1"). On
expectation, we then observe the transmission flows
![n\_{ij}= z\_{ij} \* s\_i \* s\_j](https://latex.codecogs.com/png.latex?n_%7Bij%7D%3D%20z_%7Bij%7D%20%2A%20s_i%20%2A%20s_j "n_{ij}= z_{ij} * s_i * s_j"),
if we assume that sampling within and between the two groups is
independent of each other, and if we ignore the fact that in reality
sampling is without replacement. It won t matter much in reasonably
large populations anyway, considering the extent of typical sampling
differences. So what we observe is on expectation

![
(n\_{11}=18, n\_{12}=6, n\_{21}=12, n\_{22}=20),
](https://latex.codecogs.com/png.latex?%0A%28n_%7B11%7D%3D18%2C%20n_%7B12%7D%3D6%2C%20n_%7B21%7D%3D12%2C%20n_%7B22%7D%3D20%29%2C%0A "
(n_{11}=18, n_{12}=6, n_{21}=12, n_{22}=20),
")

If we don’t adjust for the sampling, a crude estimate of the true
proportion of transmission flows would be

![
(\\hat{\\pi}\_{11}= 18/56= 0.321, \\hat{\\pi}\_{12}=0.107, \\hat{\\pi}\_{21}=0.214, \\hat{\\pi}\_{22}=0.357)
](https://latex.codecogs.com/png.latex?%0A%28%5Chat%7B%5Cpi%7D_%7B11%7D%3D%2018%2F56%3D%200.321%2C%20%5Chat%7B%5Cpi%7D_%7B12%7D%3D0.107%2C%20%5Chat%7B%5Cpi%7D_%7B21%7D%3D0.214%2C%20%5Chat%7B%5Cpi%7D_%7B22%7D%3D0.357%29%0A "
(\hat{\pi}_{11}= 18/56= 0.321, \hat{\pi}_{12}=0.107, \hat{\pi}_{21}=0.214, \hat{\pi}_{22}=0.357)
")

and the worst case error is
![50\\%-32.1\\%=17.9\\%](https://latex.codecogs.com/png.latex?50%5C%25-32.1%5C%25%3D17.9%5C%25 "50\%-32.1\%=17.9\%").
That is pretty bad. **The job of phyloflow is to return better
estimates**, with a worst case error lower than 2-3%.

Our solution
------------

One very quick solution to the sampling problem would be to infer the
actual transmission flows from what we know, i.e.

![
\\hat{z}\_{ij}= n\_{ij}/(s\_i\*s\_j),
](https://latex.codecogs.com/png.latex?%0A%5Chat%7Bz%7D_%7Bij%7D%3D%20n_%7Bij%7D%2F%28s_i%2As_j%29%2C%0A "
\hat{z}_{ij}= n_{ij}/(s_i*s_j),
")

and then to calculate

![
\\hat{\\pi}\_{ij}= \\hat{z}\_{ij} / \\sum\_{kl} \\hat{z}\_{kl}.
](https://latex.codecogs.com/png.latex?%0A%5Chat%7B%5Cpi%7D_%7Bij%7D%3D%20%5Chat%7Bz%7D_%7Bij%7D%20%2F%20%5Csum_%7Bkl%7D%20%5Chat%7Bz%7D_%7Bkl%7D.%0A "
\hat{\pi}_{ij}= \hat{z}_{ij} / \sum_{kl} \hat{z}_{kl}.
")

But this does not address the further problems that the sampling groups
may be different to the population groups between whom we want to infer
transmission flows; that the sampling probabilities are usually not
known themselves; and that we also want interpretable uncertainty
estimates for
![\\hat{\\pi}](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cpi%7D "\hat{\pi}").
So we use a Bayesian approach, and write down the likelihood of the
complete data,

![
p(z\|Z,\\pi)= Multinomial(z;Z,\\pi),
](https://latex.codecogs.com/png.latex?%0Ap%28z%7CZ%2C%5Cpi%29%3D%20Multinomial%28z%3BZ%2C%5Cpi%29%2C%0A "
p(z|Z,\pi)= Multinomial(z;Z,\pi),
")

where ![Z](https://latex.codecogs.com/png.latex?Z "Z") is the total
number of transmissions,
![Z=\\sum\_{kl} z\_{kl}](https://latex.codecogs.com/png.latex?Z%3D%5Csum_%7Bkl%7D%20z_%7Bkl%7D "Z=\sum_{kl} z_{kl}").
Next, we add the likelihood of the sampling process,

![
p(n\_{ij}\|z\_{ij},s\_i,s\_j)= Binomial(n\_{ij};z\_{ij},s\_i\*s\_j).
](https://latex.codecogs.com/png.latex?%0Ap%28n_%7Bij%7D%7Cz_%7Bij%7D%2Cs_i%2Cs_j%29%3D%20Binomial%28n_%7Bij%7D%3Bz_%7Bij%7D%2Cs_i%2As_j%29.%0A "
p(n_{ij}|z_{ij},s_i,s_j)= Binomial(n_{ij};z_{ij},s_i*s_j).
")

Finally, we complete the model with prior distributions on the
proportions,
![p(\\pi)](https://latex.codecogs.com/png.latex?p%28%5Cpi%29 "p(\pi)"),
which we generally choose in an uninformative way; prior distributions
on the sampling probabilities,
![p(s\_i)](https://latex.codecogs.com/png.latex?p%28s_i%29 "p(s_i)"),
which we generally choose based on other availabe information; and a
prior distribution on the total number of transmissions,
![p(Z)](https://latex.codecogs.com/png.latex?p%28Z%29 "p(Z)"), which we
also choose based on available information. And then we use Bayes
Theorem.

This looks pretty straightforward. The issue is that in real-world data
analyses, the total number of unknown parameters is 1,000 or even more.
And so we also need a nifty computational inference engine.
**phyloflow** uses a special type of Markov Chain Markov Carlo
algorithms for this job. The main feature is that it has zero tuning
variables, so the blame is 100% on us if it does not work. Just get in
touch then.

Examples
--------

1.  [Simulating data.](articles/01_simulating_data.html)

2.  [A very first example. Always start with
    me.](articles/02_basic_example.html)

3.  [Performing MCMC diagnostic checks](articles/03_diagnostics.html)

4.  [Aggregating MCMC output to a new set of
    variables](articles/04_aggregating.html)
