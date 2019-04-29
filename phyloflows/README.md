# phyloflows

**phyloflows** provides an MCMC algorithm to infer transmission flows within and between population groups from observed flow counts. 

**phyloflows** was primarily written to estimate transmission flows from counts of phylogenetically reconstructed source-recipient pairs. These input data can be obtained with [phyloscanner](https://github.com/BDI-pathogens/phyloscanner), though different types of input data can also be used, so long as they represent observed counts of flow events between two discrete groups. The primary feature of **phyloflows** is that the underlying statistical model allows adjusting for sampling heterogeneity across population groups while estimating transmission flows, a common problem in molecular epidemiologic analyses.


### Getting started

If you are just getting started with **phyloflows**, why don t you take a look at the tutorial vignettes.

### Installation

Install the latest development version from **GitHub**:

```{r}
install_github("BDI-pathogens/phyloscanner/phyloflows", dependencies=TRUE, build_vignettes=FALSE)
require(phyloflows)
```

If you have issues with installation/running of phyloflows, please raise an Issue on the GitHub page and we will get back to you.

## Quick examples

**phyloflows** expects input data in a very simple format. First, a data.frame of observed transmission counts within and between population groups, which we call `dobs`. Second a data.frame that summarises prior information on how population groups were sampled, which we call `dprior`. 

To get you started, **phyloflows** comes with a small simulated example data set of transmission counts and sampling information between two population groups, denoted by "1" and "2":
```
data(twoGroupFlows1)
dobs <- twoGroupFlows1$dobs
dprior <- twoGroupFlows1$dprior
```