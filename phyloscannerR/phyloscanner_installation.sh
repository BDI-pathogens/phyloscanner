#!/bin/bash


# Load anaconda3/personal module and install if not installed already
module load anaconda3/personal
if [ -d $HOME/anaconda3 ]; then
    echo -e "\nanaconda3/personal already installed\n\n";
else
    echo -e "\ninstalling anaconda3/personal\n\n";
    anaconda-setup
fi


# Create new conda environment called "phylo"
if [ -d $HOME/anaconda3/envs/phylo ]; then
    echo "###############################################"
    echo -e "\nphylo conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n phylo --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment for R packages: phylo"
    conda create -n phylo -y
    source activate phylo
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi

# Install initial dependencies
echo -e "\nInstalling Dependencies: R packages via conda"
conda install cython matplotlib pysam biopython mafft raxml curl -y

conda install r r-base r-devtools r-ggplot2 r-network r-scales r-sna r-igraph r-intergraph r-knitr r-rcolorbrewer r-testthat r-argparse r-ape r-extradistr r-ff r-ggally r-gtable r-kimisc r-pegas r-phangorn r-phytools r-prodlim r-reshape2 r-tidyverse bioconductor-treeio r-viridis r-digest r-gtable r-lazyeval r-mass r-mgcv r-matrix r-lattice r-nlme r-rlang r-r6 r-rcpp r-cli r-tibble r-assertthat r-fansi r-pillar r-vctrs r-backports r-ellipsis r-glue r-zeallot r-vctrs r-findpython r-jsonlite r-colorspace r-data.table r-callr r-processx r-git2r r-httr r-mime r-pkgbuild r-pkgload r-rstudioapi r-rcmdcheck r-roxygen2 r-brew r-commonmark r-purrr r-stringi r-stringr r-xml2 r-evaluate r-dplyr r-bh r-tidyselect r-dtplyr r-highr r-markdown r-xfun r-quadprog r-numderiv r-reshape r-rmarkdown r-tinytex r-units r-sf -y

conda install bioconductor-rsamtools bioconductor-rbgl bioconductor-ggtree bioconductor-genomeinfodb bioconductor-genomicranges bioconductor-genomeinfodbdata -y

# Install additional dependencies not available via anaconda
echo -e "\nInstalling Additional Dependencies not available via anaconda"
R -e 'options(unzip = "internal");library(devtools);devtools::install_github("briatte/ggnet")'

# Installing phyloscanner
echo -e "\nInstalling phyloscanner"
cd $HOME
git clone https://github.com/BDI-pathogens/phyloscanner.git
cd phyloscanner/phyloscannerR
R CMD INSTALL .

# Testing
echo -e "\nTesting phyloscanner library loads correctly"
R -e 'library(phyloscannerR)'

# Installing phyloscanner.R.utilities
echo -e "\nInstalling Phyloscanner.R.utilities."
R -e 'options(unzip = "internal");library(devtools);devtools::install_github("olli0601/Phyloscanner.R.utilities")'


# Testing
echo -e "\nTesting phyloscanner.R.utilities"
R -e 'library(Phyloscanner.R.utilities)'

echo "========================================="
echo "INSTALLATION COMPLETE"


echo -e "#############################################################"
echo -e "#############################################################\n\n"

# INSTRUCTIONS OF USE:
echo -e "\t\tINSTRUCTIONS FOR USE"
echo -e "\n#############################################################\n\n"
echo -e "Ensure you add the following lines in your jobscript:\n\n"
echo -e "\tmodule load anaconda3/personal"
echo -e "\tsource activate phylo"
echo -e "\t"'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH'
echo -e "\n\n"
echo -e "\n#############################################################\n\n"


