# Phyloscanner's analysis of phylogenies with within- and between-host diversity is achived with the
# phyloscanner_analyse_trees.R command. Having an up-to-date installation of the R language can
# avoid some runtime errors; however BEWARE: updating the R language itself can cause problems with
# packages you have in your current R installation - they may need to be reinstalled. 

############################################################
# INSTALLING UP-TO-DATE R ON UBUNTU
# Following the steps at https://www.r-bloggers.com/updating-r-on-ubuntu/
# worked for me.
# However, the command 
# sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
# may fail, for example due to firewall issues. In that case, consult the section
# 'Secure APT' at https://cran.r-project.org/bin/linux/ubuntu/README.html.
############################################################

############################################################
# INSTALLING UP-TO-DATE R ON MAC OS
???? TODO
############################################################

# Next, install the packages needed for phyloscanner_analyse_trees.R
# Change directory to the 'phyloscannerR' subdirectory of the main phyloscanner code directory,
# e.g. cd ~/phyloscanner/phyloscannerR/, then run
sudo R
# Inside the interactive R console, run
install.packages("devtools")
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
install(".", dependencies = T)
