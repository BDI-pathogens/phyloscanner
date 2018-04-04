# Phyloscanner's analysis of phylogenies with within- and between-host diversity is achived with the
# phyloscanner_analyse_trees.R command. Having an up-to-date installation of the R language can
# avoid some runtime errors; however BEWARE: updating the R language itself can cause problems with
# packages you have in your current R installation - they may need to be reinstalled. 

############################################################
# INSTALLING UP-TO-DATE R ON UBUNTU
# The comands below *should* do it.

# Add an R repository to the list of things Ubuntu should keep up to date.
# The command below uses the R studio mirror (you can use others), and specifies 'trusty' as the
# codename for the Ubuntu version you're running (you should replace this appropriately if running
# a different version of Ubuntu). See https://cran.r-project.org/bin/linux/ubuntu/README.html for
# details.
sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'

# Gather information about what needs updating:
sudo apt-get update
# That command may fail with an error message saying "GPG error",
# "the public key is not available", "The repository... is not signed".
# In that case, try adding the key thus:
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
# That may fail, for example due to firewall issues. In that case, consult the section
# 'Secure APT' at https://cran.r-project.org/bin/linux/ubuntu/README.html.
# Once you've managed to add the key, rerun the 'sudo apt-get update' command. Then,
sudo apt-get install r-base
sudo apt-get install r-base-dev
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
