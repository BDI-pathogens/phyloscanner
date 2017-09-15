# Phyloscanner's analysis of phylogenies with within- and between-host diversity, achived with the
# phyloscanner_analyse_trees.R command, requires first of all an up-to-date installation of the R
# language.
# BEWARE: unfortunately, updating the R language itself can cause problems with
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
# Once you've managed to add the key, rerun the sudo apt-get update command. Then,
sudo apt-get install r-base
sudo apt-get install r-base-dev
############################################################

############################################################
# INSTALLING UP-TO-DATE R ON MAC OS
???? TODO
############################################################

# Next, install the packages needed for phyloscanner_analyse_trees.R, by running the installation
# script.
# If you don't have sudo priviledges on your system (e.g. you're using a computing cluster), 
# instead of executing the command below, try entering an interactive R session by typing 'R' at
# the command line, then running the commands contained inside the package_install.R file one
# by one. (R should then prompt you to install things in your home directory, where you have
# permission.)
sudo Rscript ~/phyloscanner/tools/package_install.R
# If ggtree fails to install, this could be because of an out of date BiocInstaller package.
# Removing it worked for me: entering an interactive R session and running
remove.packages("BiocInstaller")
# exiting the R session (Ctrl+D), then rerunning
sudo Rscript ~/phyloscanner/tools/package_install.R

# Finally, ensuring all installed R packages are up-to-date can help: enter an interactive R
# session and run
library(rvcheck)
update_all()
