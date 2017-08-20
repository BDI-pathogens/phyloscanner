# Phyloscanner's analysis of phylogenies with within- and between-host diversity, achived with the
# phyloscanner_analyse_trees.R command, requires first of all an up-to-date installation of the R
# language.

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
# That command may fail, and complain about a key. In that case, try adding the key thus:
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

# Next, install the packages needed by phyloscanner_analyse_trees.R.
sudo Rscript ~/phyloscanner/tools/package_install.R
# If ggtree fails to install, this could be because of an out of date BiocInstaller package.
# Removing it worked for me: entering an interactive R session by typing 'R' at the command
# line, running the single command
remove.packages("BiocInstaller")
# exiting the R session (Ctrl+D), then rerunning
sudo Rscript ~/phyloscanner/tools/package_install.R

# Finally, ensure all installed R packages are up-to-date by entering an interactive R session
# and running
library(rvcheck)
update_all()
