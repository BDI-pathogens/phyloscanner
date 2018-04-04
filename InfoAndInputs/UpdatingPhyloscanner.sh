# To update your phyloscanner code after you've installed it: from the
# command line change directory to the main phyloscanner code directory, e.g.
cd ~/phyloscanner/
# if that is where the code lives on your machine, then run
git pull
# If you get the message "Already up-to-date." you're done.
# If any of the R code has been updated (which you can check by seeing if any
# of the names of altered files that were just printed to the screen end in
# '.R', or if you prefer to keep it simple you can just assume that the R code 
# was updated) you then need to
cd phyloscannerR
sudo R
# then inside the interactive R console, run
library(devtools)
install(".", dependencies = T)
