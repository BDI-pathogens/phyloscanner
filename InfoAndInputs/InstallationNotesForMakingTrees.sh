# INSTALLING DEPENDENCIES OF phyloscanner_make_trees.py

# Note that phyloscanner_make_trees.py is written in python 2 (not 3), so biopython and pysam
# need to be installed as part of your python 2 installation.

################################################################################
# ON LINUX:

# RAxML can be installed in different ways from https://github.com/stamatak/standard-RAxML
# Here follows one way to install it directly from the command line.
# The three different 'make' commands try to compile faster versions
# suitable for more recent processors where possible.
sudo apt install git
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML/
make -f Makefile.AVX.gcc; rm *.o
make -f Makefile.SSE3.gcc; rm *.o
make -f Makefile.gcc; rm *.o
cd ..
# Optionally, add RAxML to your PATH environment variable:
echo 'PATH=$PATH:~/standard-RAxML/' >> ~/.bashrc
source ~/.bashrc

# MAFFT can be installed in different ways from http://mafft.cbrc.jp/alignment/software/linux.html
# Here follows one way to install version 7.310 directly from the command line.
wget http://mafft.cbrc.jp/alignment/software/mafft-7.310-without-extensions-src.tgz
tar -xzf mafft-7.310-without-extensions-src.tgz
cd mafft-7.310-without-extensions/core/
make clean
make
sudo make install
cd ../..

# samtools
sudo apt-get install samtools
# Test it works by running the command 'samtools'

# Get pip
sudo apt-get install python-pip
pip install --upgrade pip

# biopython
sudo pip install biopython
# Test it works by running the command 'python' to start an interactive python
# session, then typing 'import Bio'.

# cython
sudo apt-get install cython

# We need pysam version at least 0.8.1. It can be tricky to install.
# Installing these packages may help subsequently:
sudo apt-get install libz-dev zlib1g-dev libxml2-dev libxslt1-dev libbz2-dev liblzma-dev
# Then
git clone https://github.com/pysam-developers/pysam.git
cd pysam/
python setup.py build
sudo python setup.py install
# Test it works by running the command 'python' to start an interactive python
# session, then typing 'import pysam'. If you get an error, try closing your terminal,
# reopening it and trying again: bizzarely this has worked for me for one error.

# Optional: the python module matplotlib. If installed, it is used by the helper
# script tools/EstimateReadCountPerWindow.py to plot its output data.
sudo apt-get install python-matplotlib
# Test it works by running the command 'python' to start an interactive python
# session, then typing 'import matplotlib.pyplot'.

# This completes the installation of the dependencies of phyloscanner_make_trees.py.
# Of course to run phyloscanner you'll need to download the phyloscanner code.
# If you haven't already done that, do this:
git clone https://github.com/BDI-pathogens/phyloscanner.git

# Optionally, add phyloscanner to your PATH environment variable: if you downloaded
# the phyloscanner code directly into your home directory, that's achieved with
echo 'PATH=$PATH:~/phyloscanner/' >> ~/.bashrc
# otherwise it's the above command with the directory modified accordingly. Then
source ~/.bashrc
################################################################################

################################################################################
# ON MAC OS:

# Get xcode
xcode-select --install

# RAxML can be installed in different ways from https://github.com/stamatak/standard-RAxML
# Here follows one way to install it directly from the command line.
# The three different 'make' commands try to compile faster versions
# suitable for more recent processors where possible.
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML/
make -f Makefile.AVX.gcc; rm *.o
make -f Makefile.SSE3.gcc; rm *.o
make -f Makefile.gcc; rm *.o
cd ..
# Optionally, add RAxML to your PATH environment variable:
echo 'PATH=$PATH:~/standard-RAxML/' >> ~/.bashrc
source ~/.bashrc

# MAFFT can be installed in different ways from http://mafft.cbrc.jp/alignment/software/linux.html
# Here follows one way to install version 7.310 directly from the command line.
curl http://mafft.cbrc.jp/alignment/software/mafft-7.310-without-extensions-src.tgz > mafft-7.310-without-extensions-src.tgz 
tar -xzf mafft-7.310-without-extensions-src.tgz
cd mafft-7.310-without-extensions/core/
make clean
make
sudo make install
cd ../..

# home brew: 
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# samtools
brew tap homebrew/science
brew install samtools
# Test it works by running the command 'samtools'

# Install up-to-date python
brew install python

# pip
sudo easy_install pip

# biopython
pip install biopython
# Test it works by running the command 'python' to start an interactive python
# session, then typing 'import Bio'.

# pysam (we need version 0.8.1 or later)
pip install pysam --upgrade
# Test it works by running the command 'python' to start an interactive python
# session, then typing 'import pysam'. If you get an error, try closing your terminal,
# reopening it and trying again: bizzarely this has worked for me for one error.
################################################################################

# Note that the python modules below are also required, however unlike pysam and Biopython,
# these would normally be included in a standard installation of the python language:
# os, collections, itertools, subprocess, sys, re, copy, shutil, glob, time, argparse
