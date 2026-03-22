#!/usr/bin/env bash

# set download locations ======================================================
cluster=arc78@mccleary.ycrc.yale.edu:/home/arc78/csf1/bulk/BafA1.MEF/data/02_kallisto
computer=~/Dropbox/research/projects/csf1/bulk/BafA1.MEF/data/

# download ====================================================================

rsync -avz --exclude "*.h5" $cluster $computer
