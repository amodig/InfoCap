#!/bin/sh

logpath='./log/batcher.log'
# remove old output
#rm ${logpath}

r_exec='sudo -H -u amodig /usr/bin/Rscript'

# argument for R (script)
cmdR="${1}"

# Note: will concatenate log file
date >> ${logpath}

# execute R script
echo 'Running Rscript:'
echo ${cmdR}
(${r_exec} ${cmdR} 2>&1) | tee -a ${logpath} 

#rm matlab_command.m
