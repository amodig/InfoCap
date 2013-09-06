#!/bin/sh

logpath='./log/batcher.log'
matlab_exec='sudo -H -u amodig /usr/local/MATLAB/R2012a/bin/matlab'
r_exec='sudo -H -u amodig /usr/bin/Rscript'

# first argument for MATLAB (CTW)
cmdM1="${1}"
# second argument for MATLAB (GP)
cmdM2="${2}"
# third argument for R (TPs)
cmdR="${3}"

# NOTE: will concatenate all outputs into a single file
'---------------------------------------------------' >> output.log
date >> output.log

echo 'Running MATLAB command:'
echo ${cmdM1} > matlab_command1.m
cat matlab_command1.m
# execute MATLAB script (CTW)
# stdout and stderr to stdout and append to log  
(${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command1.m 2>&1) | tee -a ${logpath} 

echo 'Running MATLAB command:'
echo ${cmdM2} > matlab_command2.m
cat matlab_command2.m
# execute MATLAB script (GP)
(${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command2.m 2>&1) | tee -a ${logpath} 

# execute R script (TP)
echo 'Running Rscript:'
echo ${cmdR}
(${r_exec} ${cmdR} 2>&1) | tee -a ${logpath}
