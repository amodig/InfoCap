#!/bin/sh

logpath='./log/batcher.log'
# remove old output
#rm ${logpath}

matlab_exec='sudo -H -u amodig /usr/local/MATLAB/R2012a/bin/matlab'
r_exec='sudo -H -u amodig /usr/bin/Rscript'

# first argument for MATLAB (function)
cmdMATLAB="${1}"
# second argument for R (script)
cmdR="${2}"

# Note: will concatenate log file
'---------------------------------------------------' >> output.log
date >> ${logpath}

# execute MATLAB script
echo 'Running MATLAB command:'
echo ${cmdMATLAB} > matlab_command.m
cat matlab_command.m
# stdout and stderr to stdout and append to log  
(${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command.m 2>&1) | tee -a ${logpath} 

# execute R script
echo 'Running Rscript:'
echo ${cmdR}
(${r_exec} ${cmdR} 2>&1) | tee -a ${logpath} 

#rm matlab_command.m
