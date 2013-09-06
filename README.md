Information Capacity
====================

MATLAB and R code for Information Capacity Computation.

Requirements
------------
1) Canonical Time Warping (CTW) code
http://www.f-zhou.com/ta_code.html

2) Var-GP-LVM software including the other Matlab toolboxes from N.Lawrence
http://staffwww.dcs.shef.ac.uk/people/N.Lawrence/vargplvm/

Extract CTW code to "ctw" subfolder (you should run "make" according to the CTW installation instructions) and all the GP-related toolboxes to "gp" subfolder. The MATLAB scripts will load the paths then automatically.

You'll also probably need to install R packages called "methods" and "kernlab".


Brief How-To-Use
----------------
The server-side uses shell scripts to execute both R and Matlab scripts. There's 6 main Matlab functions: 4 CTW aligning functions (align_pairdir, align_dir, align_subdir, align_pair) and 2 GP function (evalGPdir, evalVarGP). The R functions are controlled with Rscript frontend "frontend.R".

There's shell scripts for:

1) Just for R code (to use without alignment):
"batcher-R.sh"
- 1 argument for Rscript

2) Matlab and R with CTW:
"batcher.sh"
- first arg for Matlab (CTW function)
- second arg for Rscript

3) Matlab and R with CTW and GP:
"batcher-gp.sh"
- first arg for Matlab (CTW function)
- second arg for Matlab (GP function)
- third arg for Rscript

To use them, you have to change Matlab and Rscript exec paths.


Examples
--------
Parameter-blown examples for shell:

1) With 2 files + CTW + PCA:
$ ./batcher.sh "align_pair('./path/to/dir/testSeq1.txt', './path/to/dir/testSeq2.txt','./path/to/dir/alignment/')" "./frontend.R --fps 120 --method pairdir --maindir ./path/to/dir/ --dim_reduction pca --calculate_residuals TRUE --save_residuals FALSE --residualdir residuals/ --resultfile results.txt --align TRUE"

2) With 2 files + CTW + Var-GP:
$ ./batcher_gp.sh "align_pair('./path/to/dir/9.txt', './path/to/dir/10.txt','./path/to/dir/alignment/')" "evalGPdir('./path/to/dir/', './path/to/dir/gp')" "./frontend.R --fps 200 --method pairdir --maindir ./path/to/dir/gp --dim_reduction gp --calculate_residuals TRUE --save_residuals FALSE --residualdir residuals/ --resultfile results.txt --align TRUE --alidir ../alignment"
