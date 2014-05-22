#!/usr/bin/Rscript --vanilla --verbose

## This R script is a front-end for running the infocapacity functions from command line

# The MIT License (MIT)
# Copyright (c) 2013 Arttu Modig
# Version 1.03

cat("Starting frontend.R\n")
cat("Arguments passed to Rscript:\n")
cat(commandArgs(TRUE))
cat("\n")

library("getopt")

# Rscript does not, by default, load the “methods” package.
require(methods)
# Print traceback when there is an error.
options(error = quote({traceback(); q()}))

# Read functions
source('eval.R')
# sink error messages to stdout
sink(stdout(), type="message")

# Output redirection
# Send messages and output to separate files
#msg <- file("Rmessage.log", open="wt")
#out <- file("Routput.log", open="wt")
#sink(msg, type="message")
#sink(out, type="output")

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c( # short flags in use: a b c d e f g h i j k l m n o p q r s t u v w x y z
    'resultfile', 'o', 2, "character",
    'fps', 'f', 2, "integer",
    'method', 'm', 2, "character",
    'maindir', 'd', 2, "character",
    'seqfile_a', 'a', 2, "character",
    'seqfile_b', 'b', 2, "character",
    'seqnum_a', 'x', 2, "integer",
    'seqnum_b', 'y', 2, "integer",
    'dim_reduction', 'i', 2, "character",
    'calculate_residuals', 'c', 2, "logical",
    'save_residuals', 's', 2, "logical",
    'save_prediction', 'j', 2, "logical",
    'save_pca', 'p', 2, "logical",
    'normalize', 'n', 2, "logical",
    'residualdir', 'r', 2, "character",
    'prediction_dir', 'k', 2, "character",
    'remove_duplicates', 'e', 2, "logical",
    'feature_throughputs', 't', 2, "logical",
    'select_features', 'l', 2, "character",
    'username', 'u', 2, "character",
    'append_results', 'v', 2, "logical",
    'align', 'g', 2, "logical",
    'write_results', 'w', 2, "logical",
    'alidir', 'z', 2, "character",
    'pca_dir', 'q', 2, "character",
	  'gp_dir', 'h', 2, "character"), ncol=4,byrow=TRUE);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$resultfile ) ) { opt$resultfile = "results.txt" }
if ( is.null(opt$fps ) ) { opt$fps = 120 }
if ( is.null(opt$method ) ) { opt$method = "subdir" }
if ( is.null(opt$maindir ) ) { opt$maindir = "." }
if ( is.null(opt$seqfile_a ) ) { opt$seqfile_a = "01.txt" }
if ( is.null(opt$seqfile_b ) ) { opt$seqfile_b = "02.txt" }
if ( is.null(opt$seqnum_a ) ) { opt$seqnum_a = 1 }
if ( is.null(opt$seqnum_b ) ) { opt$seqnum_b = 2 }
if ( is.null(opt$dim_reduction ) ) { opt$dim_reduction = "pca" }
if ( is.null(opt$calculate_residuals ) ) { opt$calculate_residuals = TRUE }
if ( is.null(opt$save_residuals ) ) { opt$save_residuals = FALSE }
if ( is.null(opt$save_prediction ) ) { opt$save_prediction = FALSE }
if ( is.null(opt$save_pca ) ) { opt$save_pca = FALSE }
if ( is.null(opt$normalize ) ) { opt$normalize = FALSE }
if ( is.null(opt$residualdir ) ) { opt$residualdir = "residuals" }
if ( is.null(opt$prediction_dir ) ) { opt$prediction_dir = "prediction" }
if ( is.null(opt$alidir ) ) { opt$alidir = "alignment" }
if ( is.null(opt$pca_dir ) ) { opt$pca_dir = "pca" }
if ( is.null(opt$gp_dir ) ) { opt$gp_dir = "gp" }
if ( is.null(opt$remove_duplicates ) ) { opt$remove_duplicates = TRUE }
if ( is.null(opt$feature_throughputs ) ) { opt$feature_throughputs = FALSE }
if ( is.null(opt$select_features ) ) { opt$select_features = "0" }
if ( is.null(opt$append_results ) ) { opt$append_results = FALSE }
if ( is.null(opt$align ) ) { opt$align = TRUE }
if ( is.null(opt$write_results ) ) { opt$write_results = TRUE }

cat("Use method: ")
cat(opt$method)
cat("\n")

# current working dir
cat("Current working dir: ")
cat(getwd())
cat("\n")
cat("Change to dir: ")
cat(opt$maindir)
cat("\n")

if (opt$method == "subdir") {
    # CD to directory, default "."
    try(setwd(opt$maindir))
    output = subdir_residual_complexity(
                fps = opt$fps,
                dim_reduction = opt$dim_reduction,
                calculate_residuals = opt$calculate_residuals,
                save_residuals = opt$save_residuals,
                save_prediction = opt$save_prediction,
                save_pca = opt$save_pca,
                normalize = opt$normalize,
                residualdir = opt$residualdir,
                prediction_dir = opt$prediction_dir,
                remove_duplicates = opt$remove_duplicates,
                feature_throughputs = opt$feature_throughputs,
                align = opt$align,
                write_results = opt$write_results,
                alignment_dir = opt$alidir,
                pca_dir = opt$pca_dir,
                gp_dir = opt$gp_dir,
                results_file = opt$resultfile)
} else if (opt$method == "pairdir") {
    try(setwd(opt$maindir))
    output = pairdir_residual_complexity(
                fps = opt$fps,
                dim_reduction = opt$dim_reduction,
                calculate_residuals = opt$calculate_residuals,
                save_residuals = opt$save_residuals,
                save_prediction = opt$save_prediction,
                save_pca = opt$save_pca,
                normalize = opt$normalize,
                residualdir = opt$residualdir,
                prediction_dir = opt$prediction_dir,
                remove_duplicates = opt$remove_duplicates,
                feature_throughputs = opt$feature_throughputs,
                align = opt$align,
                alignment_dir = opt$alidir,
                pca_dir = opt$pca_dir,
                gp_dir = opt$gp_dir)
} else if (opt$method == "singledir") {
    try(setwd(opt$maindir))
    output = singledir_residual_complexity(
                fps = opt$fps,
                dim_reduction = opt$dim_reduction,
                calculate_residuals = opt$calculate_residuals,
                save_residuals = opt$save_residuals,
                save_prediction = opt$save_prediction,
                save_pca = opt$save_pca,
                normalize = opt$normalize,
                residualdir = opt$residualdir,
                prediction_dir = opt$prediction_dir,
                remove_duplicates = opt$remove_duplicates,
                feature_throughputs = opt$feature_throughputs,
                align = opt$align,
                alignment_dir = opt$alidir,
                pca_dir = opt$pca_dir,
                gp_dir = opt$gp_dir)
} else if (opt$method == "pair") {
    try(setwd(opt$maindir))
    output = pair_residual_complexity(
                opt$seqfile_a,
                opt$seqfile_b,
                opt$seqnum_a,
                opt$seqnum_b,
                fps = opt$fps,
                dim_reduction = opt$dim_reduction,
                feature_throughputs = opt$feature_throughputs,
                features = as.numeric(opt$select_features),
                calculate_residuals = opt$calculate_residuals,
                save_residuals = opt$save_residuals,
                save_prediction = opt$save_prediction,
                save_pca = opt$save_pca,
                normalize = opt$normalize,
                residualdir = opt$residualdir,
                prediction_dir = opt$prediction_dir,
                remove_duplicates = opt$remove_duplicates,
                align = opt$align,
                alignment_dir = opt$alidir,
                pca_dir = opt$pca_dir,
                gp_dir = opt$gp_dir)
} else {
    stop("Unknown method. Please re-check arguments.\n")
    q(status=10)
}

# print output file
#cat(paste(output))
#cat("\n")

# write results in a separate file
if (opt$write_results) {
    if (opt$method != "subdir") { # subdirs are handled differently (results are appended)
        write.table(output$results, file = opt$resultfile, sep = '\t')
    }
}

if (opt$feature_throughputs) {
  # reduced features
  write.table(output$feature_tp, file="feature-TP.txt")
  cat("Wrote feature TPs into: \n")
  cat("feature-TP.txt\n")
  # original features
  feature_tp_orig <- rbind(get_orig_feature_throughputs(output$results_a_b),
                           get_orig_feature_throughputs(output$results_b_a),
                           deparse.level = 0)
  write.table(feature_tp_orig, file="feature-TP-orig.txt")
  cat("Wrote original feature TPs into: \n")
  cat("feature-TP-orig.txt\n")
}

if (opt$append_results) {
    try(setwd("..")) # go back up!
    score <- max(output$results[,1]) # choose maximum score of the pair
    cat("Name: ")
    cat(opt$username)
    cat("\n")
    cat("Pair: ")
    cat(output$results[,1])
    cat("\n")
    scoreline <- sprintf("%s\t%.2f", opt$username, score)
    cat(scoreline)
    cat("\n")
    write(scoreline, file="scores.txt", append=TRUE)
    cat("Appended TP-score into: \n")
    cat("scores.txt\n")
}

# signal success and exit.
cat("Successfully wrote results into: \n")
cat(opt$resultfile)
cat("\n")
q(status=0)
