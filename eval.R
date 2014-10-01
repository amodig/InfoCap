# General R software for evaluating movement throughput for data sequences

# The MIT License (MIT)
# Copyright (c) 2013 Arttu Modig
# Version 1.03

source("infocapacity.R")
source("data_handling.R")
library("methods")
#options(error = quote({traceback(); q()}))

# subdir_residual_complexity
#
# Calculate residual complexity assuming that each sequence type is in
# its own subdirectory. The original sequence files for different
# sequences should be in subdirectories with zero-padding: 01*, 02*, 03* etc.
# These subdirectories in turn should all have a subdirectory called
# aligneddata that contains all sequence repetitions aligned with
# each other.
#
# Parameters:
# - fps                     frames per second in the sequences
# - dim_reduction           which dimension reduction method to use
# - calculate_residuals     whether calculate residuals (default: TRUE)
# - save_residuals          whether save residuals (default: FALSE)
# - save_prediction         whether save predictions (default: FALSE)
# - save_pca                whether save PCA components
# - normalize               whether normalize features before analysis
# - residual_dir            name of the residual data folder
# - prediction_dir          name of the prediction data folder (if saving predictions)
# - pca_dir                 name of the PCA folder (if saving PCA components)
# - gp_dir                  name of the GP-LVM folder (if using the components)
# - remove_duplicates       whether remove duplicates after alignment (default: TRUE)
# - feature_throughputs     whether save feature throughputs
# - align                   whether time-aligning data or not
# - write_results           whether write results in to an ASCII file
# - alignment_dir           name of the alignment index folder
# - dimred_dir              (not in use)
# - results_file            name of the TSV-ASCII results file
#
# Example of the directory tree:
# - working directory
#   |
#   -- 01 -- contains three repetitions of sequence #1
#   |  |
#   |  * 01.txt
#   |  * 02.txt
#   |  * 03.txt
#   |  -- alignment
#   |  |  |
#   |  |  * 1_ali_2.txt
#   |  |  * 1_ali_3.txt
#   |  |  * ...
#   |  -- gp (with GP-LVM method)
#   |     |
#   |     * 01.txt
#   |     * 02.txt
#   |     * 03.txt
#   -- 02 -- contains two repetitions of sequence #2
#      |
#      * 01.txt
#      * 02.txt
#      -- alignment
#      |  |
#      |  * 1_ali_2.txt
#      |  * 2_ali_1.txt
#      -- gp (with GP-LVM method)
#         |
#         * 01.txt
#         * 02.txt
#
subdir_residual_complexity <- function(fps = 120, dim_reduction = "pca",
                            calculate_residuals = TRUE, save_residuals = FALSE,
                            save_prediction = FALSE, save_pca = FALSE,
                            normalize = FALSE, residual_dir = "residuals",
                            prediction_dir = "prediction",
                            pca_dir = "pca", gp_dir = "gp", remove_duplicates = TRUE,
                            feature_throughputs = FALSE, align = TRUE,
                            write_results = FALSE, alignment_dir = "alignment",
                            dimred_dir = ".", results_file = "results.txt") {
    # get maindir path
    maindir <- getwd()
    
    #subdirs <- dir(".", "^[[:digit:]][[:digit:]]+$")
    subdirs <- dirs() # in data_handling.R
    results <- list()
    for (i in 1:length(subdirs)) {
        print(sprintf("------- Evaluating subdir: %s -------", subdirs[i]))
        # sets the path to subdir
        setwd(subdirs[i])
        # sets the path to dimension reduced data folder
        setwd(dimred_dir)
        # count the number of txt-type (ASCII) data
        #sequences <- length(dir(".", "^[[:digit:]]+.txt$"))
        sequences <- dir(".", "*.(txt|tsv)$")
        N <- length(sequences)
        
        rownames <- c()
        combos <- combinations(N, 2) * 2
        col_names <- c("TP", "shared", "RSS", "RSS_cond", "quotient")
        all_results <- matrix(nrow = combos, ncol = 5,
                              dimnames = list(1:combos, col_names))
        rows <- c(1, 2)
        
        for (j in 1:(N-1)) {
            for (k in (j+1):N) {
                rownames = append(rownames, sprintf("(%s, %s)",
                                                    sequences[j], sequences[k]))
                rownames = append(rownames, sprintf("(%s, %s)",
                                                    sequences[k], sequences[j]))
                M <- pair_residual_complexity(
                    seqfile_a = sequences[j], seqfile_b = sequences[k],
                    seqnum_a = j, seqnum_b = k, fps = fps,
                    dim_reduction = dim_reduction,
                    feature_throughputs = feature_throughputs, features = 0,
                    calculate_residuals = calculate_residuals,
                    save_residuals = save_residuals,
                    save_prediction = save_prediction,
                    save_pca = save_pca, normalize = normalize,
                    residual_dir = residual_dir,
                    prediction_dir = prediction_dir, pca_dir = pca_dir,
                    remove_duplicates = remove_duplicates,
                    align = align, alignment_dir = alignment_dir)
                all_results[rows[1]:rows[2],] <- M$results
                rows <- rows + 2
                if (feature_throughputs) {
                    all_features <- M$features
                }
            }
        }
        rownames(all_results) <- rownames
        if (feature_throughputs)
            results[[i]] <- list(dir = subdirs[i], results = all_results,
                                feature_results = all_features)
        else
            results[[i]] <- list(dir = subdirs[i], results = all_results)
        
        # set path back to maindir
        setwd(maindir)
        if (write_results) {
            # write dir name
            write(sprintf("Folder: %s", subdirs[i]), file=results_file,
                            append = TRUE, ncolumns=1)
            # write dir results
            result_frame <- results[[i]]$results
            write.table(result_frame, file=results_file, sep='\t',
                            append=TRUE, col.names=TRUE, row.names=TRUE)
        }
    }
    print("Finished evaluating subdirs.")
    return(results)
}

# Calculate the residual complexity of all sequences cross-wise
# in the current directory (with the directory structure specified
# in the readme file). Sequences are considered to be repetitions of
# a same motion/gesture.
singledir_residual_complexity <- function(fps = 120, dim_reduction = "pca",
                                calculate_residuals = TRUE,
                                save_residuals = FALSE,
                                save_prediction = FALSE,
                                save_pca = FALSE,
                                normalize = FALSE,
                                residual_dir = "residuals",
                                prediction_dir = "prediction",
                                pca_dir = "pca",
                                gp_dir = "gp",
                                remove_duplicates = TRUE,
                                feature_throughputs = FALSE,
                                align = TRUE,
                                alignment_dir = "alignment",
                                dimred_dir = ".") {
    # sets the path to dimension reduced data folder
    setwd(dimred_dir)
    #filenames <- dir("alignment", "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    #sequences <- length(filenames)
    sequences <- dir(".", "*.(txt|tsv)$")
    N <- length(sequences)
    
    rownames <- c()
    combos <- combinations(N, 2) * 2
    col_names <- c("TP", "shared", "RSS", "RSS_resid", "quotient")
    TP_results <- matrix(nrow = combos, ncol = 5,
        dimnames = list(1:combos, col_names))
    rows <- c(1, 2)

    for (j in 1:(N-1)) {
        for (k in (j+1):N) {
            rownames = append(rownames, sprintf("(%s, %s)",
                                                sequences[j], sequences[k]))
            rownames = append(rownames, sprintf("(%s, %s)",
                                                sequences[k], sequences[j]))
            M <- pair_residual_complexity( seqfile_a = sequences[j],
                seqfile_b = sequences[k], seqnum_a = j, seqnum_b = k,
                fps = fps, dim_reduction = dim_reduction,
                feature_throughputs = feature_throughputs,
                features = 0, calculate_residuals = calculate_residuals,
                save_residuals = save_residuals,
                save_prediction = save_prediction,
                save_pca = save_pca, normalize = normalize,
                residual_dir = residual_dir,
                prediction_dir = prediction_dir, pca_dir = pca_dir,
                remove_duplicates = remove_duplicates,
                align = align, alignment_dir = alignment_dir)
            TP_results[rows[1]:rows[2],] <- M$results
            rows <- rows + 2
            if (feature_throughputs) {
                all_features <- M$features
            }
        }
    }
    rownames(TP_results) <- rownames
    if (feature_throughputs)
        all_results <- list(results = TP_results, feature_results = all_features)
    else 
        all_results <- list(results = TP_results)
        
    setwd('..')
    print("Done")
    return(all_results)
}

# Calculate the residual complexity of all pairs of sequences
# in the current directory (with the directory structure specified
# in the readme file). Sequences with consecutive numbers are
# assumed to be pairs (e.g. 1 and 2, 3 and 4 and so on).
pairdir_residual_complexity <- function(fps = 120, dim_reduction = "pca",
                                calculate_residuals = TRUE, save_residuals = FALSE,
                                save_prediction = FALSE, save_pca = FALSE,
                                normalize = FALSE, residual_dir = "residuals",
                                prediction_dir = "prediction", pca_dir = "pca",
                                gp_dir = "gp", remove_duplicates = TRUE,
                                feature_throughputs = FALSE, align = TRUE,
                                alignment_dir = "alignment", dimred_dir = ".") {
    dir <- getwd()
    print(sprintf("Working dir: %s", dir))
    print(sprintf("Alignment dir: %s", alignment_dir))
    
    setwd(dimred_dir)
    sequences <- dir(".", "*.(txt|tsv)$")
    N <- length(sequences)
    
    TP_results <- matrix(nrow = N, ncol = 5,
        dimnames = list(1:N, c("TP", "shared", "RSS", "RSS_resid", "quotient")))
    all_features <- list()

    for (i in 1:(N/2)) {
        j = 2 * i - 1
        k = j + 1      
        if (feature_throughputs) {
            M <- pair_residual_complexity(seqfile_a = sequences[j],
                seqfile_b = sequences[k], seqnum_a = j, seqnum_b = k, fps = fps,
                dim_reduction = dim_reduction, feature_throughputs = TRUE,
                features = 0, calculate_residuals = calculate_residuals,
                save_residuals = save_residuals, save_prediction = save_prediction,
                save_pca = save_pca, normalize = normalize, residual_dir = residual_dir,
                prediction_dir = prediction_dir, pca_dir = pca_dir, gp_dir = gp_dir,
                remove_duplicates = remove_duplicates, align = align,
                alignment_dir = alignment_dir)
            TP_results[j:k,] <- M$results
            all_features[[i]] <- M$features
        }
        else {
            TP_results[j:k,] <- pair_residual_complexity(seqfile_a = sequences[j],
                seqfile_b = sequences[k], seqnum_a = j, seqnum_b = k, fps = fps,
                dim_reduction = dim_reduction, feature_throughputs = FALSE,
                features = 0, calculate_residuals = calculate_residuals,
                save_residuals = save_residuals, save_prediction = save_prediction,
                save_pca = save_pca, normalize = normalize, residual_dir = residual_dir,
                prediction_dir = prediction_dir, pca_dir = pca_dir, gp_dir = gp_dir,
                remove_duplicates = remove_duplicates, align = align,
                alignment_dir = alignment_dir
                )$results
        }
    }
    
    if (feature_throughputs)
        all_results <- list(results = TP_results, feature_results = all_features)
    else 
        all_results <- list(results = TP_results)
    
    setwd(dir)
    print("Done")
    return(all_results)
}

# (Every other batch function boils down to this.)
# Calculate the residual complexity of two sequences numbered
# seqnum_a and seqnum_b. Throughputs per feature can also be printed.
#
# Parameters:
# - seqfile_a   name of the first sequence file
# - seqfile_b   name of the second sequence file
# - seqnum_a    number of first sequence
# - seqnum_b    number of second sequence
# - fps         data frames per second
# - dim_reduction method for dimension reduction
# - feature_throughputs   save throughput per each feature
# - features    select certain features (all if 0)
# - calculate_residuals calculate residuals; if false, use pre-calculated residuals 
# - save_residuals  saves the residuals
# - save_prediction saves the predicted data
# - save_pca        save PCA data
# - normalize       whether normalize sequences
# - residual_dir     folder name for residual data
# - prediction_dir  folder name for predicted data
# - pca_dir         folder name for PCA data
# - gp_dir          fodler name for GP-LVM data
# - remove_duplicates   whether duplicates from CTW are removed
# - align               whether align sequences
# - alignment_dir       folder name for alignment indeces
pair_residual_complexity <- function(seqfile_a, seqfile_b, seqnum_a = 1, seqnum_b = 2,
                                fps = 120, dim_reduction = "pca",
                                feature_throughputs = FALSE, features = 0,
                                calculate_residuals = TRUE, save_residuals = FALSE,
                                save_prediction = FALSE, save_pca = FALSE,
                                normalize = FALSE, residual_dir = "residuals",
                                prediction_dir = "prediction", pca_dir = "pca",
                                gp_dir = "gp", remove_duplicates = TRUE,
                                align = TRUE, alignment_dir = "alignment") {
    if (save_residuals) calculate_residuals <- TRUE
    if (save_prediction) calculate_residuals <- TRUE
    
    print(sprintf("Evaluating throughput between %s and %s", seqfile_a, seqfile_b))
    
    # create evaluator object
    pair_a_b <- new("residualComplexityEvaluator", seqfile_a = seqfile_a,
                    seqfile_b = seqfile_b, seqnum_a = seqnum_a, seqnum_b = seqnum_b,
                    fps = fps, dim_reduction = dim_reduction, features = features,
                    save_pca = save_pca, calculate_residuals = calculate_residuals,
                    normalize = normalize, residual_dir = residual_dir, pca_dir = pca_dir,
                    gp_dir = gp_dir, remove_duplicates = remove_duplicates, align = align,
                    alignment_dir = alignment_dir)
    # compute complexity
    results_a_b <- evaluate_complexity(pair_a_b)
    
    # create evaluator object
    pair_b_a <- new("residualComplexityEvaluator", seqfile_a = seqfile_b,
                    seqfile_b = seqfile_a, seqnum_a = seqnum_b, seqnum_b = seqnum_a,
                    fps = fps, dim_reduction = dim_reduction, features = features,
                    save_pca = save_pca, calculate_residuals = calculate_residuals,
                    normalize = normalize, residual_dir = residual_dir, pca_dir = pca_dir,
                    gp_dir = gp_dir, remove_duplicates = remove_duplicates, align = align,
                    alignment_dir = alignment_dir)
    # compute complexity
    results_b_a <- evaluate_complexity(pair_b_a)

    resultvectors <- rbind(construct_result_vector(results_a_b),
           construct_result_vector(results_b_a), deparse.level = 0)
    
    if (feature_throughputs) {
        feature_tp <- rbind(get_feature_throughputs(results_a_b),
                        get_feature_throughputs(results_b_a),
                        deparse.level = 0)   
        resultlist <- list("results" = resultvectors, "feature_tp" = feature_tp,
                           "results_a_b" = results_a_b, "results_b_a" = results_b_a)
    } else
      resultlist <- list("results" = resultvectors,
                         "results_a_b" = results_a_b, "results_b_a" = results_b_a)
        
    if (save_residuals) {
      write_data(dir=residual_dir, sequence_file=seqfile_a, data=results_a_b@res_a)
      write_data(dir=residual_dir, sequence_file=seqfile_b, data=results_b_a@res_b)
    }
    
    if (save_prediction) {
      write_data(dir=prediction_dir, sequence_file=seqfile_a, data=results_a_b@pred_a)
      write_data(dir=prediction_dir, sequence_file=seqfile_b, data=results_b_a@pred_b)
    }
    
    cat("Pair TP: \n")
    cat(sprintf("(%.2f, %.2f)", resultvectors[1,1], resultvectors[2,1]))
    cat("\n")
    
    return(resultlist)
    #return(rbind(construct_result_vector(results_a_b),
    #       construct_result_vector(results_b_a), deparse.level = 0))
}

# Create a union class for data files which can be either vectors or matrices
setClassUnion("data", c("vector","matrix"))

# A class for performing the required residual complexity evaluations,
# including the loading of files, calculation of residuals etc.
setClass("residualComplexityEvaluator", representation(
  seqfile_a = "character", # name of file A
  seqfile_b			= "character", # name of file B
  seqnum_a			= "numeric",	# number of file A
  seqnum_b			= "numeric",  # number of file B
  fps						= "numeric",  # framerate
  dim_reduction = "character", # choose reduction method
  features			= "numeric",	# choose features
  save_residuals = "logical",	# option to save residual data
  save_prediction = "logical",    # option to save prediction data
  save_pca			= "logical", # option to save pca-reduced data
  calculate_residuals = "logical", # option to calculate residuals
  normalize			= "logical", # option to normalize sequences
  residual_dir		= "character", # dir of residuals
  alignment_dir	= "character", # dir of alignment indeces
  prediction_dir = "character", # dir of prediction data
  pca_dir				= "character", # dir of pca-reduced sequences
  gp_dir				= "character", # dir where gplvm models are saved
  remove_duplicates	= "logical", # whether to remove duplicate frames
  align					= "logical", # whether align residuals
  seq_a         = "data", # original seq a data
  seq_b         = "data", # original seq b data
  res_orig_a		= "data",	# original residuals a
  res_orig_b		= "data",	# original residuals b
  res_a					= "data",	# processed residual sequence a
  res_b					= "data",	# processed residual sequence b
  pred_a				= "data",	# prediction a
  pred_b				= "data",	# prediction b
  res_pca_a 	  = "data",
  res_pca_b			= "data",
  eig						= "NULL",
  pcv_a					= "matrix",
  pcv_b					= "matrix",
  results				= "list"))	# results after evaluation is done

# Evaluates complexity and returns a new residualComplexityEvaluator
# with the 'results' field set to the result list returned by the
# evaluation function. See construct_result_vector and get_feature_throughputs
# for result formatting and throughput calculation.
setGeneric("evaluate_complexity",
           function(this) standardGeneric("evaluate_complexity"))
setMethod("evaluate_complexity", "residualComplexityEvaluator",
           function(this) {
               # load residuals (and prediction)
               this <- load_residuals(this)
               # check length
               if (length(this@res_a)==0) {
                 stop("Empty residual matrix!")
               }
               # check for dimensions
               if (any(dim(this@res_a) != dim(this@res_b)))
                    stop("Cut and/or aligned residual dimensions don't match!")
               
                if (any(is.na(this@res_a)) || any(is.na(this@res_b)))
                     stop("Residual sequences contain NaNs!")
               
               this <- do_reduction(this)
               if (this@dim_reduction == "pca") {
                 this@results <- evaluate_residual_shared_information(this@res_pca_a,
                                                                      this@res_pca_b)
               } else {
                 this@results <- evaluate_residual_shared_information(this@res_a,
                                                                      this@res_b)
               }
               this
           })

# Loads and aligns the residuals indicated by the two sequence numbers.
setGeneric("load_residuals",
           function (this) standardGeneric("load_residuals"))
setMethod("load_residuals", "residualComplexityEvaluator",
          function(this) {
            if (this@align) {
              data <- load_aligned_residuals(
                this@seqfile_a,
                this@seqfile_b,
                this@seqnum_a,
                this@seqnum_b,
                this@features,
                this@calculate_residuals,
                this@normalize,
                this@remove_duplicates,
                this@alignment_dir,
                this@residual_dir)
            } else {
              data <- load_unaligned_residuals(
                this@seqfile_a,
                this@seqfile_b,
                this@seqnum_a,
                this@seqnum_b,
                this@features,
                this@calculate_residuals,
                this@normalize,
                this@remove_duplicates,
                this@residual_dir)
            }
            # original sequences
            this@seq_a <- data$sequence_a
            this@seq_b <- data$sequence_b
            # original residuals
            this@res_orig_a <- data$residuals_orig_a
            this@res_orig_b <- data$residuals_orig_b
            # aligned residuals
            this@res_a <- data$residuals_alig_a
            this@res_b <- data$residuals_alig_b
            # predicted model
            this@pred_a <- data$prediction_a
            this@pred_b <- data$prediction_b
            # return
            this
          })

# Performs dimension reduction on the two loaded sequences.
setGeneric("do_reduction",
           function(this) standardGeneric("do_reduction"))
setMethod("do_reduction", "residualComplexityEvaluator",
          function(this) {
              if (this@dim_reduction == "pca") {
                  eigenvectors <- pca(this@res_a)
                  this@res_pca_a <- normalize_features(this@res_a) %*% eigenvectors
                  this@res_pca_b <- normalize_features(this@res_b) %*% eigenvectors
              } else if (this@dim_reduction == "pca2") {
                  data_a <- normalize_features(this@res_a)
                  data_b <- normalize_features(this@res_b)
                  data <- rbind(data_a,data_b)
                  prc <- prcomp(data, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  this@res_pca_a <- prc$x[1:nrow(data_a),]
                  this@res_pca_b <- prc$x[(nrow(data_a)+1):(2*nrow(data_a)),]
              } else if (this@dim_reduction == "pca3") {
                  pc_a <- prcomp(this@res_a, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  pc_b <- prcomp(this@res_b, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  this@res_pca_a <- pc_a$x
                  this@res_pca_b <- pc_b$x
              } else if (this@dim_reduction == "gplvm") {
                  # Just use precomputed models.
                  this@res_pca_a <- load_model(this@seqnum_a, this@gp_dir)
                  this@res_pca_b <- load_model(this@seqnum_b, this@gp_dir)
			        }							
              # optionally write PCA data
              if (this@save_pca) {
                  write_data(this@pca_dir, this@seqfile_a, this@res_pca_a)
                  write_data(this@pca_dir, this@seqfile_b, this@res_pca_b)
              }
              this
          })

# Generates a vector with throughput, shared information
# in bits, total sum of residuals, total conditional sum of
# residuals and the residual quotient. (Vector length is 5.)
setGeneric("construct_result_vector",
           function(this) standardGeneric("construct_result_vector"))
setMethod("construct_result_vector", "residualComplexityEvaluator",
          function(this) {
              quotient <- this@results$RSS / this@results$RSS_conditional
              throughput <- calculate_throughput(this, this@results$total_shared)
              shared_information_bits <- this@results$total_shared / log(2.0)

              c(throughput, shared_information_bits, this@results$RSS,
                this@results$RSS_conditional, quotient)
          })

# Generates a vector of feature throughputs. Can only be
# run after evaluation.
setGeneric("get_feature_throughputs",
           function(this) standardGeneric("get_feature_throughputs"))
setMethod("get_feature_throughputs", "residualComplexityEvaluator",
          function(this) {
              shared <- this@results$feature_shared # this is a list
              unlist(lapply(shared, function(x) calculate_throughput(this, x)))
          })

# Generates a vector of original feature throughputs.
setGeneric("get_orig_feature_throughputs",
           function(this) standardGeneric("get_orig_feature_throughputs"))
setMethod("get_orig_feature_throughputs", "residualComplexityEvaluator",
          function(this) {
            results <- evaluate_residual_shared_information(this@res_a, this@res_b)
            shared <- results$feature_shared # this is a list
            unlist(lapply(shared, function(x) calculate_throughput(this, x)))
          })

# Given a shared information value, calculate throughput.
setGeneric("calculate_throughput",
           function(this, shared) standardGeneric("calculate_throughput"))
setMethod("calculate_throughput", "residualComplexityEvaluator",
          function(this, shared) {
              shared / NROW(this@res_a) * this@fps / log(2.0)
          })

# Calculate the number of k-combinations in a set with n elements.
combinations <- function(n, k) {
    return(factorial(n) / (factorial(k) * factorial(n - k)))
}

# Calculates the averages of all individual features in a subdir results list.
# Note that the number of dimensions in each subdir component must be equal.
average_feature_throughputs <- function(subdirlist=NULL, featurelist=NULL) {
    n <- 0
    sum_ftp <- 0
    average_ftp <- 0
    
    if (!is.null(subdirlist)) {
        n <- length(subdirlist)
        sum_ftp <- array(0, dim = dim(subdirlist[[1]]$feature_results))
        for (i in 1:n) {
            sum_ftp <- sum_ftp + subdirlist[[i]]$feature_results
        }
    }
    else if (!is.null(featurelist)) {
        n <- length(featurelist)
        sum_ftp <- array(0, dim = dim(featurelist[[1]]))
        for (i in 1:n) {
            sum_ftp <- sum_ftp + featurelist[[i]]
        }
    }
    average_ftp <- sum_ftp / n
    return(average_ftp)
}
