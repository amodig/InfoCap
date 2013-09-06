# The MIT License (MIT)
# Copyright (c) 2013 Arttu Modig

source("infocapacity.R")
source("data_handling.R")
library("methods")
library("kernlab")

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
# - save_pca                whether save PCA components
# - normalize               whether normalize features before analysis
# - residualdir             name of the residual folder (if saving residuals)
# - pca_dir                 name of the PCA folder (if saving PCA components)
# - remove_duplicates       whether remove duplicates after time-alignment (default: TRUE)
# - feature_throughputs     whether save feature throughputs
# - align                   whether time-aligning data or not
# - write_results           whether write results in to an ASCII file
# - alignment_dir           name of the alignment index folder
# - dimred_dir              (not in use)
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
subdir_residual_complexity <- function(fps = 120, dim_reduction = "pca", calculate_residuals = TRUE,
                            save_residuals = FALSE, save_pca = FALSE, normalize = FALSE,
                            residualdir = "residuals", pca_dir = "pca", remove_duplicates = TRUE,
                            feature_throughputs = FALSE, align = TRUE, write_results = FALSE,
                            alignment_dir = "alignment", dimred_dir = ".") {
    if (save_residuals) calculate_residuals <- TRUE
    
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
                rownames = append(rownames, sprintf("(%s, %s)", sequences[j], sequences[k]))
                rownames = append(rownames, sprintf("(%s, %s)", sequences[k], sequences[j]))
                M <- pair_residual_complexity(seqfile_a = sequences[j], seqfile_b = sequences[k],
                    seqnum_a = j, seqnum_b = k, fps = fps, dim_reduction = dim_reduction,
                    feature_throughputs = feature_throughputs, features = 0,
                    calculate_residuals = calculate_residuals, save_residuals = save_residuals,
                    save_pca = save_pca, normalize = normalize, residualdir = residualdir, pca_dir = pca_dir,
                    remove_duplicates = remove_duplicates, align = align, alignment_dir = alignment_dir)
                all_results[rows[1]:rows[2],] <- M$results
                rows <- rows + 2
                if (feature_throughputs) {
                    all_features <- M$features
                }
            }
        }
        rownames(all_results) <- rownames
        if (feature_throughputs)
            results[[i]] <- list(dir = subdirs[i], results = all_results, feature_results = all_features)
        else
            results[[i]] <- list(dir = subdirs[i], results = all_results)
        
        # set path back to maindir
        setwd(maindir)
        if (write_results) {
            filename <- "results.txt"
            # write dir name
            write(sprintf("Folder: %s", subdirs[i]), file = filename, append = TRUE, ncolumns=1)
            # write dir results
            result_frame <- results[[i]]$results
            write.table(result_frame, file=filename, sep='\t', append=TRUE, col.names=TRUE, row.names=TRUE)
        }
    }

    return(results)
}

# Calculate the residual complexity of all sequences cross-wise
# in the current directory (with the directory structure specified
# in the readme file). Sequences are considered to be repetitions of
# a same motion/gesture.
singledir_residual_complexity <- function(fps = 120, dim_reduction = "pca", calculate_residuals = TRUE,
                                save_residuals = FALSE, save_pca = FALSE, normalize = FALSE,
                                residualdir = "residuals", pca_dir = "pca",
                                remove_duplicates = TRUE, feature_throughputs = FALSE, align = TRUE,
                                alignment_dir = "alignment", dimred_dir = ".") {
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
            rownames = append(rownames, sprintf("(%s, %s)", sequences[j], sequences[k]))
            rownames = append(rownames, sprintf("(%s, %s)", sequences[k], sequences[j]))
            M <- pair_residual_complexity(seqfile_a = sequences[j], seqfile_b = sequences[k],
                seqnum_a = j, seqnum_b = k, fps = fps, dim_reduction = dim_reduction,
                feature_throughputs = feature_throughputs, features = 0,
                calculate_residuals = calculate_residuals, save_residuals = save_residuals,
                save_pca = save_pca, normalize = normalize, residualdir = residualdir,
                pca_dir = pca_dir, remove_duplicates = remove_duplicates, align = align,
                alignment_dir = alignment_dir)
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
pairdir_residual_complexity <- function(fps = 120, dim_reduction = "pca", calculate_residuals = TRUE,
                                save_residuals = FALSE, save_pca = FALSE, normalize = FALSE,
                                residualdir = "residuals", pca_dir = "pca", remove_duplicates = TRUE,
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
            M <- pair_residual_complexity(seqfile_a = sequences[j], seqfile_b = sequences[k],
                seqnum_a = j, seqnum_b = k, fps = fps, dim_reduction = dim_reduction,
                feature_throughputs = TRUE, features = 0, calculate_residuals = calculate_residuals,
                save_residuals = save_residuals, save_pca = save_pca, normalize = normalize,
                residualdir = residualdir, pca_dir = pca_dir, remove_duplicates = remove_duplicates,
                align = align, alignment_dir = alignment_dir)
            TP_results[j:k,] <- M$results
            all_features[[i]] <- M$features
        }
        else
            TP_results[j:k,] <- pair_residual_complexity(seqfile_a = sequences[j], seqfile_b = sequences[k],
                seqnum_a = j, seqnum_b = k, fps = fps, dim_reduction = dim_reduction,
                feature_throughputs = FALSE, features = 0, calculate_residuals = calculate_residuals,
                save_residuals = save_residuals, save_pca = save_pca, normalize = normalize,
                residualdir = residualdir, pca_dir = pca_dir, remove_duplicates = remove_duplicates,
                align = align, alignment_dir = alignment_dir)$results
    }
    
    if (feature_throughputs)
        all_results <- list(results = TP_results, feature_results = all_features)
    else 
        all_results <- list(results = TP_results)
    
    setwd(dir)
    print("Done")
    return(all_results)
}

# Calculate the residual complexity of two sequences numbered
# seqnum_a and seqnum_b. Throughputs per feature can also be printed.
#
# Parameters:
# - seqfile_a
# - seqfile_b
# - seqnum_a    number of first sequence
# - seqnum_b    number of second sequence
# - fps         data frames per second
# - method      dimension reduction method
# - print_feature_throughputs   prints throughputs per each feature
# - features    select certain features (all if 0)
# - save_residuals  saves the residuals
# - calculate_residuals calculate residuals; if false, use pre-calculated residuals 
pair_residual_complexity <- function(seqfile_a, seqfile_b, seqnum_a = 1, seqnum_b = 2, fps = 120,
                                     dim_reduction = "pca", feature_throughputs = FALSE, features = 0,
                                     calculate_residuals = TRUE, save_residuals = FALSE,
                                     save_pca = FALSE, normalize = FALSE, residualdir = "residuals",
                                     pca_dir = "pca", remove_duplicates = TRUE, align = TRUE,
                                     alignment_dir = "alignment") {
    print(sprintf("Evaluating throughput between %s and %s", seqfile_a, seqfile_b))
                                         
    pair_a_b <- new("residualComplexityEvaluator", seqfile_a = seqfile_a, seqfile_b = seqfile_b,
    seqnum_a = seqnum_a, seqnum_b = seqnum_b, fps = fps, dim_reduction = dim_reduction, features = features,
    save_residuals = save_residuals, save_pca = save_pca, calculate_residuals = calculate_residuals,
    normalize = normalize, residualdir = residualdir, pca_dir = pca_dir,
    remove_duplicates = remove_duplicates, align = align, alignment_dir = alignment_dir)
    
    results_a_b <- evaluate_complexity(pair_a_b)
    
    pair_b_a <- new("residualComplexityEvaluator", seqfile_a = seqfile_b, seqfile_b = seqfile_a,
    seqnum_a = seqnum_b, seqnum_b = seqnum_a, fps = fps, dim_reduction = dim_reduction, features = features,
    save_residuals = save_residuals, save_pca = save_pca, calculate_residuals = calculate_residuals,
    normalize = normalize, residualdir = residualdir, pca_dir = pca_dir,
    remove_duplicates = remove_duplicates, align = align, alignment_dir = alignment_dir)
    
    results_b_a <- evaluate_complexity(pair_b_a)

    #if (save_residuals) {
    #    # in data_handling.R
    #    write_residuals(residualdir, seqnum_a, results_a_b@seq_a)
    #    write_residuals(residualdir, seqnum_b, results_a_b@seq_b)
    #} # changed for clipdata

    resultvectors <- rbind(construct_result_vector(results_a_b),
           construct_result_vector(results_b_a), deparse.level = 0)
    
    if (feature_throughputs) {
        featurevectors <- rbind(get_feature_throughputs(results_a_b),
                        get_feature_throughputs(results_b_a),
                        deparse.level = 0)   
        resultlist <- list("results"=resultvectors, "features"=featurevectors, results_a_b,
                            results_b_a)
    }
    else
        resultlist <- list("results"=resultvectors, results_a_b, results_b_a)
    
    cat("Pair TP: \n")
    cat(sprintf("(%.3f, %.3f)", resultvectors[1,1], resultvectors[2,1]))
    cat("\n")
    
    return(resultlist)
    #return(rbind(construct_result_vector(results_a_b),
    #       construct_result_vector(results_b_a), deparse.level = 0))
}

# Create a union class for data files which can be either vectors or matrices
setClassUnion("data", c("vector","matrix"))
# A class for performing the required residual complexity evaluations,
# including the loading of files, calculation of residuals etc.
setClass("residualComplexityEvaluator",
         representation(seqfile_a = "character", # name of file A
                        seqfile_b = "character", # name of file B
                        seqnum_a = "numeric",	# number of file A
                        seqnum_b = "numeric",   # number of file B
                        fps     = "numeric",    # framerate
                        dim_reduction = "character", # choose reduction method
                        features = "numeric",	# choose features
						save_residuals = "logical",	# option to save residual data
						save_pca = "logical", # option to save pca-reduced data
						calculate_residuals = "logical", # option to calculate residuals
						normalize = "logical", # option to normalize sequences before residual calculation
						residualdir = "character", # choose the dir of residuals
						alignment_dir = "character", # choose the dir of alignment indeces
						pca_dir = "character", # choose the dir of pca-reduced sequences
						remove_duplicates = "logical", # choose whether remove duplicate frames
						align = "logical", # whether align residuals
                        res_a = "data",	# the original residuals a
						res_b = "data",	# the original residuals b
						seq_a   = "data",	# the processed residual sequence a
                        seq_b   = "data",	# the processed residual sequence b
                        pc_a    = "list",
                        pc_b    = "list",
                        kpc     = "kpca",
                        eig   = "NULL",
                        pcv_a   = "matrix",
                        pcv_b   = "matrix",
                        results = "list"))	# results after evaluation is done

# Evaluates complexity and returns a new residualComplexityEvaluator
# with the 'results' field set to the result list returned by the
# evaluation function. See construct_result_vector and get_feature_throughputs
# for result formatting and throughput calculation.
setGeneric("evaluate_complexity",
           function(this) standardGeneric("evaluate_complexity"))
setMethod("evaluate_complexity", "residualComplexityEvaluator",
           function(this) {
               this <- load_residuals(this)
               
               # check for dimensions
               if (any(dim(this@seq_a) != dim(this@seq_b)))
                    stop("Cut and/or aligned residual dimensions don't match!")
               
                if (any(is.na(this@seq_a)) || any(is.na(this@seq_b)))
                     stop("Residual sequences contain NaNs!")
               
               this <- do_reduction(this)
               this@results <- evaluate_residual_shared_information(this@seq_a, this@seq_b)
               this
           })

# Loads and aligns the residuals indicated by the two sequence numbers.
setGeneric("load_residuals",
           function (this) standardGeneric("load_residuals"))
setMethod("load_residuals", "residualComplexityEvaluator",
		function(this) {
		    if (this@align) {
			    data <- load_aligned_residuals(this@seqfile_a, this@seqfile_b,
			                                this@seqnum_a, this@seqnum_b, this@features,
			                                this@calculate_residuals, this@normalize,
			                                this@save_residuals, this@residualdir,
			                                this@remove_duplicates, this@alignment_dir)
		    }
		    else {
		        data <- load_unaligned_residuals(this@seqfile_a, this@seqfile_b,
			                                this@seqnum_a, this@seqnum_b, this@features,
			                                this@calculate_residuals, this@normalize,
			                                this@save_residuals, this@residualdir,
			                                this@remove_duplicates, this@alignment_dir)
		    }
			# original residuals
			this@res_a <- data[[1]]
            this@res_b <- data[[2]]
            # aligned residuals
            this@seq_a <- data[[3]]
            this@seq_b <- data[[4]]
            this
        })

# Performs PCA on the two loaded sequences.
setGeneric("do_reduction",
           function(this) standardGeneric("do_reduction"))
setMethod("do_reduction", "residualComplexityEvaluator",
          function(this) {
              if (this@dim_reduction == "pca") {
                  eigenvectors <- pca(this@seq_a)
                  this@seq_a <- normalize_features(this@seq_a) %*% eigenvectors
                  this@seq_b <- normalize_features(this@seq_b) %*% eigenvectors
              }
              else if (this@dim_reduction == "pca2") {
                  data_a <- normalize_features(this@seq_a)
                  data_b <- normalize_features(this@seq_b)
                  data <- rbind(data_a,data_b)
                  prc <- prcomp(data, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  this@seq_a <- prc$x[1:nrow(data_a),]
                  this@seq_b <- prc$x[(nrow(data_a)+1):(2*nrow(data_a)),]
              }
              else if (this@dim_reduction == "pca3") {
                  pc_a <- prcomp(this@seq_a, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  pc_b <- prcomp(this@seq_b, retx=TRUE, center=TRUE, scale=TRUE, tol=0.1)
                  this@seq_a <- pc_a$x
                  this@seq_b <- pc_b$x
              }
              else if (this@dim_reduction == "kpca") {
                  eigenvectors <- kernelpca(this@seq_a)
                  this@seq_a <- normalize_features(this@seq_a) %*% eigenvectors
                  this@seq_b <- normalize_features(this@seq_b) %*% eigenvectors
              }
              else if (this@dim_reduction == "kpca2") {
                  data_a <- normalize_features(this@seq_a)
                  data_b <- normalize_features(this@seq_b)
                  data <- rbind(data_a,data_b)
                  this@kpc <- kpca(data, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  # treshold won't work, why?
                  #kpc_b <- kpca(data_b, kernel="rbfdot", kpar = list(sigma = 1), th = 0.1)
                  this@eig <- this@kpc@eig
                  this@pcv_a <- this@kpc@pcv[1:nrow(data_a),]
                  this@pcv_b <- this@kpc@pcv[(nrow(data_a)+1):(2*nrow(data_a)),]
                  this@seq_a <- this@kpc@rotated[1:nrow(data_a),]
                  this@seq_b <- this@kpc@rotated[(nrow(data_a)+1):(2*nrow(data_a)),]
              }
              else if (this@dim_reduction == "kpca3") {
                  data_a <- normalize_features(this@seq_a)
                  data_b <- normalize_features(this@seq_b)
                  kpc_a <- kpca(data_a, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  kpc_b <- kpca(data_a, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  this@seq_a <- kpc_a@rotated
                  this@seq_b <- kpc_b@rotated
              }
              else if (this@dim_reduction == "kfa") {
                  data_a <- this@seq_a
                  data_b <- this@seq_b
                  data <- rbind(data_a,data_b)
                  princ <- kfa(data, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  this@seq_a <- princ@xmatrix[1:nrow(data_a),]
                  this@seq_b <- princ@xmatrix[(nrow(data_a)+1):(2*nrow(data_a)),]
              }
              else if (this@dim_reduction == "kfa2") {
                  data_a <- normalize_features(this@seq_a)
                  data_b <- normalize_features(this@seq_b)
                  princ_a <- kfa(data_a, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  princ_b <- kfa(data_b, kernel="rbfdot", kpar = list(sigma = 0.01),na.action=na.omit)
                  this@seq_a <- princ_a@xmatrix
                  this@seq_b <- princ_b@xmatrix
              }
              
              if (this@save_pca) {
                  write_pca(this@pca_dir, this@seqfile_a, this@seq_a)
                  write_pca(this@pca_dir, this@seqfile_b, this@seq_b)
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

# Given a shared information value, calculate throughput.
setGeneric("calculate_throughput",
           function(this, shared) standardGeneric("calculate_throughput"))
setMethod("calculate_throughput", "residualComplexityEvaluator",
          function(this, shared) {
              shared / NROW(this@seq_a) * this@fps / log(2.0)
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
