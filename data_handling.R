# Functions for handling the data.

# The MIT License (MIT)
# Copyright (c) 2013 Arttu Modig
# Version 1.02

source("infocapacity.R")
data_handling_debug <- FALSE

# Loads the alignment files corresponding to the given
# sequence numbers and calculates residuals.
# Returns a list with:
# 1) residuals_orig_a  first original residuals
# 2) residuals_orig_b  second original residuals
# 3) residuals_alig_a  first aligned residuals
# 4) residuals_alig_b  second aligned residuals
load_aligned_residuals <- function(seqfile_a, seqfile_b, seqnum_a, seqnum_b, features,
                                    calculate_residuals, normalize, remove_duplicates,
                                   alignment_dir, residual_dir) {
    # get residuals
    orig_a <- get_residuals(seqfile_a, features, calculate_residuals, normalize, residual_dir)
    orig_b <- get_residuals(seqfile_b, features, calculate_residuals, normalize, residual_dir)
    
    sequence_a <- orig_a$sequence
    sequence_b <- orig_b$sequence
    
    residuals_orig_a <- orig_a$residuals
    residuals_orig_b <- orig_b$residuals
    
    prediction_a <- orig_a$prediction
    prediction_b <- orig_b$prediction
    
    # get alignment
    print(sprintf("Load alignment indeces from folder: %s", alignment_dir))            
    ali_a <- load_aligned(seqnum_a, seqnum_b, alignment_dir)
    ali_b <- load_aligned(seqnum_b, seqnum_a, alignment_dir)
    
    # warning
    if (length(ali_a)!=length(ali_b))
        stop("STOP: Alignment vectors have different length!")
        
    # remove duplicate frames
    if (remove_duplicates) {
        new_alis <- remove_duplicate_frames(ali_a, ali_b)
        ali_a <- new_alis[[1]]
        ali_b <- new_alis[[2]]
    }
    
    # DIRTY FIX!
    # remove NaNs that the method above creates
    ali_a <- ali_a[!is.na(ali_a)]
    ali_b <- ali_b[!is.na(ali_b)]
        
    #### align residuals ####
    
    #residuals_a <- align_residuals(ali_a, residuals_orig_a)
    #residuals_b <- align_residuals(ali_b, residuals_orig_b)
    new_residuals <- align_residuals2(ali_a, ali_b, residuals_orig_a, residuals_orig_b)
    residuals_a <- new_residuals[[1]]
    residuals_b <- new_residuals[[2]]
    
    # warning                    
    if (any(is.na(residuals_a[[1]])) || any(is.na(residuals_b[[1]])))
        stop("load_aligned_residuals: 1. There is NaNs!")
	
	  # cut aligned residuals to equal length
    residuals_alig <- cut_to_equal_length(residuals_a, residuals_b)
    residuals_alig_a <- residuals_alig[[1]]
    residuals_alig_b <- residuals_alig[[2]]
    
    # warning
    if (any(is.na(residuals_alig_a[[1]])) || any(is.na(residuals_alig_b[[1]]) ||
            is.na(residuals_alig_a[[2]])) || any(is.na(residuals_alig_b[[2]])))
        stop("load_aligned_residuals: 2. There is NaNs!")
        
    # check for dimensions
    if (any(dim(residuals_alig_a) != dim(residuals_alig_b)))
        stop("load_aligned_residuals: Cut and/or aligned residual dimensions don't match!")

    return(list("sequence_a"=sequence_a,
                "sequence_b"=sequence_b,
                "residuals_orig_a"=residuals_orig_a,
                "residuals_orig_b"=residuals_orig_b,
                "residuals_alig_a"=residuals_alig_a,
                "residuals_alig_b"=residuals_alig_b,
                "prediction_a"=prediction_a,
                "prediction_b"=prediction_b))
}

# Loads sequence files and calculates unaligned residuals.
# Returns a list with:
# 1) residuals_orig_a  first original residuals
# 2) residuals_orig_b  second original residuals
# 3) residuals_alig_a  first aligned residuals
# 4) residuals_alig_b  second aligned residuals
load_unaligned_residuals <- function(seqfile_a, seqfile_b, seqnum_a, seqnum_b, features,
                                    calculate_residuals, normalize, remove_duplicates,
                                    residual_dir) {
    # get residuals
    orig_a <- get_residuals(seqfile_a, features, calculate_residuals, normalize, residual_dir)
    orig_b <- get_residuals(seqfile_b, features, calculate_residuals, normalize, residual_dir)
    
    sequence_a <- orig_a$sequence
    sequence_b <- orig_b$sequence
    
    residuals_orig_a <- orig_a$residuals
    residuals_orig_b <- orig_b$residuals
    
    prediction_a <- orig_a$prediction
    prediction_b <- orig_b$prediction
    
    # cut unaligned residuals to equal length
    residuals <- cut_to_equal_length(residuals_orig_a, residuals_orig_b)
    residuals_a <- residuals[[1]]
    residuals_b <- residuals[[2]]
               
    # check for dimensions
    if (any(dim(residuals_a) != dim(residuals_b)))
        stop("load_unaligned_residuals: Cut and/or aligned residual dimensions don't match!")
        
    return(list("sequence_a"=sequence_a,
                "sequence_b"=sequence_b,
                "residuals_orig_a"=residuals_orig_a,
                "residuals_orig_b"=residuals_orig_b,
                "residuals_alig_a"=residuals_a,
                "residuals_alig_b"=residuals_b,
                "prediction_a"=prediction_a,
                "prediction_b"=prediction_b))
}

# Returns the residuals for a given sequence, aligned according
# to a pre-made residual alignment.
get_residuals <- function(sequence_file, features, calculate_residuals, normalize, residual_dir) {
    # empty prediction array
    prediction <- array()
  
    # load sequence
    sequence <- load_sequence(sequence_file, features)
    # calculate residuals
    if (calculate_residuals) {
        model <- model(sequence) # See in: infocapacity.R
        residuals <- model$residuals
        prediction <- model$prediction
    }
    # optionally use pre-calculated residuals
    else residuals <- load_residual(sequence_file, features, residual_dir)
    
    # error
    if (any(is.na(residuals)))
        stop("get_residuals: There is NaNs!")
    
    # normalize residual vectors
    if (normalize) residuals <- normalize_features(residuals)	  
        
    return(list("sequence"=data.matrix(sequence), "residuals"=residuals, "prediction"=prediction))
}

# Loads and returns the sequence denoted by the given
# sequence number (assumed to be in the current working directory)
load_sequence <- function(sequence_file, features) {
    print(sprintf("Reading sequence %s", sequence_file))
    seq <- read.table(sprintf("%s", sequence_file), sep="\t")
    if (features[1]) {
        print("Choosing features:")
        print(features)
        seq <- seq[,features,drop=FALSE]
    }
    # remove zero feature columns
    seq <- seq[, colSums(abs(seq)) != 0]
}

# Residual files must have the same name as parent files
load_residual <- function(sequence_file, features, residual_dir) {
    print(sprintf("Reading residual %s/%s", residual_dir, sequence_file))
    res <- read.table(sprintf("%s/%s", residual_dir, sequence_file), sep="\t")
    if (features[1]) { #  if  first element exists i.e. any elements exist
        print("Choosing features:") #  pick only the features indexed in vector features
        print(features)
        res <- res[,features]
    }
    return(data.matrix(res))
}

load_model <- function(sequence_file, dir) {
    print(sprintf("Reading residual %s", sequence_file))
    data <- read.table(sprintf("./%s/%s", dir, sequence_file), sep="\t")
    return(data.matrix(data))
}

# Removes duplicated rows from alignment index sequences a and b.
#
# Returns the altered sequences as a list with two elements.
remove_duplicate_frames <- function(a, b) {
    skip <- c(FALSE, a[2:NROW(a),]-a[1:(NROW(a)-1),] == 0)

    return(list(a[!skip,], b[!skip,]))
}

# Loads the aligned data file "seqnum1_ali_seqnum2.txt" or the corresponding
# alignment file as a data frame and returns it. (Note: the alignment file
# has only one column, so it is practically a vector.)
#
# Parameters:
# seqnum_a       number of the first sequence
# seqnum_b       number of the second sequence
# folder         alignment index folder (default: "alignment")
load_aligned <- function(seqnum_a, seqnum_b, folder) {
    print(sprintf("Reading alignment vector: %s/%d_ali_%d.txt", folder, seqnum_a, seqnum_b))
    ali <- read.table(sprintf("%s/%d_ali_%d.txt", folder, seqnum_a, seqnum_b), sep="\t")
    
    if (data_handling_debug)
        print(sprintf("Vector length: %d", NROW(ali)))
    
    return(ali)
}

# Cuts two sequences to equal length by removing frames from the
# beginning or the end.
#
# Parameters:
# - a, b    residual sequences
# - end     choose to cut from end
cut_to_equal_length <- function(a, b, end = TRUE) {
    diff <- NROW(a) - NROW(b)
    if (end) {
        if (diff > 0)
            a <- a[1:(NROW(a)-diff),]
        else if (diff < 0)
            b <- b[1:(NROW(b)-abs(diff)),]
        return(list(a, b))
    }
    if (diff < 0)
        b <- b[(abs(diff)+1):NROW(b),]
    else if (diff > 0)
        a <- a[(diff+1):NROW(a),]    
    return(list(a, b))
}

# Deprecated! Please see "align_residuals2"
# Aligns the given residuals according to the given pre-aligned coordinate
# sequence and returns it in data frame format.
# Parameters:
# - alignment   the alignment file, indicating frame duplications
# - residuals   residuals for the original (non-aligned) sequence
align_residuals <- function(alignment, residuals) {
    # errors (remember that residuals are 2 frames shorter than original data)
    if (tail(alignment, n=1) > length(residuals)+2)
        stop("align_residuals: alignment indeces outmatch data length!")
    # The first two frames (numbered 1 and 2) are cut out because
    # no residuals are calculated for those frames.
    alignment <- alignment[alignment > 2] - 2
    aligned_residuals <- residuals[alignment,]
    
    return(aligned_residuals)
}

# Aligns the given residuals according to the given pre-aligned coordinate
# sequence and returns it in data frame format.
# Parameters:
# - ali_a           the first alignment file, indicating frame duplications
# - ali_b           the second alignment file, indicating frame duplications
# - residuals_a     residuals for the first original (non-aligned) sequence
# - residuals_b     residuals for the second original (non-aligned) sequence
align_residuals2 <- function(ali_a, ali_b, residuals_a, residuals_b) {
    # errors (remember that residuals are 2 frames shorter than original data)
    if ((tail(ali_a, n=1) > length(residuals_a)+2) || (tail(ali_b, n=1) > length(residuals_b)+2))
        stop("align_residuals2: Alignment indeces outmatch data length!")
    
    # debug
    if (data_handling_debug) {
        print(sprintf("alignment A length: %d", length(ali_a)))
        print(sprintf("alignment B length: %d", length(ali_b)))
        print(sprintf("ali A tail: %d", tail(ali_a, n=1)))
        print(sprintf("ali B tail: %d", tail(ali_b, n=1)))
        print(sprintf("residuals A frames: %d", NROW(residuals_a)))
        print(sprintf("residuals B frames: %d", NROW(residuals_b)))
    }
    
    # length of the frames (1,2) block at the start
    len_a <- length(ali_a[ali_a < 3])
    len_b <- length(ali_b[ali_b < 3])
    
    if (len_a > len_b) {
        ali_a <- ali_a[-(1:len_a)] - 2
        ali_b <- ali_b[-(1:len_a)] - 2
    } else {
        ali_a <- ali_a[-(1:len_b)] - 2
        ali_b <- ali_b[-(1:len_b)] - 2
    }
    
    if (data_handling_debug) {
        print("ALIS CUT")
        print(sprintf("alignment A length: %d", length(ali_a)))
        print(sprintf("alignment B length: %d", length(ali_b)))
        print(sprintf("ali A tail: %d", tail(ali_a, n=1)))
        print(sprintf("ali B tail: %d", tail(ali_b, n=1)))
    }
    
    aligned_residuals_a <- residuals_a[ali_a,]
    aligned_residuals_b <- residuals_b[ali_b,]
    
    return(list(aligned_residuals_a, aligned_residuals_b))
}

# Writes data in a given folder.
# Parameters:
# - dir             the folder for residuals
# - sequence_file   sequence filename
# - data            residual data
write_data <- function(dir, sequence_file, data) {
	if (!file.exists(dir)) dir.create(dir)
  print(sprintf("Writing data: %s/%s ...", dir, sequence_file))
	write.table(data, sprintf("%s/%s", dir, sequence_file),
	sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Lists all directories in the current working path.
# Based on function by Joshua Ulrich (2011)
dirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=FALSE,
    ignore.case=FALSE) {
    all <- list.files(path, pattern, all.dirs,
           full.names, recursive=FALSE, ignore.case)
    # leave only dirs
    all <- all[file.info(all)$isdir]
    # remove stupid dirs
    remove1 <- list.files(path, pattern="__MACOSX", all.dirs,
            full.names, recursive=FALSE, ignore.case)
    remove2 <- list.files(path, pattern="alignment", all.dirs,
                      full.names, recursive=FALSE, ignore.case)
    all <- all[!all %in% remove1]
    all <- all[!all %in% remove2]
}
