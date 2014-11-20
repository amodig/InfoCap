# Algorithms for evaluating movement throughput

# The MIT License (MIT)
# Copyright (c) 2013 Arttu Modig
# Version 1.1

# parameters
step <- 2
pca_variance_ratio <- 0.9

# Return the AR(2) model as a list of:
# 1) residuals: each element is the residual vector for the corresponding feature,
# given the observed sequence.
# 2) prediction: each element is the predicted vector for the corresponding feature,
# given the observed sequence.
model <- function(sequence) {
    n <- NROW(sequence) - step
    if (NCOL(sequence)==1) {
        sequence_predictors <- cbind(sequence[1:n], sequence[1:n+step-1])
        observed_sequence <- sequence[1:n+step]
    } else {
        sequence_predictors <- cbind(sequence[1:n,], sequence[1:n+step-1,])    
        observed_sequence <- sequence[1:n+step,]
    }

    frames <- NROW(observed_sequence)
    num_features <- NCOL(observed_sequence)
    residuals = matrix(nrow = frames, ncol = num_features)
    prediction = matrix(nrow = frames, ncol = num_features)
    
    for (feature_idx in 1:num_features) {
        QR <- get_qr_for_feature(feature_idx, sequence_predictors, observed_sequence)
        residuals[, feature_idx] <- get_residuals_for_feature(QR, observed_sequence,
                                                              feature_idx)
        prediction[, feature_idx] <- get_prediction_for_feature(QR, observed_sequence,
                                                                feature_idx)
    }
    return(list("residuals"=residuals, "prediction"=prediction))
}

get_qr_for_feature <- function(feat, predictors, observed) {
    frames <- NROW(observed)
    # extract regressors depending on input size
    if (NCOL(predictors) == NCOL(observed)) {
        QR <- qr(cbind(matrix(1,frames,1), predictors[, feat]))
    } else {
        QR <- qr(cbind(matrix(1,frames,1), predictors[, feat],
            predictors[, NCOL(observed) + feat]))
    }
}

# Calculates the residuals for a given feature in the regression
# model.
#
# Parameters:
# - QR              the QR decomposition for a sequence
# - observed        the observed sequence
# - feature_idx     feature index
get_residuals_for_feature <- function(QR, observed, feature_idx) {
    if (NCOL(observed)==1) {# needs special handling in R
        return(qr.resid(QR, observed))
    } else {
        return(qr.resid(QR, observed[, feature_idx]))
    }   
}

# Calculates the prediction for a given feature in the regression
# model.
#
# Parameters:
# - QR              the QR decomposition for a sequence
# - observed        the observed sequence
# - feature_idx     feature index
get_prediction_for_feature <- function(QR, observed, feature_idx) {
    if (NCOL(observed)==1) {# needs special handling in R
        return(qr.fitted(QR, observed))
    } else {
        return(qr.fitted(QR, observed[, feature_idx]))
    }   
}

# Evaluates shared information according to residuals of residuals.
#
# Parameters:
# - residuals_a, residuals_b    residual sequences (data frames or matrices)
evaluate_residual_shared_information <- function(residuals_a, residuals_b) {
    total_shared <- 0
    total_RSS <- 0
    total_RSS_residual <- 0
    feature_shared <- array(0, NCOL(residuals_a))
    n <- NROW(residuals_a)

    for (i in 1:NCOL(residuals_a)) { # calculate shared information by feature
        if (NCOL(residuals_a)==1) {
            RSS <- sum(residuals_a^2)
        } else {
            RSS <- sum(residuals_a[, i]^2)
        }
        total_RSS <- total_RSS + RSS
        
        if (NCOL(residuals_a)==1) {
            RSS_residual <- sum_squared_residuals(residuals_a, residuals_b)
        } else {
            RSS_residual <- sum_squared_residuals(residuals_a[, i], residuals_b[, i])
        }
        total_RSS_residual <- total_RSS_residual + RSS_residual

        feature_shared[i] <- (n/2*log(RSS/RSS_residual) - 1/2) # fixed! (24.4.2013)
        #feature_shared[i] <- (n/2*log(RSS/RSS_residual) - log(n)/2) # old version
        
        # find possible error
        if (is.na(feature_shared[i])) {
            stop("Feature shared NA!")
        }
        
        total_shared <- total_shared + feature_shared[i]
    }

    return(list(total_shared = total_shared, feature_shared = feature_shared,
        RSS = total_RSS, RSS_conditional = total_RSS_residual))
}

# THIS CALCULATES MI ACCORDING TO THE SP_KNN method  
evaluate_residual_mi <- function(residuals_a, residuals_b, k, feature_mi = TRUE) {
  if (is.null(k)) {
    k <- 2
  }
  
  # for single feature vectors
  if (!is.matrix(residuals_a)) {
    m <- cbind(residuals_a,residuals_b);                  
    
    neighbours <- get_all_knn(m,k);
      
    y_points <- count_points(neighbours[,1],m[,1]);
    x_points <- count_points(neighbours[,2],m[,2]);
    
    # calculate mutual information
    MI <- ( -1*(mean(digamma(x_points)) + mean(digamma(y_points))) - 1/(k-1)
            + digamma(k-1) + digamma(NROW(m)) );  
    
    return (list(total_shared = MI, feature_shared = 1, RSS = 1, RSS_conditional = 1));
  }
  
  # this is for calculating without separating features, assuming sequences are
  # of the same dimensionality
  m <- cbind(residuals_a,residuals_b);  
  neighbours <- get_all_knn(m,k);
  l = NCOL(residuals_a);
  y_points <- count_points(neighbours[,1:l],m[,1:l]);
  x_points <- count_points(neighbours[,(l+1):(NCOL(m))],m[,(l+1):(NCOL(m))]);
  MI = ( -1*(mean(digamma(x_points)) + mean(digamma(y_points))) - 1/(k-1) +
         digamma(k-1) + digamma(NROW(m)) );  
  f_shared = 1;
  
  if (feature_mi && is.matrix(residuals_a)) { 
    registerDoParallel();
    data <- foreach( i = 0:((NCOL(residuals_a)/3)-1), .combine='cbind') %dopar%  {  
      lim_l <-  i*3+1;
      lim_r <- (i+1)*3;
      m <- cbind(residuals_a[,lim_l:lim_r],residuals_b);                        
      neighbours <- get_all_knn(m,k,3);
      
      y_points <- count_points(neighbours[,1:3],m[,1:3]);
      x_points <- count_points(neighbours[,4:NCOL(m)],m[,4:NCOL(m)]);
      
      -1*(digamma(x_points)+digamma(y_points)) - 1/(k-1) +
        digamma(k-1) + digamma(NROW(m));  
    }
    f_shared = apply(data,2,mean);
    #browser();
    #write.table(data,file=sprintf("%s/feature_mi_sequence.txt",getwd()));
  }
  
  return(list(total_shared = MI, feature_shared = f_shared,
              RSS = 1, RSS_conditional = 1));  
}

# Fit a linear model seq_a ~ seq_b and return the sum of squared
# residuals from the model.
sum_squared_residuals <- function(seq_a, seq_b) {
    lmfit = lm(seq_a ~ seq_b)
    residuals = resid(lmfit)
    return(sum(residuals^2))
}

# Performs PCA on data and returns a matrix where each column is an eigenvector
# The number of eigenvectors is chosen so that 90% of variance
# is covered.
pca <- function(data) {
    data <- normalize_features(data)

    principal_components <- prcomp(data)
    num_eigenvecs <- choose_number_of_eigenvectors(principal_components$sdev, pca_variance_ratio)

    eigenvectors <- principal_components$rotation[,1:num_eigenvecs]
    return(eigenvectors)
}

# Given a list of standard deviations (from the prcomp function
# output), choose a number of eigenvectors in such a way that
# e.g. 90% (ratio = 0.9) of variance is covered.
choose_number_of_eigenvectors <- function(sdevs, ratio) {
    treshold <- ratio * sum(sdevs**2)
    sum <- 0
    for (k in 1:length(sdevs)) {
        sum <- sum + sdevs[k]**2
        if (sum >= treshold)
            return(k)
    }
}

# Performs normalization on the features of sequence a.
# Returns the altered sequence.
normalize_features <- function(a) {
    # apply over matrix margins
    a <- t(apply(a,1,'-',apply(a,2,mean))) # zero mean
    a <- t(apply(a,1,'/',apply(a,2,var))) # unit variance
    a[is.nan(a)] <- 0 # NaNs to zero
    a
}
