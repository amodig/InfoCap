# Algorithm for calulating mutual information based on k-nn search.
# MI value in nats.
# Based on Kraskov et. al.
# Santeri Räisänen, 2014
#
#
#
#

library("iterators",lib.loc="~/R/Rpackages");
library("foreach",lib.loc="~/R/Rpackages");
library("doParallel",lib.loc="~/R/Rpackages");

construct_sp_tree <- function(data, METRIC, b_size)  {      
  sp_tree <- list(data)
  i <- 1;
  
  while(1==1) {
    p <- part(sp_tree[[i]],METRIC);
    
    p <- list(p);
    if (is.matrix(data)) {
      
      pr <- projections(p[[1]],sp_tree[[i]]);
      
      SU <- summary(pr);
      sp_tree[[2*i]] <- rbind(sp_tree[[i]][which(pr>=SU[3]),]);
      sp_tree[[(2*i)+1]] <- rbind(sp_tree[[i]][which(pr<SU[3]),]);
      
      p[[2]] <- SU[3];
      p[[3]] <- SU[2];
      p[[4]] <- SU[5];
      
      
      sp_tree[[i]] <- p;
      
      if ((log2(2*(i+1)))==floor(log2(2*(i+1))) && (NROW(sp_tree[[2*i]])<=b_size)){
        
        return(sp_tree);;                                                                      
      }
      i<-i+1;
    }
  }
}


part <- function(data, METRIC) {
  
  
  
  A <- euclid_d(data,data[1,]);
  a <- data[which.max(A),];
  A <- euclid_d(data,a);
  b <- data[which.max(A),];
  
  a <- a-b;
  
  
  
  return(a);       #Retuns vector which defines the partition of the set of vectors
}



projections <- function(pr_v, data) {
  k <- sqrt(sum(pr_v*pr_v));
  if (is.matrix(data)) {
    return(apply(t(t(data)*pr_v),1,sum)/k);
  }
  return(sum(data*pr_v)/k);
}

evaluate_mi <- function(seq1,seq2,k=2) {                                
  d_size <- NCOL(seq1);                                                  
  data <- cbind(seq1,seq2);
  radii = get_all_knn(data,k,d_size);
  print("all neighbours found!")
  y_points <- count_points(radii[,1],seq1);
  x_points <- count_points(radii[,2],seq2)
  
  return(-1*(mean(digamma(x_points)+digamma(y_points))) - (1/(k-1)) + digamma(k-1) + digamma(NROW(seq1)));
}


get_all_knn <- function(data,k, partition) {                 #this calulates the k-nearest neighbour for each row in matrix 
  
  m_part <<- partition;
  
  sp_tree <- construct_sp_tree(data,max_d,NROW(data)/25);                   #with regards to the metric-function given.
  neighbours <- c();
  registerDoParallel();
  limit <- log2(length(sp_tree)+1)-1
  radii <- foreach(i=1:NROW(data), .combine='rbind') %dopar% {
    level <- 1;
    a <- recursive_search(sp_tree,data[i,],k,max_d,level,limit)[,2:(NCOL(data)+1)];
    d_x <- euclid_d(a[,1:m_part],data[i,1:m_part]);
    d_y <- euclid_d(a[,(m_part+1):NCOL(data)],data[i,(m_part+1):NCOL(data)]);
    return(c(max(d_x),max(d_y)));
  }
  
  return(radii);
}



recursive_search <- function(sp_tree,Vector,k,METRIC,level,limit)  {
  if (level>=(2^limit)){
    a <- METRIC(sp_tree[[level]],Vector);
    
    ord <- order(a);
    copy <- sp_tree[[level]][ord,];
    
    return(cbind(a[ord[1:k]],copy[1:k,]));
  }
  
  r = projections(sp_tree[[level]][[1]],Vector)
  a <- c(-1);
  b <- c(-1);
  
  if (r>sp_tree[[level]][[3]]) {
    a <- recursive_search(sp_tree,Vector,k,METRIC,2*level,limit);
  }
  if (r<sp_tree[[level]][[4]]) {
    b <- recursive_search(sp_tree,Vector,k,METRIC,2*level+1,limit);
  }
  
  if (!is.matrix(a)) {
    return(b);
  }
  if (!is.matrix(b)) {
    return(a);
  }
  else {
    ret <- rbind(a,b);
    ret <- ret[order(ret[,1]),];
    ret <- ret[1:k,];
    return(ret);
  }
}



max_d <- function(a, b){
  l <- m_part;                   
  r <- NCOL(a);
  
  n1 <- euclid_d(a[,1:l],b[1:l]);
  n2 <- euclid_d(a[,(l+1):r],b[(l+1):r]);
  
  return(apply(cbind(n1,n2),1,max));
}




euclid_d <- function(a, b) {        
  if (NROW(a)==NROW(b)&&!is.matrix(a)&&!is.matrix(b)) {
    
    return(sqrt(sum((a-b)^2)));
  }
  else {
    
    return(sqrt(apply(t(t(a)-b)^2,1,sum)));
  }
}


recursive_count <- function(sp_tree,Vector,radius,METRIC,level,limit)  {
  if (level>=(2^limit)){
    a <- METRIC(sp_tree[[level]],Vector);
    return(NROW(which(a<=radius)));
  }
  
  r = projections(sp_tree[[level]][[1]],Vector)
  a <- 0;
  b <- 0;
  if (r>sp_tree[[level]][[3]]) {
    a <- recursive_count(sp_tree,Vector,radius,METRIC,2*level,limit);
  }
  if (r<sp_tree[[level]][[4]]) {
    b <- recursive_count(sp_tree,Vector,radius,METRIC,2*level+1,limit);
  }
  
  return(a+b);
}



count_points <- function(radii, data) {
  if (!is.matrix(data)) {
    return(count_points_s(radii,data));
  } 
  
  ret <- c();
  
  
  sp_tree <- construct_sp_tree(data,euclid_d,NROW(data)/25); 
  
  limit <- log2(length(sp_tree)+1)-1
  ret <- foreach (j = 1:NROW(data),.combine="rbind") %dopar% {
    level<-1
    recursive_count(sp_tree,data[j,],radii[j],euclid_d,level,limit);
  }
  return(ret-1);
}


count_points_s <- function(radii, data) {
  ret <- matrix(nrow = NROW(data), ncol = 1);
  n_data <- data[order(data)];
  for (i in 1:NROW(data)) {
    counter <- 0;
    ind <- which(n_data==data[i]);
    radius = radii[i];
    while (ind<=NROW(data)&&euclid_d(n_data[ind],data[i])<=radius) {
      counter = counter + 1;
      ind = ind + 1;
    }
    ind <- which(n_data==data[i])-1;
    while (ind > 0 && euclid_d(n_data[ind],data[i])<=radius) {
      counter <- counter + 1;
      ind = ind - 1;
    }
    ret[i] <- counter;
  }
  return(ret-1);
}
