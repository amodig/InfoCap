# Algorithm for calulating k-nearest neighbors using sp-trees.
# Based on Kraskov et. al.
# Santeri Räisänen, 2014
#
#
#
#


construct_sp_tree <- function(data, METRIC = max_d, size = 800)  {      
  sp_tree <- list(data)
  i <- 1;
  
  while(1==1) {
    p <- part(sp_tree[[i]],METRIC);
    
    p <- list(p[1,],p[2,]);
    if (is.matrix(data)) {
      pr <- projections(p,sp_tree[[i]]);
      
      SU <- summary(pr);
      sp_tree[[2*i]] <- rbind(sp_tree[[i]][which(pr>SU[2]),]);
      sp_tree[[(2*i)+1]] <- rbind(sp_tree[[i]][which(pr<SU[5]),]);
      
      
      p[[3]] <- SU[3];
  
      sp_tree[[i]] <- p;
    
      if ((log2(2*(i+1)))==floor(log2(2*(i+1))) && (NROW(sp_tree[[2*i]])<=size)){  ## the last inequality controls how large the final bins are
      break;                                                                      ## larger value means slower algorithm but more accurate results.
    }
    i<-i+1;
  }
  
  return(sp_tree);
}
}


part <- function(data, METRIC) {
  if (is.matrix(data)) {  
    m_data <- apply(data,2,mean);
  }
  else {
    m_data <- mean(data);
  }
  
  n_data <- scale(data, scale = FALSE);
  
  
  A <- METRIC(n_data,rep(0,NCOL(data)));
  a <- n_data[which.max(A),];
  
  
  return(rbind(a,m_data));       #Retuns vector which defines the partition of the set of vectors
}



projections <- function(p, data) {
  pr_v <- p[[1]]
  if (is.matrix(data)) {
    d <- t(t(data) - p[[2]])
  }
  else {
    d <-data-p[[2]];    
  }
  
  k <- sqrt(sum(pr_v*pr_v));
  if (is.matrix(data)) {
    return(apply(t(t(d)*pr_v),1,sum)/k);
  }
  return(sum(d*pr_v)/k);
}

evaluate_mi <- function(seq1,seq2,k=2) {                                # parallellization might cause problems
  d_size <- NCOL(seq1);                                                 # Calculate for all K at the same time!?
  data <- cbind(seq1,seq2);
  neighbours = get_all_knn(data,k,d_size);
  y_points <- count_points(seq1,neighbours[,1:d_size]);
  x_points <- count_points(seq2,neighbours[,(d_size+1):NCOL(data)])
  
  return(-1*(mean(digamma(x_points))+mean(digamma(y_points))) - 1/(k-1) + digamma(k-1) + digamma(NROW(data)));
}



get_all_knn <- function(data,k, partition,metric = max_d) {                 #this calulates the k-nearest neighbour for each row in matrix 
  
  m_part <<- partition;
  
  t <- construct_sp_tree(data,max_d);                   #with regards to the metric-function given.
  neighbors <- matrix(0,nrow=NROW(data), ncol=NCOL(data));
  
  for (i in 1:NROW(data)) {
    neighbors[i,] <- get_knn(t,data[i,],k,max_d);
  }
  return(neighbors);
}



get_knn <- function(sp_tree,Vector,k,METRIC) {                       # get k-nearest neighbor for a single vector
  level <- 1;
  limit <- log2(length(sp_tree)+1)-1
  for (i in 1:limit) {
    r = projections(sp_tree[[level]],Vector)
    if (!is.logical(r>sp_tree[[level]][[3]])) {
      browser();
    }
    if (r>sp_tree[[level]][[3]]) {
      level <- 2*level;
    }
    else {
      level <- 2*level+1;
    }
  }
  
  a <- METRIC(sp_tree[[level]],Vector);
  if (min(a)!=0) {
    return(NULL);
  }
  copy <- sp_tree[[level]][order(a),];
  return(copy[k,]);
}




max_d <- function(a, b){         #metric for the R^6 vector space, but this could be any metric which fits the vector-space in question
  l <- m_part;
  r <- NCOL(a);
  if (is.matrix(a)) {
    if (is.matrix(b)) {
      
      n1 <- sqrt(apply(t(t(a[,1:l])-b[,1:l])^2,1,sum));
      n2 <- sqrt(apply(t(t(a[,(l+1):r])-b[,(l+1):r])^2,1,sum));
      m <- cbind(n1,n2);
     
      return(apply(m,1,max));
    }
    else {
      n1 <- sqrt(apply(t(t(a[,1:l])-b[1:l])^2,1,sum));
      
      n2 <- sqrt(apply(t(t(a[,(l+1):r])-b[(l+1):r])^2,1,sum));
      m <- cbind(n1,n2);
      
      return(apply(m,1,max));
    }
  }
  else {
    r <- NROW(b);
    n1 <- sqrt(sum(t(t(a[1:l])-b[1:l])^2));
    n2 <- sqrt(sum(t(t(a[(l+1):r])-b[(l+1):r])^2));
    m <- cbind(n1,n2);
    return(max(m));
  }
}




euclid_d <- function(a, b) {        
  if (is.matrix(a)) {
    
    return(sqrt(apply(t(t(a)-b)^2,1,sum)));
  }
  else {
    
    return(sqrt(sum((a-b)^2)));
  }
}





count_points <- function(neighbours, data) {
  if (!is.matrix(neighbours)) {
    return(count_points_s(neighbours,data));
  } 
  
  ret <- matrix(nrow = NROW(data), ncol = 1);
  
  
  sp_tree <- construct_sp_tree(data,euclid_d,1500); 
  
  limit <- log2(length(sp_tree)+1)-1
  for (j in 1:NROW(data)) {
    level <- 1;
    for (i in 1:limit) {
      
      r <- projections(sp_tree[[level]],data[j,]);
      
      
      if (r>sp_tree[[level]][[3]]) {
        level <- 2*level;
      }
      else {
        level <- (2*level)+1;
      }
      
    }
    
    a <- euclid_d(sp_tree[[level]],data[j,]);
    radius <- euclid_d(data[j,],neighbours[j,]);
    
    
    ret[j] = NROW(which(a<=radius));
    if (ret[j]==0) {
      browser();
    }
    
  }
  return(ret)
}


count_points_s <- function(neighbours, data) {
  ret <- matrix(nrow = NROW(data), ncol = 1);
  n_data <- data[order(data)];
  for (i in 1:NROW(data)) {
    counter <- 0;
     ind <- which(n_data==data[i]);
     radius = abs(data[i]-neighbours[i]);
     while (ind<=NROW(data)&&abs(n_data[ind]-data[i])<radius) {
       counter = counter + 1;
       ind = ind + 1;
     }
    ind <- which(n_data==data[i])-1;
    if (!is.logical(ind > 0 && abs(n_data[ind]-data[i])<radius)) {
      browser();
    }
    while (ind > 0 && abs(n_data[ind]-data[i])<radius) {
      counter <- counter + 1;
      ind = ind - 1;
    }
    ret[i] <- counter;
  }
  return(ret);
}
