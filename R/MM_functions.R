gen.data <- function(Y, nN, nD, nR, ordinal = FALSE){
  # IMPORTANT NOTES: 
  #   1. Within each reader, first nN subjects should be healthy and following nD subjects should be diseased.
  #      Calculations will be incorrect if diseased subjects are located above healthy subject in data matrix.
  #   2. Last column of Y MUST BE readers.
  
  # ordinal: a logical. If TRUE, ties (i.e Y_healthy = Y_diseased) set to 1 instead 0.5.
  
  output <- matrix(nrow = nR * nN * nD, ncol = 5);
  
  # Make the data set
  # fill in reader, normal case, and diseased case numbers  
  for (i in 1:nR){    # go through the readers
    count <- nN * nD * (i - 1)     # for indexing
    subset <- Y[Y[ ,5] == i, ]    # data just for reader i    
    
    # reader number in 5th column:
    output[(count + 1):(count + nN * nD), 5] = i
    
    for (j in 1:nN){
      count.2 <- (j - 1) * nD
      
      # case number s (healthy subject) in 3rd column, k (diseased subject) in 4th column:
      output[(count + count.2 + 1):(count + count.2 + nD), 3] <- j
      output[(count + count.2 + 1):(count + count.2 + nD), 4] <- (nN + 1):(nD + nN)
      
      # fill in indicators
      if (!ordinal){
        output[(count + count.2 + 1:nD), 1] <- (subset[j, 1] < subset[(nN + 1:nD), 1]) + .5 * (subset[j, 1] == subset[(nN + 1:nD), 1])
        output[(count + count.2 + 1:nD), 2] <- (subset[j, 2] < subset[(nN + 1:nD), 2]) + .5 * (subset[j, 2] == subset[(nN + 1:nD), 2])  
      } else {
        output[(count + count.2 + 1:nD), 1] <- (subset[j, 1] <= subset[(nN + 1:nD), 1])
        output[(count + count.2 + 1:nD), 2] <- (subset[j, 2] <= subset[(nN + 1:nD), 2])
      }
    }
  }
  
  colnames(output) <- c("Test1", "Test2", "subjHealthy", "subjDiseased", "Reader")
  return(output)
}


getINF <-function(data, xbeta){
  #dU<- matrix(nrow = nR*nN*nD, ncol = 5);
  dU <- matrix(nrow = dim(data)[1], ncol = 5)
  temp1 <- matrix(nrow = dim(data)[1], ncol = 2)
  temp2 <- matrix(nrow = dim(data)[1], ncol = 2)
  
  # dU[ ,3:5] = data[ ,3:5]
  # column 3-> normal case, 4-> diseased case, 5 -> reader
  
  dnorm.1 = dnorm(xbeta[1])
  dnorm.2 = dnorm(xbeta[2])
  phi.1 = pnorm(xbeta[1])
  phi.2 = pnorm(xbeta[2])
  
  INF <- matrix(nrow = 2, ncol = 2)
  temp1[ ,1] <- (data[ ,1] - phi.1) / (phi.1 * (1 - phi.1))
  temp1[ ,2] <- (data[ ,2] - phi.2) / (phi.2 * (1 - phi.2))
  
  temp2[ ,1] <- (-1) * (temp1[ ,1] ** 2)
  temp2[ ,2] <- (-1) * (temp1[ ,2] ** 2)
  
  dU[ ,1] <- temp2[ ,1]*(dnorm.1 ** 2) - temp1[ ,1] * xbeta[1] * dnorm.1
  dU[ ,2] <- temp2[ ,2]*(dnorm.2 ** 2) - temp1[ ,2] * xbeta[2] * dnorm.2
  
  INF[1, 1] = sum(dU[ ,1], na.rm = TRUE) + sum(dU[ ,2], na.rm = TRUE)
  INF[1, 2] = sum(dU[ ,2], na.rm = TRUE)
  INF[2, 1] = INF[1, 2]
  INF[2, 2] = sum(dU[ ,2], na.rm = TRUE)
  
  return(INF)
}


getB <- function(data, xbeta, nN, nD, nR){
  B <- matrix(c(0, 0, 0, 0), 2, 2)
  u <- matrix(NA, nrow = dim(data)[1], ncol = 5)
  u[ ,3:5] <- data[ ,3:5]
  # column 3-> normal case, 4-> diseased case, 5 -> reader
  
  phi.1 <- pnorm(xbeta[1])
  phi.2 <- pnorm(xbeta[2])
  dnorm.1 <- dnorm(xbeta[1])
  dnorm.2 <- dnorm(xbeta[2])
  
  u[ ,1] <- ((data[ ,1] - phi.1) * dnorm.1) / (phi.1 * (1 - phi.1))
  u[ ,2] <- ((data[ ,2] - phi.2) * dnorm.2) / (phi.2 * (1 - phi.2))
  
  #U<-matrix(nrow = nR*nN*nN*2, ncol =2);
  #U[1:(nN*nR*nD),1]<-u[,1];
  #U[(nN*nR*nD=1):(nN*nR*nD*2), 1]<-u[,2];
  #U[1:(nN*nR*nD),2]<-0;
  #U[(nN*nR*nD=1):(nN*nR*nD*2), 2]<-u[,2];
  
  # healthy subjects
  for (s in 1:nN){
    temp <- u[u[ ,3] == s, ]
    len <- dim(temp)[1]
    subset <- matrix(c(temp[ ,1], temp[ ,2], rep(0, len), temp[ ,2]), 2 * len, 2)
    SI <- kronecker(subset, subset)
    
    si = matrix(c(sum(SI[ ,1], na.rm = TRUE), sum(SI[ ,3], na.rm = TRUE), sum(SI[ ,2], na.rm = TRUE), sum(SI[ ,4], na.rm = TRUE)), 2, 2)
    B = B + si
  }
  
  # diseased subjects
  for (k in (nN + 1):(nN + nD)){
    temp <- u[u[ ,4] == k,]
    len <- dim(temp)[1]
    subset <- matrix(c(temp[ ,1], temp[ ,2], rep(0, len), temp[ ,2]), 2 * len, 2)
    SI = kronecker(subset, subset)
    si = matrix(c(sum(SI[ ,1], na.rm = TRUE), sum(SI[ ,3], na.rm = TRUE), sum(SI[ ,2], na.rm = TRUE), sum(SI[ ,4], na.rm = TRUE)), 2, 2)
    B = B + si
  }
  
  for (s in 1:nN){
    for (k in (nN + 1):(nN + nD)){
      temp <- u[(u[ ,4] == k & u[ ,3] == s), ]
      len <- dim(temp)[1]
      subset <-matrix(c(temp[ ,1], temp[ ,2], rep(0, len), temp[ ,2]), 2 * len, 2)
      SI = kronecker(subset, subset)
      si = matrix(c(sum(SI[ ,1], na.rm = TRUE), sum(SI[ ,3], na.rm = TRUE), sum(SI[ ,2], na.rm = TRUE), sum(SI[ ,4], na.rm = TRUE)), 2, 2)
      B = B + si
    }
  }
  return(B)
}


transformData_zhou <- function(x, reader = "Reader", test = "Test", 
                               score = "Score", case = "Cases", status = "Status"){
  # IMPORTANT NOTES: 
  #   1. Within each reader, first nN subjects should be healthy and following nD subjects should be diseased.
  #      Calculations will be incorrect if diseased subjects are located above healthy subject in data matrix.
  #   2. Last column of Y MUST BE readers.
  
  library(dplyr)
  
  tmp <- split(x, x[ ,test])
  tmp1 <- tmp[[1]]
  tmp2 <- cbind(tmp1, tmp[[2]][ ,score])
  
  colnames(tmp2)[ncol(tmp2)] <- paste0(score, ".2")
  
  tmp2 <- tmp2[ ,c(score, paste0(score, ".2"), status, case, reader)]
  
  tmp2 <- eval(
    substitute(arrange(tmp2, reader, status, case), 
               list(reader = as.name(reader),
                    status = as.name(status),
                    case = as.name(case)))
  )
  
  colnames(tmp2) <- c("Test1", "Test2", "Status", "Cases", "Reader")
  return(tmp2)
}