
findCov2 <- function(data, variable = NULL, marker = NULL, status = NULL, ...){
  require(pROC)
  # Covariances among readers.
  nReader <- length(unique(data[ ,variable]))
  covMat <- matrix(NA, nrow = nReader, ncol = nReader)
  
  for (r in 1:(nReader - 1)){
    for (j in (r + 1):nReader){
      roc1 <- roc(data[data[ ,variable] == r ,status] ~ data[data[ ,variable] == r ,marker])
      roc2 <- roc(data[data[ ,variable] == j ,status] ~ data[data[ ,variable] == j ,marker])
      
      covEstimate <- pROC::cov(roc1, roc2, ...)
      covMat[r, j] <- covMat[j, r] <- covEstimate
    }
  }
  
  return(covMat)
}


# data = dataCov
# variable = "Reader"
# marker = "Score"
# status = "Status"
# test = "Test"

findCov3 <- function(data, variable = NULL, marker = NULL, status = NULL, test = NULL, ...){
  require(pROC)
  
  Readers <- as.numeric(sort(unique(data[ ,variable])))
  Tests <- as.numeric(sort(unique(data[ ,test])))
  # 
  # nReader <- length(Readers)
  # nTest <- length(Tests)
  
  covkk <- NULL
  for (r in Readers[-length(Readers)]){
    for (t in Tests){
      for (k in Readers[!(Readers %in% r) & Readers > r]){
        roc1 <- roc(data[data[ ,variable] == r & data[ ,test] == t, status] ~ data[data[ ,variable] == r & data[ ,test] == t, marker])
        roc2 <- roc(data[data[ ,variable] == k & data[ ,test] != t, status] ~ data[data[ ,variable] == k & data[ ,test] != t, marker])

        cov.k <- data.frame(Reader = r, Test = t, Reader2 = k, Covariance = cov(roc1, roc2))
        covkk <- rbind(covkk, cov.k)
      }
    }
  }
  
  return(covkk)
}












