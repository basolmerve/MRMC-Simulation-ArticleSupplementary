# Estimation Funciton for AUC (multiple tests, multiple readers)
# Reader'lar üzerinden bootstrap yapılmayacak ise (reader fixed) boot.reader = FALSE olarak bırakılacak.
estimationFunBWC <- function(x, data, nR = 3, nT = 2, nC = NULL, boot.reader = FALSE, return.diff = FALSE, auc.parametric = FALSE, return.all = FALSE){
  # Estimation Funciton for AUC (single test, single reader)
  aucEst <- function(x = NULL, marker, status, parametric = FALSE){
    # Args:
    #   x: indices of the selected cases
    if (is.null(x)){
      x <- 1:length(marker)
    }
    
    status.x <- status[x]
    nHD = as.numeric(table(status.x))
    
    nH <- nHD[1]
    nD <- nHD[2]
    marker0 <- marker[x][status.x == 0]
    marker1 <- marker[x][status.x != 0]
    
    if (!parametric){
      W <- wilcox.test(marker1, marker0, exact = FALSE)$statistic
      AUC <- W / (nD * nH)
    } else {
      mu1 <- mean(marker1)
      mu0 <- mean(marker0)
      sd1 <- sd(marker1)
      sd0 <- sd(marker0)
      
      a <- (mu1 - mu0) / sd1
      b <- sd0 / sd1
      
      AUC <- pnorm(a / sqrt(1 + b^2))
    }
    
    return(AUC)
  }
  
  ##
  ## Readerlar üzerinden bootstrap çekimine ait kodlamalar buraya eklendi
  ##
  if (boot.reader){
    r.boot <- sample(nR, nR, TRUE)
    data.tmp <- plyr:::ldply(lapply(r.boot, function(r){
      data[data[ ,"Reader"] == r, ]
    }), rbind)
    
    readerIDs <- rep(c(1:nR), rep(nC*nT, nR))
    data.tmp$Reader <- readerIDs
    data <- data.tmp
  }
  ##
  
  ## WAY I:
  data.split <- lapply(split(data, data[ ,"Test"]), function(y){
    split(y, y[ ,"Reader"])
  })
  
  ests <- unlist(lapply(data.split, function(y){
    unlist(lapply(y, function(z){
      aucEst(x, marker = z[ ,"Score"], status = z[ ,"Status"], parametric = auc.parametric)
    }))
  }))
  auc1 <- mean(ests[1:nR])
  auc2 <- mean(ests[(nR + 1):(length(ests))])
  
  
  # # WAY II:  (WAY I is faster)
  # data2 <- transformData_zhou(data)
  # data.split <- split(data2, data2[ ,"Reader"])
  # 
  # ests <- rowMeans(data.frame(lapply(data.split, function(y){
  #   auc1.tmp <- aucEst(x, marker = y[ ,"Test1"], status = y[ ,"Status"])
  #   auc2.tmp <- aucEst(x, marker = y[ ,"Test2"], status = y[ ,"Status"])
  #   c(auc1.tmp, auc2.tmp)
  # })))
  # 
  # auc1 <- ests[1]
  # auc2 <- ests[2]
  
  
  # return(auc1 - auc2)
  if (return.all){
    auc.all <- expand.grid(1:nR, 1:nT)
    colnames(auc.all) <- c("Reader", "Test")
    auc.all[ ,"AUCs"] <- ests
    auc.all <- dplyr::arrange(auc.all, Reader, Test)
    
    if (!return.diff){
      return(list(AUCs_all = auc.all, AUC_tests = c(auc1, auc2)))
    } else {
      return(list(AUCs_all = auc.all, AUC_tests = auc1 - auc2))
    }
  } else {
    if (!return.diff){
      return(c(auc1, auc2))
    } else {
      return(auc1 - auc2)
    }
  }
}