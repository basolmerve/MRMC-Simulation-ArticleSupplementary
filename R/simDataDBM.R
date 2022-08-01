# prmtrs <- modelParameters

simDataDBM <- function(prmtrs = NULL, varcomp.unequal = FALSE){
  
  if (!varcomp.unequal){
    nh <- nd <- prmtrs$nhd
    n <- nh + nd
    nr <- prmtrs$nr
    nt <- prmtrs$nt
    
    var_case = prmtrs$var_case
    var_reader = prmtrs$var_reader
    var_readercase = prmtrs$var_readercase
    var_casetest = prmtrs$var_casetest
    var_readertest = prmtrs$var_readertest
    var_error = prmtrs$var_error
    
    mu <- prmtrs$mu
    CaseIDs <- 1:(nh + nd)
    status <- rep(c(0, 1), c(nh, nd))
    
    DesignMatrix <- expand.grid(Status = status, Reader = 1:nr, Test = 1:nt)
    DesignMatrix <- data.frame(Cases = CaseIDs, DesignMatrix)
    
    d <- rnorm(n = nd, mean = 0, sd = sqrt(var_case))
    h <- rnorm(n = nh, mean = 0, sd = sqrt(var_case))
    caseEffects <- c(h, d)
    
    # 1. Ana etkiler
    # Reader Effect
    r <- rnorm(n = nr, mean = 0, sd = sqrt(var_reader))
    
    # 2. İkili Etkileşimler
    # Reader-Case Etkileşimi
    rc <- rnorm(n = nr * (nh + nd), mean = 0, sd = sqrt(var_readercase))
    rc <- rep(rc, nt)
    
    # Case-Test etkileşimi
    ct <- c(rep(rnorm(n = (nd + nh), mean = 0, sd = sqrt(var_casetest)), nr), rep(rnorm(n = (nd + nh), mean = 0, sd = sqrt(var_casetest)), nr))
    
    # Reader-Test etkileşimi
    rt <- rnorm(n = nr * nt, mean = 0, sd = sqrt(var_readertest))
    rt <- rep(rt, times = rep(n, length(rt)))
    
    # 3. Rastgele hata
    error <- rnorm(n = nrow(DesignMatrix), mean = 0, sd = sqrt(var_error))
    
    scores <- NULL
    for (i in 1:nrow(DesignMatrix)){
      scores[i] <- mu * DesignMatrix[i, "Status"] + caseEffects[DesignMatrix[i, "Cases"]] + r[DesignMatrix[i, "Reader"]] + rc[i] + ct[i] + rt[i] + error[i]
    }
    
    DesignMatrix[ ,"Score"] <- scores
    
  } else {
    ### Unequal variance components.
    nh <- nd <- prmtrs$nhd
    n <- nh + nd
    nr <- prmtrs$nr
    nt <- prmtrs$nt
    
    var_readertest = prmtrs$var_readertest
    var_reader = prmtrs$var_reader
    
    var_case0 = prmtrs$var_case0
    var_readercase0 = prmtrs$var_readercase0
    var_casetest0 = prmtrs$var_casetest0
    var_error0 = prmtrs$var_error0
    
    var_case1 = prmtrs$var_case1
    var_readercase1 = prmtrs$var_readercase1
    var_casetest1 = prmtrs$var_casetest1
    var_error1 = prmtrs$var_error1
    
    mu <- prmtrs$mu
    CaseIDs <- 1:(nh + nd)
    status <- rep(c(0, 1), c(nh, nd))
    
    DesignMatrix <- expand.grid(Status = status, Reader = 1:nr, Test = 1:nt)
    DesignMatrix <- data.frame(Cases = CaseIDs, DesignMatrix)
    
    d <- rnorm(n = nd, mean = 0, sd = sqrt(var_case1))
    h <- rnorm(n = nh, mean = 0, sd = sqrt(var_case0))
    caseEffects <- c(h, d)
    
    # 1. Ana etkiler
    # Reader Effect
    r <- rnorm(n = nr, mean = 0, sd = sqrt(var_reader))
    
    # 2. İkili Etkileşimler
    # Reader-Case Etkileşimi
    # rc <- rnorm(n = nr*(nh + nd), mean = 0, sd = sqrt(var_readercase))
    # rc <- rep(rc, nt)
    
    rc <- NULL
    for (i in 1:nr){
      # Her reader için hasta ve sağlıklıların etkileri hesaplanıyor.
      # Bütün readerlar için hasta ve sağlıklılar ayrı ayrı hesaplandıktan sonra bu değerler bütün testler için aynı şekilde kullanılıyor.
      rc0 <- rnorm(n = nh, mean = 0, sd = sqrt(var_readercase0))
      rc1 <- rnorm(n = nd, mean = 0, sd = sqrt(var_readercase1))
      
      rc <- c(rc, rc0, rc1)
    }
    rc <- rep(rc, nt)
    
    # Case-Test etkileşimi
    # ct <- c(rep(rnorm(n = (nd + nh), mean = 0, sd = sqrt(var_casetest)), nr), rep(rnorm(n = (nd + nh), mean = 0, sd = sqrt(var_casetest)), nr))
    ct <- NULL
    for (i in 1:nt){
      # Her reader için hasta ve sağlıklıların etkileri hesaplanıyor.
      # Bütün readerlar için hasta ve sağlıklılar ayrı ayrı hesaplandıktan sonra bu değerler bütün testler için aynı şekilde kullanılıyor.
      ct0 <- rnorm(n = nh, mean = 0, sd = sqrt(var_readercase0))
      ct1 <- rnorm(n = nd, mean = 0, sd = sqrt(var_readercase1))
      
      ct.tmp <- c(ct0, ct1)
      ct.tmp <- rep(ct.tmp, nr)
      
      ct <- c(ct, ct.tmp)
    }
    
    # Reader-Test etkileşimi
    rt <- rnorm(n = nr*nt, mean = 0, sd = sqrt(var_readertest))
    rt <- rep(rt, times = rep(n, length(rt)))
    
    # 3. Rastgele hata
    error <- numeric(nrow(DesignMatrix))
    n00 <- sum(DesignMatrix$Status == 0)
    n11 <- sum(DesignMatrix$Status == 1)
    
    error0 <- rnorm(n = n00, mean = 0, sd = sqrt(var_error0))
    error1 <- rnorm(n = n11, mean = 0, sd = sqrt(var_error1))
    
    error[DesignMatrix$Status == 0] <- error0
    error[DesignMatrix$Status != 0] <- error1
    
    scores <- NULL
    for (i in 1:nrow(DesignMatrix)){
      scores[i] <- mu*DesignMatrix[i, "Status"] + caseEffects[DesignMatrix[i, "Cases"]] + 
                   r[DesignMatrix[i, "Reader"]] + rc[i] + ct[i] + rt[i] + error[i]
    }
    
    DesignMatrix[ ,"Score"] <- scores
  }
  
  return(DesignMatrix)
}