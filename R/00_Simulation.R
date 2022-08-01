#### Load libraries
library(parallel)
library(pROC)
library(dplyr)

# Load source codes
source("R/rocdata.R")
source("R/pseudo_values.R")
source("R/simDataDBM.R")
source("R/covOR.R")
source("R/estimationFunBWC.R")
source("R/bcanon2.R")
source("R/boott.R")
source("R/MM_functions.R")
source("R/helper_functions.R")

# Activate parallel CPU cores. "spec" is the number of CPUs.
cl <- makeCluster(spec = 120, type = "PSOCK")

# Methods to be used in the MRMC study.
methods <- c("DBM", "OR", "MM", "BCa_Percentile")   # c("DBM", "OR", "MM", "BCa_Percentile", "boott", "stan")

## Which variance components to use? Equal or Unequal. If FALSE, equal variance components are used.
varcomp.unequal <- FALSE    

# Number of generated datasets in each scenario. Note that this is not the number of bootstrap samples.
B <- 1000 

# Number of bootstraps for the BWC method (resampling-based methods)
nboot.bca <- 2000

# Some extra parameters for the boot-t method. This method is not included in the simulation study.
# Therefore, below parameters has no effect on the results. Keep it unchanged.
nbootsd = 25
nboott = 200
v.nbootg = 100
v.nbootsd = 25
v.nboott = 200
VS = TRUE

# Read the varaince components from local files.
if (varcomp.unequal){
  simCombs <- read.table("data/simCombs_unequal.txt", header = TRUE, stringsAsFactors = FALSE)
} else {
  simCombs <- read.table("data/simCombs_equal.txt", header = TRUE, stringsAsFactors = FALSE)
}


simCombsResults <- NULL
elapsedTimes <- NULL

# Start the simulation here ...
start_all <- Sys.time()
for (i in 1:nrow(simCombs)){
  simCombs.i <- simCombs[i, ]
  modelParameters <- as.list(simCombs.i)
  scen <- simCombs.i$scen
  
  clearConsole()
  printSimulationHeader(step = i, total = nrow(simCombs), equal.var = !varcomp.unequal, 
                        seed = simCombs.i$seed, id = simCombs.i$simID)
  printSimulationInfo(.info = simCombs.i, equal.var = !varcomp.unequal, .ndatasets = B, .nbootstrap = nboot.bca)
  printVarianceComponents(.info = simCombs.i, equal.var = !varcomp.unequal)
  cat("\n")
  
  MU <- if (modelParameters$mu == 0.75){
    "075"
  } else if (modelParameters$mu == 1.5){
    "15"
  } else if (modelParameters$mu == 2.5){
    "25"
  } else if (modelParameters$mu == 0.821){
    "821"
  } else if (modelParameters$mu == 1.831){
    "1831"
  } else {
    "3661"
  }
  
  if (!varcomp.unequal){
    # Check if required folders exist.
    if (!dir.exists(paths = "results/CIs_Equal")){
      dir.create("results/CIs_Equal", recursive = TRUE)
    }
    fpName <- file.path(getwd(), "results/CIs_Equal", paste("H", modelParameters$nhd, "D", modelParameters$nhd, "_",
                                              "R", modelParameters$nr, "_", "T", modelParameters$nt, "_",
                                              "MU", MU, "_", scen, ".txt", sep = ""))
  } else {
    # Check if required folders exist.
    if (!dir.exists(paths = "results/CIs_Unequal")){
      dir.create("results/CIs_Unequal", recursive = TRUE)
    }
    fpName <- file.path(getwd(), "results/CIs_Unequal", paste("H", modelParameters$nhd, "D", modelParameters$nhd, "_",
                                              "R", modelParameters$nr, "_", "T", modelParameters$nt, "_",
                                              "MU", MU, "_", scen, ".txt", sep = ""))
  }
  
  clusterSetRNGStream(cl = cl, iseed = modelParameters$seed)  ## set seed for parallel clusters (child nodes).
  set.seed(modelParameters$seed, kind = "L'Ecuyer-CMRG")  # set seed for master node.
  
  ## Analizler için kullanılan bütün datalar:
  cat(rep(".", 6), " Generating datasets...", sep = "")
  datasets <- lapply(1:B, function(x){
    simDataDBM(modelParameters, varcomp.unequal)
  })
  cat(" [OK]", "\n")
  
  ### If boot-t method is not selected:
  if (!("boott" %in% methods)){
    nbootsd <- nboott <- v.nbootg <- v.nbootsd <- v.nboott <- NA
    VS <- TRUE
  }
  
  # Paralel çekirdeklerde kullanılacak elemanların çekirdeklere gönderilmesi.
  clusterExport(cl, c("modelParameters", "pseudo_values", "rocdata", "simDataDBM", "fpName",
                      "findCov2", "findCov3", "nboot.bca", "bcanon2", "estimationFunBWC",
                      "boott2", "nbootsd", "nboott", "v.nbootg", "v.nbootsd", "v.nboott", "VS",
                      "gen.data", "getB", "getINF", "transformData_zhou"))
  
  DBM_SimRes <- OR_SimRes <- MM_SimRes <- BootT_SimRes <- NULL
  BCA_Percentile_SimRes <- NULL
  Coverage_DBM <- Coverage_OR <- Coverage_MM <- NA
  Coverage_BCa <- Coverage_Percentile <- Coverage_Boott <- Coverage_Stan <- NA
  elapsed_DBM <- elapsed_OR <- elapsed_MM <- elapsed_BCaPerc <- elapsed_boott <- NA
    
  for (m in methods){
    if (m == "DBM"){
      cat(rep(".", 6), " Fitting model (DBM)...", sep = "")
      ####### 1. DBM ######
      start_DBM <- Sys.time()
      DBM_Results <- parLapply(cl = cl, X = datasets, fun = function(simulatedData = X, filePath = fpName, 
                                                                     modelParams = modelParameters){
        
        nR <- modelParams$nr
        
        # We used equal sample sizes for healthy and diseased subjects.
        # One may easily modify this line to set unequal sample size scenarios.
        nD <- nH <- modelParams$nhd   
        n <-  nH + nD
        nT <- modelParams$nt
        
        #dat <- simulatedData[[x]]
        dat <- simulatedData
        
        dat$Cases <- as.factor(dat$Cases)
        dat$Reader <- as.factor(dat$Reader)
        dat$Test <- as.factor(dat$Test)
        dat$Status <- as.factor(dat$Status)
        
        data_split <- split(dat, dat$Reader:dat$Test)
        
        pseudos <- lapply(data_split, function(x){
          pseudo_values(object = x, statusVariable = "Status", marker = "Score")
        })
        
        dat <- plyr::ldply(pseudos, rbind)[ ,-1]

        DBM_Data <- estimationFunBWC(x = 1:n, data = dat, nR = nR, nT = nT, nC = n, boot.reader = FALSE, 
                                     return.diff = FALSE, auc.parametric = FALSE, return.all = TRUE)
      
        # AUCs:
        AUC1 <- as.numeric(DBM_Data$AUC_tests[1])
        AUC2 <- as.numeric(DBM_Data$AUC_tests[2])
        auc.diff <- AUC1 - AUC2
        
        ### ANOVA Model DBM
        model <- anova(aov(pseudo ~ Cases + Test + Reader + Cases:Test + Reader:Test + Cases:Reader + Cases:Reader:Test, data = dat))
        
        Residual <- model[which(rownames(model) == "Test:Reader"), "Mean Sq"] + max(0, (model[which(rownames(model) == "Cases:Test"), "Mean Sq"] - model[which(rownames(model) == "Cases:Test:Reader"), "Mean Sq"]))
        FStat <- model[which(rownames(model) == "Test"), "Mean Sq"] / Residual
        df <- (Residual^2) / ((((model[which(rownames(model) == "Test:Reader"), "Mean Sq"])^2) / (model[which(rownames(model) == "Test:Reader"), "Df"])))
        pValue <- 1 - pf(q = FStat, df1 = 1, df2 = df)
        
        
        # CI for AUC differences. (Random, DBM)
        lowerLimit <- upperLimit <- NA
        lowerLimit <- as.numeric(auc.diff - qt(p = 0.975, df = df) * sqrt(Residual * (2 / (nR * n))))
        upperLimit <- as.numeric(auc.diff + qt(p = 0.975, df = df) * sqrt(Residual * (2 / (nR * n))))
        
        result <- data.frame(Healthy = nH, Diseased = nD,
                             Reader = nR, Test = modelParams$nt, Model = "DBM",
                             AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                             pValue = pValue, LowerLimit = lowerLimit, UpperLimit = upperLimit)
        
        return(result)
      })
      cat("[OK]", "\n")
      
      # Sonuçların tek bir data.frame olarak birleştirilmesi
      DBM_SimRes <- bind_rows(DBM_Results)
      
      DBM_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "pValue", "LowerLimit", "UpperLimit")] <- round(DBM_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "pValue", "LowerLimit", "UpperLimit")], 6)
      Coverage_DBM <- sum(DBM_SimRes$pValue > 0.05) / length(DBM_SimRes$pValue)
      
      end_DBM <- Sys.time()
      elapsed_DBM <- as.numeric(difftime(end_DBM, start_DBM, units = "mins"))
      
    } else  if (m == "OR"){
      cat(rep(".", 6), " Fitting model (OR)...", sep = "")
      ####### 2. OR ######
      start_OR <- Sys.time()
      OR_Results <- parLapply(cl = cl, X = datasets, fun = function(simulatedData = X, filePath = fpName, modelParams = modelParameters){
        
        dat <- dataCov <- simulatedData
        
        nR <- modelParams$nr
        nD <- nH <- modelParams$nhd
        n <- nH + nD
        nT <- modelParams$nt
        
        ## OR Model
        # Split data via treatment:reader combinations
        # reader ve test'lerin factor hale getirilmesi
        dat[ ,"Reader"] <- as.factor(dat[ ,"Reader"])
        dat[ ,"Test"] <- as.factor(dat[ ,"Test"])
        dat[ ,"Cases"] <- as.factor(dat[ ,"Cases"])
        dat[ ,"Status"] <- as.factor(dat[ ,"Status"])
        
        OR_Data <- estimationFunBWC(x = 1:n, data = dat, nR = nR, nT = nT, nC = n, boot.reader = FALSE, 
                                    return.diff = FALSE, auc.parametric = FALSE, return.all = TRUE)
        
        # AUCs:
        AUC1 <- as.numeric(OR_Data$AUC_tests[1])
        AUC2 <- as.numeric(OR_Data$AUC_tests[2])
        auc.diff <- AUC1 - AUC2
        
        OR_Data$AUCs_all$Reader <- as.factor(OR_Data$AUCs_all$Reader)
        OR_Data$AUCs_all$Test <- as.factor(OR_Data$AUCs_all$Test)
        
        OR_ModelData <- OR_Data$AUCs_all
        
        ### ANOVA MODEL OR
        modelOR <- anova(aov(AUCs ~ Reader + Test + Test:Reader, data = OR_ModelData))
        
        ## Covariance Estimates:
        dataCov_split <- split(dataCov, dataCov$Test)
        cov2Estimates <- lapply(dataCov_split, findCov2, variable = "Reader", marker = "Score", status = "Status", method = "delong")
        
        summ <- nn <- NULL
        for (i in 1:length(cov2Estimates)){
          tmp.i <- cov2Estimates[[i]]
          sum.i <- sum(tmp.i[upper.tri(tmp.i)])
          
          nn <- c(nn, length(tmp.i[upper.tri(tmp.i)]))
          summ <- c(sum.i, summ)
        }
        
        cov2Est <- sum(summ) / sum(nn)
        
        cov3Estimates <- findCov3(data = dataCov, variable = "Reader", marker = "Score", status = "Status", test = "Test")
        cov3Est <- mean(cov3Estimates$Covariance)
        
        ResidualOR <- modelOR[which(rownames(modelOR) == "Reader:Test"), "Mean Sq"] + max(0, (cov2Est - cov3Est)) * (modelOR[which(rownames(modelOR) == "Reader"), "Df"] + 1)
        FStatOR <- modelOR[which(rownames(modelOR) == "Test"), "Mean Sq"] / ResidualOR
        dfOR <- ResidualOR^2 / ((modelOR[which(rownames(modelOR) == "Reader:Test"), "Mean Sq"])^2 / (modelOR[which(rownames(modelOR) == "Reader:Test"), "Df"]))
        pValueOR <- 1 - pf(q = FStatOR, df1 = 1, df2 = dfOR)
        
        # CI for AUC differences. (OR, Random)
        lowerLimitOR <- upperLimitOR <- NA
        StdErr <- sqrt((2 / (modelOR[which(rownames(modelOR) == "Reader"), "Df"] + 1)) * ResidualOR)
        lowerLimitOR <- as.numeric(auc.diff - qt(p = 0.975, df = dfOR) * StdErr)
        upperLimitOR <- as.numeric(auc.diff + qt(p = 0.975, df = dfOR) * StdErr)
        
        result <- data.frame(Healthy = nH, Diseased = nD,
                             Reader = nR, Test = modelParams$nt, Model = "OR", 
                             AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                             pValue = pValueOR, LowerLimit = lowerLimitOR, 
                             UpperLimit = upperLimitOR)
        
        return(result)
      })
      cat("[OK]", "\n")
      
      # Sonuçların tek bir data.frame olarak birleştirilmesi
      OR_SimRes <- bind_rows(OR_Results)
      
      OR_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "pValue", "LowerLimit", "UpperLimit")] <- round(OR_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "pValue", "LowerLimit", "UpperLimit")], 6)
      Coverage_OR <- sum(OR_SimRes$pValue > 0.05) / length(OR_SimRes$pValue)
      
      end_OR <- Sys.time()
      elapsed_OR <- as.numeric(difftime(end_OR, start_OR, units = "mins"))
      
    } else if (m == "MM"){
      cat(rep(".", 6), " Fitting model (MM)...", sep = "")
      ####### 3. MM ######
      start_MM <- Sys.time()
      MM_Results <- parLapply(cl = cl, X = datasets, fun = function(simulatedData = X, filePath = fpName, modelParams = modelParameters){
        
        # dat <- simulatedData[[x]]
        dat <- simulatedData
        
        nR <- modelParams$nr
        nD <- nH <- modelParams$nhd
        n <- nH + nD
        nT <- modelParams$nt
        
        Y <- transformData_zhou(dat, reader = "Reader", test = "Test", 
                                score = "Score", case = "Cases", status = "Status")
        
        BIdata <- gen.data(Y, nH, nD, nR, ordinal = FALSE)
        
        response <- c(BIdata[ ,1], BIdata[ ,2])
        len <- length(BIdata[ ,1])
        predictor.1 <- c(rep(0, len), rep(1, len))
        lin.mod.1 <- glm(response ~ predictor.1, family = binomial(link = "probit"))
        beta0 <- lin.mod.1$coefficients[1]
        beta1 <- lin.mod.1$coefficients[2]
        beta <- c(beta0, beta1, beta0 + beta1)
        xbeta <- c(beta0, beta0 + beta1)
        
        AUC <- c(pnorm(beta0), pnorm(beta0 + beta1), 1)
        AUC[3] <- AUC[1] - AUC[2]
        se <- matrix(c(0, 0, 0), 3, 1)
        se_AUC <- matrix(c(0, 0, 0), 3, 1)
        temp1 <- dnorm(beta0)
        temp2 <- dnorm(beta0 + beta1)
        
        #compute A
        A = getINF(BIdata, xbeta)
        A = solve(A)
        
        # compute B
        B = getB(BIdata, xbeta, nH, nD, nR)
        V = A %*% B %*% t(A)
        D <- matrix(c(1, 1), 1, 2)
        V_beta2 = D %*% V %*% t(D)
        
        D <- matrix(c(temp1, temp2, temp1 - temp2, 0, temp2, (-1) * temp2), 3, 2)
        V_AUC <- D %*% V %*% t(D)
        
        se = c(sqrt(abs(diag(V))), sqrt(abs(V_beta2)))
        se_AUC = sqrt(abs(diag(V_AUC)))
        
        AUC_low = c(pnorm(beta[1] - qnorm(0.975) * se[1]), pnorm(beta[3] - qnorm(0.975) * se[3]), (AUC[3] - qnorm(0.975) * se_AUC[3]))
        AUC_up = c(pnorm(beta[1] + qnorm(0.975) * se[1]), pnorm(beta[3] + qnorm(0.975) * se[3]), (AUC[3] + qnorm(0.975) * se_AUC[3]))
        
        AUC1 <- AUC[1]
        AUC2 <- AUC[2]
        auc.diff <- AUC[3]
        
        lowerLimitMM <- upperLimitMM <- NA
        if (!is.null(AUC_low) & !is.null(AUC_up)){
          lowerLimitMM <- AUC_low[3]
          upperLimitMM <- AUC_up[3]
        }
        
        result <- data.frame(Healthy = nH, Diseased = nD,
                             Reader = nR, Test = modelParams$nt, Model = "MM",
                             AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                             pValue = NA, LowerLimit = lowerLimitMM, UpperLimit = upperLimitMM)
        
        return(result)
      })
      cat("[OK]", "\n")
      
      # Sonuçların tek bir data.frame olarak birleştirilmesi
      MM_SimRes <- bind_rows(MM_Results)
      
      MM_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")] <- round(MM_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")], 6)
      
      Coverage_MM <- apply(as.matrix(MM_SimRes[ ,c("LowerLimit", "UpperLimit")]), 1, function(x){
        x["LowerLimit"] <= 0 & x["UpperLimit"] >= 0
      })
      
      Coverage_MM <- Coverage_MM[complete.cases(Coverage_MM)]
      
      if (length(Coverage_MM) > 0){
        Coverage_MM <- sum(Coverage_MM)/length(Coverage_MM)
      } else {
        Coverage_MM <- NA
      }
      
      end_MM <- Sys.time()
      elapsed_MM <- as.numeric(difftime(end_MM, start_MM, units = "mins"))
      
    } else if (m == "BCa_Percentile"){
      cat(rep(".", 6), " Fitting model (BCa and Percentile Bootstrap)...", sep = "")
      ####### 4. BCa and Percentile ######
      start_BCaPerc <- Sys.time()
      BCA_Percentile_Results <- parLapply(cl = cl, X = datasets, fun = function(simulatedData = X, filePath = fpName, modelParams = modelParameters){
        library(dplyr)
        library(bootstrap)
        
        nR <- modelParams$nr
        nD <- nH <- modelParams$nhd
        n <- nH + nD
        nT <- modelParams$nt
        
        # dat <- simulatedData[[x]]
        dat <- simulatedData
        
        # status <- dat$Status
        nCase <- c(modelParams$nhd, modelParams$nhd)
        
        ## BCa güven aralığı hesaplamaları.
        bcaRes <- try({bcanon2(x = 1:sum(nCase), theta2 = estimationFunBWC, nboot = nboot.bca, alpha = c(0.025, 0.975), nC = sum(nCase),
                               data = dat, nR = modelParams$nr, nT = modelParams$nt, boot.reader = TRUE)})
        if (class(bcaRes) == "try-error"){
          bcaRes <- try({bcanon2(x = 1:sum(nCase), theta2 = estimationFunBWC, nboot = nboot.bca, alpha = c(0.025, 0.975), nC = sum(nCase),
                                 data = dat, nR = modelParams$nr, nT = modelParams$nt, boot.reader = TRUE)})
          if (class(bcaRes) == "try-error"){
            bcaRes <- NULL
          }
        }
        
        est <- try({estimationFunBWC(x = 1:sum(nCase), data = dat, nR = modelParams$nr, nT = modelParams$nt,
                                      nC = sum(nCase), boot.reader = FALSE, return.diff = FALSE)})
        
        if (class(est) != "try-error"){
          AUC1 <- est[1]
          AUC2 <- est[2]
        } else {
          AUC1 <- AUC2 <- NA
        }
        
        auc.diff <- AUC1 - AUC2

        lowerPerc <- upperPerc <- lowerBCa <- upperBCa <- NA
        if (!is.null(bcaRes)){
          lowerPerc = bcaRes$confpoints[1, 2]
          upperPerc = bcaRes$confpoints[2, 2]
          lowerBCa = bcaRes$confpoints[1, 3]
          upperBCa = bcaRes$confpoints[2, 3]
        }
        
        result_bca <- data.frame(Healthy = nH, Diseased = nD,
                             Reader = nR, Test = modelParams$nt, Model = "BCa",
                             AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                             pValue = NA, LowerLimit = lowerBCa, UpperLimit = upperBCa)
        
        result_perc <- data.frame(Healthy = nH, Diseased = nD,
                                 Reader = nR, Test = modelParams$nt, Model = "Percentile",
                                 AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                                 pValue = NA, LowerLimit = lowerPerc, UpperLimit = upperPerc)

        
        result <- rbind(result_bca, result_perc)
        return(result)
        
      })
      cat("[OK]", "\n")
      
      # Sonuçların tek bir data.frame olarak birleştirilmesi
      BCA_Percentile_SimRes <- bind_rows(BCA_Percentile_Results)
      
      BCA_Percentile_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")] <- round(BCA_Percentile_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")], 6)
      
      Coverage_BCa <- apply(as.matrix(BCA_Percentile_SimRes[BCA_Percentile_SimRes$Model == "BCa" ,c("LowerLimit", "UpperLimit")]), 1, function(x){
        x["LowerLimit"] <= 0 & x["UpperLimit"] >= 0
      })
      
      Coverage_BCa <- Coverage_BCa[complete.cases(Coverage_BCa)]
      
      if (length(Coverage_BCa) > 0){
        Coverage_BCa <- sum(Coverage_BCa)/length(Coverage_BCa)
      } else {
        Coverage_BCa <- NA
      }
      
      Coverage_Percentile <- apply(as.matrix(BCA_Percentile_SimRes[BCA_Percentile_SimRes$Model == "Percentile" ,c("LowerLimit", "UpperLimit")]), 1, function(x){
        x["LowerLimit"] <= 0 & x["UpperLimit"] >= 0
      })
      
      Coverage_Percentile <- Coverage_Percentile[complete.cases(Coverage_Percentile)]
      
      if (length(Coverage_Percentile) > 0){
        Coverage_Percentile <- sum(Coverage_Percentile)/length(Coverage_Percentile)
      } else {
        Coverage_Percentile <- NA
      }
      
      end_BCaPerc <- Sys.time()
      elapsed_BCaPerc <- as.numeric(difftime(end_BCaPerc, start_BCaPerc, units = "mins"))
      
    } else if (m == "boott"){
      cat(rep(".", 6), " Fitting model (boott)...", sep = "")
      ####### 5. Boot-t ######
      start_boott <- Sys.time()
      BootT_Results <- parLapply(cl = cl, X = datasets, fun = function(simulatedData = X, filePath = fpName, modelParams = modelParameters){
        library(dplyr)
        library(bootstrap)
        
        nR <- modelParams$nr
        nD <- nH <- modelParams$nhd
        n <- nH + nD
        nT <- modelParams$nt
        
        # dat <- simulatedData[[x]]
        dat <- simulatedData
        
        # status <- dat$Status
        nCase <- c(modelParams$nhd, modelParams$nhd)
        
        ## bootT güven aralığı hesaplamaları.
        boottRes <- try({boott2(x = 1:sum(nCase), theta = estimationFunBWC, nbootsd = nbootsd, nboott = nboott, v.nbootg = v.nbootg, 
                              v.nbootsd = v.nbootsd, v.nboott = v.nboott, VS = VS, perc = c(0.025, 0.975),
                              data = dat, nR = modelParams$nr, nT = modelParams$nt, nC = sum(nCase), boot.reader = TRUE, return.diff = TRUE)})
        if (class(boottRes) == "try-error"){
          boottRes <- try({boott2(x = 1:sum(nCase), theta = estimationFunBWC, nbootsd = nbootsd, nboott = nboott, v.nbootg = v.nbootg, 
                                v.nbootsd = v.nbootsd, v.nboott = v.nboott, VS = VS, perc = c(0.025, 0.975),
                                data = dat, nR = modelParams$nr, nT = modelParams$nt, nC = sum(nCase), boot.reader = TRUE, return.diff = TRUE)})
          if (class(boottRes) == "try-error"){
            boottRes <- NULL
          }
        }
        
        est <- try({estimationFunBWC(x = 1:sum(nCase), data = dat, nR = modelParams$nr, nT = modelParams$nt,
                                     nC = sum(nCase), boot.reader = FALSE, return.diff = FALSE)})
        
        if (class(est) != "try-error"){
          AUC1 <- est[1]
          AUC2 <- est[2]
        } else {
          AUC1 <- AUC2 <- NA
        }
        
        auc.diff <- AUC1 - AUC2
        
        lowerBootT <- upperBootT <- NA
        if (!is.null(boottRes)){
          lowerBootT = boottRes$confpoints[1, 1]
          upperBootT = boottRes$confpoints[1, 2]
        }
        
        result <- data.frame(Healthy = nH, Diseased = nD,
                             Reader = nR, Test = modelParams$nt, Model = "bootT",
                             AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                             pValue = NA, LowerLimit = lowerBootT, UpperLimit = upperBootT)
      
        return(result)
      })
      cat("[OK]", "\n")
      
      # Sonuçların tek bir data.frame olarak birleştirilmesi
      BootT_SimRes <- bind_rows(BootT_Results)
      
      BootT_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")] <- round(BootT_SimRes[ ,c("AUC1", "AUC2", "AUC.Diff.", "LowerLimit", "UpperLimit")], 6)
      
      Coverage_Boott <- apply(as.matrix(BootT_SimRes[ ,c("LowerLimit", "UpperLimit")]), 1, function(x){
        x["LowerLimit"] <= 0 & x["UpperLimit"] >= 0
      })
      Coverage_Boott <- Coverage_Boott[complete.cases(Coverage_Boott)]
      if (length(Coverage_Boott) > 0){
        Coverage_Boott <- sum(Coverage_Boott)/length(Coverage_Boott)
      } else {
        Coverage_Boott <- NA
      }
      
      end_boott <- Sys.time()
      elapsed_boott <- as.numeric(difftime(end_boott, start_boott, units = "mins"))
      
    } else if (m == "stan"){
      ####### 6. Stan ######
      
    }
  }
  
  # Create "results" folder to print simulation results
  if (!dir.exists("results")){
    dir.create("results")
  }
  
  CI_Res <- bind_rows(DBM_SimRes, OR_SimRes, MM_SimRes, BCA_Percentile_SimRes, BootT_SimRes)
  write.table(CI_Res, file = fpName, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  simCombsResults.i <- data.frame(simCombs.i, Coverage_DBM, Coverage_OR, Coverage_MM, Coverage_Boott,
                                  Coverage_Stan, Coverage_Percentile, Coverage_BCa)
  
  ## Computation time:
  elapsedTimes.i <- data.frame(simCombs.i, DBM = elapsed_DBM, OR = elapsed_OR, MM = elapsed_MM, BootT = elapsed_boott,
                               Stan = NA, Percentile = NA, BCa = elapsed_BCaPerc)
  
  if (varcomp.unequal){
    resultFileName <- "simCombsResults_Unequal.txt"
    elapsedFileName <- "elapsedTimes_Unequal.txt"
  } else {
    resultFileName <- "simCombsResults_Equal.txt"
    elapsedFileName <- "elapsedTimes_Equal.txt"
  }
  
  if (file.exists(file.path("results", resultFileName))){
    tmp1 <- read.table(file = file.path("results", resultFileName), header = TRUE, sep = "\t")
    tmp2 <- read.table(file = file.path("results", elapsedFileName), header = TRUE, sep = "\t")
    
    if (simCombsResults.i$simID %in% tmp1$simID){
      rowIdx <- which(tmp1$simID == simCombsResults.i$simID)
      tmp1[rowIdx, ] <- simCombsResults.i
    } else {
      tmp1 <- bind_rows(tmp1, simCombsResults.i)
    }
    
    write.table(tmp1, file = file.path("results", resultFileName), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    
    if (elapsedTimes.i$simID %in% tmp2$simID){
      rowIdx <- which(tmp2$simID == elapsedTimes.i$simID)
      tmp2[rowIdx, ] <- elapsedTimes.i
    } else {
      tmp2 <- bind_rows(tmp2, elapsedTimes.i)
    }
    
    write.table(tmp2, file = file.path("results", elapsedFileName), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    cat(rep(".", 6), " Results saved... [OK]", sep = "", "\n")
  } else {
    write.table(simCombsResults.i, file = file.path("results", resultFileName), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(elapsedTimes.i, file = file.path("results", elapsedFileName), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(rep(".", 6), " Results saved... [OK]", sep = "", "\n")
  }
  
  simCombsResults <- rbind(simCombsResults, simCombsResults.i)
  elapsedTimes <- rbind(elapsedTimes, elapsedTimes.i)
}

end_all <- Sys.time()
elapsed <- end_all - start_all

difftime(end_all, start_all, units = "hours")

