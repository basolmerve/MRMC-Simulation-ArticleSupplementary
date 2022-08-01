# Clear the cache.
rm(list = ls())

# Load the libraries and the source codes.
library(parallel)
library(pROC)
library(dplyr)
source("R/rocdata.R")
source("R/pseudo_values.R")
source("R/simDataDBM.R")
source("R/covOR.R")
source("R/estimationFunBWC.R")
source("R/bcanon2.R")
source("R/boott.R")
source("R/MM_functions.R")


## Read data
# data <- read.table("datasets/prostate_2.txt", header = TRUE, sep = ",")
data <- read.table("data/mammogram.txt", header = TRUE, sep = "\t")

## Variables MUST be in {"Cases", "Status", "Reader", "Test", "Score"} order.
## Otherwise, change the order of columns of ROC data as above.
head(data)

length(unique(data$Reader))
length(unique(data$Test))

tmp <- data[data$Reader == 1 & data$Test == 1, ]
table(tmp$Status)

# Define data properties
nR <- 4  # Number of readers
# nD <- nH <- 100  # Number of healthy and diseased subjects.
nD <- 100
nH <- 100
n <- nH + nD  # Total sample size
nT <- 2   # Number of tests


# Fit Models
####### 1. DBM ######
dat <- data  # We use a copy data to keep original one unchanged.

# Split data via treatment:reader combinations
# reader ve test'lerin factor hale getirilmesi
dat[ ,"Reader"] <- as.factor(dat[ ,"Reader"])
dat[ ,"Test"] <- as.factor(dat[ ,"Test"])
dat[ ,"Cases"] <- as.factor(dat[ ,"Cases"])
dat[ ,"Status"] <- as.factor(dat[ ,"Status"])

data_split <- split(dat, dat[ ,"Reader"]:dat[ ,"Test"])

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

result_DBM <- data.frame(Healthy = nH, Diseased = nD,
                     Reader = nR, Test = nT,
                     Model = "DBM", AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                     pValue = pValue, LowerLimit = lowerLimit, UpperLimit = upperLimit)

####### 2. OR ######
dat <- dataCov <- data

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

result_OR <- data.frame(Healthy = nH, Diseased = nD,
                     Reader = nR, Test = nT,
                     Model = "OR",
                     AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                     pValue = pValueOR, LowerLimit = lowerLimitOR, UpperLimit = upperLimitOR)

####### 3. MM ######
# For all these analyses, the data should be in the same format. There are 5 columns. The first has the ratings or scores 
# from the first treatment, the second is the second treatment.
# The third column is disease status, fourht is case, and 5th is reader
# Disease status, case, and reader should all be sorted low to high
# (eg, start with ratings from the first reader, non-diseased patients, case 1 to number non-diseased etc).
# Also, cases must be numbered starting from 1.
# The columns are:
#   
#   Treat1	Treat2	Disease	Case	Reader
dat <- data

# This function is used to rearrange variables (columns of data frame) to make MM data appropriate to use
# within 'gen.data' function. 
# 'gen.data' requires data matrix with columns in {"ScoresForTest1", "ScoresForTest2", "Status", "Cases", "Reader"} order.
Y <- transformData_zhou(dat, reader = "Reader", test = "Test", 
                        score = "Score", case = "Cases", status = "Status")    

BIdata <- gen.data(Y, nH, nD, nR, ordinal = TRUE)  # Generate appropriate data for MM model.

response <- c(BIdata[ ,1], BIdata[ ,2])
len <- length(BIdata[ ,1])
predictor.1 <- c(rep(0, len), rep(1, len))
lin.mod.1 <- glm(response ~ predictor.1, family = binomial(link = "probit"))
beta0 = lin.mod.1$coefficients[1]
beta1 = lin.mod.1$coefficients[2]
beta = c(beta0, beta1, beta0 + beta1)
xbeta <- c(beta0, beta0 + beta1)

AUC = c(pnorm(beta0), pnorm(beta0 + beta1), 1)
AUC[3] = AUC[1] - AUC[2]
se <- matrix(c(0, 0, 0), 3, 1)
se_AUC <- matrix(c(0, 0, 0), 3, 1)
temp1 = dnorm(beta0)
temp2 = dnorm(beta0 + beta1)

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

result_MM <- data.frame(Healthy = nH, Diseased = nD,
                     Reader = nR, Test = nT,
                     Model = "MM",
                     AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                     pValue = NA, LowerLimit = lowerLimitMM, UpperLimit = upperLimitMM)

####### 4. BWC ######
nboot.bca <- 2000   # number of bootstrap samples
dat <- data

set.seed(1)  # For reproducble bootstrap results.
## BCa güven aralığı hesaplamaları.
bcaRes <- try({bcanon2(x = 1:n, theta2 = estimationFunBWC, nboot = nboot.bca, alpha = c(0.025, 0.975), nC = n,
                       data = dat, nR = nR, nT = nT, boot.reader = TRUE)})

if (class(bcaRes) == "try-error"){
  bcaRes <- try({bcanon2(x = 1:n, theta2 = estimationFunBWC, nboot = nboot.bca, alpha = c(0.025, 0.975), nC = n,
                         data = dat, nR = modelParams$nr, nT = modelParams$nt, boot.reader = TRUE)})
  if (class(bcaRes) == "try-error"){
    bcaRes <- NULL
  }
}

est <- try({estimationFunBWC(x = 1:n, data = dat, nR = nR, nT = nT,
                             nC = n, boot.reader = FALSE, return.diff = FALSE)})

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
                         Reader = nR, Test = nT,
                         Model = "BCa", 
                         AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                         pValue = NA, LowerLimit = lowerBCa, UpperLimit = upperBCa)

result_perc <- data.frame(Healthy = nH, Diseased = nD,
                          Reader = nR, Test = nT,
                          Model = "Percentile",
                          AUC1 = AUC1, AUC2 = AUC2, AUC.Diff. = auc.diff,
                          pValue = NA, LowerLimit = lowerPerc, UpperLimit = upperPerc)

result_BWC <- rbind(result_bca, result_perc)

# Combine results from all models into single matrix.
results <- rbind(result_DBM, result_OR, result_MM, result_BWC)
results[ ,"LowerLimit"] <- round(results[ ,"LowerLimit"], 3)
results[ ,"UpperLimit"] <- round(results[ ,"UpperLimit"], 3)
rownames(results) <- NULL
results


## Plot Codes
library(ggplot2)

plotData <- mutate(results, methodType = factor(rep(c("parametric", "bwc"), c(3, 2))))

plotData <- plotData %>% 
  mutate(
    Model = factor(Model, levels = c("DBM", "OR", "MM", "BCa", "Percentile"))
  )

p <- ggplot(plotData, aes(x = Model, ymin = LowerLimit, ymax = UpperLimit, colour = methodType, linetype = methodType)) + 
  geom_point(aes(y = LowerLimit), size = 2) + 
  geom_point(aes(y = UpperLimit), size = 2) + 
  geom_errorbar(width = 0, size = 0.5) + 
  theme_bw(base_size = 12) + 
  xlab("Methods") + 
  ylab("95% Confidence interval of \n difference in diagnostic accuracies") + 
  ylim(c(-0.05, 0.10)) + 
  geom_hline(yintercept = 0, linetype = 3, col = "gray50") + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face = "plain", margin = margin(0, 0, 5, 0), hjust = 0.5),
        axis.title.x = element_text(margin = margin(5,0,5,0), face = "plain"),
        axis.title.y = element_text(margin = margin(0,5,0,5), face = "plain"),
        axis.text.x = element_text(margin = margin(5,0,0,0), face = "plain"),
        axis.text.y = element_text(margin = margin(0,5,0,0), face = "plain"),
        legend.position = "right") + 
  scale_colour_manual(values = c("gray50", "gray30")) + 
  scale_linetype_manual(values = c(1, 2), labels = c("BWC", "Model-based")) + 
  guides(linetype = guide_legend(title = "", keywidth = 2), colour = "none")

print(p)

# Save Figures
# cairo_pdf(file = "figure/CIs.pdf", height = 4.5, width = 6.5, pointsize = 12)
# p
# dev.off()
# 
# ggsave(filename = "figure/CI_figure.eps", plot = p, width = 6.5, height = 4.5, pointsize = 12)
# 
# png(filename = "figure/CIs.png", width = 3250, height = 2250, res = 500, pointsize = 12)
# p
# dev.off()










