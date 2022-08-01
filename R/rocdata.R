rocdata <- function(status, marker, event, higherValuesDiseased){
  
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: status - labels classifying subject status
  #            marker - values of each observation
  #            event  - disease value
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval

    
  if (length(marker) != length(status)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(unique(status)) != 2) {
    stop("There must only be 2 values for the status")
  }
  
  cut <- unique(marker)  # possible cutoff points
  mrk.nd = marker[status != event]  ## markers for non-diseased
  mrk.d = marker[status == event]  ## markers for diseased
  
  if (higherValuesDiseased){
    tp <- sapply(cut, function(x) length(which(marker >= x & status == event)))
    fn <- sapply(cut, function(x) length(which(marker < x & status == event)))
    fp <- sapply(cut, function(x) length(which(marker >= x & status != event)))
    tn <- sapply(cut, function(x) length(which(marker < x & status != event)))
  }                                            
  
  if (!higherValuesDiseased){
    tp <- sapply(cut, function(x) length(which(marker <= x & status == event)))
    fn <- sapply(cut, function(x) length(which(marker > x & status == event)))
    fp <- sapply(cut, function(x) length(which(marker <= x & status != event)))
    tn <- sapply(cut, function(x) length(which(marker > x & status != event)))
  }
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)

  roc = data.frame(CutOff = cut, FPR = fpr, TPR = tpr)
  roc <- roc[order(roc$FPR, roc$TPR),]
  
  roc = rbind(c(roc[1,"Marker"], Inf, 0, 0), roc)
  roc = rbind(roc, c(roc[1,"Marker"], -Inf, 1, 1))
  
  i <- 2:nrow(roc)
  auc <- ((roc$FPR[i] - roc$FPR[i - 1]) %*% (roc$TPR[i] + roc$TPR[i - 1]))/2
  
  return(auc)
}
