
# x = 1:50
# theta2 = estimationFun
# nboot = nboot.bca
# alpha = c(0.025, 0.975)
# nC = sum(nCase)
# data = dat
# nR = modelParams$sampleSizes$nr
# nT = modelParams$sampleSizes$nt

bcanon2 <- function (x, nboot, theta2, data = NULL, nR = NULL, nT = NULL, nC = NULL,
                     alpha = c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975), boot.reader = TRUE, ...){
  if (!all(alpha < 1) || !all(alpha > 0))
    stop("All elements of alpha must be in (0,1)")
  alpha_sorted <- sort(alpha)
  if (nboot <= 1/min(alpha_sorted[1], 1 - alpha_sorted[length(alpha_sorted)]))
    warning("nboot is not large enough to estimate your chosen alpha.")
  
  call <- match.call()
  n <- length(x)
  
  trueAUCs <- theta2(x, data = data, nR = nR, nT = nT, nC = nC, boot.reader = FALSE, ...)  ## TRUE AUC estiamtes.
  thetahat <- trueAUCs[1] - trueAUCs[2]

  bootsam <- matrix(sample(x, size = n * nboot, replace = TRUE), nrow = nboot)
  thetastar <- apply(bootsam, 1, theta2, data = data, nR = nR, nT = nT, nC = nC, boot.reader = boot.reader, return.diff = TRUE, ...)

  confpoints.perc <- quantile(thetastar, alpha)

  z0 <- qnorm(sum(thetastar < thetahat) / nboot)
  u <- rep(0, n)
  for (i in 1:n) {
    u[i] <- theta2(x[-i], data = data, nR = nR, nT = nT, nC = nC, boot.reader = boot.reader, return.diff = TRUE, ...)
  }
  uu <- mean(u) - u
  acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)

  zalpha <- qnorm(alpha)
  tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))

  confpoints <- quantile(x = thetastar, probs = tt)
  names(confpoints) <- NULL
  confpoints <- cbind(alpha, confpoints.perc, confpoints)
  dimnames(confpoints)[[2]] <- c("alpha", "percentile", "bca point")
  return(list(confpoints = confpoints, z0 = z0, acc = acc,
              u = u, call = call))
}
