clearConsole <- function(){
  # cat(rep("\n", 100))
  cat("\014")
}

wsPre <- function(n = 2, ...){
  paste0(rep(" ", n), collapse = "")
}

printSimulationHeader <- function(step = 1, total = 10, equal.var = TRUE, seed = NULL, id = NULL, ...){
  cat("\n", " Running simulations...", "\n", sep = "")
  cat(wsPre(...), "- Step:", step, "/", total, paste0("(", round(100 * (step - 1) / total, 2),"% completed.)", collapse = ""), "\n")
  cat(wsPre(...), "- Simulation ID:", id, "\n")
  cat(wsPre(...), "- Equal variance components:", as.character(equal.var), "\n")
  if (!is.null(seed) & is.numeric(seed)){
    cat(wsPre(...), "- Seed (Random number generation):", seed, "\n")
  }
}

printSimulationInfo <- function(.info, equal.var = TRUE, .ndatasets = NULL, .nbootstrap = NULL, ...){
  cat("\n", " Simulation scenario...", "\n", sep = "")
  cat(wsPre(...), "- # mu:", .info$mu, "\n")
  cat(wsPre(...), "- # samples:", .info$nhd * 2, "(i.e.,", .info$nhd, "healthy and", .info$nhd, "diseased.)", "\n")
  cat(wsPre(...), "- # readers:", .info$nr, "\n")
  cat(wsPre(...), "- # tests:", .info$nt, "\n")
  if (!is.null(.ndatasets) & is.numeric(.ndatasets)){
    cat(wsPre(...), "- # generated data sets:", .ndatasets, "\n")
  }
  if (!is.null(.nbootstrap) & is.numeric(.nbootstrap)){
    cat(wsPre(...), "- # bootstrap samples:", .nbootstrap, "\n")
  }
}

printVarianceComponents <- function(.info, equal.var = TRUE, ...){
  dots <- list(...)
  n <- dots$n
  if (is.null(n)){
    n <- as.list(args(wsPre))[["n"]]
  }
  cat("\n", " Variance components...", "\n")
  if (equal.var){
    cat(wsPre(n), "- Correlation structure:", .info$scen, "\n")
    cat(wsPre(n), "- Var(Case):", .info$var_case, "\n")
    cat(wsPre(n), "- Var(Case*Test):", .info$var_casetest, "\n")
    cat(wsPre(n), "- Var(Case*Reader):", .info$var_readercase, "\n")
    cat(wsPre(n), "- Var(Reader):", .info$var_reader, "\n")
    cat(wsPre(n), "- Var(Reader*Test):", .info$var_readertest, "\n")
    cat(wsPre(n), "- Var(Error):", .info$var_error, "\n")
  } else {
    cat(wsPre(n + 1), "(1) Healthy subjects", "\n")
    cat(wsPre(n + 3), "- Correlation structure:", .info$scen, "\n")
    cat(wsPre(n + 3), "- Var(Case):", .info$var_case0, "\n")
    cat(wsPre(n + 3), "- Var(Case*Test):", .info$var_casetest0, "\n")
    cat(wsPre(n + 3), "- Var(Case*Reader):", .info$var_readercase0, "\n")
    cat(wsPre(n + 3), "- Var(Reader):", .info$var_reader, "\n")
    cat(wsPre(n + 3), "- Var(Reader*Test):", .info$var_readertest, "\n")
    cat(wsPre(n + 3), "- Var(Error):", .info$var_error0, "\n\n")
    
    cat(wsPre(n + 2), "(2) Diseased subjects", "\n")
    cat(wsPre(n + 4), "- Correlation structure:", .info$scen, "\n")
    cat(wsPre(n + 4), "- Var(Case):", .info$var_case1, "\n")
    cat(wsPre(n + 4), "- Var(Case*Test):", .info$var_casetest1, "\n")
    cat(wsPre(n + 4), "- Var(Case*Reader):", .info$var_readercase1, "\n")
    cat(wsPre(n + 4), "- Var(Reader):", .info$var_reader, "\n")
    cat(wsPre(n + 4), "- Var(Reader*Test):", .info$var_readertest, "\n")
    cat(wsPre(n + 4), "- Var(Error):", .info$var_error1, "\n")
  }
}

