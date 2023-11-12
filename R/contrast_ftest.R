
#' Multiple Comparison F Test
#'
#' @param y Vector of Sample
#' @param gr Factor of Group
#' @param coef
#' @param alpha Level of significance; defaults to 0.05
#' @param conf.int If TRUE, the (1-alpha)% Confidence Interval. Default is FALSE
#'
#'
#' @return A summary table from the result
#' @export
#'
#' @examples contrast.ftest(y = Treatment,gr = Stimulant,coef, alpha = 0.05,conf.int = TRUE)
#' @examples

contrast_ftest <- function(y, gr, coef, alpha = 0.05, conf.int = FALSE) {
  ## Mean of Each Group
  rbm <- aggregate(y, by = list(Category = gr), FUN = mean)

  ## Length of each group
  lpg <- aggregate(y, by = list(Category = gr), FUN = length)

  ## Fit ANOVA Model
  fit <- aov(y ~ gr)

  ## Other Parameters
  a <- length(levels(gr))
  N <- length(y)
  n <- lpg[, 2]
  dfg <- a - 1
  dfy <- N - a
  dft <- N - 1

  ## Mean Square Error
  mse <- sum((fit$residuals)^2) / dfy

  ## Mean Estimates and Numerator
  sum.mean <- sum(coef * rbm[, 2])

  ## Mean of Squares of Coefficients and denominator
  sum.sqc <- sum((coef^2) / n)

  ## Test Statistic
  fval <- sum.mean^2 / (mse * sum.sqc)

  ## Computed p-value (from f-distribution)
  pvalue <- pf(fval, 1, dfy, lower.tail = FALSE)

  ## Evaluate Confidence Interval
  if (conf.int) {
    ## (1-alpha)% Confidence Interval
    t_val <- qt(1 - alpha / 2, dfy)
    margin <- t_val * sqrt(mse * sum.sqc)

    lower <- sum.mean - margin
    upper <- sum.mean + margin
    CI <- c(Lowerlvl = lower, Upperlvl = upper)

    ## Summary
    test.summary <- c(
      Estimate = sum.mean,
      MSE = mse,
      n = n,
      means = rbm[, 2],
      N = N,
      a = a,
      dfy = dfy,
      dfg = dfg,
      Fstatistic = fval,
      pvalue = pvalue,
      CI
    )
    r <- data.frame(t(test.summary))
  } else {
    test.summary <- c(
      Estimate = sum.mean,
      MSE = mse,
      n = n,
      means = rbm[, 2],
      N = N,
      a = a,
      dfy = dfy,
      dfg = dfg,
      Fstatistic = fval,
      pvalue = pvalue
    )
    r <- data.frame(t(test.summary))
  }

  return(r)
}





#######################################################################


#' @export

print.ftest <- function(x,...){
  df1 <- x$dfg
  df2 <- x$dfy

  cat("Degrees of freedom by group:  ", df1, "\n")
  cat("Degrees of freedom of residuals:  ", df2, "\n")
}

