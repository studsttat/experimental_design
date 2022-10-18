
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

contrast_ftest <- function(y,gr,coef,alpha = 0.05,
                           conf.int = TRUE){



  if(!conf.int == TRUE){
    # With confidence interval

    ## Mean of Each Group
    rbm <- aggregate(y,
                     by = list(Category = gr),
                     FUN = mean
    )

    ## Length of each group
    lpg <- aggregate(y,
                     by = list(Category = gr),
                     FUN = length
    )

    ## ANOVA Model

    fit = aov(y~gr)

    ## Size Parameters

    a = length(lpg$x)
    N = length(y)
    n = lpg[,2]
    dfg = a-1
    dfy = N-a
    dft = N-1

    ## Mean Square Error
    mse = sum((fit$residuals)^2)/dfy

    ## Mean Estimates and Numerator
    sum.mean = round(sum(coef*rbm[,2]),digits = 2)
    sq.sum = (sum.mean)^2

    ## Mean of Squares of Coefficients and denominator
    sum.sqc = sum((coef^2)/n)

    ## Test Statistic

    fval = round(sq.sum/(mse*sum.sqc),digits = 2)

    ## Computed p-value (f-distribution)
    pvalue = pf(fval,1,dfy,lower.tail = F)

    ## (1-alpha)% Confidence Interval

    lower = sum.mean - (qt((alpha)/2,N-a,lower.tail = F))*(mse*sum.sqc)
    upper = sum.mean + (qt((alpha)/2,N-a,lower.tail = F))*(mse*sum.sqc)
    CI = c(Lowerlvl = lower, Upperlvl = upper)

    ## Summary

    test.summary = round(c(Estimate = sum.mean, MSE = mse,
                           n = n,
                           N = N, a = a, dfy = dfy, dfg = dfg,
                           Fstatistic = fval, pvalue=pvalue, CI),
                         digits = 5)
    r = rbind(test.summary)
    return(r)


  } else
  {
    # Without confidence interval


    ## ANOVA Model

    fit = aov(y~gr)

    ## Mean of each group

    rbm <- aggregate(y,
                     by = list(Category = gr),
                     FUN = mean
    )

    ## Length of each group

    lpg <- aggregate(y,
                     by = list(Category = gr),
                     FUN = length
    )

    ## Size Parameters

    a = length(lpg$x)
    N = length(y)
    n = lpg[,2]
    dfg = a-1
    dfy = N-a
    dft = N-1

    ## MSE

    mse = sum((fit$residuals)^2)/dfy

    ## Mean Estimates and Numerator

    sum.mean = round(sum(coef*rbm[,2]),digits = 2)
    sq.sum = (sum.mean)^2

    ## Mean of Squares of Coefficients and denominator

    sum.sqc = sum((coef^2)/n)

    ## Test Statistic

    fval = round(sq.sum/(mse*sum.sqc),digits = 2)

    ## Computed p-value

    pvalue = pf(fval,1,dfy,lower.tail = F)

    ## Summary

    test.summary = round(c(Estimate = sum.mean, MSE = mse,
                           n=n,
                           N = N, a = a, dfy = dfy, dfg = dfg,
                           Fstatistic = fval, pvalue=pvalue),
                         digits = 5)
    r = rbind(test.summary)
    return(r)

  }

}





#######################################################################



