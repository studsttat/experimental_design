\name{contrast_ftest}
\alias{contrast_ftest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post - Hoc F Test for multiple treatments
%%  ~~function to do ... ~~
}


\description{Multiple Comparison Test of ANOVA with Linear Combinations
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}



\usage{
contrast_ftest(y,gr,coef = c(2,1,1),alpha = 0.05,
                           conf.int = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{response}
  \item{gr}{group factors}
  \item{coef}{Linear combination coefficients}
  \item{alpha}{level of significance}
  \item{conf.int}{If TRUE, return confidence interval}
}

