\name{scheffe_contrast}
\alias{scheffe_contrast}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{contrast_ftest
%%  ~~function to do ... ~~
}

\description{Multiple Comparison Test of ANOVA by Scheffe Method with Linear Combinations
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}


\usage{
scheffe_contrast(y,gr,coef1 = c(2,1,1), coef2 = c(1,1,1),
                             alpha = 0.05, conf.int = TRUE)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{response}
  \item{gr}{group factors}
  \item{coef1}{Second Linear combination coefficients}
  \item{coef2}{Second Linear combination coefficients}
  \item{alpha}{level of significance}
  \item{conf.int}{If TRUE, return confidence interval}
}

