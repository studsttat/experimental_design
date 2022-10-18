\name{scheffe_contrast}
\alias{scheffe_contrast}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post - Hoc Test with Scheffe Method
%%  ~~function to do ... ~~
}

\description{Multiple Comparison Test of ANOVA by Scheffe Method with Linear Combinations
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}


\usage{
scheffe_contrast(y,gr,coef1 = c(2,1,1), coef2 = NULL,
                             alpha = 0.05, conf.int = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ response}
  \item{gr}{ group factors}
  \item{coef1}{ First Linear combination coefficients}
  \item{coef2}{ Second Linear combination coefficients. Defaults at NULL}
  \item{alpha}{ level of significance}
  \item{conf.int}{ If TRUE, return confidence interval}

}

\examples{
st <- c(1.53,1.61,3.75,2.89,3.26,
        3.89,3.68,5.70,5.62,5.79,
        8.18,5.64,7.36,8.82,7.10)
gr <- factor(rep(c("A", "B", "C"), c(5,5,5)))
data_rabbit <- data.frame(st,gr)
colnames(data_rabbit) <- c("Treatment", "Stimulant")
scheffe_contrast(data_rabbit$Treatment,data_rabbit$Stimulant,coef1 = c(2,1,1))
}
