\name{fspe}
\alias{fspe}

\title{
Estimating number of factors with out-of-sample prediction error
}
\description{
Estimates the number of factors in Exploratory Factor Analysis (EFA) by approximating the out-of-sample prediction errors using a cross-validation scheme.
}
\usage{
fspe(data, maxK, nfold = 10, rep = 1, method = "PE",
     rotate = "oblimin", pbar = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A n x p matrix or data.frame.
}
  \item{maxK}{
The largest number of factors considered. For example, if \code{maxK = 8}, the factors \code{1, 2, ..., 8} are considered.
}
  \item{nfold}{
The number of folds in the cross-validation scheme. Defaults to \code{nfold = 10}.
}
  \item{rep}{
The number of repetitions of the cross-validation scheme. The repetiations are aggregated by averaging the prediction errors for each number of factor and subsequently selecting the number of factors that minimizes this average. Defaults to \code{rep = 1}.
}
  \item{method}{
If \code{method = "PE"}, the factor estimation is based on the prediction error on the variables; If \code{method = "Cov"} the error is computed on the covariance matrix. Defaults to \code{method = "PE"} since it has performed better in simulations (see Haslbeck & van Bork, 2021 below.)
}
  \item{rotate}{
  The totation being used in the underlying factor analysis models. The options are the same as in the \code{fa()} function in the \code{psych} package, which is called by \code{fpse}. Defaults to \code{rotate = "oblimin"}.
}
  \item{pbar}{
If \code{pbar = TRUE}, a progress bar is shown.
}
  \item{...}{
Additional arguments passed to the \code{fa()} function in the \code{psych} package.
}

}
\details{
The function returns:
}
\value{
  \item{nfactor}{An integer indicating the estimated number of factors.}
  \item{nfactor}{A 4-dimensional (for the PE method: variables, candidate models, folds, repetitions) or a 3-dimensional (for the Cov method: candidate models, folds, repetitions) array containing all prediction errors.}

}
\references{
Haslbeck & van Bork (2021). Estimating the number of factors in Exploratory Factor Analysis using out-of-sample prediction errors. PsyArXiv Preprint: https://psyarxiv.com/qktsd
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{
\donttest{

data(holzinger19)

fspe_out <- fspe(holzinger19,
                 maxK = 10,
                 nfold = 10,
                 method = "PE")

fspe_out$nfactor # estimated factors = 4

# Plot OoS PE path
plot(fspe_out$PEs, xlab="No. Factors", ylab="OoS PE")

}
}

