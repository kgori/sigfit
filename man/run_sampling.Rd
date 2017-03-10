\name{run_sampling}
\alias{run_sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run_sampling
}
\description{
Generates samples using Stan
}
\usage{
run_sampling(counts, signatures, prior = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{signatures}{
%%     ~~Describe \code{signatures} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (counts, signatures, prior = NULL, ...) 
{
    if (is.null(prior)) {
        prior = rep(1, ncol(signatures))
    }
    stopifnot(length(counts) == nrow(signatures))
    stopifnot(length(prior) == ncol(signatures))
    dat = list(C = length(counts), S = ncol(signatures), counts = counts, 
        signatures = signatures, alpha = prior)
    sampling(model = stan_models$multinomial, data = dat, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
