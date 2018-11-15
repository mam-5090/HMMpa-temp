\name{dgenpois}
\alias{dgenpois}
\alias{pgenpois}
\alias{rgenpois}

\title{
The Generalized Poisson Distribution
}
\description{
Density, distribution function and random generation function for the generalized Poisson distribution. 
}
\usage{
dgenpois(x, lambda1, lambda2)
pgenpois(q, lambda1, lambda2)
rgenpois(n, lambda1, lambda2)
}

\arguments{
  \item{x}{a vector object of (non-negative integer) quantiles.}
  \item{q}{a numeric value.}
  \item{n}{number of random values to return.}
  \item{lambda1}{a single numeric value for parameter \code{lambda1} with \eqn{lambda1 > 0}.}
  \item{lambda2}{a single numeric value for parameter \code{lambda2} with \eqn{0 \le lamdba2 < 1}.  When \code{lambda2=0}, the generalized Poisson distribution reduces to the Poisson distribution.}
}
\details{
The generalized Poisson distribution has the density

\deqn{ p(x) = \lambda_1 (\lambda_1 + \lambda_2 \cdot x)^{x-1} \frac{ \exp(-\lambda_1-\lambda_2 \cdot x) )}{x!}}{%
       p(x) = lambda1 (lambda1 + lambda2 x)^(x-1)  exp(-lambda1-lambda2 x) )/x!}
for \eqn{x = 0,1,2,\ldots},b

with \eqn{\mbox{E}(X)=\frac{\lambda_1}{1-\lambda_2}}{E(x)=lambda1/(1-lambda2)} and variance \eqn{\mbox{var}(X)=\frac{\lambda_1}{(1-\lambda_2)^3}}{var(x)=lambda1/(1-lambda2)^3}.

}
\value{
\code{dgenpois} gives the density, \code{pgenpois} gives the distribution function and \code{rgenpois} generates random deviates.
}

\references{
Joe, H., Zhu, R. (2005). Generalized poisson distribution: the property of mixture of poisson and comparison with negative binomial distribution. Biometrical Journal \bold{47}(2):219--229. 
}

\author{
Based on Joe and Zhu (2005). Implementation by  Vitali Witowski (2013).
}

\seealso{
\link{Distributions} for other standard distributions, including \code{\link{dpois}} for the Poisson distribution.
}

\examples{
dgenpois(x = seq(0,20), lambda1 = 10, lambda2 = 0.5) 

pgenpois(q = 5, lambda1 = 10, lambda2 = 0.5) 

hist(rgenpois(n = 1000, lambda1 = 10, lambda2 = 0.5) )

}
\keyword{ distribution }
