\name{initial_parameter_training}
\alias{initial_parameter_training}

\title{
Algorithm to Find Plausible Starting Values for Parameter Estimation
}
\description{
The function computes plausible starting values for both the Baum-Welch algorithm and the algorithm for directly maximizing the log-Likelihood.  Plausible starting values can potentially diminish problems of (i) numerical instability and (ii) not finding the global optimum.
}
\usage{
initial_parameter_training(x, m, distribution_class, n = 100, 
                           discr_logL = FALSE, discr_logL_eps = 0.5)
}

\arguments{
  \item{x}{%%
a vector object containing the time-series of observations that are assumed to be realizations of the (hidden Markov state dependent) observation process of the model.}
  \item{m}{a (finite) number of states in the hidden Markov chain.}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the \code{m} observation distributions of the Markov dependent observation process.  The following distributions are supported by this algorithm: Poisson (\code{pois}); generalized Poisson (\code{genpois}); normal (\code{norm}); geometric (\code{geom}).}
  \item{n}{%%
a single numerical value specifying the number of samples to find the best starting value for the training algorithm.  Default value is \code{100}.}
  \item{discr_logL}{%%
a logical object.  \code{TRUE}, if the discrete log-likelihood shall be calculated (for \code{distribution_class="norm"} instead of the general log-likelihood.  Default is \code{FALSE}.}
  \item{discr_logL_eps}{%%
discrete log-likelihood for a hidden Markov model based on nomal distributions (for \code{distribution_class="norm"}).  The default value is \code{0.5}.}
}


\details{
From our experience, parameter estimation for long time-series of observations (\code{T>1000}) or observation values \code{>1500} tend to be numerical instable and does not necessarily find a global maximum.  Both problems can eventually be diminished with plausible starting values.  Basically, the idea behind \code{initial_parameter_training} is to sample randomly \code{n} sets of \code{m} observations from the time-series \code{x}, as means (\code{E}) of the state-dependent distributions. This \code{n} samplings of \code{E}, therefore induce \code{n} sets of parameters (\code{distribution_theta}) for the HMM without running a (slow) parameter estimation algorithm. Furthermore, \code{initial_parameter_training} calculates the log-Likelihood for all those \code{n} sets of parameters.  The set of parameters with the best Likelihood are outputted as plausible starting values.
(Additionally to the \code{n} sets of randomly chosen observations as means, the \code{m} quantiles of the observations are also checked as plausible means within this algorithm.)


}

\value{
\code{ initial_parameter_training } returns a list containing the following components:
  \item{m}{input number of states in the hidden Markov chain.}
  \item{k}{a single numerical value representing the number of parameters of the defined distribution class of the observation process.}
  \item{logL}{logarithmized likelihood of the model evaluated at the HMM with given starting values (\code{delta, gamma, distribution theta}) induced by \code{E}.}
  \item{E}{randomly choosen means of the observation time-series \code{x}, used for the observation distributions, for which the induced parameters 
  
  (\code{delta, gamma, distribution theta}) produce the largest Likelihood.}
  \item{distribution_theta}{a list object containing the plausible starting values for the parameters of the \code{m} observation distributions that are dependent on the hidden Markov state.}
  \item{delta}{a vector object containing plausible starting values for the marginal probability distribution of the \code{m} states of the Markov chain at the time point \code{t=1}.  At the moment:\cr \code{delta = rep(1/m, times=m)}.}
  \item{gamma}{a matrix (\code{nrow=ncol=m}) containing the plausible starting values for the transition matrix of the hidden Markov chain.  At the moment:\cr \code{gamma = 0.8 * diag(m) + rep(0.2/m, times=m)}.} 
}

\author{
Vitali Witowski (2013).
}

\seealso{
\code{\link{Baum_Welch_algorithm}}
\code{\link{direct_numerical_maximization}}
\code{\link{HMM_training}}
}

\examples{

################################################################
### Fictitious observations ####################################
################################################################

x <- c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,
  14,16,13,11,11,10,12,19,23,25,24,23,20,21,22,22,18,7,
  5,3,4,3,2,3,4,5,4,2,1,3,4,5,4,5,3,5,6,4,3,6,4,8,9,12,
  9,14,17,15,25,23,25,35,29,36,34,36,29,41,42,39,40,43,
  37,36,20,20,21,22,23,26,27,28,25,28,24,21,25,21,20,21,
  11,18,19,20,21,13,19,18,20,7,18,8,15,17,16,13,10,4,9,
  7,8,10,9,11,9,11,10,12,12,5,13,4,6,6,13,8,9,10,13,13,
  11,10,5,3,3,4,9,6,8,3,5,3,2,2,1,3,5,11,2,3,5,6,9,8,5,
  2,5,3,4,6,4,8,15,12,16,20,18,23,18,19,24,23,24,21,26,
  36,38,37,39,45,42,41,37,38,38,35,37,35,31,32,30,20,39,
  40,33,32,35,34,36,34,32,33,27,28,25,22,17,18,16,10,9,
  5,12,7,8,8,9,19,21,24,20,23,19,17,18,17,22,11,12,3,9,
  10,4,5,13,3,5,6,3,5,4,2,5,1,2,4,4,3,2,1) 


### Finding plausibel starting values for the parameter estimation 
### for a generealized-Pois-HMM with m=4 states
m <- 4 

## TIME CONSUMING:
\dontrun{
plausible_starting_values <- 
   initial_parameter_training(x = x, 
     m = m, 
     distribution_class = "genpois", 
     n=100)

print(plausible_starting_values)
}
}
\keyword{ iteration }
