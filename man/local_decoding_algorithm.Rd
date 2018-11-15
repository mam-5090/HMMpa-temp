\name{local_decoding_algorithm}
\alias{local_decoding_algorithm}

\title{
Algorithm for Decoding Hidden Markov Models (local)
}
\description{
The function decodes a hidden Markov model into a most likely sequence of hidden states. Different to the \code{\link{Viterbi_algorithm}}, this algorithm determines the most likely hidden state for each time point seperately.}

\usage{
local_decoding_algorithm(x, m, delta, gamma, distribution_class, 
      distribution_theta, discr_logL = FALSE, discr_logL_eps = 0.5)
}


\arguments{
  \item{x}{a vector object containing the time-series of observations that are assumed to be realizations of the (hidden Markov state dependent) observation process of the model.}
  \item{m}{a (finite) number of states in the hidden Markov chain.}
  \item{delta}{%%
a vector object containing values for the marginal probability distribution of the \code{m} states of the Markov chain at the time point \code{t=1}.}
  \item{gamma}{a matrix (\code{ncol=nrow=m}) containing values for the transition matrix of the hidden Markov chain.}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the \code{m} observation distributions of the Markov dependent observation process.  The following distributions are supported by this algorithm: Poisson (\code{pois}); generalized Poisson (\code{genpois}); normal (\code{norm}); geometric (\code{geom}).}
  \item{distribution_theta}{%%
a list object containing the parameter values for the \code{m} observation distributions that are dependent on the hidden Markov state.}
  \item{discr_logL}{%%
a logical object. It is \code{TRUE} if the discrete log-likelihood shall be calculated (for \code{distribution_class="norm"} instead of the general log-likelihood).  Default is \code{FALSE}.}
  \item{discr_logL_eps}{%%
a single numerical value to approximately determine the discrete log-likelihood for a hidden Markov model based on nomal distributions (for \code{"norm"}).  The default value is \code{0.5}.}
}

\value{
\code{ local_decoding_algorithm } returns a list containing the following two components:
  \item{state_probabilities}{%%
a (T,m)-matrix (when T indicates the length/size of the observation time-series and m the number of states of the HMM) containing probabilities (conditional probability of a state i=1,...,m at a time point t=1,...,T given all observations x) calculated by the algorithm. See MacDonald & Zucchini (2009, Paragraph 5.3.1) for further details.}
  \item{decoding}{%%
a numerical vector containing the locally most likely sequence of hidden states as decoded by the local_decoding_algorithm.}

}
\references{
MacDonald, I. L.,  Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
}
\author{
The basic algorithm for a Poisson-HMM can be found in MacDonald & Zucchini (2009, Paragraph A.2.6).  Extension and implementation by Vitali Witowski (2013).
}



\seealso{
\code{\link{Viterbi_algorithm}},
\code{\link{HMM_decoding}}		
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

### Train hidden Markov model for m=4
\donttest{
m_trained_HMM <- 
    HMM_training(x = x, 
      min_m = 4, 
      max_m = 4, 
      distribution_class = "pois")$trained_HMM_with_selected_m

### Decode the trained HMM using the local-decoding algorithm 
### to get the locally most likely sequence of hidden states 
### for the time-series of observations
local_decoding <- 
    local_decoding_algorithm(x = x, 
       m = m_trained_HMM$m, 
       delta = m_trained_HMM$delta, 
       gamma = m_trained_HMM$gamma, 
       distribution_class = m_trained_HMM$distribution_class, 
       distribution_theta = m_trained_HMM$distribution_theta)

### Most likely sequence of hidden states
print(local_decoding$decoding)
plot(local_decoding$decoding)
}
}

