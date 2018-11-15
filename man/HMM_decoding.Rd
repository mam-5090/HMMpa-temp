\name{HMM_decoding}
\alias{HMM_decoding}

\title{
Algorithm for Decoding Hidden Markov Models (local or global)
}
\description{
The function decodes a hidden Markov model into a most likely sequence of hidden states. Furthermore this function provides estimated observation values along the most likely sequence of hidden states. See Details for more information.
} 
\usage{
		HMM_decoding(x, m, delta, gamma, distribution_class, 
		     distribution_theta, decoding_method = "global", 
		     discr_logL = FALSE, discr_logL_eps = 0.5)
}

\arguments{
  \item{x}{%%
  	a vector object containing the time-series of observations that are assumed to be realizations of the (hidden Markov state dependent) observation process of the model.}
  \item{m}{%%
  	a (finite) number of states in the hidden Markov chain.}
  \item{delta}{%%
a vector object containing values for the marginal probability distribution of the \code{m} states of the Markov chain at the time point \code{t=1}.}
  \item{gamma}{%%
  	a matrix (\code{ncol=nrow=m}) containing values for the transition matrix of the hidden Markov chain.}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the \code{m} observation distributions of the Markov dependent observation process.  The following distributions are supported by this algorithm: Poisson (\code{pois}); generalized Poisson (\code{genpois}); normal (\code{norm}); geometric (\code{geom}).}
  \item{distribution_theta}{%%
a list object containing the parameter values for the \code{m} observation distributions that are dependent on the hidden Markov state.}
  \item{decoding_method}{%%
a string object to choose the applied decoding-method to decode the HMM given the time-series of observations \code{x}.  Possible values are \code{"global"} (for the use of the \code{Viterbi_algorithm}) and \code{"local"} (for the use of the 

\code{local_decoding_algorithm}). Default value is \code{"global"}.}
  \item{discr_logL}{%%
a logical object. It is \code{TRUE} if the discrete log-likelihood shall be calculated (for \code{distribution_class="norm"} instead of the general log-likelihood).  Default is \code{FALSE}.}
  \item{discr_logL_eps}{%%
a single numerical value to approximately determine the discrete log-likelihood for a hidden Markov model based on nomal distributions (for \code{"norm"}).  The default value is \code{0.5}.}
}
\details{
More precisely, the function works as follows:

\bold{Step 1:}
In a first step, the algorithm decodes a HMM into the most likely sequence of hidden states, given a time-series of observations. The user can choose between a global and a local approch.\cr
If \code{decoding_method="global"} is applied, the function calls \code{\link{Viterbi_algorithm}} to determine the sequence of most likely hidden states for all time points simultaneously.\cr 
If \code{decoding_method="local"} is applied, the function calls \code{\link{local_decoding_algorithm}} to determine the most likely hidden state for each time point seperately. 
 
\bold{Step 2:}
In a second step, this function links each observation to the mean of the distribution, that corresponds to the decoded state at this point in time. 
}
\value{

\code{ HMM_decoding } returns a list containing the following two components:
  \item{decoding_method}{%%
a string object indicating the applied decoding method.
}
  \item{decoding}{%%
a numerical vector containing the most likely sequence of hidden states as decoded by the \code{\link{Viterbi_algorithm}} (if \code{"global"} was applied) or by the 

\code{\link{local_decoding_algorithm}} (if \code{"local"} was applied).
}
  \item{decoding_distr_means}{%%
a numerical vector of estimated oberservation values along the most likely seuquence of hidden states (see \code{decoding} and Step 2).}
}

\references{
MacDonald, I. L.,  Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
}

\author{
Vitali Witowski (2013).
}

\seealso{
\code{\link{local_decoding_algorithm}},		
\code{\link{Viterbi_algorithm}}	
}
\examples{


################################################################
### i) HMM-training  ###########################################
################################################################

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

### Train hidden Markov model for different number of 
### states m=2,...,6 and select the optimal model
\donttest{
m_trained_HMM <- 
    HMM_training(x = x, 
       min_m = 2, 
       max_m = 6, 
       distribution_class="pois")$trained_HMM_with_selected_m
         
                 
################################################################
################################################################
### ii) Global decoding ########################################
################################################################
################################################################

### Decode the trained HMM using the Viterbi algorithm to get 
### the estimated sequence of hidden physical activity levels
global_decoding <- 
    HMM_decoding(x = x, 
        m = m_trained_HMM$m, 
        delta = m_trained_HMM$delta, 
        gamma = m_trained_HMM$gamma, 
        distribution_class = m_trained_HMM$distribution_class, 
        distribution_theta = m_trained_HMM$distribution_theta,
        decoding_method = "global")
          
### Globally most likely sequence of hidden states, 
### i.e. in this case sequence of activity levels
global_decoding$decoding

par(mfrow = c(1,1))
plot(global_decoding$decoding)


### Plot the observed impulse counts and the most likely 
### sequence (green) according to the Viterbi algorithm that 
### generated these observations
plot(x)
lines(global_decoding$decoding_distr_means, col = "green")

################################################################
################################################################ 
### iii) Local decoding ########################################
################################################################ 
################################################################

### Decode the trained HMM using the local decoding algorithm 
### to get the estimated sequence of hidden physical activity 
### levels
local_decoding <- 
   HMM_decoding(x = x, 
   m = m_trained_HMM$m, 
   delta = m_trained_HMM$delta, 
   gamma = m_trained_HMM$gamma, 
   distribution_class = m_trained_HMM$distribution_class, 
   distribution_theta = m_trained_HMM$distribution_theta,
   decoding_method = "local")
        

### Locally most likely sequence of hidden states, 
### i.e. in this case sequence of activity levels
local_decoding$decoding

par(mfrow=c(1,1))
plot(local_decoding$decoding)


### Plot the observed impulse counts and the most likely 
### sequence (green) according to the local decoding algorithm 
### that generated these observations
plot(x)
lines(local_decoding$decoding_distr_means, col = "red")

################################################################
################################################################
### iv) Comparison of global and local decoding ################
################################################################
################################################################

### Comparison of global decoding (green), local decoding (red) 
### and the connection to the closest mean (blue)
print(global_decoding$decoding)  
print(local_decoding$decoding)

### Plot comparison 
par(mfrow = c(2,2))
plot(global_decoding$decoding[seq(230,260)], col = "green", 
  ylab = "global decoding", main = "(zooming)")
  plot(x[seq(230,260)], ylab = "global decoding", 
  main = "(zooming x[seq(230,260)])")
lines(global_decoding$decoding_distr_means[seq(230,260)], 
  col = "green")
plot(local_decoding$decoding[seq(230,260)], col = "red", 
  ylab = "local decoding", main = "(zooming)")
plot(x[seq(230,260)], ylab = "local decoding", 
  main = "(zooming x[seq(230,260)])")
lines(local_decoding$decoding_distr_means[seq(230,260)], 
  col = "red")
par(mfrow = c(1,1))
}
}

