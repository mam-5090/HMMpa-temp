\name{HMM_training}
\alias{HMM_training}

\title{
Training of Hidden Markov Models
}
\description{
Function to estimate the model specific parameters (\code{delta, gamma, distribution_theta}) for a hidden Markov model, given a time-series and a user-defined distribution class. Can also be used for model selection (selecting the optimal number of states \code{m}). See Details for more information. 
}
\usage{
HMM_training (x, distribution_class, min_m = 2, max_m = 6, 
              n = 100, training_method = "EM", discr_logL = FALSE, 
              discr_logL_eps = 0.5, Mstep_numerical = FALSE, 
              dynamical_selection = TRUE, BW_max_iter = 50, 
              BW_limit_accuracy = 0.001, BW_print = TRUE,
              DNM_max_iter = 50, DNM_limit_accuracy = 0.001, 
              DNM_print = 2)
}

\arguments{
  \item{x}{%%
a vector object of length \code{T} containing observations of a time-series \code{x}, which are assumed to be realizations of the (hidden Markov state dependent) observation process of the HMM.}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the $m$ observation distributions of the Markov dependent observation process.  The following distributions are supported:  Poisson (\code{pois}); generalized Poisson (\code{genpois}, only available for \code{training_method="numerical"}); normal (\code{norm})).}
  \item{min_m}{minimum number of hidden states in the hidden Markov chain.  Default value is \code{2}.}
  \item{max_m}{maximum number of hidden states in the hidden Markov chain.  Default value is \code{6}.}
  \item{n}{%%
a single numerical value specifying the number of samples to find the best starting values for the training algorithm.  Default value is \code{n=100}.}
  \item{training_method}{%%
a logical value indicating whether the Baum-Welch algorithm (\code{"EM"}) or the method of direct numerical maximization (\code{"numerical"}) should be applied for estimating the model specific parameters. See \code{\link{Baum_Welch_algorithm}} and \code{\link{direct_numerical_maximization}} for further details.  

Default is \code{training_method="EM"}.}
  \item{discr_logL}{%%
a logical object.  Default is \code{FALSE} for the general log-likelihood, \code{TRUE} for the discrete log-likelihood (for \code{distribution_class="norm"}).
}
  \item{discr_logL_eps}{%%
a single numerical value, used to approximate the discrete log-likelihood for a hidden Markov model based on nomal distributions (for \code{"norm"}).  The default value is \code{0.5}.
}
  \item{Mstep_numerical}{%%
a logical object indicating whether the Maximization Step of the Baum-Welch algorithm should be performed by numerical maximization. Default is \code{FALSE}.}
  \item{dynamical_selection}{%%
a logical value indicating whether the method of dynamical initial parameter selection should be applied (see Details).  Default is \code{TRUE}.}
  \item{BW_max_iter}{%%
a single numerical value representing the maximum number of iterations in the Baum-Welch algorithm.  Default value is \code{50}.}
  \item{BW_limit_accuracy}{%%
a single numerical value representing the convergence criterion of the Baum-Welch algorithm. Default value is is \code{0.001}.}
  \item{BW_print}{%%
a logical object indicating whether the log-likelihood at each iteration-step shall be printed. Default is \code{TRUE}.
}
  \item{DNM_max_iter}{%%
a single numerical value representing the maximum number of iterations of the numerical maximization using the nlm-function (used to perform the Maximization Step of the Baum-Welch-algorithm).  Default value is \code{50}.}
  \item{DNM_limit_accuracy}{%%
a single numerical value representing the convergence criterion of the numerical maximization algorithm using the \link[stats]{nlm} function (used to perform the Maximization Step of the Baum-Welch- algorithm).  Default value is \code{0.001}.} 
  \item{DNM_print}{%%
a single numerical value to determine the level of printing of the \code{nlm}-function.  See \code{nlm}-function for further informations. The value \code{0} suppresses, that no printing will be outputted. Default value is \code{2} for full printing.}
}


\details{
More precisely, the function works as follows:

\bold{Step 1:}
In a first step, the algorithm estimates the model specific parameters for different values of \code{m} (indeed for \code{min_m,...,max_m}) using either the function \code{\link{Baum_Welch_algorithm}} or 

\code{\link{direct_numerical_maximization}}. Therefore, the function first searches for plausible starting values by using the function \code{\link{initial_parameter_training}}. 

\bold{Step 2:}
In a second step, this function evaluates the AIC and BIC values for each HMM (built in Step 1) using the functions \code{\link{AIC_HMM}} and \code{\link{BIC_HMM}}. Then, based on that values, this function decides for the most plausible number of states \code{m} (respectively for the most appropriate HMM for the given time-series of observations). In case when AIC and BIC claim for a different \code{m}, the algorithm decides for the smaller value for \code{m} (with the background to have a more simplistic model). 
If the user is intereseted in having a HMM with a fixed number for \code{m}, \code{min_m} and \code{max_m} have to be chosen equally (for instance \code{min_m=4=max_m} for a HMM with \code{m=4} hidden states).

To speed up the parameter estimation for each \eqn{m > m_min}, the user can choose the method of dynamical initial parameter selection.

If the method of dynamical intial parameter selection \bold{is not applied}, the function 

\code{\link{initial_parameter_training}} will be called to find plausible starting values for each state 
\eqn{ m \in \{min_m, \ldots, max_m\}}{ min_m \le m \le max_m}. \cr

If the method of dynamical intial parameter selection \bold{is applied}, then starting parameter values using the function \code{\link{initial_parameter_training}} will be found only for the first HMM (respectively the HMM with \code{m_min} states). The further starting parameter values for the next HMM (with \code{m+1} states and so on) are retained from the trained parameter values of the last HMM (with \code{m} states and so on). 

}
\value{
\code{HMM_training} returns a list containing the following components:

  \item{trained_HMM_with_selected_m}{%%
a list object containg the key data of the optimal trained HMM (HMM with selected \code{m}) -- summarized output of the \code{\link{Baum_Welch_algorithm}} or 

\code{\link{direct_numerical_maximization}} algorithm, respectively.}
  \item{list_of_all_initial_parameters}{%%
a list object containing the plausible starting values for all HMMs (one for each state \code{m}).}
  \item{list_of_all_trained_HMMs}{%%
a list object containing all trained m-state-HMMs. See \code{\link{Baum_Welch_algorithm}} or \code{\link{direct_numerical_maximization}} for \code{training_method="EM"} or 

\code{training_method="numerical"}, respectively.}
  \item{list_of_all_logLs_for_each_HMM_with_m_states}{%%
a list object containing all logarithmized Likelihoods of each trained HMM.}
  \item{list_of_all_AICs_for_each_HMM_with_m_states}{%%
a list object containing the AIC values of all trained HMMs.}
  \item{list_of_all_BICs_for_each_HMM_with_m_states}{%%
a list object containing the BIC values of all trained HMMs.}
  \item{model_selection_over_AIC}{is logical.  \code{TRUE}, if model selection was based on AIC and \code{FALSE}, if model selection was based on BIC.}
}

\references{
MacDonald, I. L.,  Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
}

\author{
Vitali Witowski (2013).
}

\seealso{
  \code{\link{initial_parameter_training}},
  \code{\link{Baum_Welch_algorithm}},
  \code{\link{direct_numerical_maximization}},
  \code{\link{AIC_HMM}},
  \code{\link{BIC_HMM}}
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


## Train a poisson hidden Markov model using the Baum-Welch 
## algorithm for different number of states m=2,...,6

## TIME CONSUMING:
\dontrun{
trained_HMMs <- 
    HMM_training(x = x, 
      distribution_class = "pois", 
      min_m = 2, 
      max_m = 6, 
      training_method = "EM")


## Various output values for the HMM
names(trained_HMMs)

## Print details of the most plausible HMM for the given 
## time-series of observations
print(trained_HMMs$trained_HMM_with_selected_m)

## Print details of all trained HMMs (by this function) 
## for the given time-series of observations
print(trained_HMMs$list_of_all_trained_HMMs)

## Print the BIC-values of all trained HMMs for the given 
## time-series of observations  
print(trained_HMMs$list_of_all_BICs_for_each_HMM_with_m_states)

## Print the logL-values of all trained HMMs for the 
## given time-series of observations  
print(trained_HMMs$list_of_all_logLs_for_each_HMM_with_m_states)
}
}
\keyword{ iteration }
