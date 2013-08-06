\name{HMM_based_method}
\alias{HMM_based_method}

\title{
Hidden Markov Method for Predicting Physical Activity Patterns
}
\description{
This function assigns a physical activity range to each observation of a time-series (such as a sequence of impulse counts recorded by an accelerometer) using hidden Markov models (HMM). The activity ranges are defined by thresholds called cut-off points.  Basically, this function combines \code{\link{HMM_training}}, \code{\link{HMM_decoding}} and \code{\link{cut_off_point_method}}. See Details for further information.
}
\usage{
HMM_based_method(x, cut_points, distribution_class, 
                 min_m = 2, max_m = 6, n = 100,
                 max_scaled_x = NA, names_activity_ranges = NA,  
                 discr_logL = FALSE, discr_logL_eps = 0.5, 
                 dynamical_selection = TRUE, training_method = "EM", 
                 Mstep_numerical = FALSE, BW_max_iter = 50, 
                 BW_limit_accuracy = 0.001, BW_print = TRUE,
                 DNM_max_iter = 50, DNM_limit_accuracy = 0.001, 
                 DNM_print = 2, decoding_method = 'global',
                 bout_lengths = NULL, plotting = 0)
}

\arguments{
  \item{x}{%%
a vector object of length \code{T} containing non-negative observations of a time-series, such as a sequence of accelerometer impulse counts, which are assumed to be realizations of the (hidden Markov state dependent) observation process of a HMM.}
  \item{cut_points}{%%
  a vector object containing cut-off points to separate activity ranges.  For instance, the vector \code{c(7,15,23)} separates the four activity ranges [0,7), [7,15), [15,23) and [23,Inf).}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the \code{m} observation distributions of the Markov dependent observation process. The following distributions are supported: Poisson (\code{pois}); generalized Poisson (\code{genpois}); normal (\code{norm})).}
  \item{min_m}{miminum number of hidden states in the hidden Markov chain. Default value is \code{2}.}
  \item{max_m}{maximum number of hidden states in the hidden Markov chain. Default value is \code{6}.}
  \item{n}{%%
  a single numerical value specifying the number of samples.  Default value is \code{100}.}
  \item{max_scaled_x}{%%
  an optional numerical value,  to be used to scale the observations of the time-series \code{x} before the hidden Markov model is trained and decoded (see Details). Default value is \code{NA}.}
  \item{names_activity_ranges}{%%
  an optional character string vector to name the activity ranges induced by the cut-points. This vector must contain one element more than the vector \code{cut_points}.}
  \item{discr_logL}{%%
a logical object indicating whether the discrete log-likelihood should be used (for \code{"norm"}) for estimating the model specific parameters instead
of the general log-likelihood. See MacDonald & Zucchini (2009, Paragraph 1.2.3) for further details.  Default is \code{FALSE}.}
  \item{discr_logL_eps}{%%
  a single numerical value to approximate the discrete log-likelihood for a hidden Markov model based on nomal distributions (for \code{distribution_class="norm"}).  The default value is \code{0.5}.}
  \item{dynamical_selection}{%%
  a logical value indicating whether the method of dynamical initial parameter selection should be applied (see \code{\link{HMM_training}} for details).  Default is \code{TRUE}.}
  \item{training_method}{%%
a logical value indicating whether the Baum-Welch algorithm (\code{"EM"}) or the method of direct numerical maximization (\code{"numerical"}) should be applied for estimating the model specific parameters of the HMM. See \code{\link{Baum_Welch_algorithm}} and \code{\link{direct_numerical_maximization}} for further details.  Default is 

\code{training_method="EM"}.}
  \item{Mstep_numerical}{%%
a logical object indicating whether the Maximization Step of the Baum-Welch algorithm shall be performed by numerical maximization.  Default is FALSE.}
  \item{BW_max_iter}{%%
a single numerical value representing the maximum number of iterations in the Baum-Welch algorithm. Default value is \code{50}.
}
  \item{BW_limit_accuracy}{%%
a single numerical value representing the convergence criterion of the 
Baum-Welch algorithm. Default value is \code{0.001}.}
  \item{BW_print}{%%
a logical object indicating whether the log-likelihood at each iteration-step shall be printed. Default is \code{TRUE}.
}
  \item{DNM_max_iter}{%%
a single numerical value representing the maximum number of iterations
of the numerical maximization using the nlm-function (used to perform the M-step of the 
Baum-Welch-algorithm). Default value is \code{50}.}
  \item{DNM_limit_accuracy}{%%
a single numerical value representing the convergence criterion of the numerical
maximization algorithm using the \link[stats]{nlm} function (used to perform the M-step of 
the Baum-Welch-algorithm). Default value is \code{0.001}.  
}
  \item{DNM_print}{%%
a single numerical value to determine the level of printing of the \code{nlm}-function.  See \code{nlm}-function for further informations. The value \code{0} suppresses, that no printing will be outputted. Default value is \code{2} for full printing.}
 \item{decoding_method}{%%
a string object to choose the applied decoding-method to decode the HMM given the time-series of observations \code{x}.  Possible values are \code{"global"} (for the use of the \code{\link{Viterbi_algorithm}}) and \code{"local"} (for the use of the \code{\link{local_decoding_algorithm}}). Default value is \code{"global"}.}
  \item{bout_lengths}{%%
a vector object (with even number of elemets) to define the range of the bout intervals (see Details for the definition of bouts).  
For instance, 

\code{bout_lengths=c(1,1,2,2,3,10,11,20,1,20)} defines the five bout intervals 
[1,1] (1 count); [2,2] (2 counts); [3,10] (3-10 counts); [11,20] (11-20 counts); [1,20] 
(1-20 counts - overlapping with other bout intervalls is possible). Default value is \code{bout_lengths=NULL}. 
 }
  \item{plotting}{%%
a numeric value between 0 and 5 (generates different outputs). NA suppresses graphical output. Default value is \code{0}.\cr
\code{0}: output 1-5 \cr
\code{1}: summary of all results \cr
\code{2}: time series of activity counts, classified into activity ranges  \cr
\code{3}: time series of bouts (and, if available, the sequence of the estimated hidden physical activity levels, extracted by decoding a trained HMM, in green colour) \cr
\code{4}: barplots of absolute and relative frequencies of time spent in different activity ranges  \cr
\code{5}: barplots of relative frequencies of the lenghts of bout intervals (overall and by activity ranges )}
}

\details{
The function combines \code{\link{HMM_training}}, \code{\link{HMM_decoding}} and \code{\link{cut_off_point_method}} as follows: \cr

\bold{Step 1:} \code{\link{HMM_training}} trains the most likely HMM for a given time-series of accelerometer counts. \cr
\bold{Step 2:} \code{\link{HMM_decoding}} decodes the trained HMM (Step 1) into the most likely sequence of hidden states corresponding to the given time-series of observations (respectively the most likely sequence of physical activity levels corresponding to the time-series of accelerometer counts). \cr
\bold{Step 3}. \code{\link{cut_off_point_method}} assigns an activity range to each accelerometer count by its hidden physical activity level (extracted in Step 2).
}

\value{
\code{HMM_based_method} returns a list containing the output of the trained hidden Markov model, including the selected number of states \code{m} (i.e., number of physical activities) and plots key figures.
  \item{trained_HMM_with_selected_m}{%%
a list object containing the trained hidden Markov model including the selected number of states \code{m} (see \code{\link{HMM_training}} for further details).}
  \item{decoding}{%%
a list object containing the output of the decoding (see \code{\link{HMM_decoding}} for further details)}.
  \item{extendend_cut_off_point_method}{%%
a list object containing the output of the cut-off point method. 
The counts \code{ x } are classified into the activity ranges by the corresponding sequence of hidden PA-levels, which were decoded by the HMM (see \code{\link{cut_off_point_method}} for further details).}
}
\note{
The parameter \code{ max_scaled_x } can be applied to scale the values of the observations. This might prevent the alogrithm from numerical instabilities.  At the end, the results are internaly rescaled to the original scale.  For instance, a value of \code{ max_scaled_x=200 } shrinks the count values of the complete time-series \code{ x } to a maximum of 200. Training and decoding of the HMM is carried out using the scaled time-series. \cr
From our experience, especially time-series with observations values \code{ >1500}, or where \code{T > 1000}, show numerical instabilities. We then advice to make use of \code{ max_scaled_x }.
 
The extention of the cut-off point method using a Poisson based HMM has been provided and evaluated successfully on simulated data firstly by Barbara Brachmann in her diploma thesis (see References).
}

\references{
Brachmann, B. (2011). Hidden-Markov-Modelle fuer Akzelerometerdaten. Diploma Thesis, University Bremen - Bremen Institute for Prevention Research and Social Medicine (BIPS).

MacDonald, I. L., Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.

Witowski, V. (2013). Hidden-Markov-Modelle fuer Zeitreihen - Analyse von Akzelerometerdaten der IDEFICS-Studie. Diploma Thesis, University Bremen - Leibniz Institute for Preventions Research and Epidemiology - BIPS GmbH
}
\author{
Vitali Witowski (2013). 
}

\seealso{
\code{\link{initial_parameter_training}},
\code{\link{Baum_Welch_algorithm}},
\code{\link{direct_numerical_maximization}},
\code{\link{AIC_HMM}},
\code{\link{BIC_HMM}},
\code{\link{HMM_training}},
\code{\link{Viterbi_algorithm}},
\code{\link{local_decoding_algorithm}},
\code{\link{cut_off_point_method}}
}
\examples{

################################################################
### Fictitious activity counts #################################
################################################################

x <- 100 * c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,
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

	   
### Fictitious cut-off points that produce four so-called 
### activity ranges "sedentary", "light", "moderate", 
### and "vigorous".

cut_points <- 100 * c(7,15,23)
names_activity_ranges <- c("SED","LIG","MOD","VIG")


### Plot fictitious activity counts

plot(x, main = "counts with high values", 
     xlab = "time/epoch", ylab = "counts")
abline(h = cut_points, col = "grey50", lty = "dashed")


################################################################
### Comparing the results of the traditional ################### 
### cut-off point method and the new HMM-based method ##########
################################################################

### Apply the traditional cut-off point method to assign 
### physical activity ranges to each observed count

solution_of_tradtionional_cut_off_point_method <-
   cut_off_point_method(x = x, 
       hidden_PA_levels = NA, 
       cut_points = cut_points, 
       names_activity_ranges = names_activity_ranges, 
       bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12, 
       13,40,41,265,1,265), 
	     plotting = 1)

### Apply the HMM-based method to assign physical activity 
### ranges to the hidden physical activity level of each count

## TIME CONSUMING:
\dontrun{
solution_of_HMM_based_method <- 
    HMM_based_method(x = x, 
      max_scaled_x = 50, 
      cut_points  =cut_points, 
    	min_m = 2, 
    	max_m = 6, 
    	names_activity_ranges = names_activity_ranges, 
      distribution_class = "pois", 
      training_method = "EM", 
      decoding_method = "global", 
      bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12,
      13,40,41,265,1,265), 
      plotting = 1)

		
### Print details of the traditional cut-off point method 
### and the new HMM-based method

print(solution_of_tradtionional_cut_off_point_method)
print(solution_of_HMM_based_method)
}
}

