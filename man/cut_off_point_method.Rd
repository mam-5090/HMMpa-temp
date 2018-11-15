\name{cut_off_point_method}
\alias{cut_off_point_method}

\title{
Cut-Off Point Method for Assigning Physical Activity Patterns
}
\description{
This function assigns an activity range to each observation of a time-series, such as for a sequence of impulse counts recorded by an accelerometer.  The activity ranges are defined by thresholds called \dQuote{cut-off points}.  Furthermore, bout periods are analysed (see Details for further informations). 
}
\usage{
cut_off_point_method (x, cut_points, 
        names_activity_ranges = NA, hidden_PA_levels = NA, 
        bout_lengths = NULL, plotting = 0)
}

\arguments{
  \item{x}{%%
a vector object of length \code{T}  containing non-negative observations of a time-series, such as a sequence of accelerometer impulse counts.}
  \item{cut_points}{%%
a vector object containing cut-off points to separate activity ranges.  For instance, the vector c(7,15,23) separates the four activity ranges [0,7);[7,15);[15,23);[23,Inf).}
  \item{names_activity_ranges}{%%
an optional character string vector to name the activity ranges induced by the cut-points.  This vector must contain one element more than the vector \code{cut_points}.}
  \item{bout_lengths}{%%
a vector object (with even number of elemets) to define the range of the bout intervals (see Details for the definition of bouts).  
For instance, 

\code{bout_lengths=c(1,1,2,2,3,10,11,20,1,20)} defines the five bout intervals 
[1,1] (1 count); [2,2] (2 counts); [3,10] (3-10 counts); [11,20] (11-20 counts); [1,20] 
(1-20 counts - overlapping with other bout intervalls is possible). Default value is \code{bout_lengths=NULL}. 
 }
  \item{hidden_PA_levels}{%%
an optional vector object of length \code{T} containing a sequence of the estimated hidden physical activity levels (i.e. means) underlying the time-series of accelerometer counts.  Such a sequence can be extracted by decoding a trained hidden Markov model.  The cut-point method classifies then each count by its level in the hidden Markov chain that generates the physical activity counts, and does not use the observed count value (see \code{\link{HMM_based_method}} for further details).  Default is NA (for the traditional cut-point method).}
  \item{plotting}{%%
a numeric value between \code{0} and \code{5} (generates different outputs). NA suppresses graphical output. Default value is \code{0}.\cr
\code{0}: output 1-5 \cr
\code{1}: summary of all results \cr
\code{2}: time series of activity counts, classified into activity ranges \cr
\code{3}: time series of bouts (and, if available, the sequence of the estimated hidden physical activity levels, extracted by decoding a trained HMM, in green colour)\cr
\code{4}: barplots of absolute and relative frequencies of time spent in different activity ranges\cr
\code{5}: barplots of absolute frequencies of different bout intervals (overall and by activity ranges )}
}

\details{
A bout is defined as a period of time spending a defined intensity of physical activities in a specified physical activity range, without switching to activity intensities in a different activity range.
}

\value{
\code{ cut_off_point_method } returns a list containing the extracted sequence of activity ranges  and plots key figures. 
  \item{activity_ranges}{an array object containing the cut-off intervals that indicate the activity ranges.}
  \item{classification}{an integer vector containing the sequence of activity ranges  that were assigned to the observed time-series of accelerometer counts.  If \code{hidden_PA_levels=NA}, then \code{classification} is the output of the traditional cut-point method, meaning that an activity range has been assigned to each accelerometer count over its observed value actual position.  In case when \code{hidden_PA_levels} is available, \code{classification} is the output of the extendend cut-point method using hidden Markov models (see \code{\link{HMM_based_method}} for further details). }
  \item{classification_per_activity_range}{a pairlist object containing the classification of the observed counts by the assigned activity range.}
  \item{freq_acitvity_range}{table object containing the absolute frequencies of classifications into activity ranges.}
  \item{rel_freq_acitvity_range}{table object containing the relative frequencies of classifications into activity ranges.}
  \item{quantity_of_bouts}{overall number of bouts.}
  \item{bout_periods}{an array including the bout length assigned to acitiy ranges.}
  \item{abs_freq_bouts_el}{a pairlist object containing the absolute frequency of bout length per epoch length (aggregated).}
}

\author{
Vitali Witowski (2013).
}

\seealso{
\code{\link{HMM_based_method}}
}


\examples{


################################################################
### Fictitious activity counts #################################
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


################################################################
### i.) Tradionional cut_point method ##########################
################################################################

### Assigning activity ranges to activity counts using 
### fictitious cut-off points that produce the four activity 
### ranges "sedentary"", "light"", "moderate"", and "vigorous". 
\donttest{
solution_of_traditional_cut_off_point_method <- 
    cut_off_point_method(x = x, 
      cut_points = c(5,15,23), 
      names_activity_ranges = c("SED","LIG","MOD","VIG"), 
      bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12,
      13,40,41,265,1,265), 
      plotting = 0)

print(solution_of_traditional_cut_off_point_method)
}

###############################################################
### ii.) Extension of the tradionional cut_point method #######
###      using HMMs      ######################################
###############################################################

## The following three steps define an extension of the 
## traditional cut-off method by first extracting the hidden 
## physical activity pattern behind the accelerometer counts 
## using a HMM (those three steps are basically combined in 
## the function HMM_based_method, see HMM_based_method for 
## further details and references): 


### Step 1 ##################################################### 

## Train hidden Markov model for different number of 
## states m=2,...,6 and select the model with the most 
## plausible m
\donttest{
m_trained_HMM <- 
    HMM_training(x = x, 
      min_m = 2, 
      max_m = 6, BW_print=FALSE,
      distribution_class = "pois")$trained_HMM_with_selected_m
}

### Step 2 ##################################################### 

## Decode the trained HMM (by using the 
## Viterbi algorithm (global decoding)) to get the estimated 
## sequence of hidden physical activity levels 
## underlying the the accelerometer counts 

## You have to compute 'm_trained_HMM' first (see Step 1)
\donttest{
global_decoding <- 
    HMM_decoding(x = x, 
      m = m_trained_HMM$m, 
      delta = m_trained_HMM$delta, 
      gamma = m_trained_HMM$gamma, 
      distribution_class = m_trained_HMM$distribution_class, 
      distribution_theta = m_trained_HMM$distribution_theta,
      decoding_method = "global")
        
hidden_PA_levels <- 
    global_decoding$decoding_distr_means
}

### Step 3 #####################################################

## Assigning activity ranges to activity counts using the 
## information extracted by decoding the HMM for the counts 
## (PA-levels) and fictitious cut-off points that produce 
## four so-called activity ranges:"sedentary", "light", 
## "moderate" and "vigorous":

## You have to compute 'm_trained_HMM' and 
## 'hidden_PA_levels' first (see above)
\donttest{
solution_of_HMM_based_cut_off_point_method <- 
 cut_off_point_method(x = x, 
  hidden_PA_levels = hidden_PA_levels, 
  cut_points = c(5,15,23), 
  names_activity_ranges = c("SED","LIG","MOD","VIG"), 
  bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12,13,40,41,265,1,265), 
  plotting=1)
}
}
\keyword{ ts } 
