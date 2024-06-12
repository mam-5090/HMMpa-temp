# library(HMMpa)
?Baum_Welch_algorithm
?dnbinom
?log
?optim
?neg_log_liklihood
Baum_Welch_algorithm <- 
function(x, m, delta, gamma, distribution_class, distribution_theta, discr_logL = FALSE, 
    discr_logL_eps = 0.5, BW_max_iter = 50, BW_limit_accuracy = 0.001, BW_print=TRUE,
    Mstep_numerical = FALSE, DNM_limit_accuracy = 0.001, DNM_max_iter = 50, DNM_print = 2) 
{
 
################################################################################
### Needed variables and functions #############################################
################################################################################

 
  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k=1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois" |
      distribution_class == "bivariate_pois" | distribution_class == "nbin") 
  {
    k=2
  }	
  
  size <- length(x)	
  oldlogL <- -Inf
  reached_limit_of_accuracy <- FALSE
  underflow_error_1 <- FALSE
  underflow_error_2 <- FALSE
  pos_logL_error <- FALSE
  worse_logL_error <- FALSE

  
  DNM_log_n2w <- function(np) 
  {
    wp <- log(np)
    return(wp)
  }
  
  DNM_exp_w2n <- function(wp) 
  {
    np <- exp(wp)
    return(np)
  }
 
  DNM_logit_n2w <- function(np) 
  {
    wp = log(np / (1 - np) )
    return(wp)
  }
 
  DNM_invlogit_w2n <- function(wp)
  {
    np = exp(wp) /(1 + exp(wp))
    return(np)
  }
  
  # WIP
  neg_log_likelihood <- function(params, x) 
  {
    mu <- exp(params[1])
    size <- exp(params[2])
    -sum(dnbinom(x, size))
  }
  
  DNM_MLE <- function(xt)
  {
    result <- optim(par = c(log(mean(xt)), log(var(xt) / mean(xt) - 1)),  # initial guesses for log(mu) and log(size)
                    fn = neg_log_likelihood,
                    x = xt,
                    method = "BFGS")
    
    # Extract the estimated parameters
    mu <- exp(result$par[1])
    size <- exp(result$par[2])
    prob <- size / (size + mu)
    return(list(mu = mu, size = size, prob = prob))
  }
   
  negterm3 <- function(x, p, distribution_class, m, list_eta_zeta)
  { 
    size <- length(x)
    
    # WIP done-ish
    if (distribution_class == "nbin") 
    {
      
      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          est <- DNM_MLE(x[tt])
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log(dnbinom(x=x[tt], size=est$size,
                                                                  prob = est$prob ) ) 
        }
      }    
      
      if (is.na(term3) | term3 == Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }
      negterm3 <- -term3
    }
    
    if (distribution_class == "pois") 
    {

      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dpois(x=x[tt], lambda=DNM_exp_w2n(p[i]))) 
        }
      }    
      
      if (is.na(term3) | term3 == Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }
      negterm3 <- -term3
    }
    
    
    if (distribution_class == "norm")
    {
      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dnorm(x[tt], mean=p[i], sd=DNM_exp_w2n(p[(m + i)])))
        }
      }
      if(is.na(term3) | term3 == Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }
      negterm3 <- -term3
    }
    
    if (distribution_class == "genpois") 
    {
      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dgenpois(x[tt], lambda1=DNM_exp_w2n(p[i]), lambda2=DNM_invlogit_w2n(p[(m + i)])))
        }
      } 
      if(is.na(term3) | term3==Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }  
      negterm3 <- -term3
    }	
    return(negterm3)
  }
  

################################################################################
######################### The training        ##################################
################################################################################
     
  for (l in 1:BW_max_iter)
  {
    
################################################################################
############### Estimation Step (E-Step)            ############################
################################################################################
   fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, 
            distribution_class = distribution_class, 
            distribution_theta = distribution_theta, 
            discr_logL = discr_logL, discr_logL_eps = discr_logL_eps), 
            silent = FALSE)    

    if (inherits(fb, "try-error")) 
    {  
       underflow_error_1 <- TRUE
       l <- l - 1
       logL <- oldlogL
       break
    }
    
    if (any(is.na(fb$log_alpha)) | any(is.na(fb$log_beta)))
    { 
      underflow_error_2 <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }
    
    
    logL <- fb$logL
    
    
    if (BW_print == TRUE) 
    { 
      print(list("##############################", 
            m = paste("train (EM) HMM with m =", toString(m)), 
            distribution = distribution_class, iteration = l, logL = logL, 
            "##############################"))
    }
    
    if (logL >= 0 ) 
    { 
      pos_logL_error <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }
    
    if (oldlogL > logL )     
    { 
      worse_logL_error <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }   

    # WIP done-ish
    if (distribution_class == "nbin")
    {

      zeta <- matrix(c(0), ncol = m, nrow=size)  
      zeta <- exp(fb$log_alpha + fb$log_beta - logL) 
      
      eta <- array(NA, dim = c((size - 1), m, m))
      for (j in 1:m)
      { 
        for (i in 1:m)
        { 
          for (t in 2:size)
          { 
            eta[t - 1, i, j] <- exp( fb$log_alpha[(t - 1),i] + log(gamma[i, j]) 
                                     + log(dnbinom(x=x[t], size=size, prob = DNM_MLE(x[t])) 
                                     + fb$log_beta[t, j] - logL))          
          }
        }
      }    
    }
    
    if (distribution_class == "pois")
    { 	
        zeta <- matrix(c(0), ncol=m, nrow = size)  
        zeta <- exp(fb$log_alpha + fb$log_beta - logL) 

        eta <- array(NA, dim = c((size-1), m, m))
        for (j in 1:m)
        { 
        	for (i in 1:m)
            { 
            	for (t in 2:size)
            	{ 
          			eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
          			    log( dpois(x[t], distribution_theta$lambda[j]) ) + fb$log_beta[t,j] - logL)
      	  		}
      		}
      	}    
      	      	 
    }
    
    
    
    if (distribution_class == "geom") 
    { 
        zeta <- matrix(c(0), ncol = m, nrow = size)  
        zeta <- exp( fb$log_alpha + fb$log_beta - logL) 

      eta <- array(NA, dim = c( (size - 1), m, m))
      for (j in 1:m) 
      	{ 
      		for (i in 1:m) 
      		{ 
      			for(t in 2:size)
      			{ 
      				eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
      				    log( dgeom(x[t], distribution_theta$prob[j]) ) + fb$log_beta[t,j] - logL)
      			}
      		}
      	}
      	
    }	
    
    
        
    if (distribution_class == "genpois")
    { 
      zeta <- matrix(c(0), ncol = m, nrow = size)  
      zeta <- exp( fb$log_alpha + fb$log_beta - logL) 
     
     	eta <- array(NA, dim = c((size - 1), m, m))
      	for (j in 1:m) 
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) 
      				    + log(dgenpois(x[t], distribution_theta$lambda1[j], distribution_theta$lambda2[j]) ) 
      				    + fb$log_beta[t,j] - logL)
      			}
      		}
      	}
      	
    }
    
    
    if (distribution_class == "norm")
    { 
      zeta <- matrix(c(0), ncol = m, nrow = size )  
      zeta <- exp( fb$log_alpha + fb$log_beta - logL) 
      
    	eta <- array(NA, dim = c((size-1), m, m))
      	for (j in 1:m)
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t-1,i , j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
      				    log(dnorm(x[t], distribution_theta$mean[j], distribution_theta$sd[j]) ) 
      				    + fb$log_beta[t,j]  - logL)
      			}
      		}
      	}
      	
    }
    
    
    
    if (distribution_class == "bivariate_pois")
    { 	
      size=length(x[,1])
       
      zeta <- matrix(c(0), ncol = m, nrow=size)  
      zeta <- exp(fb$log_alpha + fb$log_beta - logL) 
      
      eta <- array(NA, dim = c((size - 1), m, m))
      	for (j in 1:m)
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t - 1, i, j] <- exp( fb$log_alpha[(t - 1),i] + log(gamma[i, j]) 
      				  + log(dpois(x[t, 1], distribution_theta$lambda_1[j]) ) 
      				  + log(dpois(x[t, 2], distribution_theta$lambda_2[j])) 
      				  + fb$log_beta[t, j] - logL)          
      			}
      		}
      	}    	
    }


  list_eta_zeta = list(zeta = zeta, eta = eta)
          
################################################################################
############### Maximization Step (M-Step)  ####################################
################################################################################
    
 delta <- list_eta_zeta$zeta[1,]
    
 gamma <- matrix(0, ncol=m, nrow=m)
    for (i in 1:m)
    { 
    	for (j in 1:m)
    	{ 
    		sum_numerator <- 0
    		sum_numerator <- sum(list_eta_zeta$eta[, i, j]) 
    		sum_denominator <- 0
      		for (k in 1:m) 
      		{ 
        		sum_denominator <- sum_denominator + sum(list_eta_zeta$eta[, i, k]) 
      		}
      
      	gamma[i, j] <-  sum_numerator / sum_denominator      
    	}
    }
    
    
    # WIP
    if ( distribution_class == "nbin" & Mstep_numerical == FALSE) 
    {
      mu <- rep(0, times = m)
      size <- rep(0, times = m)
      prob <- rep(0, times = m)
      
      for (i in 1:m) 
      {
        sum_numerator <- 0
        sum_denominator <- 0
        
        for (tt in 1:size) 
        {
          sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt, i] * x[tt]
          sum_denominator <- sum_denominator + list_eta_zeta$zeta[tt, i]
        }
        
        mu[i] <- sum_numerator / sum_denominator
        size[i] <- var(x) / mu[i]
        prob[i] <- size[i] / (size[i] + mu[i])
      }
      
      estimated_mean_values <- mu
      distribution_theta <- list(mu = mu, size = size, prob = prob)
    }
 
 
    if (distribution_class == "pois" & Mstep_numerical == FALSE) 
    { 
      
      lambda <- rep(0, times = m)
      
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[, i] * x)
        sum_denominator <- 0
        sum_denominator <- sum(list_eta_zeta$zeta[, i])
                
        lambda[i] <- sum_numerator / sum_denominator
      }
      
      estimated_mean_values <- lambda
      distribution_theta <- list(lambda = lambda)
    }
    
    
    if (distribution_class == "norm" & Mstep_numerical == FALSE)
    { 	
      mean <- rep(0, times = m)
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[, i] * x)
          
        sum_denominator <- 0
        sum_denominator <- sum(list_eta_zeta$zeta[, i])
     
      mean[i] <- sum_numerator / sum_denominator
      }
      
      sd <- rep(0, times = m)
      for (i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt,i] * ((x[tt] - mean[i])^2)
        }
                  
       sum_denominator <- 0
       sum_denominator <- sum(list_eta_zeta$zeta[, i])
       sd[i] <- sum_numerator / sum_denominator
       sd[i] <- sqrt(sd[i])
      }
      
      estimated_mean_values <- mean
      distribution_theta <- list(mean = mean, sd = sd)
    }
    
    
    if(distribution_class == "geom") 
    { 
      prob <- rep(0, times = m)
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[,i] * x)
        

        sum_denominator <- 0
        for (tt in 1:size)
        {
          sum_denominator <- sum_denominator + list_eta_zeta$zeta[tt, i] * (x[tt] - 1)
        }
        
        prob[i] <- sum_numerator / sum_denominator
      }
      
      estimated_mean_values <- 1 / prob
      distribution_theta <- list(prob = prob)
    }
    
    
    if(distribution_class == "bivariate_pois")
    { 
      size <- length(x[, 1])	
      lambda_1 <- rep(0, times = m)
      lambda_2 <- rep(0, times = m)
      
      for(i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt, i] * x[tt, 1]
        }
    
        sum_denominator <- 0
        for (tt in 1:size)
        { 
          sum_denominator + list_eta_zeta$zeta[tt, i]
        }
        
        lambda_1[i] <- sum_numerator / sum_denominator
      }
      
      for (i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt, i] * x[tt, 2]
        }
        
        sum_denominator <- 0
        for (tt in 1:size)
        { 
          sum_denominator <- sum_denominator + list_eta_zeta$zeta[tt, i]
        }
        
        lambda_2[i] <- sum_numerator / sum_denominator
      }
      
      
      estimated_mean_values <- lambda_1
      
      distribution_theta <- list(lambda_1 = lambda_1, lambda_2 = lambda_2)
    }
    
 
    # WIP
    if(distribution_class == "nbin" & Mstep_numerical == TRUE)
    {
      # x = distribution_theta$x, size = distribution_theta$size, 
      distribution_theta <- list(x = , size = , prob = prob)
    }
 
 
    if(distribution_class == "pois" & Mstep_numerical == TRUE) 
    {
      trans_lambda <- DNM_log_n2w(distribution_theta$lambda)
      
      vector_of_parameters <- c(trans_lambda)   
      
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, print.level = DNM_print, 
                      gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
            
      est_lambda <- DNM_exp_w2n(minterm3$estimate)
      
      estimated_mean_values <- est_lambda 
      
      estimated_var <- est_lambda
      estimated_sd <- sqrt(estimated_var)
      
      distribution_theta <- list(lambda=est_lambda)
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))  
    }
    
    
    if (distribution_class == "norm" & Mstep_numerical == TRUE)
    {	
      trans_mean <- distribution_theta$mean
      
      trans_sd <- DNM_log_n2w(distribution_theta$sd)   
      
      vector_of_parameters <- c(trans_mean, trans_sd) 
        
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, print.level = DNM_print, 
                      gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
      
      estimated_mean <- minterm3$estimate[1:m]
      
      sum_denominator <- minterm3$estimate[1:m]
      
      estimated_sd <- DNM_exp_w2n(minterm3$estimate[(m + 1):(2 * m)])
      
      estimated_mean_values <- estimated_mean
      
      estimated_var <- estimated_sd^2
      
      estimated_sd <- estimated_sd      
      
      distribution_theta <- list(mean = estimated_mean, sd = estimated_sd) 
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))
    }
    
    if (distribution_class == "genpois" & Mstep_numerical == TRUE) 
    {	
      trans_lambda1 <- DNM_log_n2w(distribution_theta$lambda1)
      
      trans_lambda2 <- DNM_logit_n2w(distribution_theta$lambda2)
      
      vector_of_parameters <- c(trans_lambda1, trans_lambda2)   
      
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, 
                      print.level = DNM_print, gradtol = DNM_limit_accuracy, 
                      iterlim = DNM_max_iter)
      
      estimated_lambda1 <- DNM_exp_w2n(minterm3$estimate[1:m])
      
      estimated_lambda2 <- DNM_invlogit_w2n(minterm3$estimate[(m + 1):(m + m)])
      
      estimated_mean_values <- estimated_lambda1 / (1 - estimated_lambda2)
      
      estimated_var <- estimated_lambda1 / ((1-estimated_lambda2)^3)
      
      estimated_sd <- sqrt(estimated_var)
      
      distribution_theta <- list(lambda1 = estimated_lambda1, lambda2 = estimated_lambda2)
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))
    }
    
    difference_old_logL_and_new_logL = abs(oldlogL - logL)
    
    if (difference_old_logL_and_new_logL < BW_limit_accuracy) 
    { 
      reached_limit_of_accuracy = TRUE
      break
    }
    
    oldlogL <- logL
  }


################################################################################
######################### Accessing AIC and BIC for the trained HMM ############
################################################################################
  
  AIC <- AIC_HMM(logL = logL, m = m, k = k) 
  BIC <- BIC_HMM(size = size, logL = logL, m = m, k = k) 
  
################################################################################
######################### Return results #######################################
################################################################################
  
  return(list(x = x,
              m = m,
              zeta = zeta,
              eta = eta,
              iter = l,
              logL = logL,
              AIC = AIC,
              BIC = BIC,
              delta = delta,
              gamma = gamma,
              distribution_class = distribution_class,
              distribution_theta = distribution_theta,
              estimated_mean_values = estimated_mean_values,
              reached_limit_of_accuracy = reached_limit_of_accuracy,
              underflow_error_1 = underflow_error_1,
              underflow_error_2 = underflow_error_2,
              pos_logL_error = pos_logL_error,
              worse_logL_error = worse_logL_error))
}
