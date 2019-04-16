## Functions for summarizing results

SummariseResults  <- function(samples.Y,samples.fitted,samples.beta, samples.re, samples.tau2, samples.sigma2, samples.loglike, X, X.standardised, X.indicator, X.mean, X.sd, Y, trials,offset,formula, n.miss,n.keep, p, accept.all) {
  accept.beta <- 100 * accept.all[1] / accept.all[2]
  accept.phi <- 100 * accept.all[3] / accept.all[4]
  accept.theta <- 100 * accept.all[5] / accept.all[6]
  accept.tau2 <- 100
  accept.sigma2 <- 100
  accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2)
  names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2")
 

#### Compute the fitted deviance
  mean.beta <- apply(samples.beta, 2, mean)
  mean.re <- apply(samples.re, 2, mean)
  mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.re + offset    
  mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
  fitted.mean <- trials * mean.prob
  deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE), na.rm=TRUE)
  
  
  #### Model fit criteria
  ##### returns some estimators given deviance and loglikelihood
  modelfit <- common.modelfit(samples.loglike, deviance.fitted)
  
  
  #### transform the parameters back to the origianl covariate scale.
  samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
  
  
  #### Create a summary object
  samples.beta.orig <- mcmc(samples.beta.orig)
  summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  rownames(summary.beta) <- colnames(X)
  colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  
  #### Summary of hyperparameters
  #### tau2, sigma2, 
  summary.hyper <- array(NA, c(2,7))
  summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
  summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
  summary.hyper[2, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
  summary.hyper[2, 4:7] <- c(n.keep, accept.sigma2, effectiveSize(samples.sigma2), geweke.diag(samples.sigma2)$z)
  
  summary.results <- rbind(summary.beta, summary.hyper)
  rownames(summary.results)[(nrow(summary.results)-1):(nrow(summary.results))] <- c("tau2", "sigma2")
  summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
 
  
  fitted.values <- apply(samples.fitted, 2, mean)
  response.residuals <- as.numeric(Y) - fitted.values

  pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
  residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)


  ### Compile and return the results
  model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - BYM CAR\n")
  if(n.miss==0) samples.Y = NA

  samples <- list(beta=samples.beta.orig, psi=mcmc(samples.re), tau2=mcmc(samples.tau2), sigma2=mcmc(samples.sigma2), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y),loglike=mcmc(samples.loglike))


  results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)



#  exports  <- list("results" = summary.results, 
#                   "mean.prob" = mean.prob,
#                   "samples.beta.orig" = samples.beta.orig,
#                   "modelfit" = modelfit,
#                   "accept.final" = accept.final) 
#
#
  return(results)
}
