# This is a script for loading all the necessary data for processing 

source("single.binomial.bymCAR.R")

joint.binomial.bymCARModel2  <- function(formula_S, formula_P, data=DataFrame,n.sample,burnin, postburnin,thin,verbose){

### Prepare common frame for both

frame.sample.results <- common.frame(formula_S, DataFrame, "binomial")
K.sample <- frame.sample.results$n
p.sample <- frame.sample.results$p
X.sample <- frame.sample.results$X
X.standardised.sample <- frame.sample.results$X.standardised
X.sd.sample <- frame.sample.results$X.sd
X.mean.sample <- frame.sample.results$X.mean
X.indicator.sample <- frame.sample.results$X.indicator 
offset.sample <- frame.sample.results$offset
Y.sample <- frame.sample.results$Y
which.miss.sample <- frame.sample.results$which.miss
n.miss.sample <- frame.sample.results$n.miss  
Y.DA.sample <- Y.sample


frame.presence.results <- common.frame(formula_P, DataFrame, "binomial")
K.presence <- frame.presence.results$n
p.presence <- frame.presence.results$p
X.presence <- frame.presence.results$X


X.standardised.presence <- frame.presence.results$X.standardised
X.sd.presence <- frame.presence.results$X.sd
X.mean.presence <- frame.presence.results$X.mean
X.indicator.presence <- frame.presence.results$X.indicator 
offset.presence <- frame.presence.results$offset
Y.presence <- frame.presence.results$Y
which.miss.presence <- frame.presence.results$which.miss
n.miss.presence <- frame.presence.results$n.miss  
Y.DA.presence <- Y.presence




#### Matrices to store samples
#n.keep <- floor((n.sample - burnin)/thin)

n.keep <- floor((n.sample)/thin)
samples.beta.sample <- array(NA, c(n.keep, p.sample))
samples.beta.presence <- array(NA, c(n.keep, p.sample))
## The same because they share the GMRF
samples.re <- array(NA, c(n.keep, K.sample))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.sigma2 <- array(NA, c(n.keep, 1))

samples.loglike.presence <- array(NA, c(n.keep, K.presence))
samples.fitted.presence <- array(NA, c(n.keep, K.presence))
if(n.miss.presence>0) samples.Y.presence <- array(NA, c(n.keep, n.miss.presence))

samples.loglike.sample <- array(NA, c(n.keep, K.sample))
samples.fitted.sample <- array(NA, c(n.keep, K.sample))
if(n.miss.sample>0) samples.Y.sample <- array(NA, c(n.keep, n.miss.sample))




model.sample  <- single.binomial.bymCAR(formula=formula_S, name='Sample Effort Model',data=DataFrame, trials=trials, W=M_bis, burnin=burnin, n.sample=postburnin, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


model.presence  <- single.binomial.bymCAR(formula=formula_P, name='Presence model', data=DataFrame, trials=trials, W=M_bis, burnin=burnin, n.sample=postburnin, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


presence.state  <-  model.presence$state
sample.state  <- model.sample$state
#### Init data sandboxing the carsampler


#### Start timer
    if(verbose)
    {
    cat("Sampling from the joint model", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }

CarSampler.sample <- model.sample$f_sampler 
CarSampler.presence  <- model.presence$f_sampler
###[J] Trace compilation
    for(j in 1:n.sample)
    {
    ## Iteration of the CarSampler
    ## What happen with the hyperparameters ??  
      sample.state <-  CarSampler.sample(Y.DA = sample.state$Y.DA,
                      beta =  sample.state$beta,
                      phi = presence.state$phi,
                      tau2 =  presence.state$tau2,
                      theta =  presence.state$theta,
                      sigma2 = presence.state$sigma2,
                      prob = sample.state$prob,
                      proposal.sd.theta = presence.state$proposal.sd.theta,
                      proposal.sd.phi = presence.state$proposal.sd.phi,
                      proposal.sd.beta = sample.state$proposal.sd.beta,
                      accept.all = sample.state$accept.all,
                      accept = sample.state$accept,
                      iter_index = j) 


      presence.state <-  CarSampler.presence(Y.DA = presence.state$Y.DA,
                      beta =  presence.state$beta,
                      phi = sample.state$phi,
                      tau2 =  sample.state$tau2,
                      theta = sample.state$theta,
                      sigma2 = sample.state$sigma2,
                      prob = presence.state$prob,
                      proposal.sd.theta = sample.state$proposal.sd.theta,
                      proposal.sd.phi = sample.state$proposal.sd.phi,
                      proposal.sd.beta = presence.state$proposal.sd.beta,
                      accept.all = presence.state$accept.all,
                      accept = presence.state$accept,
                      iter_index = j) 


    ###################
    ## Save the results
    ###################
        if(j %%thin==0)
        {

        ele <- (j) / thin
        #samples.beta[ele, ] <- model.sample$beta
        samples.beta.sample[ele, ]  <- sample.state$beta
        samples.beta.presence[ele, ]  <- presence.state$beta
        ## Here it's the same because they are sharing the same GMRF
        samples.re[ele, ] <- presence.state$phi + presence.state$theta
        samples.tau2[ele, ] <- presence.state$tau2
        samples.sigma2[ele, ] <- presence.state$sigma2
        samples.loglike.presence[ele, ] <- presence.state$loglike
        samples.fitted.presence[ele, ] <- presence.state$fitted
            if(n.miss.presence>0) samples.Y.presence[ele, ] <- presence.state$Y.DA[which.miss.presence==0]
        samples.loglike.sample[ele, ] <- sample.state$loglike
        samples.fitted.sample[ele, ] <- sample.state$fitted
            if(n.miss.sample>0) samples.Y.sample[ele, ] <- sample.state$Y.DA[which.miss.sample==0]
        
               }else
        {
        }

    ################################       
    ## print progress to the console
    ################################
          if(j %in% percentage.points & verbose)
          {
          setTxtProgressBar(progressBar, j/n.sample)
          }
     }

##### end timer
    if(verbose)
    {

    close(progressBar)
    }else
    {}


print("Compiling summary for S process")
results.sample  <- SummariseResults(samples.Y = samples.Y.sample, 
                             samples.fitted = samples.fitted.sample, 
                             samples.beta = samples.beta.sample,  
                             samples.re = samples.re, 
                             samples.tau2 = samples.tau2, 
                             samples.sigma2 =  samples.sigma2, 
                             samples.loglike = samples.loglike.sample, 
                             X = X.sample, 
                             X.standardised = X.standardised.sample, 
                             X.indicator = X.indicator.sample, 
                             X.mean = X.mean.sample, 
                             X.sd = X.sd.sample, 
                             Y = Y.sample, 
                             trials = trials, 
                             p = p.sample, 
                             offset = offset.sample, 
                             accept.all = model.sample$state$accept.all,
                             formula = formula_S, 
                             n.miss = n.miss.sample,
                             n.keep = n.keep)

print("Compiling summary for the P process")
results.presence  <- SummariseResults(samples.Y = samples.Y.presence, 
                             samples.fitted = samples.fitted.presence, 
                             samples.beta = samples.beta.presence,  
                             samples.re = samples.re, 
                             samples.tau2 = samples.tau2, 
                             samples.sigma2 =  samples.sigma2, 
                             samples.loglike = samples.loglike.presence, 
                             X = X.presence, 
                             X.standardised = X.standardised.presence, 
                             X.indicator = X.indicator.presence, 
                             X.mean = X.mean.presence, 
                             X.sd = X.sd.presence, 
                             Y = Y.presence, 
                             trials = trials, 
                             p = p.presence, 
                             offset = offset.presence, 
                             accept.all = model.presence$state$accept.all,
                             formula = formula_P, 
                             n.miss = n.miss.presence,
                             n.keep = n.keep)


exports = list("S"=results.sample, "P"=results.presence)

return(exports)
} 
