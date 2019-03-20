# This is a script for loading all the necessary data for processing 
# i.e. it is a preprocessing script.
# Change it appropriately

#my.binomial.bymCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


#my.binomial.bymCAR(formula=formula_sample,W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE)


#model.sample <-S.CARbym(formula=formula_sample,family="binomial",W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE)

### Order of the code to be executed
# 1. Load init_data
# 2. Load imports.R
# 3. Load as function my_CARbym in S.CARbym.R
# 4. Load as function my.binomial.bymCAR in my.binomial.bymCAR.R
# 5. Run the above line

#model.sample <-my_CARbym(formula=formula_sample,family="binomial",W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE,custom_function=my.binomial.bymCAR) 


#binomial.bymCAR <- getFromNamespace('binomial.bymCAR','CARBayes')
#
#model_std <- binomial.bymCAR(formula=formula_sample, data=DataFrame, trials=trials, W=M_bis, burnin=10000, n.sample=15000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)
#
#


rm(list=ls())
source("init_data.R")
# load the building function
source("samplerCarFunction.R")
source("binomial.jointCARModel.R")

model.sample  <- binomial.bymCAR2(formula=formula_sample, data=DataFrame, trials=trials, W=M_bis, burnin=10000, n.sample=15000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


model.pines  <- binomial.bymCAR2(formula=formula_presence, data=DataFrame, trials=trials, W=M_bis, burnin=10000, n.sample=15000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


pines.state  <-  model.pines$state
sample.state  <- model.sample$state
#### Init data sandboxing the carsampler
n.sample = 5000

### Remove and leave only the necesary

#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.re <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
#if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))



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


###[J] Trace compilation
    for(j in 1:n.sample)
    {
    ## Iteration of the CarSampler
    ## What happen with the hyperparameters ??  
      sample.state <-  CarSampler(Y.DA = sample.state$Y.DA,
                      beta =  sample.state$beta,
                      phi = pines.state$phi,
                      tau2 =  pines.state$tau2,
                      theta =  pines.state$theta,
                      sigma2 = pines.state$sigma2,
                      prob = sample.state$prob,
                      proposal.sd.theta = pines.state$proposal.sd.theta,
                      proposal.sd.phi = pines.state$proposal.sd.phi,
                      proposal.sd.beta = sample.state$proposal.sd.beta,
                      accept.all = sample.state$accept.all,
                      accept = sample.state$accept,
                      iter_index = j) 


      pines.state <-  CarSampler(Y.DA = pines.state$Y.DA,
                      beta =  pines.state$beta,
                      phi = sample.state$phi,
                      tau2 =  sample.state$tau2,
                      theta = sample.state$theta,
                      sigma2 = sample.state$sigma2,
                      prob = pines.state$prob,
                      proposal.sd.theta = sample.state$proposal.sd.theta,
                      proposal.sd.phi = sample.state$proposal.sd.phi,
                      proposal.sd.beta = pines.state$proposal.sd.beta,
                      accept.all = pines.state$accept.all,
                      accept = pines.state$accept,
                      iter_index = j) 






    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {

        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- model.sample$beta
        samples.re[ele, ] <- model.sample$phi + model.sample$theta
        samples.tau2[ele, ] <- model.sample$tau2
        samples.sigma2[ele, ] <- model.sample$sigma2
        samples.loglike[ele, ] <- model.sample$loglike
        samples.fitted[ele, ] <- model.sample$fitted
            if(n.miss>0) samples.Y[ele, ] <- model.sample$Y.DA[which.miss==0]
        
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
    cat("\nSummarising new results.")
    close(progressBar)
    }else
    {}








