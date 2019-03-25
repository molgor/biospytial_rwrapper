source("imports.R")
source("samplerCarFunction.R")
source("summary.results.R")


single.binomial.bymCAR <- function(formula, data=NULL,name="<no name>", trials, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)
{

##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "binomial")
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.DA <- Y

    
#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check and format the trials argument
#####[J] This are only checks for assessing data 
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- K-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
failures.DA <- trials - Y.DA
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(1, 0.01)

###[J] This functions stop the execution
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)
common.prior.var.check(prior.sigma2)


## Compute the blocking structure for beta     
###[J] Structure for the fixed effects
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
    for(r in 1:n.beta.block)
    {
    list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
    list.block[[r+n.beta.block]] <- length(list.block[[r]])
    }


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)



#############################
#### Initial parameter values
#############################
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
### Extract the sd from the covariance matrix of the beta estimators, 
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
## generate random betas 
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

## Init values for the theta parameter (i.e \phi = \psi + \theta 
theta.hat <- Y / trials
## Avoid infinity values when applying map

theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
tau2 <- var(phi) / 10
theta <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
sigma2 <- var(theta) / 10
lp <- as.numeric(X.standardised %*% beta) + phi + theta + offset
prob <- exp(lp)  / (1 + exp(lp))
##
######

###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.re <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

  
#### Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.theta <- 0.1
tau2.posterior.shape <- prior.tau2[1] + 0.5 * (K-1)
sigma2.posterior.shape <- prior.sigma2[1] + 0.5 * K 
 


##################################
#### Set up the spatial quantities
##################################
####[J] CAR quantities
W.quants <- common.Wcheckformat(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
tau2.posterior.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   



###########################
#### Run the Bayesian model
###########################
###[J] Partial function of the sampler, used only to iterate over the parameters. 
CarSampler = partial(CompleteCarSampler,
                            X.standardised = X.standardised, 
                            K = K, 
                            p = p,
                            trials = trials,
                            prior.mean.beta = prior.mean.beta,
                            prior.var.beta = prior.var.beta,
                            prior.tau2 = prior.tau2,
                            prior.sigma2 = prior.sigma2,
                            n.beta.block = n.beta.block,
                            list.block = list.block,
                            W.triplet=W.triplet,
                            W.begfin = W.begfin,
                            W.triplet.sum = W.triplet.sum,
                            n.miss = n.miss,
                            n.triplet = n.triplet,
                            which.miss = which.miss,
                            MALA = MALA,
                            tau2.posterior.shape = tau2.posterior.shape,
                            sigma2.posterior.shape = sigma2.posterior.shape,
                            Y = Y,
                            offset)



###[J] Initialize sampler 
    model.sample <-  CarSampler(Y.DA = Y.DA,
                    beta = beta,
                    phi = phi,
                    tau2 = tau2,
                    theta = theta,
                    sigma2 = sigma2,
                    prob = prob,
                    proposal.sd.theta = proposal.sd.theta,
                    proposal.sd.phi = proposal.sd.phi,
                    proposal.sd.beta = proposal.sd.beta,
                    accept.all = accept.all,
                    accept = accept,
                    iter_index = 1) 
 



#### Start timer
    if(verbose)
    {
    cat("Performing burnin period for",name,"\n", sep = " ")
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
     model.sample <-  CarSampler(Y.DA = model.sample$Y.DA,
                      beta = model.sample$beta,
                      phi = model.sample$phi,
                      tau2 =  model.sample$tau2,
                      theta =  model.sample$theta,
                      sigma2 = model.sample$sigma2,
                      prob = model.sample$prob,
                      proposal.sd.theta = model.sample$proposal.sd.theta,
                      proposal.sd.phi = model.sample$proposal.sd.phi,
                      proposal.sd.beta = model.sample$proposal.sd.beta,
                      accept.all = model.sample$accept.all,
                      accept = model.sample$accept,
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
    cat("\nSummarising results.")
    close(progressBar)
    }else
    {}



###################################
#### Summarise and save the results 
###################################

accept.all  <- model.sample$accept.all
results  <- SummariseResults(samples.Y = samples.Y, samples.fitted = samples.fitted, samples.beta = samples.beta,  samples.re = samples.re, samples.tau2 = samples.tau2, samples.sigma2 =  samples.sigma2, samples.loglike = samples.loglike, X = X, X.standardised = X.standardised, X.indicator = X.indicator, X.mean = X.mean, X.sd = X.sd, Y = Y, trials = trials, p = p, offset = offset, accept.all = accept.all,formula = formula, n.miss = n.miss,n.keep = n.keep)



class(results) <- "CARBayes"


#### Finish by stating the time taken    
    if(verbose)
    {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    }else
    {}
exports = list('model.results'=results, 'state'=model.sample,'f_sampler'=CarSampler)    
return(exports)
}

