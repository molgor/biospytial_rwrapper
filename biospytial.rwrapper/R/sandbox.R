# This is a script for loading all the necessary data for processing 

rm(list=ls())
source("init_data.R")
# load the building function
source("joint.binomial.bymCAR.R")

# The results are partial.
n.sample = 10000                                                                                  
burnin=10000                                                                                      
postburnin = burnin +1000                                                                         
thin = 100                                                                                          
verbose = TRUE                                                                                  
rr  <- na.omit(DataFrame)

results  <- joint.binomial.bymCARModel2(formula_S = formula_sample, 
                                        formula_P = formula_presence,
                                        n.sample=n.sample,
                                        burnin=burnin,
                                        postburnin=postburnin,
                                        thin=thin,
                                        verbose=TRUE,
                                       prior.tau2=c(0.01,1),
                                       prior.sigma2=c(0.01,1),
                                       data = rr) 




