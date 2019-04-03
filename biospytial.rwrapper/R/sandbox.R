# This is a script for loading all the necessary data for processing 

rm(list=ls())
source("init_data.R")
# load the building function
source("joint.binomial.bymCAR.R")

# The results are partial.
n.sample = 20000                                                                                  
burnin=5000                                                                                      
postburnin = burnin +10                                                                         
thin = 1
verbose = TRUE
c <- 100

N = 2

loglike  <- data.frame(matrix(ncol=34,nrow=0))
names(loglike)  <-  c("llS", "llP","DICS","DICP",
                   "dist2rd_med", "dist2rd_25","dist2rd_975","dist2rd_accept",
                    "pop_med","pop_25","pop_975","pop_accept",
                    "elev_med", "elev_25","elev_975","elev_accept",
                    "mtemp_m", "mtemp_25","mtemp_975","mtemp_accept",
                    "tau2_m","tau2_25","tau2_975",
                    "sigma2_m", "sigma2_25","sigma2_975",
                    "interceptS_med", "interceptS_25","interceptS_975","interceptS_accept",
                    "interceptP_med", "interceptP_25","interceptP_975","interceptP_accept"
                    )

for (c in seq(1,1000000,10000)){

results  <- joint.binomial.bymCARModel2(formula_S = formula_sample, 
                                        formula_P = formula_presence,
                                        n.sample=n.sample,
                                        burnin=burnin,
                                        postburnin=postburnin,
                                        thin=thin,
                                        verbose=TRUE,
                                       prior.tau2=c(0.1,0.01),
                                       prior.sigma2=c(0.1,0.01),
                                        tau2_sigma2_denominator=c,
                                       #prior.tau2=c(0.1,1),
                                       #prior.sigma2=c(0.1,1),
                                       data = rr)

loglike[nrow(loglike) + 1, ]  <- list("llS"=results$S$modelfit[6],
                                      "llP"=results$P$modelfit[6],
                                      "DICS"= results$S$modelfit[1],
                                      "DICP"= results$P$modelfit[1],
                                      "dist2rd_med" = results$S$summary.results[2,1],
                                      "dist2rd_25" = results$S$summary.results[2,2],
                                      "dist2rd_975" = results$S$summary.results[2,3],
                                      "dist2rd_accept" = results$S$summary.results[2,5],
                                      "pop_med" =  results$S$summary.results[3,1], 
                                       "pop_25" = results$S$summary.results[3,2], 
                                      "pop_975"= results$S$summary.results[3,3],
                                      "pop_accept"= results$S$summary.results[3,5],
                                      "elev_med" = results$P$summary.results[2,1],
                                      "elev_25"= results$P$summary.results[2,2],
                                      "elev_975"= results$P$summary.results[2,3],
                                      "elev_accept"= results$P$summary.results[2,5],
                                      "mtemp_m"= results$P$summary.results[3,1],
                                      "mtemp_25" = results$P$summary.results[3,2],
                                      "mtemp_975" = results$P$summary.results[3,3],
                                      "mtemp_accept" = results$P$summary.results[3,5],
                                      "tau2_m" = results$P$summary.results[4,1],
                                      "tau2_25" = results$P$summary.results[4,2],
                                      "tau2_975"= results$P$summary.results[4,3],
                                      "sigma2_m"= results$P$summary.results[5,1],
                                      "sigma2_25" = results$P$summary.results[5,2],
                                      "sigma2_975" = results$P$summary.results[5,3],
                                      "interceptS_med" = results$S$summary.results[1,1],
                                      "interceptS_25"=results$S$summary.results[1,2],
                                      "interceptS_975"=results$S$summary.results[1,3],
                                      "interceptS_accept"=results$S$summary.results[1,5],
                                      "interceptP_med"=results$P$summary.results[1,1],
                                      "interceptP_25"=results$P$summary.results[1,2],
                                      "interceptP_975"=results$P$summary.results[1,3],
                                      "interceptP_accept"=results$P$summary.results[1,5]
                                      )
#
print(paste("Finish iteration for",c))
                    }
                                      

