# This is a script for loading all the necessary data for processing 

rm(list=ls())
source("init_data.R")
# load the building function
source("joint.binomial.bymCAR.R")


n.sample = 10000
burnin=10000
postburnin = burnin +1000 
thin = 1
#p = 3
#aK = 4061
verbose = TRUE
### Remove and leave only the necesary

results  <- joint.binomial.bymCARModel2(formula_S = formula_sample, formula_P = formula_presence,n.sample=n.sample,burnin=burnin,postburnin=postburnin,thin=thin,verbose=TRUE)

# The results are partial.




