# This is the script that implements the Model 2 (common GMRF)

## Set Working directory
## Import code:
setwd('/apps/external_plugins/biospytial_rwrapper/biospytial.rwrapper/R/')
## remove any previous object
rm(list=ls())                                                                                     

## Init data source with variables
print('Load data source and preprocess')
source("init_data.R")                                                                             
# load the building function                                                                      
#source("joint.binomial.bymCAR.R")

## Build dataframes, S<- Sample, P<- Presence
S <- model.frame(formula_sample, DataFrame,na.action='na.pass')
P <- model.frame(formula_presence, DataFrame,na.action='na.pass')

## Split lhs and rhs from the design matrix
SX <- select(S, -c(1))
PX <- select(P, -c(1))
Sy <- select(S, c(1))
Py <- select(P, c(1))

names(Sy) <- 'response'
names(Py) <- names(Sy)
Y = rbind(Sy,Py)

### Let's build covariance matrix
T1 <- matrix(rep(0,4), ncol = 2)
T2 <- matrix(rep(0,4), ncol = 2)
T1[1,1] <- 1
T2[2,2] <- 1

## Perform Kronnecker with different covariates (Block diagonal)
X <- data.frame((T1 %x% as.matrix(SX)) + (T2 %x% as.matrix(PX)))
names(X) <- c(names(SX),names(PX))

DD <- cbind(Y,X)

nK <- dim(M_bis)[1]
## make sequence vector for id.area
ida <- data.frame(seq(nK))
idarea <- unlist(rbind(ida,ida))
## make sequence vector for correlation 
corx <- rep(x = 1,times = nK)
cory <- rep(x = 2,times = nK)

indre <- c(corx,cory)
## A general formula for all the covariates
formula <- response ~ Disttoroadm + Populationm + Elevationm + MeanTempm

###### Runnning the model
## now, assuming that the order in M_bis is the same as in cellids (OOOORDEEER, not value)
## Run the model
trials = rep(1,2 * nK)
burnin = 50000
n.sample = 100000
thin = 50
#model2 <- S.CARmultilevel(formula,family = 'binomial',
#                          trials=trials, 
#                          W=M_bis, 
#                          ind.area = idarea,
#                         ind.re=factor(idarea),
#                          rho = 1,
#                          burnin = burnin,
#                          n.sample = n.sample,
#                          data = DD
#                         )
#


## Cross validation
library(pROC)
library('caret')
#trains = createFolds(y = DataFrame$species, k=7, returnTrain = TRUE)
validate = createFolds(y = DataFrame$species, k=7, returnTrain = FALSE)

#DataFrame$presences <- DataFrame$species
#model2 <- S.CARmultilevel(formula,family = 'binomial',
#                          trials=trials, 
#                          W=M_bis, 
#                          ind.area = idarea,
#                         ind.re=factor(idarea),
#                          rho = 1,
#                          burnin = burnin,
#                          n.sample = n.sample,
#                          data = DD
#                         )
#
#

