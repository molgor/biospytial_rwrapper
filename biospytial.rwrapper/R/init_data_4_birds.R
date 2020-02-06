# This is a script for loading all the necessary data for processing 
# i.e. it is a preprocessing script.
# Change it appropriately

#file = '/outputs/presence_only_models/predictors/dataset100x100-puebla-p9/0-pred.csv'
#PDF = read.csv(file)
## REad adjancency matrix

## Load Stuff
library(CARBayes)
library(dplyr)
library(purrr)
library(reticulate)
library(biospytial.rwrapper)

## figure out a way to export the functions here. 
# It works for presence_absence_naive. I dont know why not with the others
source("SpeciesModels.R")
# Import adjancency matrix generated from region
mat_filename = "/outputs/training_data_sample_puebla_p9_abies_pinophyta_adjmat.npy"
# Use numpy functions
np <- import("numpy")
M <- np$load(mat_filename)

## Training data.frame
#TDF = read.csv("/outputs/training_data_sample_puebla_p9_abies_pinophyta.csv")
TDF = read.csv("/outputs/training_data_sample_puebla_p9_tyrannidae_birds.csv")

## Order it according to the id of the cell
### This is important because the adjancy matrix rows need to be the same
TDF = TDF[order(TDF$cell_ids),]
# Convert to numeric
TDF = mutate_at(TDF,vars(Dist.to.road_m,Elevation_m,
                         MaxTemp_m,MeanTemp_m,
                         MinTemp_m,Population_m,
                         Precipitation_m,
                         SolarRadiation_m,
                         VaporPres_m,
                         WindSp_m),as.numeric)
# Remove unnecessary symbols in variable names
names(TDF) = lapply(names(TDF),function(x) gsub("_","",x))
names(TDF) = lapply(names(TDF),function(x) gsub("\\.","",x))

##### PREPROCESS
###
# Preprocess for generating pseudo absences
# Change the name of a column that for some reason is called the same
names(TDF)[23] <- 'covid2'


## missing values only in X
#DataFrame = TDF %>% rowwise() %>% 
#            mutate(sample=pseudo_absence_naive(Aves,LUCA),species=pseudo_absence_trivial(Tyrannidae,Aves))

## missing values in X and Y
#DataFrame = TDF %>% rowwise() %>% 
#            mutate(sample=pseudo_absence_naive(Plantae,LUCA),species=pseudo_absence_naive(Pinophyta,Plantae))
DataFrame = TDF %>% rowwise() %>% 
            mutate(sample=pseudo_absence_naive(Aves,LUCA),species=pseudo_absence_naive(Tyrannidae,Aves))

D = apply(M,MARGIN = 1,sum)
### get index with 0 neighbours
idx = which(D == 0)

### select cells with no neighbours
cell_with_no_neighbour = TDF$cellids[idx]
## Remove island for TDF
TDF <- TDF[-c(idx),]
## Erase idx for M and for TDF (Or maybe only for M)
M_bis = M[-c(idx),-c(idx)]

## remove row corresponding to no neighbour (island)
DataFrame <- DataFrame[-c(idx),]

n <- dim(M_bis)[1]
trials <- rep(1,n)

## reformat dataframe
#rr  <- rr[-c(idx),]

# Formula definition
formula_sample= sample~Disttoroadm+Populationm
#+factor(tipos)
#formula_presence= species~Elevationm+MeanTempm
formula_presence= species~Elevationm+Precipitationm



