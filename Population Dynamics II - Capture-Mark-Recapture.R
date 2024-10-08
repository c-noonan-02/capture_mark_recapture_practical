# Capture-Mark-Recapture Practical

# clear current environment
rm(list=ls())

# load required packages
library(dplyr)
library(tidyr)
library(marked)
library(ggplot2)
library(R2ucare)

# check working directory
getwd()


# 1. GENERATING CAPTURE HISTORIES

# load the sparrow recapture dataset
sparrow_data_long <- read.table('./Data/sparrowrecap.txt', header = TRUE, sep = '\t')
head(sparrow_data_long)

# data exploration
length(unique(sparrow_data_long$id)) # the number of unique individuals in the dataframe
table(sparrow_data_long$sex) # equal number of observations of males and females 
table(sparrow_data_long$year) # captures from 1998-2007
table(sparrow_data_long$island) # at 4 different island locations

# the marked package requires the data to be in wide format
# use tidyr and dplyr packages to generate the correct format

temp_data <- sparrow_data_long[,1:2] # take the first two columns, id and year and put into a temporary dataframe
temp_data$detect <- 1 # add column for detection (all 1s because these represent captures) 

temp_data <- temp_data %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an sampling event
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and year where individuals were not observed
  spread(year, detect, fill = 0) %>% 
  # for every individual....
  group_by(id) %>%
  # paste together 0's and 1's using unite()
  # here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ('tail(names(.),1)')
  # use sep='' so there are no characters separating 0's and 1's
  unite('ch', 2:tail(names(.),1), sep = '')

sparrow_data_wide <- as.data.frame(temp_data) # save to new dataframe
head(sparrow_data_wide)

# we will now add back some information based on the individual IDs, using the match function
sparrow_data_wide$island <- sparrow_data_long$island[match(sparrow_data_wide$id, sparrow_data_long$id)] 
# this creates a new column called island in the sparrow df...
# using the entry from the island column in the sparrow_data_long df... 
# where id in the sparrow df matches the id in the sparrow_data_long df
sparrow_data_wide$sex <- as.factor(sparrow_data_long$sex[match(sparrow_data_wide$id, sparrow_data_long$id)])
sparrow_data_wide <- droplevels(subset(sparrow_data_wide, select = -id)) # remove id column so capture histories appear in first column
head(sparrow_data_wide)


# 2. SIMPLE CORMACK-JOLLY-SEBER MODEL

# CJS model estimates apparent survival and detection probability for open populations. Each is a linear model on a logit scale. 
# Same 'link function' as used in binomial Generalised Linear Model (GLM)
# Uses info from capture histories to estimate detection probailities. 
# Most basic form of the model estimates constant survival and detection probabilities. 

# build basic cmr model
sparrow_mod1 <- crm(sparrow_data_wide)
# examine model and coefficient estimates
sparrow_mod1 
# refit model with precision estimates
sparrow_mod1 <- cjs.hessian(sparrow_mod1)
# results of this are on the logit scale, so must transform them back to the data scale
# using plogis or predict functions
sparrow_mod1$results$reals
plogis(sparrow_mod1$results$beta$Phi)

predict(sparrow_mod1, newdata = data.frame(sex=c('Female', 'Male')), se=T)
# but no covariates in data so new.data argument is not used

