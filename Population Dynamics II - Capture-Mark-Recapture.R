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


# 3. UNEQUAL SAMPLING INTERVALS

# this model assumed an equal time between each capture event
# we can relax this assumption by including a vector of time intervals

sparrow_mod2 <- crm(sparrow_data_wide, time.intervals = c(1,2,1,1,1,1,1,3,4))
sparrow_mod2$results$reals

# question: These models assumed constant survival rates and detection
#           probabilities. Is this a realistic assumption for this system?
#           What term might you include in the model next, and why?

#           These are not realistic. Could modify the model to consider each
#           individual's survival and/or detection probability. 


# 4. INCLUDING STATIC COVARIATES

# We can test if islands within the meta-population differ in the probability of capturing an individual
# In the marked package we can do this in three stages: data processing, building the design matrix, and setting up / executing candidate models

# look into the functions we will be using here
?process.data
?make.design.data

# built in function for data processing
sparrow_proc <- process.data(sparrow_data_wide)
str(sparrow_proc)
head(sparrow_proc[[1]])
head(sparrow_proc$data)
sparrow_matrix <- make.design.data(sparrow_proc) # built in function for building design matrix 
str(sparrow_matrix)
head(sparrow_matrix[[1]])
head(sparrow_matrix$Phi)
# we can then specify the model formulation for the detection probability model
# make p dependent on island
p_island <- list(formula=~island) 
# run the cmr model using newly processed data, desifn matrix, and model specification
# specify model formulation: capture probability depends on island
sparrow_mod3 <- crm(sparrow_proc, 
            sparrow_matrix, 
            model.parameters = list(p = p_island), 
            accumulate=FALSE, hessian = TRUE)
sparrow_mod3$results$reals

# question: Does it look like detection probability varies among islands? Which
#           island has the lowest detection probability? Why might this be? How
#           might you compare this model with our simpler model above to see if#
#           it is a better fit?

#           Comparison of AIC values. Lower = simpler
(sparrow_mod1$results$AIC)
(sparrow_mod3$results$AIC)
# more complex model better fits our data here. 

# task: Test whether survival probabilities differ between the islands
island_phi <- list(formula=~island) # survival probability depends on island
island_p <- list(formula=~island) # capture probability depends on island

sparrow_mod4 <- crm(sparrow_proc, sparrow_matrix,
                    model.parameters = list(Phi = island_phi, p = island_p),
                    accumulate=FALSE, hessian = TRUE)
sparrow_mod4$results$reals

# compare AIC with current model
sparrow_mod3$results$AIC
sparrow_mod4$results$AIC
# similar AIC, model 4 is not much of an improvement so survival rates don't vary


# Do probabilities vary between sexes?
sparrow_proc <- process.data(sparrow_data_wide)
sparrow_matrix <- make.design.data(sparrow_proc)

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.sex <- list(formula=~sex) # survival differs between sexes
  Phi.island <- list(formula=~island) # survival differs between islands
  Phi.sex.island <- list(formula=~sex+island) # survival differs between sexes and islands
  p.dot <- list(formula=~1) # constant detection
  p.sex <- list(formula=~sex) # detection probability differs between sexes
  p.island <- list(formula=~island) # detection probability differs between islands
  p.sex.island <- list(formula=~sex+island) # detection probability differs between sexes and islands
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow_proc, ddl=sparrow_matrix,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}

# run function
sparrow_models <- fit.models()

# display model table
sparrow_models
# top model (lowest AIC) includes a constant model for survival probability,
# but with detection probabilities differing among islands. 
# differences not meaningfully different - so accept simplest model (fewest parameters)

# extract and plot the detection probabilities
sparrow_mod5 <- sparrow_models[[2]]

ggplot(sparrow_mod5$results$reals$p, aes(island, estimate, ymin=lcl, ymax=ucl)) +
  geom_errorbar(width=0.2, colour = "maroon") + geom_point(colour = "maroon") + ylim(0,1) +
  xlab("Island") + ylab("Estimated detection probabilities")

# can also extract sex and island differences in survival and plot those to
# confirm they aren't big differences

sparrow_mod6 <- sparrow_models[[10]]

ggplot(sparrow_mod6$results$reals$Phi, aes(sex, estimate, ymin=lcl, ymax=ucl)) +
  geom_errorbar(width=0.2, colour = "maroon") + geom_point(colour = "maroon") + ylim(0,1) +
  xlab("Sex") + ylab("Estimated survival probabilities")

sparrow_mod7 <- sparrow_models[[6]]

ggplot(sparrow_mod7$results$reals$Phi, aes(island, estimate, ymin=lcl, ymax=ucl)) +
  geom_errorbar(width=0.2, colour = "maroon") + geom_point(colour = "maroon") + ylim(0,1) +
  xlab("Island") + ylab("Estimated survival probabilities")


# 5. INCLUDING TIME-VARYING COVARIATES

# might want to test whether survival or detection probabilities differ#
# depending on a factor thst varies over time - e.g. weather

# add variable 'cold' to the design data, then test if this affects survival
# add new column
sparrow_matrix$Phi$cold <- "Cold"
# add very cold winters between capture events 2 and 3, 5 and 6, and 8 and 9
sparrow_matrix$Phi$cold[sparrow_matrix$Phi$time==2 | sparrow_matrix$Phi$time==5 | sparrow_matrix$Phi$time==8] <- "VeryCold"
head(sparrow_matrix$Phi)

Phi.cold <- list(formula=~cold) 
p.island <- list(formula=~island) 

sparrow_mod8 <- crm(sparrow_proc, 
            sparrow_matrix, 
            model.parameters = list(Phi = Phi.cold, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)

sparrow_mod8$results$reals

# compare AICs
sparrow_mod5$results$AIC
sparrow_mod8$results$AIC

# plot this
ggplot(sparrow_mod8$results$reals$Phi, aes(cold, estimate, ymin=lcl, ymax=ucl)) +
  geom_errorbar(width=0.2, colour = "maroon") + geom_point(colour = "maroon") + ylim(0,1) +
  xlab("Weather Conditions") + ylab("Estimated survival probabilities")

# no evidence that survival depended on this made up variable



