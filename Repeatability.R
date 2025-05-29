#Anthony Gagliano
#Shearing Repeatability of Covaris tubes using Generalized linear models 

# Clear the list that have in the environment
rm(list=ls())

#### Load the packages ####
library(lme4)
library(lmerTest)

#set working directory
setwd("C:/Users")

## Reading in the data
mydata<-read.csv(file.choose(), header = T, stringsAsFactors = T)

# Verify if the structure of the data frame is okay
str(mydata)

#fitting the data into the model
z <- lmer(y ~ 1 + (1|x), data = mydata)
fitted(z)
summary(z)

#Extracting the variance components
VarCorr(z)

#repeatability using formula and values from VarCorr(z)
29.8700^2/(29.8700^2+7.8243^2)


#Verifying the repeatability of our shearing collectively using code
var_comp <- as.data.frame(VarCorr(z))
varamong <- var_comp[1, "vcov"]  # among-group variance (random effect)
varwithin <- attr(VarCorr(z), "sc")^2  # residual variance (within-group)
repeatability <- varamong / (varamong + varwithin)
repeatability


################################################################################
#Checking the repeatability of each individual:
# Overall variance (variance of the response variable y)
overall_var <- var(mydata$y)

# Get unique individuals
individuals <- unique(mydata$x)
repeatabilities <- numeric(length(individuals))

for (i in seq_along(individuals)) {
  # Subset data for the individual
  individual_data <- subset(mydata, x == individuals[i])
  
  # Calculate within-group variance (residual variance)
  varwithin_ind <- var(individual_data$y)
  
  # Calculate repeatability for the individual
  repeatabilities[i] <- (overall_var - varwithin_ind) / overall_var
}

# Combine results into a data frame
individual_repeatabilities <- data.frame(
  Individual = individuals,
  Repeatability = repeatabilities
)

# View the repeatabilities for each individual
print(individual_repeatabilities)


