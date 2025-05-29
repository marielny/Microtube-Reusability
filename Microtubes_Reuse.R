# written by DR July 2024

# Clearing instances and resetting environment
rm(list=ls())

## Loading needed libraries
{
  library(ggplot2)
  library(tidyverse)
  library(parallel)
  library(dplyr)
  library(vegan)
  library(tidyr)
  library(emmeans)
  library(stats)
  library(car)
  library(rstatix)
  library(conover.test)
  library(ggpubr)
}

# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
cores <- detectCores(logical = T)

# Setting working directory
setwd("/Users/denis/Library/CloudStorage/Dropbox/McGill/students/PhD students/Marie Launay/CovarisMS/")

# Reading the raw diet data - looking for the file "diet3-20240712.csv"
covdat <- read.csv(file.choose(), header = T, stringsAsFactors = T)

# Viewing and attaching the data
head(covdat)
attach(covdat)

# Using and ifelse to set the -ve values to positive - but whether this
# is correct or not is debatable
covdat$rezo <- ifelse(covdat$zeroed < 0, 0.0001, covdat$zeroed)

# Setting the order of the treatments for consideration
covdat$treatment <- factor(covdat$treatment, levels = c('sheared_DNA', 
                                                        'Milli_Q',
                                                        'Soap', 
                                                        'DNA_Zap', 
                                                        'UV'))

# This will take the log base 10 of the rezo column
covdat$lrz <- log10(covdat$rezo)

# Below uses linear modeling to test for differences but not appropriate
# as data not normal
lmcov <- lm(zeroed ~ treatment, data = covdat)

# plot(lmcov)

# Shapiro-Wilk's demonstrates data do not fit normal distribution
shapiro_test(residuals(lmcov))

# can also plot qqplot to show lack of fit of residuals
ggqqplot(residuals(lmcov))

# Attempts to use BoxCox transformation to get normal distribution 
# but still fails
#library(MASS)
# Boxcox transformation
#bc_res <- boxcox(lmcov, plotit = T)
#bc_rez <- boxcox(covdat$rezo ~ 1, plotit = T)
#bestlam <- bc_rez$x[which.max(bc_rez$y)]
#covdat$tr_rezo <- (covdat$rezo^bestlam -1)/bestlam
#covdat$str_rezo <- sqrt(covdat$tr_rezo)

#hist(covdat$str_rezo)

# Calculate means for each treatment
meansl <- covdat %>%
  group_by(treatment) %>%
  summarise(mean_rezo = mean(lrz, na.rm = TRUE), .groups = 'drop')

meansn <- covdat %>%
  group_by(treatment) %>%
  summarise(mean_o = mean(zeroed, na.rm = TRUE), .groups = 'drop')

# Calculate the stdev for each treatment
csdevl <- covdat %>%
  group_by(treatment) %>%
  summarise(sd_rezo = sd(lrz, na.rm = TRUE), .groups = 'drop')

csdevn <- covdat %>%
  group_by(treatment) %>%
  summarise(sd_o = sd(zeroed, na.rm = TRUE), .groups = 'drop')

# from the above calculate the lower and upper 95% CI range
lwcil <- meansl$mean_rezo - 1.96*(csdevl$sd_rezo/sqrt(20))
upcil <- meansl$mean_rezo + 1.96*(csdevl$sd_rezo/sqrt(20))

lwcin <- meansn$mean_o - 1.96*(csdevn$sd_o/sqrt(20))
upcin <- meansn$mean_o + 1.96*(csdevn$sd_o/sqrt(20))

# Put means, sds, lwci and upci into new df for plotting
cleansuml <- cbind.data.frame(treatment = meansl$treatment, mean = meansl$mean_rezo,
                              lower.CL = lwcil, upper.CL = upcil)

cleansumn <- cbind.data.frame(treatment = meansn$treatment, mean = meansn$mean_o,
                              lower.CL = lwcin, upper.CL = upcin)

# view df
cleansuml
cleansumn


# Custom labels for the x-axis
custom_x_labels <- c(
  "sheared_DNA" = "Sheared DNA",
  "Milli_Q" = "Milli-Q®",
  "DNA_Zap" = "DNA Zap®",
  "Soap" = "Soap",
  "UV" = "UV Light"
)


##############################################################################
# Use ggplot to visualise the data 
##############################################################################

# Use ggplot to plot the data 
ggplot(cleansuml, aes(x = treatment, y = mean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1, colour = "grey49") +
  geom_point(shape = 19, size = 6, colour = "grey49") +
  geom_jitter(covdat, mapping = aes(x = treatment, y = lrz, color = treatment), 
              position = position_jitter(0.3), size = 6, stroke = 1, alpha = 0.6) +
  geom_line(data = cleansuml, aes(x = treatment, y = mean, group = 1), 
            color = "grey49", linewidth = 1, linetype = "solid", alpha = 0.6) + 
  labs(x = "Cleaning Step", y = expression(log[10](DNA~concentration~ng/µl)), color = "Treatment") +
  scale_color_manual(values = c("springgreen3", "skyblue3", "darkorange2","darkorchid4","firebrick4")) +
  scale_y_continuous(limits = c(-5, 2), breaks = seq(-5, 2, by = 1)) +
  scale_x_discrete(labels = custom_x_labels) +  # Adding custom x-axis labels
  theme(aspect.ratio = 0.80) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.position = "none")

# Use ggplot to plot the data 
ggplot(cleansumn, aes(x = treatment, y = mean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1, colour = "grey49") +
  geom_point(shape = 19, size = 6, colour = "grey49") +
  geom_jitter(covdat, mapping = aes(x = treatment, y = zeroed, color = treatment), 
              position = position_jitter(0.3), size = 6, stroke = 1, alpha = 0.6) +
  geom_line(data = cleansumn, aes(x = treatment, y = mean, group = 1), 
            color = "grey49", linewidth = 1, linetype = "solid", alpha = 0.6) + 
  labs(x = "Cleaning Step", y = "DNA concn. (ng/µl)", color = "treatment") +
  scale_color_manual(values = c("springgreen3", "skyblue3", "darkorange2","darkorchid4","firebrick4")) +
  scale_y_continuous(limits = c(-5, 50), breaks = seq(-5, 50, by = 5)) +
  scale_x_discrete(labels = custom_x_labels) +  # Adding custom x-axis labels
  theme(aspect.ratio = 0.80) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.position = "none")


##############################################################################
# Statistical Analyses 
##############################################################################

#Lets try some statistical tests now

## Boundary exceeded in many cases, so data is deviating from normal
## We can therefore use a Kruskall-Walis test for non-parametric 
## equivalent of testing for differences among groups.

# Check the structure of covdat
str(covdat)

kruskal.test(lrz ~ treatment, data = covdat)

kruskal.test(zeroed ~ treatment, data = covdat)


# Perform the Conover-Iman T-test corrected for multiple comparisons.
conover.test(covdat$lrz, covdat$treatment, method="by", list=TRUE)

conover.test(covdat$zeroed, covdat$treatment, method="by", list=TRUE)

# Perform Dunn's test
result <- dunn_test(lrz ~ treatment, data = covdat, p.adjust.method = "BY")

result2 <- dunn_test(zeroed ~ treatment, data = covdat, p.adjust.method = "BY")

# Print the result
print(result)
print(result2)
