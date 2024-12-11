
#set working directory and load packages

library(emmeans)
library(drc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

pamdata<-read.csv2("C:/Users/Christian W/Desktop/Seafile/Meine Bibliothek/KAUST/Fieldwork CBASS/CBASS data/CBASSJune_R_new.csv",header=T)
str(pamdata)

pamdata <- na.omit(pamdata)

pamdata$Geno<-as.factor(pamdata$Geno)
pamdata$ReefStatus<-as.factor(pamdata$ReefStatus)
pamdata$Reef<-as.factor(pamdata$Reef)
pamdata$Location<-as.factor(pamdata$Location)
pamdata$Species<-as.factor(pamdata$Species)
pamdata$PAM8New<-as.numeric(pamdata$PAM8)



# Subset for testing
offshore <- subset(pamdata, ReefStatus == "Offshore")
nearshore <- subset(pamdata, ReefStatus == "Near")

AlLith <- subset(pamdata, Location == "AlLith")
KAUST <- subset(pamdata, Location == "KAUST")

AlLithoff <- subset(AlLith, ReefStatus == "Offshore")
AlLithnear <- subset(AlLith, ReefStatus == "Near")

KAUSToff <- subset(KAUST, ReefStatus == "Offshore")
KAUSTnear <- subset(KAUST, ReefStatus == "Near")

Acro <- subset(pamdata, Species == "Acropora")
Acrooff <- subset(Acro, ReefStatus == "Offshore")
Acronear <- subset(Acro, ReefStatus == "Near")

Poci <- subset(pamdata, Species == "Pocillopora")
Pocioff <- subset(Poci, ReefStatus == "Offshore")
Pocinear <- subset(Poci, ReefStatus == "Near")

Turf <- subset(pamdata, Species == "Turf")
Turfoff <- subset(Turf, ReefStatus == "Offshore")
Turfnear <- subset(Turf, ReefStatus == "Near")

Xeni <- subset(pamdata, Species == "Xenia")
Xenioff <- subset(Xeni, ReefStatus == "Offshore")
Xeninear <- subset(Xeni, ReefStatus == "Near")





# Building Dose Response Curves (excluding limits for overall DRC) ----
#### Al Lith vs KAUST offshore comparison per species ####

#### Acropora offshore Al Lith vs KAUST ####


#### fit to each treatment individually - Acropora
DRCpamACROoffAl = drm(PAM8 ~ Temp, data = Acrooff[Acrooff$Location=="AlLith",],
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamACROoffAl)
DRCpamACROoffAl$coefficients[3]
ED(DRCpamACROoffAl, c(50))


# genotype-specific curve fits - Acropora
DRCpamACROgenooffAl = drm(PAM8 ~ Temp, data = Acrooff[Acrooff$Location=="AlLith",], curveid=Geno,
                          fct = LL.3(names = c('hill', 'max', 'ed50')),
                          upperl = c(120, 0.72, 40),
                          lowerl = c(10, 0.55, 30))
summary(DRCpamACROgenooffAl)
DRCpamACROgenooffAl$coefficients[31:45]
ED(DRCpamACROgenooffAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamACROgenooffAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


#### Acropora offshore KAUST


#### fit to each treatment individually - Acropora
DRCpamACROoffKAUST = drm(PAM8 ~ Temp, data = Acrooff[Acrooff$Location=="KAUST",],
                         fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamACROoffKAUST)
DRCpamACROoffKAUST$coefficients[3]
ED(DRCpamACROoffKAUST, c(50))


# genotype-specific curve fits - Acropora
DRCpamACROgenooffKAUST = drm(PAM8 ~ Temp, data = Acrooff[Acrooff$Location=="KAUST",], curveid=Geno,
                             fct = LL.3(names = c('hill', 'max', 'ed50')),
                             upperl = c(120, 0.72, 40),
                             lowerl = c(10, 0.55, 30))
summary(DRCpamACROgenooffKAUST)
DRCpamACROgenooffKAUST$coefficients[31:45]
ED(DRCpamACROgenooffKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamACROgenooffKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


#### Poci offshore Al Lith vs KAUST ####


#### fit to each treatment individually - Poci
DRCpamPOCIoffAl = drm(PAM8 ~ Temp, data = Pocioff[Pocioff$Location=="AlLith",],
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPOCIoffAl)
DRCpamPOCIoffAl$coefficients[3]
ED(DRCpamPOCIoffAl, c(50))


# genotype-specific curve fits - Poci
DRCpamPOCIgenooffAl = drm(PAM8 ~ Temp, data = Pocioff[Pocioff$Location=="AlLith",], curveid=Geno,
                          fct = LL.3(names = c('hill', 'max', 'ed50')),
                          upperl = c(120, 0.72, 40),
                          lowerl = c(10, 0.55, 30))
summary(DRCpamPOCIgenooffAl)
DRCpamPOCIgenooffAl$coefficients[31:45]
ED(DRCpamPOCIgenooffAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamPOCIgenooffAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


#### Pocillopora offshore KAUST


#### fit to each treatment individually - Poci
DRCpamPOCIoffKAUST = drm(PAM8 ~ Temp, data = Pocioff[Pocioff$Location=="KAUST",],
                         fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPOCIoffKAUST)
DRCpamPOCIoffKAUST$coefficients[3]
ED(DRCpamPOCIoffKAUST, c(50))


# genotype-specific curve fits - Poci
DRCpamPOCIgenooffKAUST = drm(PAM8 ~ Temp, data = Pocioff[Pocioff$Location=="KAUST",], curveid=Geno,
                             fct = LL.3(names = c('hill', 'max', 'ed50')),
                             upperl = c(120, 0.72, 40),
                             lowerl = c(10, 0.55, 30))
summary(DRCpamPOCIgenooffKAUST)
DRCpamPOCIgenooffKAUST$coefficients[31:45]
ED(DRCpamPOCIgenooffKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamPOCIgenooffKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd





#### Turf offshore Al Lith vs KAUST ####


#### fit to each treatment individually - Turf
DRCpamTURFoffAl = drm(PAM8 ~ Temp, data = Turfoff[Turfoff$Location=="AlLith",],
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamTURFoffAl)
DRCpamTURFoffAl$coefficients[3]
ED(DRCpamTURFoffAl, c(50))


# genotype-specific curve fits - Turf
DRCpamTURFgenooffAl = drm(PAM8 ~ Temp, data = Turfoff[Turfoff$Location=="AlLith",], curveid=Geno,
                          fct = LL.3(names = c('hill', 'max', 'ed50')),
                          upperl = c(120, 0.72, 40),
                          lowerl = c(10, 0.55, 30))
summary(DRCpamTURFgenooffAl)
DRCpamTURFgenooffAl$coefficients[31:45]
ED(DRCpamTURFgenooffAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamTURFgenooffAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


# test against a null model----
# Subset the data for the specific location
subset_data <- Turfoff[Turfoff$Location == "AlLith", ]
str(subset_data)
subset_data$PAM8 <- as.numeric(subset_data$PAM8)

linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())

# Perform a likelihood ratio test
anova_result <- anova(linear_model, full_model, test = "Chisq")

# Display the results
print(anova_result)


# Example: Fitting and comparing models
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())
linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
exponential_model <- drm(PAM8 ~ Temp, data = subset_data, fct = EXD.3())


# Calculate AIC for each model
aic_full <- AIC(full_model)
aic_linear <- AIC(linear_model)
aic_exponential <- AIC(exponential_model)

# Compare the AIC values
aic_comparison <- data.frame(
  Model = c("Full Model (LL.3)", "Linear Model (LL.2)", "Exponential Model (EXD.3)"),
  AIC = c(aic_full, aic_linear, aic_exponential)
)

# Display the comparison table
print(aic_comparison)

# t test to see if there are differences with temperature
# Subset the data to include only the 30°C and 39°C treatments
subset_data_temp <- subset(subset_data, Temp %in% c(30, 39))

# Perform a t-test to compare the response between 30°C and 39°C
t_test_result <- t.test(PAM8 ~ Temp, data = subset_data_temp)

# Display the results
print(t_test_result)


#### Turf Algae offshore KAUST


#### fit to each treatment individually - Turf
DRCpamTURFoffKAUST = drm(PAM8 ~ Temp, data = Turfoff[Turfoff$Location=="KAUST",],
                         fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamTURFoffKAUST)
DRCpamTURFoffKAUST$coefficients[3]
ED(DRCpamTURFoffKAUST, c(50))


# genotype-specific curve fits - Turf
DRCpamTURFgenooffKAUST = drm(PAM8 ~ Temp, data = Turfoff[Turfoff$Location=="KAUST",], curveid=Geno,
                             fct = LL.3(names = c('hill', 'max', 'ed50')),
                             upperl = c(120, 0.72, 40),
                             lowerl = c(10, 0.55, 30))
summary(DRCpamTURFgenooffKAUST)
DRCpamTURFgenooffKAUST$coefficients[31:45]
ED(DRCpamTURFgenooffKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamTURFgenooffKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


# test against a null model----
# Subset the data for the specific location
subset_data <- Turfoff[Turfoff$Location == "KAUST", ]
str(subset_data)
subset_data$PAM8 <- as.numeric(subset_data$PAM8)

linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())

# Perform a likelihood ratio test
anova_result <- anova(linear_model, full_model, test = "Chisq")

# Display the results
print(anova_result)


# Example: Fitting and comparing models
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())
linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
exponential_model <- drm(PAM8 ~ Temp, data = subset_data, fct = EXD.3())


# Calculate AIC for each model
aic_full <- AIC(full_model)
aic_linear <- AIC(linear_model)
aic_exponential <- AIC(exponential_model)

# Compare the AIC values
aic_comparison <- data.frame(
  Model = c("Full Model (LL.3)", "Linear Model (LL.2)", "Exponential Model (EXD.3)"),
  AIC = c(aic_full, aic_linear, aic_exponential)
)

# Display the comparison table
print(aic_comparison)

# t test to see if there are differences with temperature
# Subset the data to include only the 30°C and 39°C treatments
subset_data_temp <- subset(subset_data, Temp %in% c(30, 39))

# Perform a t-test to compare the response between 30°C and 39°C
t_test_result <- t.test(PAM8 ~ Temp, data = subset_data_temp)

# Display the results
print(t_test_result)





#### Xenia offshore Al Lith vs KAUST - NOT EXISTING!!!! ####
# no Xenia data from Al Lith available

#### Xenia offshore KAUST


#### fit to each treatment individually - Xenia
DRCpamXENIoffKAUST = drm(PAM8 ~ Temp, data = Xenioff[Xenioff$Location=="KAUST",],
                         fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamXENIoffKAUST)
DRCpamXENIoffKAUST$coefficients[3]
ED(DRCpamXENIoffKAUST, c(50))


# genotype-specific curve fits - Xenia
DRCpamXENIgenooffKAUST = drm(PAM8 ~ Temp, data = Xenioff[Xenioff$Location=="KAUST",], curveid=Geno,
                             fct = LL.3(names = c('hill', 'max', 'ed50')),
                             upperl = c(120, 0.72, 40),
                             lowerl = c(10, 0.55, 30))
summary(DRCpamXENIgenooffKAUST)
DRCpamXENIgenooffKAUST$coefficients[31:45]
ED(DRCpamXENIgenooffKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamXENIgenooffKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd






#### Al Lith vs KAUST nearshore comparison per species ####

#### Acropora nearshore Al Lith vs KAUST ####


#### fit to each treatment individually - Acropora
DRCpamACROnearAl = drm(PAM8 ~ Temp, data = Acronear[Acronear$Location=="AlLith",],
                       fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamACROnearAl)
DRCpamACROnearAl$coefficients[3]
ED(DRCpamACROnearAl, c(50))


# genotype-specific curve fits - Acropora
DRCpamACROgenonearAl = drm(PAM8 ~ Temp, data = Acronear[Acronear$Location=="AlLith",], curveid=Geno,
                           fct = LL.3(names = c('hill', 'max', 'ed50')),
                           upperl = c(120, 0.72, 40),
                           lowerl = c(10, 0.55, 30))
summary(DRCpamACROgenonearAl)
DRCpamACROgenonearAl$coefficients[31:45]
ED(DRCpamACROgenonearAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamACROgenonearAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


#### Acropora nearshore KAUST


#### fit to each treatment individually - Acropora
DRCpamACROnearKAUST = drm(PAM8 ~ Temp, data = Acronear[Acronear$Location=="KAUST",],
                          fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamACROnearKAUST)
DRCpamACROnearKAUST$coefficients[3]
ED(DRCpamACROnearKAUST, c(50))


# genotype-specific curve fits - Acropora
DRCpamACROgenonearKAUST = drm(PAM8 ~ Temp, data = Acronear[Acronear$Location=="KAUST",], curveid=Geno,
                              fct = LL.3(names = c('hill', 'max', 'ed50')),
                              upperl = c(120, 0.72, 40),
                              lowerl = c(10, 0.55, 30))
summary(DRCpamACROgenonearKAUST)
DRCpamACROgenonearKAUST$coefficients[31:45]
ED(DRCpamACROgenonearKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamACROgenonearKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd




#### Poci nearshore Al Lith vs KAUST ####


#### fit to each treatment individually - Poci
DRCpamPOCInearAl = drm(PAM8 ~ Temp, data = Pocinear[Pocinear$Location=="AlLith",],
                       fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPOCInearAl)
DRCpamPOCInearAl$coefficients[3]
ED(DRCpamPOCInearAl, c(50))


# genotype-specific curve fits - Poci
DRCpamPOCIgenonearAl = drm(PAM8 ~ Temp, data = Pocinear[Pocinear$Location=="AlLith",], curveid=Geno,
                           fct = LL.3(names = c('hill', 'max', 'ed50')),
                           upperl = c(120, 0.72, 40),
                           lowerl = c(10, 0.55, 30))
summary(DRCpamPOCIgenonearAl)
DRCpamPOCIgenonearAl$coefficients[31:45]
ED(DRCpamPOCIgenonearAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamPOCIgenonearAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd



#### Pocillopora nearshore KAUST


#### fit to each treatment individually - Poci
DRCpamPOCInearKAUST = drm(PAM8 ~ Temp, data = Pocinear[Pocinear$Location=="KAUST",],
                          fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPOCInearKAUST)
DRCpamPOCInearKAUST$coefficients[3]
ED(DRCpamPOCInearKAUST, c(50))


# genotype-specific curve fits - Poci
DRCpamPOCIgenonearKAUST = drm(PAM8 ~ Temp, data = Pocinear[Pocinear$Location=="KAUST",], curveid=Geno,
                              fct = LL.3(names = c('hill', 'max', 'ed50')),
                              upperl = c(120, 0.72, 40),
                              lowerl = c(10, 0.55, 30))
summary(DRCpamPOCIgenonearKAUST)
DRCpamPOCIgenonearKAUST$coefficients[31:45]
ED(DRCpamPOCIgenonearKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamPOCIgenonearKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd



#### Turf nearshore Al Lith vs KAUST ####


#### fit to each treatment individually - Turf
DRCpamTURFnearAl = drm(PAM8 ~ Temp, data = Turfnear[Turfnear$Location=="AlLith",],
                       fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamTURFnearAl)
DRCpamTURFnearAl$coefficients[3]
ED(DRCpamTURFnearAl, c(50))


# genotype-specific curve fits - Turf
DRCpamTURFgenonearAl = drm(PAM8 ~ Temp, data = Turfnear[Turfnear$Location=="AlLith",], curveid=Geno,
                           fct = LL.3(names = c('hill', 'max', 'ed50')),
                           upperl = c(120, 0.72, 40),
                           lowerl = c(10, 0.55, 30))
summary(DRCpamTURFgenonearAl)
DRCpamTURFgenonearAl$coefficients[31:45]
ED(DRCpamTURFgenonearAl, c(50))

ed50_results <- as.data.frame(ED(DRCpamTURFgenonearAl, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


# test against a null model----
# Subset the data for the specific location
subset_data <- Turfnear[Turfnear$Location == "AlLith", ]
str(subset_data)
subset_data$PAM8 <- as.numeric(subset_data$PAM8)

linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())

# Perform a likelihood ratio test
anova_result <- anova(linear_model, full_model, test = "Chisq")

# Display the results
print(anova_result)


# Example: Fitting and comparing models
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())
linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
exponential_model <- drm(PAM8 ~ Temp, data = subset_data, fct = EXD.3())


# Calculate AIC for each model
aic_full <- AIC(full_model)
aic_linear <- AIC(linear_model)
aic_exponential <- AIC(exponential_model)

# Compare the AIC values
aic_comparison <- data.frame(
  Model = c("Full Model (LL.3)", "Linear Model (LL.2)", "Exponential Model (EXD.3)"),
  AIC = c(aic_full, aic_linear, aic_exponential)
)

# Display the comparison table
print(aic_comparison)

# t test to see if there are differences with temperature
# Subset the data to include only the 30°C and 39°C treatments
subset_data_temp <- subset(subset_data, Temp %in% c(30, 39))

# Perform a t-test to compare the response between 30°C and 39°C
t_test_result <- t.test(PAM8 ~ Temp, data = subset_data_temp)

# Display the results
print(t_test_result)



#### Turf Algae nearshore KAUST


#### fit to each treatment individually - Turf
DRCpamTURFnearKAUST = drm(PAM8 ~ Temp, data = Turfnear[Turfnear$Location=="KAUST",],
                          fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamTURFnearKAUST)
DRCpamTURFnearKAUST$coefficients[3]
ED(DRCpamTURFnearKAUST, c(50))


# genotype-specific curve fits - Turf
DRCpamTURFgenonearKAUST = drm(PAM8 ~ Temp, data = Turfnear[Turfnear$Location=="KAUST",], curveid=Geno,
                              fct = LL.3(names = c('hill', 'max', 'ed50')),
                              upperl = c(120, 0.72, 40),
                              lowerl = c(10, 0.55, 30))
summary(DRCpamTURFgenonearKAUST)
DRCpamTURFgenonearKAUST$coefficients[31:45]
ED(DRCpamTURFgenonearKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamTURFgenonearKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd


# test against a null model----
# Subset the data for the specific location
subset_data <- Turfnear[Turfnear$Location == "KAUST", ]
str(subset_data)
subset_data$PAM8 <- as.numeric(subset_data$PAM8)

linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())

# Perform a likelihood ratio test
anova_result <- anova(linear_model, full_model, test = "Chisq")

# Display the results
print(anova_result)


# Example: Fitting and comparing models
full_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.3())
linear_model <- drm(PAM8 ~ Temp, data = subset_data, fct = LL.2())
exponential_model <- drm(PAM8 ~ Temp, data = subset_data, fct = EXD.3())


# Calculate AIC for each model
aic_full <- AIC(full_model)
aic_linear <- AIC(linear_model)
aic_exponential <- AIC(exponential_model)

# Compare the AIC values
aic_comparison <- data.frame(
  Model = c("Full Model (LL.3)", "Linear Model (LL.2)", "Exponential Model (EXD.3)"),
  AIC = c(aic_full, aic_linear, aic_exponential)
)

# Display the comparison table
print(aic_comparison)

# t test to see if there are differences with temperature
# Subset the data to include only the 30°C and 39°C treatments
subset_data_temp <- subset(subset_data, Temp %in% c(30, 39))

# Perform a t-test to compare the response between 30°C and 39°C
t_test_result <- t.test(PAM8 ~ Temp, data = subset_data_temp)

# Display the results
print(t_test_result)






#### Xenia nearshore Al Lith vs KAUST - NOT EXISTING!!!! ####

# Xeia was not present Al Lith

#### Xenia nearshore KAUST


#### fit to each treatment individually - Xenia
DRCpamXENInearKAUST = drm(PAM8 ~ Temp, data = Xeninear[Xeninear$Location=="KAUST",],
                          fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamXENInearKAUST)
DRCpamXENInearKAUST$coefficients[3]
ED(DRCpamXENInearKAUST, c(50))


# genotype-specific curve fits - Xenia
DRCpamXENIgenonearKAUST = drm(PAM8 ~ Temp, data = Xeninear[Xeninear$Location=="KAUST",], curveid=Geno,
                              fct = LL.3(names = c('hill', 'max', 'ed50')),
                              upperl = c(120, 0.72, 40),
                              lowerl = c(10, 0.55, 30))
summary(DRCpamXENIgenonearKAUST)
DRCpamXENIgenonearKAUST$coefficients[31:45]
ED(DRCpamXENIgenonearKAUST, c(50))

ed50_results <- as.data.frame(ED(DRCpamXENIgenonearKAUST, c(50)))
mean_sd <- c(mean = mean(ed50_results$Estimate), sd = sd(ed50_results$Estimate))
mean_sd










#### Comparison Species Offshore Separated by Location ####

#### Al Lith Offshore Species ----
#### Merging Coeffs of all species in offshore Al Lith
Coeffs<-c(DRCpamACROoffAl$coefficients[3], DRCpamPOCIoffAl$coefficients[3], DRCpamTURFoffAl$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamACROgenooffAl$coefficients[31:45], DRCpamPOCIgenooffAl$coefficients[31:45], DRCpamTURFgenooffAl$coefficients[31:45]))
GenoCoeffs$Species=c(rep("Acropora",15), rep("Pocillopora",15), rep("Turf",15))
aggregate(ED50 ~ Species, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Species, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Species, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Species, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Species, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Species"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
#write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="./CBASS_cold_ED50_genets.txt", quote=F, sep="\t", row.names=F)
#write.table(SumaryStats, file="./CBASS_cold_ED50_summarystats.txt", quote=F, sep="\t", row.names=F)

#### Plotting
temp_x<- seq(30, 40, length = 100)
#pdf("./CBASS_cold_ED50_colonies.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)


i<-1 #Acropora offshore
matplot(temp_x, predict(DRCpamACROoffAl, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col="#ffc425",lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(30,40),ylim=c(0,0.9), cex.axis=1.5, cex.lab=1.5)
with(AlLithoff[AlLithoff$Species=="Acropora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#ffc425", cex=1.5))

i<-2 #Pocillopora offshore
matpoints(temp_x, predict(DRCpamPOCIoffAl, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#fb8072",lty=c(1,3,3),lwd=line_width)
with(AlLithoff[AlLithoff$Species=="Pocillopora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#fb8072", cex=1.5))

i<-3 #Turf Algae offshore
matpoints(temp_x, predict(DRCpamTURFoffAl, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#00B4AB",lty=c(1,3,3),lwd=line_width)
with(AlLithoff[AlLithoff$Species=="Turf",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#00B4AB", cex=1.5))


legend("bottomleft",
       legend = c(expression(italic("Acropora") ~ "spp. (a)"),
                  expression(italic("P.") ~ "cf." ~ italic("verrucosa") ~ "(a)"),
                  "Turf algae (b)"),
       pch = c(18, 18, 18),
       col = c("#ffc425", "#fb8072", "#00B4AB"),
       pt.cex = 1.5, bty = "n", cex = 1.1)
title(main="Southern Red Sea - Offshore", cex.main = 1.8)
abline(v=SumaryStats$MeanED50, col=c("#fb8072","#ffc425", "#00B4AB"), lwd=3)
text(paste(SumaryStats$MeanED50), pos=c(2,4,4), c(0.9,0.8, 0.8),labels=paste(as.character(round(SumaryStats$MeanED50, digits=2)),"°C",sep=""),col=c("#fb8072","#ffc425", "#00B4AB"), cex=1.3)
mtext("C", side = 3, line = 2, cex = 2.5, at = par("usr")[1], adj = 2.5)  # Place "A" in the upper left corner outside of the plot




#### Al Lith Nearshore Species ----
#### Merging Coeffs of all species in offshore Al Lith
Coeffs<-c(DRCpamACROnearAl$coefficients[3], DRCpamPOCInearAl$coefficients[3], DRCpamTURFnearAl$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamACROgenonearAl$coefficients[31:45], DRCpamPOCIgenonearAl$coefficients[31:45], DRCpamTURFgenonearAl$coefficients[31:45]))
GenoCoeffs$Species=c(rep("Acropora",15), rep("Pocillopora",15), rep("Turf",15))
aggregate(ED50 ~ Species, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Species, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Species, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Species, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Species, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Species"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
#write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="./CBASS_cold_ED50_genets.txt", quote=F, sep="\t", row.names=F)
#write.table(SumaryStats, file="./CBASS_cold_ED50_summarystats.txt", quote=F, sep="\t", row.names=F)

#### Plotting
temp_x<- seq(30, 40, length = 100)
#pdf("./CBASS_cold_ED50_colonies.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)


i<-1 #Acropora nearshore
matplot(temp_x, predict(DRCpamACROnearAl, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col="#ffc425",lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(30,40),ylim=c(0,0.9), cex.axis=1.5, cex.lab=1.5)
with(AlLithnear[AlLithnear$Species=="Acropora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#ffc425", cex=1.5))

i<-2 #Pocillopora nearshore
matpoints(temp_x, predict(DRCpamPOCInearAl, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#fb8072",lty=c(1,3,3),lwd=line_width)
with(AlLithnear[AlLithnear$Species=="Pocillopora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#fb8072", cex=1.5))

i<-3 #Turf Algae nearshore
matpoints(temp_x, predict(DRCpamTURFnearAl, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#00B4AB",lty=c(1,3,3),lwd=line_width)
with(AlLithnear[AlLithnear$Species=="Turf",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#00B4AB", cex=1.5))


legend("bottomleft",
       legend = c(expression(italic("Acropora") ~ "spp. (a)"),
                  expression(italic("P.") ~ "cf." ~ italic("verrucosa") ~ "(a)"),
                  "Turf algae (b)"),
       pch = c(18, 18, 18),
       col = c("#ffc425", "#fb8072", "#00B4AB"),
       pt.cex = 1.5, bty = "n", cex = 1.1)
title(main="Southern Red Sea - Nearshore", cex.main = 1.8)
abline(v=SumaryStats$MeanED50, col=c("#ffc425","#fb8072", "#FFFFFF00"), lwd=3)
text(paste(SumaryStats$MeanED50), pos=c(2,4,4), c(0.9,0.8, 0.85),labels=paste(as.character(round(SumaryStats$MeanED50, digits=2)),"°C",sep=""),col=c("#ffc425","#fb8072", "#FFFFFF00"), cex=1.3)
mtext("A", side = 3, line = 2, cex = 2.5, at = par("usr")[1], adj = 2.5)





#### Thuwal Offshore Species ----
#### Merging Coeffs of all species in offshore
Coeffs<-c(DRCpamACROoffKAUST$coefficients[3], DRCpamPOCIoffKAUST$coefficients[3], DRCpamTURFoffKAUST$coefficients[3], DRCpamXENIoffKAUST$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamACROgenooffKAUST$coefficients[31:45], DRCpamPOCIgenooffKAUST$coefficients[31:45], DRCpamTURFgenooffKAUST$coefficients[31:45], DRCpamXENIgenooffKAUST$coefficients[31:45]))
GenoCoeffs$Species=c(rep("Acropora",15), rep("Pocillopora",15), rep("Turf",15), rep("Xenia",15))
aggregate(ED50 ~ Species, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Species, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Species, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Species, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Species, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Species"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
#write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="./CBASS_cold_ED50_genets.txt", quote=F, sep="\t", row.names=F)
#write.table(SumaryStats, file="./CBASS_cold_ED50_summarystats.txt", quote=F, sep="\t", row.names=F)

#### Plotting
temp_x<- seq(30, 40, length = 100)
#pdf("./CBASS_cold_ED50_colonies.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)


i<-1 #Acropora offshore
matplot(temp_x, predict(DRCpamACROoffKAUST, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col="#ffc425",lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(30,40),ylim=c(0,0.9), cex.axis=1.5, cex.lab=1.5)
with(KAUSToff[KAUSToff$Species=="Acropora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#ffc425", cex=1.5))

i<-2 #Pocillopora offshore
matpoints(temp_x, predict(DRCpamPOCIoffKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#fb8072",lty=c(1,3,3),lwd=line_width)
with(KAUSToff[KAUSToff$Species=="Pocillopora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#fb8072", cex=1.5))

i<-3 #Turf Algae offshore
matpoints(temp_x, predict(DRCpamTURFoffKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#00B4AB",lty=c(1,3,3),lwd=line_width)
with(KAUSToff[KAUSToff$Species=="Turf",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#00B4AB", cex=1.5))

i<-4 #Xenia offshore
matpoints(temp_x, predict(DRCpamXENIoffKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#4575b4",lty=c(1,3,3),lwd=line_width)
with(KAUSToff[KAUSToff$Species=="Xenia",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#4575b4", cex=1.5))


legend("bottomleft",
       legend = c(expression(italic("Acropora") ~ "spp. (a)"),
                  expression(italic("P.") ~ "cf." ~ italic("verrucosa") ~ "(b)"),
                  expression(italic("Xenia") ~ "spp. (c)"),
                  "Turf algae (d)"),
       pch = c(18, 18, 18, 18),
       col = c("#ffc425", "#fb8072", "#4575b4", "#00B4AB"),
       pt.cex = 1.5, bty = "n", cex = 1.1)
title(main="Central Red Sea - Offshore", cex.main = 1.8)
abline(v=SumaryStats$MeanED50, col=c("#ffc425","#fb8072", "#4575b4", "#FFFFFF00"), lwd=3)
text(paste(SumaryStats$MeanED50), pos=c(4,4,2,2), c(0.87,0.8, 0.87, 0.8),labels=paste(as.character(round(SumaryStats$MeanED50, digits=2)),"°C",sep=""),col=c("#ffc425","#fb8072", "#4575b4", "#FFFFFF00"), cex=1.3)
mtext("D", side = 3, line = 2, cex = 2.5, at = par("usr")[1], adj = 2.5)







#### Thuwal Nearshore Species ----
#### Merging Coeffs of all species in offshore
Coeffs<-c(DRCpamACROnearKAUST$coefficients[3], DRCpamPOCInearKAUST$coefficients[3], DRCpamTURFnearKAUST$coefficients[3], DRCpamXENInearKAUST$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamACROgenonearKAUST$coefficients[31:45], DRCpamPOCIgenonearKAUST$coefficients[31:45], DRCpamTURFgenonearKAUST$coefficients[31:45], DRCpamXENIgenonearKAUST$coefficients[31:45]))
GenoCoeffs$Species=c(rep("Acropora",15), rep("Pocillopora",15), rep("Turf",15), rep("Xenia",15))
aggregate(ED50 ~ Species, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Species, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Species, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Species, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Species, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Species"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Species, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Species, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
#write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="./CBASS_cold_ED50_genets.txt", quote=F, sep="\t", row.names=F)
#write.table(SumaryStats, file="./CBASS_cold_ED50_summarystats.txt", quote=F, sep="\t", row.names=F)

#### Plotting
temp_x<- seq(30, 40, length = 100)
#pdf("./CBASS_cold_ED50_colonies.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)


i<-1 #Acropora nearshore
matplot(temp_x, predict(DRCpamACROnearKAUST, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col="#ffc425",lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(30,40),ylim=c(0,0.9), cex.axis=1.5, cex.lab=1.5)
with(KAUSTnear[KAUSTnear$Species=="Acropora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#ffc425", cex=1.5))

i<-2 #Pocillopora nearshore
matpoints(temp_x, predict(DRCpamPOCInearKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#fb8072",lty=c(1,3,3),lwd=line_width)
with(KAUSTnear[KAUSTnear$Species=="Pocillopora",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#fb8072", cex=1.5))

i<-3 #Turf Algae nearshore
matpoints(temp_x, predict(DRCpamTURFnearKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#00B4AB",lty=c(1,3,3),lwd=line_width)
with(KAUSTnear[KAUSTnear$Species=="Turf",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#00B4AB", cex=1.5))

i<-4 #Xenia nearshore
matpoints(temp_x, predict(DRCpamXENInearKAUST, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="#4575b4",lty=c(1,3,3),lwd=line_width)
with(KAUSTnear[KAUSTnear$Species=="Xenia",],matpoints(Temp-offsets[i],PAM8,pch=18, col="#4575b4", cex=1.5))


legend("bottomleft",
       legend = c(expression(italic("Acropora") ~ "spp. (a)"),
                  expression(italic("P.") ~ "cf." ~ italic("verrucosa") ~ "(a)"),
                  expression(italic("Xenia") ~ "spp. (b)"),
                  "Turf algae (c)"),
       pch = c(18, 18, 18, 18),
       col = c("#ffc425", "#fb8072", "#4575b4", "#00B4AB"),
       pt.cex = 1.5, bty = "n", cex = 1.1)
title(main="Central Red Sea - Nearshore", cex.main = 1.8)
abline(v=SumaryStats$MeanED50, col=c("#ffc425","#fb8072", "#4575b4", "#FFFFFF00"), lwd=3)
text(paste(SumaryStats$MeanED50), pos=c(2,4,2,2), c(0.85,0.7, 0.87, 0.79),labels=paste(as.character(round(SumaryStats$MeanED50, digits=2)),"°C",sep=""),col=c("#ffc425","#fb8072", "#4575b4", "#FFFFFF00"), cex=1.3)
mtext("B", side = 3, line = 2, cex = 2.5, at = par("usr")[1], adj = 2.5)









#### Acropora Reef Location Comparison ----
#### Merging Coeffs of all species in offshore
Coeffs<-c(DRCpamACROnearKAUST$coefficients[3], DRCpamACROoffKAUST$coefficients[3], DRCpamACROnearAl$coefficients[3], DRCpamAcrooffAl$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamACROgenonearKAUST$coefficients[31:45], DRCpamACROgenooffKAUST$coefficients[31:45], DRCpamACROgenonearAl$coefficients[31:45], DRCpamACROgenooffAl$coefficients[31:45]))
GenoCoeffs$Loc=c(rep("nearKAUST",15), rep("offKAUST",15), rep("nearAl",15), rep("offAl",15))
GenoCoeffs$Reef=c(rep("nearshore",15), rep("offshore",15), rep("nearshore",15), rep("offshore",15))
GenoCoeffs$Site=c(rep("Thuwal",30), rep("Al Lith",30))
aggregate(ED50 ~ Loc, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Loc, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Loc, data=GenoCoeffs)
bartlett.test(ED50 ~ Site, data=GenoCoeffs)
bartlett.test(ED50 ~ Reef, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Loc, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Loc, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Loc"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]


AcroReef_labels <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("", "*"))
AcroReef_labels2 <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("a", "b"))
facet_labels <- c("Al Lith" = "Southern Red Sea", "Thuwal" = "Central Red Sea")

AcroPlotReef <- ggplot(GenoCoeffs, aes(y = ED50, x = Reef))+
  ## geom_point(aes(colour = Site)) +
  geom_boxplot(aes(colour = Reef, fill = Reef), coef = Inf, alpha = 0.2) +
  geom_jitter(aes(colour = Reef), position = position_jitterdodge(jitter.width = 0.4), alpha = 0.2, size = 2) +
  theme_bw()+
  facet_wrap(~ Site, ncol = 2, labeller = labeller(Site = facet_labels)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19))+
  ylab(expression(ED50 ~(degree*C))) +  
  xlab('') +
  ylim(33,40) +
  geom_text(x = 1.5, y = 39.8, aes(label = label), data = AcroReef_labels, size = 12) +
  geom_text(x = 2, y = 33.1, aes(label = label), data = AcroReef_labels2, size = 7) +
  scale_colour_manual(values = c("#00B4AB","#4575b4")) +
  scale_fill_manual(values = c("#00B4AB","#4575b4")) +
  labs(tag = "A") +
  theme(legend.position = "none")
AcroPlotReef




#### Pocillopora Reef Location Comparison ----
#### Merging Coeffs of all species in offshore
Coeffs<-c(DRCpamPOCInearKAUST$coefficients[3], DRCpamPOCIoffKAUST$coefficients[3], DRCpamPOCInearAl$coefficients[3], DRCpamPOCIoffAl$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamPOCIgenonearKAUST$coefficients[31:45], DRCpamPOCIgenooffKAUST$coefficients[31:45], DRCpamPOCIgenonearAl$coefficients[31:45], DRCpamPOCIgenooffAl$coefficients[31:45]))
GenoCoeffs$Loc=c(rep("nearKAUST",15), rep("offKAUST",15), rep("nearAl",15), rep("offAl",15))
GenoCoeffs$Reef=c(rep("nearshore",15), rep("offshore",15), rep("nearshore",15), rep("offshore",15))
GenoCoeffs$Site=c(rep("Thuwal",30), rep("Al Lith",30))
aggregate(ED50 ~ Loc, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Loc, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Loc, data=GenoCoeffs)
bartlett.test(ED50 ~ Site, data=GenoCoeffs)
bartlett.test(ED50 ~ Reef, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Loc, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Loc, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Loc"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]



PociReef_labels <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("", "*"))
facet_labels <- c("Al Lith" = "Southern Red Sea", "Thuwal" = "Central Red Sea")

PociPlotReef <- ggplot(GenoCoeffs, aes(y = ED50, x = Reef))+
  ## geom_point(aes(colour = Site)) +
  geom_boxplot(aes(colour = Reef, fill = Reef), coef = Inf, alpha = 0.2) +
  geom_jitter(aes(colour = Reef), position = position_jitterdodge(jitter.width = 0.4), alpha = 0.2, size = 2) +
  theme_bw()+
  facet_wrap(~ Site, ncol = 2, labeller = labeller(Site = facet_labels)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19))+
  ylab(expression(ED50 ~(degree*C))) +  
  xlab('') +
  ylim(33,40) +
  geom_text(x = 1.5, y = 39.8, aes(label = label), data = PociReef_labels, size = 12) +
  scale_colour_manual(values = c("#00B4AB","#4575b4")) +
  scale_fill_manual(values = c("#00B4AB","#4575b4")) +
  labs(tag = "B") +
  theme(legend.position = "none")
PociPlotReef







#### Turf Algae Reef Location Comparison ----
#### Merging Coeffs of all species in offshore
# Coeffs<-c(DRCpamTURFnearKAUST$coefficients[3], DRCpamTURFoffKAUST$coefficients[3], DRCpamTURFnearAl$coefficients[3], DRCpamTURFoffAl$coefficients[3])
# GenoCoeffs<-data.frame("ED50"=c(DRCpamTURFgenonearKAUST$coefficients[31:45], DRCpamTURFgenooffKAUST$coefficients[31:45], DRCpamTURFgenonearAl$coefficients[31:45], DRCpamTURFgenooffAl$coefficients[31:45]))
# GenoCoeffs$Loc=c(rep("nearKAUST",15), rep("offKAUST",15), rep("nearAl",15), rep("offAl",15))
# GenoCoeffs$Reef=c(rep("nearshore",15), rep("offshore",15), rep("nearshore",15), rep("offshore",15))
# GenoCoeffs$Site=c(rep("Thuwal",30), rep("Al Lith",30))
# aggregate(ED50 ~ Loc, data=GenoCoeffs, mean)
# 
# Coeffs
# 
# 
# #### ED50 ANoVA
# aggregate(ED50 ~ Loc, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
# bartlett.test(ED50 ~ Loc, data=GenoCoeffs)
# bartlett.test(ED50 ~ Site, data=GenoCoeffs)
# bartlett.test(ED50 ~ Reef, data=GenoCoeffs)
# 
# ED50.aov<-aov(ED50 ~ Loc, GenoCoeffs)
# summary(ED50.aov)
# TukeyHSD(ED50.aov)
# 
# aggregate(ED50 ~ Loc, data=GenoCoeffs, summary)
# GenoCoeffs
# SumaryStats<-data.frame("Loc"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, length)))
# SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
# 
# TurfReef_labels <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("", ""))
# TurfReef_labels2 <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("a", "b"))
# facet_labels <- c("Al Lith" = "Southern Red Sea", "Thuwal" = "Central Red Sea")
# 
# TurfPlotReef <- ggplot(GenoCoeffs, aes(y = ED50, x = Reef))+
#   ## geom_point(aes(colour = Site)) +
#   geom_boxplot(aes(colour = Reef, fill = Reef), coef = Inf, alpha = 0.2) +
#   geom_jitter(aes(colour = Reef), position = position_jitterdodge(jitter.width = 0.4), alpha = 0.2, size = 2) +
#   theme_bw()+
#   facet_wrap(~ Site, ncol = 2, labeller = labeller(Site = facet_labels)) +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         text = element_text(size = 18), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19))+
#   ylab(expression(ED50 ~(degree*C))) +  
#   xlab('') +
#   ylim(33,40) +
#   geom_text(x = 1.5, y = 39.8, aes(label = label), data = TurfReef_labels, size = 12) +
#   geom_text(x = 2, y = 33, aes(label = label), data = TurfReef_labels2, size = 8) +
#   scale_colour_manual(values = c("#00B4AB","#4575b4")) +
#   scale_fill_manual(values = c("#00B4AB","#4575b4")) +
#   labs(tag = "C") +
#   theme(legend.position = "none")
# TurfPlotReef
# 



#### Merging Coeffs of all species in offshore
# Keep only valid data for "Al Lith Offshore"
# Manually create ED50 values
ED50_values <- c(rep(NA, 45), DRCpamTURFgenooffAl$coefficients[31:45])

# Recreate GenoCoeffs dataframe with the manually created ED50 values
GenoCoeffs <- data.frame(
  ED50 = ED50_values,
  Loc = c(rep("nearKAUST", 15), rep("offKAUST", 15), rep("nearAl", 15), rep("offAl", 15)),
  Reef = c(rep("nearshore", 15), rep("offshore", 15), rep("nearshore", 15), rep("offshore", 15)),
  Site = c(rep("Thuwal", 30), rep("Al Lith", 30))
)


# Plot as before
TurfPlotReef <- ggplot(GenoCoeffs, aes(y = ED50, x = Reef)) +
  geom_boxplot(aes(colour = Reef, fill = Reef), coef = Inf, alpha = 0.2, na.rm = TRUE) +
  geom_jitter(aes(colour = Reef), position = position_jitterdodge(jitter.width = 0.4), alpha = 0.2, size = 2, na.rm = TRUE) +
  theme_bw() +
  facet_wrap(~ Site, ncol = 2, labeller = labeller(Site = facet_labels)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 18),
    axis.text.x = element_text(size = 19),
    axis.text.y = element_text(size = 19)
  ) +
  ylab(expression(ED50 ~(degree*C))) +
  xlab('') +
  ylim(33, 40) +
  scale_colour_manual(values = c("#4575b4", "#4575b4")) +
  scale_fill_manual(values = c("#4575b4", "#4575b4")) +
  labs(tag = "C") +
  theme(legend.position = "none")

# Display the plot
TurfPlotReef




#### Xenia Reef Location Comparison ----
#### Merging Coeffs of all species in offshore
Coeffs<-c(DRCpamXENInearKAUST$coefficients[3], DRCpamXENIoffKAUST$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamXENIgenonearKAUST$coefficients[31:45], DRCpamXENIgenooffKAUST$coefficients[31:45]))
GenoCoeffs$Loc=c(rep("nearKAUST",15), rep("offKAUST",15))
GenoCoeffs$Reef=c(rep("nearshore",15), rep("offshore",15))
GenoCoeffs$Site=c(rep("Thuwal",30))
aggregate(ED50 ~ Loc, data=GenoCoeffs, mean)

Coeffs


#### ED50 ANoVA
aggregate(ED50 ~ Loc, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Loc, data=GenoCoeffs)
bartlett.test(ED50 ~ Reef, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Loc, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Loc, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Loc"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Loc, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]


XeniReef_labels <- data.frame(Site = c("Al Lith", "Thuwal"), label = c("", ""))

# Define custom labels
facet_labels <- c("Al Lith" = "Southern Red Sea", "Thuwal" = "Central Red Sea")

XeniPlotReef <- ggplot(GenoCoeffs, aes(y = ED50, x = Reef)) +
  geom_boxplot(aes(colour = Reef, fill = Reef), coef = Inf, alpha = 0.2) +
  geom_jitter(aes(colour = Reef), position = position_jitterdodge(jitter.width = 0.4), alpha = 0.2, size = 2) +
  theme_bw() +
  facet_wrap(~ Site, ncol = 2, labeller = labeller(Site = facet_labels)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) +
  ylab(expression(ED50 ~(degree * C))) +  
  xlab('') +
  ylim(33, 40) +
  geom_text(x = 1.5, y = 39.8, aes(label = label), data = XeniReef_labels, size = 12) +
  scale_colour_manual(values = c("#00B4AB", "#4575b4")) +
  scale_fill_manual(values = c("#00B4AB", "#4575b4")) +
  labs(tag = "D") +
  theme(legend.position = "none")

XeniPlotReef





# Test a run to get all ED50 at the same time----

# Extract coefficients
coefficients_df <- data.frame(
  Acro_Near_Al = DRCpamACROgenonearAl$coefficients[31:45],
  Acro_Near_KAUST = DRCpamACROgenonearKAUST$coefficients[31:45],
  Acro_Off_Al = DRCpamACROgenooffAl$coefficients[31:45],
  Acro_Off_KAUST = DRCpamACROgenooffKAUST$coefficients[31:45],
  Poci_Near_Al = DRCpamPOCIgenonearAl$coefficients[31:45],
  Poci_Near_KAUST = DRCpamPOCIgenonearKAUST$coefficients[31:45],
  Poci_Off_Al = DRCpamPOCIgenooffAl$coefficients[31:45],
  Poci_Off_KAUST = DRCpamPOCIgenooffKAUST$coefficients[31:45],
  Turf_Near_Al = DRCpamTURFgenonearAl$coefficients[31:45],
  Turf_Near_KAUST = DRCpamTURFgenonearKAUST$coefficients[31:45],
  Turf_Off_Al = DRCpamTURFgenooffAl$coefficients[31:45],
  Turf_Off_KAUST = DRCpamTURFgenooffKAUST$coefficients[31:45],
  Xeni_Near_KAUST = DRCpamXENIgenonearKAUST$coefficients[31:45],
  Xeni_Off_KAUST = DRCpamXENIgenooffKAUST$coefficients[31:45]
)

# Add a column for the model names
coefficients_df$model <- rownames(coefficients_df)
# Reset rownames
rownames(coefficients_df) <- NULL
# Remove the "model" column
coefficients_df <- subset(coefficients_df, select = -c(model))
# Reshape the dataframe into long format
coefficients_df <- gather(coefficients_df, key = "Treatment", value = "ED50")
# Separate treatment column
coefficients_df <- separate(coefficients_df, Treatment, into = c("Species", "Distance", "Location"), sep = "_")
# Print the resulting dataframe
print(coefficients_df)


# run a two way ANOVA on them

# check for outliers
coefficients_df %>%
  group_by(Species, Distance, Location) %>%
  identify_outliers(ED50) # there are a few outliers but they are correct in the dataset, so we leave them in

# normality assumption
# Build the linear model
model  <- lm(ED50 ~ Species*Distance*Location,
             data = coefficients_df)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# check by group
coefficients_df %>%
  group_by(Species, Distance, Location) %>%
  shapiro_test(ED50)

# no transformations yielded a normal distribution (log, sqrt, reciprocal, logit, boxcox, etc)
# so we continue with a kruskal wallis test as a non parametric alternative

# Function to perform Kruskal-Wallis test for Distance within each Location and for Location within each Distance
kruskal_test <- function(df) {
  results <- list()
  
  # Test for Distance within each Location
  for (loc in unique(df$Location)) {
    df_loc <- df %>% filter(Location == loc)
    if (length(unique(df_loc$Distance)) > 1) {
      kruskal_distance <- kruskal.test(ED50 ~ Distance, data = df_loc)
      results <- rbind(results, data.frame(Species = unique(df$Species), Comparison = paste("Distance in", loc), p.value = kruskal_distance$p.value))
    }
  }
  
  # Test for Location within each Distance
  for (dist in unique(df$Distance)) {
    df_dist <- df %>% filter(Distance == dist)
    if (length(unique(df_dist$Location)) > 1) {
      kruskal_location <- kruskal.test(ED50 ~ Location, data = df_dist)
      results <- rbind(results, data.frame(Species = unique(df$Species), Comparison = paste("Location in", dist), p.value = kruskal_location$p.value))
    }
  }
  
  return(as.data.frame(results))
}

# Group by Species and perform Kruskal-Wallis test for Distance within each Location and for Location within each Distance
results <- coefficients_df %>%
  group_by(Species) %>%
  do(kruskal_test(.)) %>%
  ungroup()

# Adjust p-values using Benjamini-Hochberg method
results <- results %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH")) #or use BH (Benjamini-Hochberg)

print(results)
