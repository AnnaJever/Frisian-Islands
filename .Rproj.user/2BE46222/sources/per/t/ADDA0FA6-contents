# Frisian Islands MS
# Walentowitz et al.

# last changes on 22.09.2022





# load packages
library(lavaan)
library(semPlot) 
library(tidyverse)
library(gridExtra)
library(fastDummies)
library(performance)
library(MVN)
library(rms)
require(ggeffects)
require(jtools)
library(lmtest)
library(sandwich)




### Data ###

islands <- read.csv2("Frisian_Islands_data.csv")




### Analysis ###

# correlation between spec number and isolation
cor(islands[,c("totspec", "totnat", "totinv",
               "isolation_100", "isolation_80",
               "isolation_60", "isolation_40",
               "isolation_20", "isolation_10")])

# isolation_100 has the highest correlation coefficient and will be used in subsequent analysis

# check for differences in species pending country and island type

# country
boxplot(totspec ~ country, data = islands, col = c("deepskyblue4", "darkred", "gold"))
boxplot(totnat ~ country, data = islands, col = c("deepskyblue4", "darkred", "gold"))
boxplot(totinv ~ country, data = islands, col = c("deepskyblue4", "darkred", "gold"))
re_aov <- aov(totspec ~ country, data = islands)
summary(re_aov) # not significantly different

#type
boxplot(totspec ~ type, data = islands, col = c("lightgrey", "forestgreen", "deepskyblue3"))
boxplot(totnat ~ type, data = islands, col = c("lightgrey","forestgreen", "deepskyblue3"))
boxplot(totinv ~ type, data = islands, col = c("lightgrey","forestgreen", "deepskyblue3"))
summary(aov(totspec ~ type, data = islands))
TukeyHSD(aov(totspec ~ type, data = islands), conf.level=.95)
# "Halligen" have significantly less, native and alien species, but they are commonly smaller




### Exploratory 2-dimensional analysis ###

# subset for analysis of inhabited islands
islands_inhabt <- islands[which(islands$inhabited == "yes"),]

# area
plot_a <- ggplot(islands) +
  geom_point(aes(log(area),log(totnat)), col = "steelblue3") +
  geom_point(aes(log(area),log(totinv)), col = "salmon") +
  geom_smooth(aes(log(area),log(totnat)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(log(area),log(totinv)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Area, log", y = "Species richness, log")

glm_totnat_area <- glm(log(totnat) ~ log(area), data = islands)
summary(glm_totnat_area)
lrm(islands$totnat ~ islands$area)

glm_totinv_area <- glm(log(totinv) ~ log(area), data = islands)
summary(glm_totinv_area)
lrm(log(islands$totinv) ~ log(islands$area))

# isolation
plot_b <- ggplot(islands) +
  geom_point(aes(isolation_100,totnat), col = "steelblue3") +
  geom_point(aes(isolation_100,totinv), col = "salmon") +
  geom_smooth(aes(isolation_100,totnat),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(isolation_100,totinv),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x =  "Isolation", y = "Species richness")

glm_totnat_isolation<- glm(totnat ~ isolation_100, data = islands)
summary(glm_totnat_isolation)
lrm(islands$totnat ~ islands$isolation_100)

glm_totinv_isolation <- glm(totinv ~ isolation_100, data = islands)
summary(glm_totinv_isolation)
lrm(islands$totinv ~ islands$isolation_100)

# inhabitants
plot_c <- ggplot(islands[which(islands$inhabited == "yes"),]) +
  geom_point(aes(log(inhabitants),log(totnat)), col = "steelblue3") +
  geom_point(aes(log(inhabitants),log(totinv)), col = "salmon") +
  geom_smooth(aes(log(inhabitants),log(totnat)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(log(inhabitants),log(totinv)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Inhabitants, log", y = "Species richness, log")

glm_totnat_density_inh <- glm(log(totnat) ~ log(inhabitants), data = islands[which(islands$inhabited == "yes"),])
summary(glm_totnat_density_inh)
lrm(log(islands_inhabt$totnat) ~ log(islands_inhabt$inhabitants))

glm_totinv_density_inh <- glm(log(totinv) ~ log(inhabitants), data = islands[which(islands$inhabited == "yes"),])
summary(glm_totinv_density_inh)
lrm(log(islands_inhabt$totinv) ~ log(islands_inhabt$inhabitants))

# annual tourist visits
plot_d <- ggplot(islands[which(islands$inhabited == "yes"),]) +
  geom_point(aes(log(tourism),log(totnat)), col = "steelblue3") +
  geom_point(aes(log(tourism),log(totinv)), col = "salmon") +
  geom_smooth(aes(log(tourism),log(totnat)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(log(tourism),log(totinv)),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Tourists per year", y = "Species richness")

glm_totnat_tou <- glm(log(totnat) ~ log(tourism), data = islands[which(islands$inhabited == "yes"),])
summary(glm_totnat_tou)
lrm(log(islands_inhabt$totnat) ~ log(islands_inhabt$tourism))

glm_totinv_tou <- glm(log(totinv) ~ log(tourism), data = islands[which(islands$inhabited == "yes"),])
summary(glm_totinv_tou)
lrm(log(islands_inhabt$totinv) ~ log(islands_inhabt$tourism))

# land anthro
plot_e <- ggplot(islands[which(islands$inhabited == "yes"), ]) +
  geom_point(aes(land_anthro,totnat), col = "steelblue3") +
  geom_point(aes(land_anthro,totinv), col = "salmon") +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Human land use", y = "Species richness")

glm_totnat_land_anthro <- glm(totnat ~ land_anthro, data = islands[which(islands$inhabited == "yes"),])
summary(glm_totnat_land_anthro)

glm_totinv_land_anthro <- glm(totinv ~ land_anthro, data = islands[which(islands$inhabited == "yes"),])
summary(glm_totinv_land_anthro)

# habitat heterogeneity
plot_f <- ggplot(islands) +
  geom_point(aes(habitat_heter,totnat), col = "steelblue3") +
  geom_point(aes(habitat_heter, totinv), col = "salmon") +
  geom_smooth(aes(habitat_heter,totnat),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes( habitat_heter,totinv),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Habitat heterogeneity", y = "Species richness")

glm_totnat_habitat_heter <- glm(totnat ~ habitat_heter, data = islands)
summary(glm_totnat_habitat_heter)
lrm(islands$totnat ~ islands$habitat_heter)

glm_totinv_habitat_heter <- glm(totinv ~ habitat_heter, data = islands)
summary(glm_totinv_habitat_heter)
lrm(islands$totinv ~ islands$habitat_heter)



# GLMs with standardized spec numbers (to max = 1)

# area
plot_a2 <- ggplot(islands) +
  geom_point(aes(area,totnat_stdz), col = "steelblue3") +
  geom_point(aes(area,totinv_stdz), col = "salmon") +
  geom_smooth(aes(area,totnat_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(area,totinv_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Area", y = "Species richness")+
  ylim(c(0,1))


# isolation
plot_b2 <- ggplot(islands) +
  geom_point(aes(isolation_100,totnat_stdz), col = "steelblue3") +
  geom_point(aes(isolation_100,totinv_stdz), col = "salmon") +
  geom_smooth(aes(isolation_100,totnat_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(isolation_100,totinv_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x =  "Isolation", y = "Species richness")+
  ylim(c(0,1))


# population density
plot_c2 <- ggplot(islands[which(islands$inhabited == "yes"),]) +
  geom_point(aes(inhabitants,totnat_stdz), col = "steelblue3") +
  geom_point(aes(inhabitants,totinv_stdz), col = "salmon") +
  geom_smooth(aes(inhabitants,totnat_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(inhabitants,totinv_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Inhabitants", y = "Species richness")+
  ylim(c(0,1))


# annual tourist visits
plot_d2 <- ggplot(islands[which(islands$inhabited == "yes"),]) +
  geom_point(aes(tourism,totnat_stdz), col = "steelblue3") +
  geom_point(aes(tourism,totinv_stdz), col = "salmon") +
  geom_smooth(aes(tourism,totnat_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(tourism,totinv_stdz),method = "glm", se = T, 
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Tourists per year", y = "Species richness")+
  ylim(c(0,1))


# land anthro
plot_e2 <- ggplot(islands[which(islands$inhabited == "yes"),]) +
  geom_point(aes(land_anthro,totnat_stdz), col = "steelblue3") +
  geom_point(aes(land_anthro,totinv_stdz), col = "salmon") +
  geom_smooth(aes(land_anthro,totnat_stdz),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes(land_anthro,totinv_stdz),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Human land use", y = "Species richness")+
  ylim(c(0,1))


# habitat heterogeneity
plot_f2 <- ggplot(islands) +
  geom_point(aes(habitat_heter,totnat_stdz), col = "steelblue3") +
  geom_point(aes( habitat_heter, totinv_stdz), col = "salmon") +
  geom_smooth(aes( habitat_heter,totnat_stdz),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "steelblue3")+
  geom_smooth(aes( habitat_heter,totinv_stdz),method = "glm", se = T,
              method.args = list(family = "poisson"),
              colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Habitat heterogeneity", y = "Species richness")+
  ylim(c(0,1))




# Multipanel plot
x11()
grid.arrange(plot_a, plot_b,
             plot_c, plot_d,
             plot_e, plot_f,
             ncol=3)
dev.off()

x11()
grid.arrange(plot_a2, plot_b2,
             plot_c2, plot_d2,
             plot_e2, plot_f2,
             ncol=3)
dev.off()
### end 2 dimensional exploratory analysis ###





### SEM ###

#correlation matrix of all variables used in SEM

data_all<- islands[,c(22, 23, 29, 32, 13, 19, 33)]
cor(data_all)

data_nat <- islands[,c(29, 32, 13, 19, 33, 22)]

data_inv <- islands[,c(29, 32, 13, 19, 33, 23)]

result1 <- mvn(data = data_nat, mvnTest = "mardia",
               univariatePlot = "qqplot",
               multivariatePlot = "qq",
               univariateTest = "SW")
mvn(data = data_nat, mvnTest = "mardia",
    multivariatePlot = "qq")

result1$multivariateNormality[,]
result1$univariateNormality


result2 <- mvn(data = data_inv, mvnTest = "mardia",
               multivariatePlot = "qq",
               univariatePlot = "qqplot",
               univariateTest = "SW")
mvn(data = data_inv, mvnTest = "mardia",
    multivariatePlot = "qq")

result2$multivariateNormality[,]
result2$univariateNormality





# SEM native species

covmat.nat <- cor(data_nat)

# saturated model natives
island.model.nat <- "
                totnat ~ area_log
                totnat ~ tourism_sqrt
                totnat ~ inhabitants_sqrt
                totnat ~ isolation_100
                totnat ~ habitat_heter

                inhabitants_sqrt ~ area_log
                
                habitat_heter ~ area_log
                
                tourism_sqrt ~ area_log
                "

# confirmatory factor analysis
island.fit.nat <- sem(model = island.model.nat,
                      sample.cov = covmat.nat,
                      sample.nobs = 31)


summary(island.fit.nat, standardized = T, fit = T, rsq = T)
semPaths(island.fit.nat, "std", "hide")
coef(island.fit.nat)


# Calculate fit indices
fitMeasures(island.fit.nat, c("chisq", "df", "pvalue", "rmsea", "srmr", "nnfi", "cfi"))




# SEM alien

covmat.inv <- cor(data_inv)

# saturated model alien
island.model.inv <- "
                totinv ~ area_log
                totinv ~ tourism_sqrt
                totinv ~ inhabitants_sqrt
                totinv ~ isolation_100
                totinv ~ habitat_heter

                inhabitants_sqrt ~ area_log
                
                habitat_heter ~ area_log
                
                tourism_sqrt ~ area_log
                "

island.fit.inv <- sem(model = island.model.inv,
                      sample.cov = covmat.inv,
                      sample.nobs = 31)


summary(island.fit.inv, standardized = T, fit = T, rsq = T)
semPaths(island.fit.inv, "std", "hide")
coef(island.fit.inv)


# Calculate fit indices
fitMeasures(island.fit.inv, c("chisq", "df", "pvalue", "rmsea", "srmr", "nnfi", "cfi"))

### end SEM ###





### additional SEM ###

# use reduced set of islands (only inhabited)
# without tourism as this variable was not important when analyzing univariate relationships

data_nat2 <- islands[which(islands$inhabitants > 1),]
data_nat2 <- data_nat2[,c( 29, 32, 13, 19, 33, 22)]
data_inv2 <- islands[which(islands$tourism > 0),]
data_inv2 <- data_inv2[,c(29, 32, 13, 19, 33, 23)]

result1a <- mvn(data = data_nat2, mvnTest = "mardia",
                univariatePlot = "qqplot",
                multivariatePlot = "qq",
                univariateTest = "SW")
mvn(data = data_nat2, mvnTest = "mardia",
    multivariatePlot = "qq")

result1a$multivariateNormality[,]
result1a$univariateNormality

result2a <- mvn(data = data_inv2, mvnTest = "mardia",
                multivariatePlot = "qq",
                univariatePlot = "qqplot",
                univariateTest = "SW")
mvn(data = data_inv2, mvnTest = "mardia",
    multivariatePlot = "qq")

result2a$multivariateNormality[,]
result2a$univariateNormality

covmat.nat2 <- cor(data_nat2)

# saturated model natives
island.model.nat2 <- "
                totnat ~ area_log
                totnat ~ tourism_sqrt
                totnat ~ inhabitants_sqrt
                totnat ~ isolation_100
                totnat ~ habitat_heter

                inhabitants_sqrt ~ area_log
                
                habitat_heter ~ area_log
                
                tourism_sqrt ~ area_log
                "

# confirmatory factor analysis
island.fit.nat2 <- sem(model = island.model.nat2,
                       sample.cov = covmat.nat2,
                       sample.nobs = 31)


summary(island.fit.nat2, standardized = T, fit = T, rsq = T)
semPaths(island.fit.nat2, "std", "hide")
coef(island.fit.nat2)


# Calculate fit indices
fitMeasures(island.fit.nat2, c("chisq", "df", "pvalue", "rmsea", "srmr", "nnfi", "cfi"))

# SEM alien

covmat.inv2 <- cor(data_inv2)

# saturated model alien
island.model.inv2 <- "
                totinv ~ area_log
                totinv ~ tourism_sqrt
                totinv ~ inhabitants_sqrt
                totinv ~ isolation_100
                totinv ~ habitat_heter

                inhabitants_sqrt ~ area_log
                
                habitat_heter ~ area_log
                
                tourism_sqrt ~ area_log
                "

island.fit.inv2 <- sem(model = island.model.inv2,
                       sample.cov = covmat.inv2,
                       sample.nobs = 31)


summary(island.fit.inv2, standardized = T, fit = T, rsq = T)
semPaths(island.fit.inv2, "std", "hide")
coef(island.fit.inv2)


# Calculate fit indices
fitMeasures(island.fit.inv2, c("chisq", "df", "pvalue", "rmsea", "srmr", "nnfi", "cfi"))

### end additional SEM ###