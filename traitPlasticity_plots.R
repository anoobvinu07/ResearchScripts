## notes
# gw - growth
# bb - budbreak
# bs - budset
# B - lmer family blups (means)
# P - plasticity value for the trait

# packages
require(ggplot2)
require(lmerTest)
require(dplyr)
require(tidyr)
require(MCMCglmm)
require(plyr)

# setwd("Z:/Anoob/MCMCglmm")

exome_meta <- read.table("./Exome/RS_Exome_metadata.txt", sep="\t", header = T)

#---------------------------------------Height growth-----------------------------------------------#

#########################################################
## Height growth of 2019 - Means vs Plasticity  
#########################################################

# read the output
gw_19_B_gw_19_P_corr <- readRDS("./plasticity_outputs/gw_19_B_gw_19_P_corr")
plot(gw_19_B_gw_19_P_corr)

effectiveSize(gw_19_B_gw_19_P_corr$VCV)
heidel.diag(gw_19_B_gw_19_P_corr$VCV)

## setting up dataframe for correlation based on the model blups
gw_19_B_gw_19_P_coefs <- data_frame(Trait = attr(colMeans(gw_19_B_gw_19_P_corr$Sol), "names"),
                                    Value = colMeans(gw_19_B_gw_19_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitPlasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_19_B_gw_19_P_coefs <-  merge(x=gw_19_B_gw_19_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_19_B_gw_19_P_fit_slope <- gw_19_B_gw_19_P_corr$VCV[,"traitHeight_Growth:traitPlasticity.Population"]/
  gw_19_B_gw_19_P_corr$VCV[,"traitPlasticity:traitPlasticity.Population"]

# rename regions
gw_19_B_gw_19_P_coefs$Region <- revalue(gw_19_B_gw_19_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_19_B_gw_19_P_coefs$Region <- factor(gw_19_B_gw_19_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_19_B_gw_19_P_coefs, aes(x = traitPlasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_19_B_gw_19_P_fit_slope)) + 
  labs(x = "Plasticity of Height Growth in 2019 (BLUP)", y = "Height Growth in 2019 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 
  



#########################################################
## Height growth of 2020 - Means vs Plasticity  
#########################################################
gw_20_B_gw_20_P_corr <- readRDS("./plasticity_outputs/gw_20_B_gw_20_P_corr")
plot(gw_20_B_gw_20_P_corr)
effectiveSize(gw_20_B_gw_20_P_corr$VCV)
heidel.diag(gw_20_B_gw_20_P_corr$VCV)

## setting up dataframe for correlation based on the model blups
gw_20_B_gw_20_P_coefs <- data_frame(Trait = attr(colMeans(gw_20_B_gw_20_P_corr$Sol), "names"),
                                    Value = colMeans(gw_20_B_gw_20_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitPlasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_20_B_gw_20_P_coefs <-  merge(x=gw_20_B_gw_20_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_20_B_gw_20_P_fit_slope <- gw_20_B_gw_20_P_corr$VCV[,"traitHeight_Growth:traitPlasticity.Population"]/
  gw_20_B_gw_20_P_corr$VCV[,"traitPlasticity:traitPlasticity.Population"]


# rename regions
gw_20_B_gw_20_P_coefs$Region <- revalue(gw_20_B_gw_20_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_20_B_gw_20_P_coefs$Region <- factor(gw_20_B_gw_20_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_20_B_gw_20_P_coefs, aes(x = traitPlasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_20_B_gw_20_P_fit_slope)) + 
  labs(x = "Plasticity of Height Growth in 2020 (BLUP)", y = "Height Growth in 2020 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 





#----------------------------------Height growth vs budbreak-------------------------------------#

##########################################################
## Height growth of 2019 means vs budbreak 2020 Plasticity  
##########################################################

gw_19_B_bb_20_P_corr <- readRDS("./plasticity_outputs/gw_19_B_bb_20_P_corr")

plot(gw_19_B_bb_20_P_corr)

effectiveSize(gw_19_B_bb_20_P_corr$VCV)
heidel.diag(gw_19_B_bb_20_P_corr$VCV)



## setting up dataframe for correlation based on the model blups
gw_19_B_bb_20_P_coefs <- data_frame(Trait = attr(colMeans(gw_19_B_bb_20_P_corr$Sol), "names"),
                                    Value = colMeans(gw_19_B_bb_20_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitBudBreak2020_Plasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_19_B_bb_20_P_coefs <-  merge(x=gw_19_B_bb_20_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_19_B_bb_20_P_fit_slope <- gw_19_B_bb_20_P_corr$VCV[,"traitHeight_Growth:traitBudBreak2020_Plasticity.Population"]/
  gw_19_B_bb_20_P_corr$VCV[,"traitBudBreak2020_Plasticity:traitBudBreak2020_Plasticity.Population"]

# rename regions
gw_19_B_bb_20_P_coefs$Region <- revalue(gw_19_B_bb_20_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_19_B_bb_20_P_coefs$Region <- factor(gw_19_B_bb_20_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_19_B_bb_20_P_coefs, aes(x = traitBudBreak2020_Plasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_19_B_bb_20_P_fit_slope)) + 
  labs(x = "Plasticity of Bud break in 2020 (BLUP)", y = "Height Growth in 2019 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 

cor.test(gw_19_B_bb_20_P_coefs$traitBudBreak2020_Plasticity,gw_19_B_bb_20_P_coefs$traitHeight_Growth)

##########################################################
## Height growth of 2020 means vs budbreak 2020 Plasticity
##########################################################

gw_20_B_bb_20_P_corr <- readRDS("./plasticity_outputs/gw_20_B_bb_20_P_corr")
plot(gw_20_B_bb_20_P_corr)
effectiveSize(gw_20_B_bb_20_P_corr$VCV)
heidel.diag(gw_20_B_bb_20_P_corr$VCV)

## setting up dataframe for correlation based on the model blups
gw_20_B_bb_20_P_coefs <- data_frame(Trait = attr(colMeans(gw_20_B_bb_20_P_corr$Sol), "names"),
                                    Value = colMeans(gw_20_B_bb_20_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitBudBreak2020_Plasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_20_B_bb_20_P_coefs <-  merge(x=gw_20_B_bb_20_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_20_B_bb_20_P_fit_slope <- gw_20_B_bb_20_P_corr$VCV[,"traitHeight_Growth:traitBudBreak2020_Plasticity.Population"]/
  gw_20_B_bb_20_P_corr$VCV[,"traitBudBreak2020_Plasticity:traitBudBreak2020_Plasticity.Population"]

# rename regions
gw_20_B_bb_20_P_coefs$Region <- revalue(gw_20_B_bb_20_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_20_B_bb_20_P_coefs$Region <- factor(gw_20_B_bb_20_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_20_B_bb_20_P_coefs, aes(x = traitBudBreak2020_Plasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_20_B_bb_20_P_fit_slope)) + 
  labs(x = "Plasticity of Bud break in 2020 Plasticity (BLUP)", y = "Height Growth in 2020 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 

cor.test(gw_20_B_bb_20_P_coefs$traitBudBreak2020_Plasticity,gw_20_B_bb_20_P_coefs$traitHeight_Growth)
#------------------------------------Height growth vs budset---------------------------------------#

##########################################################
## Height growth of 2019 means vs budset 2019 Plasticity  
##########################################################

gw_19_B_bs_19_P_corr <- readRDS("./plasticity_outputs/gw_19_B_bs_19_P_corr")
plot(gw_19_B_bs_19_P_corr)
effectiveSize(gw_19_B_bs_19_P_corr$VCV)
heidel.diag(gw_19_B_bs_19_P_corr$VCV)

## setting up dataframe for correlation based on the model blups
gw_19_B_bs_19_P_coefs <- data_frame(Trait = attr(colMeans(gw_19_B_bs_19_P_corr$Sol), "names"),
                                    Value = colMeans(gw_19_B_bs_19_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitBudSet2019_Plasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_19_B_bs_19_P_coefs <-  merge(x=gw_19_B_bs_19_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_19_B_bs_19_P_fit_slope <- gw_19_B_bs_19_P_corr$VCV[,"traitHeight_Growth:traitBudSet2019_Plasticity.Population"]/
  gw_19_B_bs_19_P_corr$VCV[,"traitBudSet2019_Plasticity:traitBudSet2019_Plasticity.Population"]

# rename regions
gw_19_B_bs_19_P_coefs$Region <- revalue(gw_19_B_bs_19_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_19_B_bs_19_P_coefs$Region <- factor(gw_19_B_bs_19_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_19_B_bs_19_P_coefs, aes(x = traitBudSet2019_Plasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_19_B_bs_19_P_fit_slope)) + 
  labs(x = "Bud set 2019 Plasticity (BLUP)", y = "Height Growth in 2019 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 
cor.test(gw_19_B_bs_19_P_coefs$traitBudSet2019_Plasticity,gw_19_B_bs_19_P_coefs$traitHeight_Growth)

##########################################################
## Height growth of 2020 means vs budset 2020 Plasticity
##########################################################

gw_20_B_bs_20_P_corr <- readRDS("./plasticity_outputs/gw_20_B_bs_20_P_corr")

plot(gw_20_B_bs_20_P_corr)
effectiveSize(gw_20_B_bs_20_P_corr$VCV)
heidel.diag(gw_20_B_bs_20_P_corr$VCV)


## setting up dataframe for correlation based on the model blups
gw_20_B_bs_20_P_coefs <- data_frame(Trait = attr(colMeans(gw_20_B_bs_20_P_corr$Sol), "names"),
                                    Value = colMeans(gw_20_B_bs_20_P_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitHeight_Growth", "traitBudSet2020_Plasticity")) %>%
  select(-Type) %>%
  spread(Trait, Value)

gw_20_B_bs_20_P_coefs <-  merge(x=gw_20_B_bs_20_P_coefs,
                                y=exome_meta[,c("Pop","Region")],
                                all.x=T,
                                by.x="Population",
                                by.y="Pop")


# Find the regression line -
# the covariance of growth, plasticity divided by
# the variance in plasticity
gw_20_B_bs_20_P_fit_slope <- gw_20_B_bs_20_P_corr$VCV[,"traitHeight_Growth:traitBudSet2020_Plasticity.Population"]/
  gw_20_B_bs_20_P_corr$VCV[,"traitBudSet2020_Plasticity:traitBudSet2020_Plasticity.Population"]

# rename regions
gw_20_B_bs_20_P_coefs$Region <- revalue(gw_20_B_bs_20_P_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_20_B_bs_20_P_coefs$Region <- factor(gw_20_B_bs_20_P_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_20_B_bs_20_P_coefs, aes(x = traitBudSet2020_Plasticity, y = traitHeight_Growth, group = Population, color = Region)) +
  geom_point(alpha = 0.7) + 
  geom_abline(intercept = 0, slope = mean(gw_20_B_bs_20_P_fit_slope)) + 
  labs(x = "Plasticity of Bud set in 2020 (BLUP)", y = "Height Growth in 2020 (BLUP)") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) 
cor.test(gw_20_B_bs_20_P_coefs$traitBudSet2020_Plasticity,gw_20_B_bs_20_P_coefs$traitHeight_Growth)

