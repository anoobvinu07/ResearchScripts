# genetic correlation plots

## notes
# gw - growth
# bb - budbreak
# bs - budset
# B - lmer family blups (means)
# P - plasticity value for the trait

# notes from: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2017q1/025532.html  
# multiresponse mcmc models, page 92: https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf 

# packages
require(MCMCglmm)
require(lmerTest)
require(dplyr)
require(tidyr)
require(plotly)
require(plyr)
# require(grid)
# require(gridExtra)

# setwd("Z:/Anoob/MCMCglmm")

exome_meta <- read.table("./Exome/RS_Exome_metadata.txt", sep="\t", header = T)

# data 
budset_19 <- read.csv("./trait_data/BudSet_2019.csv")
budset_20 <- read.csv("./trait_data/BudSet_2020.csv")
budbreak_20 <- read.csv("./trait_data/BudBreak_2020_cGDD.csv")
heightGrowth_19 <- read.csv("./trait_data/Growth_2019.csv")
heightGrowth_20 <- read.csv("./trait_data/Growth_2020.csv")

# Phenology vs growth----------------------------------------------------------------------------------#

######################################### Budset 2019-Growth 2019 ######################################

# against current year
bs_19_gw_19_corr <- readRDS("./genetic_correlation_outputs/bs_19_gw_19_corr")

effectiveSize(bs_19_gw_19_corr$VCV)
heidel.diag(bs_19_gw_19_corr$VCV)

## setting up dataframe for correlation based on the model blups
bs_19_gw_19_coefs <- data_frame(Trait = attr(colMeans(bs_19_gw_19_corr$Sol), "names"),
                                    Value = colMeans(bs_19_gw_19_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitGrowth", "traitBudSet")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bs_19_gw_19_coefs <- merge(x=bs_19_gw_19_coefs,
                              y=exome_meta[,c("Pop","Region")],
                              all.x=T,
                              by.x="Population",
                              by.y="Pop")

# rename regions
bs_19_gw_19_coefs$Region <- revalue(bs_19_gw_19_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bs_19_gw_19_coefs$Region <- factor(bs_19_gw_19_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bs_19_gw_19_coefs, aes(x = traitBudSet, y = traitGrowth, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  # geom_abline(intercept = 0, slope = mean(gw_19_B_gw_19_P_fit_slope)) + 
  labs(x = "Budset in 2019", y = "Height Growth in 2019") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.688", x = 0.26, y = -0.9)  

cor.test(bs_19_gw_19_coefs$traitBudSet,bs_19_gw_19_coefs$traitGrowth) #0.6883882


######################################### Budset 2020-Growth 2020 ######################################

# against current year
bs_20_gw_20_corr <- readRDS("./genetic_correlation_outputs/bs_20_gw_20_corr")

effectiveSize(bs_20_gw_20_corr$VCV)
heidel.diag(bs_20_gw_20_corr$VCV)

## setting up dataframe for correlation based on the model blups
bs_20_gw_20_coefs <- data_frame(Trait = attr(colMeans(bs_20_gw_20_corr$Sol), "names"),
                                Value = colMeans(bs_20_gw_20_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitGrowth", "traitBudSet")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bs_20_gw_20_coefs <- merge(x=bs_20_gw_20_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bs_20_gw_20_coefs$Region <- revalue(bs_20_gw_20_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bs_20_gw_20_coefs$Region <- factor(bs_20_gw_20_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
bs_20_gw_20_plot <- ggplot(bs_20_gw_20_coefs, aes(x = traitBudSet, y = traitGrowth, group = Population, color = Population)) +
                    geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
                    # geom_abline(intercept = 0, slope = mean(gw_19_B_gw_19_P_fit_slope)) + 
                    labs(x = "Budset in 2020", y = "Height Growth in 2020") + 
                    theme_bw() + theme(axis.text=element_text(size=14), 
                                       axis.title=element_text(size=14,face="bold"),
                                       legend.title = element_text(color = "black", size = 14),
                                       legend.text = element_text(color = "black", size = 13)) +
                    annotate("text", label = "Correlation between traits: 0.947", x = 0.6, y = -0.9)  
bs_20_gw_20_plot
ggplotly(bs_20_gw_20_plot)

cor.test(bs_20_gw_20_coefs$traitBudSet,bs_20_gw_20_coefs$traitGrowth) #0.9474837

######################################## Budbreak 2020-Growth 2019 #####################################

# against previous year
bb_20_gw_19_corr <- readRDS("./genetic_correlation_outputs/bb_20_gw_19_corr")

effectiveSize(bb_20_gw_19_corr$VCV)
heidel.diag(bb_20_gw_19_corr$VCV)

## setting up dataframe for correlation based on the model blups
bb_20_gw_19_coefs <- data_frame(Trait = attr(colMeans(bb_20_gw_19_corr$Sol), "names"),
                                Value = colMeans(bb_20_gw_19_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitGrowth", "traitcGDD")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bb_20_gw_19_coefs <- merge(x=bb_20_gw_19_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bb_20_gw_19_coefs$Region <- revalue(bb_20_gw_19_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bb_20_gw_19_coefs$Region <- factor(bb_20_gw_19_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bb_20_gw_19_coefs, aes(x = traitcGDD, y = traitGrowth, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  # geom_abline(intercept = 0, slope = mean(gw_19_B_gw_19_P_fit_slope)) + 
  labs(x = "Budbreak in 2020", y = "Height Growth in 2019") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.68", x = 0.2, y = -0.9)  

cor.test(bb_20_gw_19_coefs$traitcGDD,bb_20_gw_19_coefs$traitGrowth) #0.6799553 

######################################## Budbreak 2020-Growth 2020 #####################################

# against current year
bb_20_gw_20_corr <- readRDS("./genetic_correlation_outputs/bb_20_gw_20_corr")

effectiveSize(bb_20_gw_20_corr$VCV)
heidel.diag(bb_20_gw_20_corr$VCV)

## setting up dataframe for correlation based on the model blups
bb_20_gw_20_coefs <- data_frame(Trait = attr(colMeans(bb_20_gw_20_corr$Sol), "names"),
                                Value = colMeans(bb_20_gw_20_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitGrowth", "traitcGDD")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bb_20_gw_20_coefs <- merge(x=bb_20_gw_20_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bb_20_gw_20_coefs$Region <- revalue(bb_20_gw_20_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bb_20_gw_20_coefs$Region <- factor(bb_20_gw_20_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bb_20_gw_20_coefs, aes(x = traitcGDD, y = traitGrowth, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Budbreak in 2020", y = "Height Growth in 2020") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: -0.234", x = 0.15, y = -0.9)

cor.test(bb_20_gw_20_coefs$traitcGDD,bb_20_gw_20_coefs$traitGrowth) #-0.2341478 

## combined dataset-----------------------------------------------------------------------------------#
bb_20_gw_19_coefs$Dataset <- "2019"
bb_20_gw_20_coefs$Dataset <- "2020"
budbreak_growth_comb <- rbind(bb_20_gw_19_coefs,bb_20_gw_20_coefs)
# A data frame with labels for each facet
# f_labels <- data.frame(Dataset = c("a) Against previous year growth", "b) Against current year growth"), 
#                        label = c("Correlation between traits: 0.68", "Correlation between traits: -0.234"))

ggplot(budbreak_growth_comb, aes(x = traitGrowth, y = traitcGDD, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Height Growth", y = "Budbreak in cGDD") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  facet_grid(Dataset~.) 
# bb_plot + geom_text(x = 0.15, y = -0.9, aes(label = label), data = f_labels)

#-----------------------------------------------------------------------------------------------------#

# Phenology vs phenology-------------------------------------------------------------------------------#

######################################## Budbreak 2020-Budset 2019 #####################################

# against previous year
bb_20_bs_19_corr <- readRDS("./genetic_correlation_outputs/bb_20_bs_19_corr")

effectiveSize(bb_20_bs_19_corr$VCV)
heidel.diag(bb_20_bs_19_corr$VCV)

## setting up dataframe for correlation based on the model blups
bb_20_bs_19_coefs <- data_frame(Trait = attr(colMeans(bb_20_bs_19_corr$Sol), "names"),
                                Value = colMeans(bb_20_bs_19_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitBudSet", "traitcGDD")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bb_20_bs_19_coefs <- merge(x=bb_20_bs_19_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bb_20_bs_19_coefs$Region <- revalue(bb_20_bs_19_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bb_20_bs_19_coefs$Region <- factor(bb_20_bs_19_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bb_20_bs_19_coefs, aes(x = traitcGDD, y = traitBudSet, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Budbreak in 2020", y = "Budset in 2019") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.925", x = 0.2, y = -0.9) 

cor.test(bb_20_bs_19_coefs$traitcGDD,bb_20_bs_19_coefs$traitBudSet) #0.9246381 

######################################## Budbreak 2020-Budset 2020 #####################################

# against current year
bb_20_bs_20_corr <- readRDS("./genetic_correlation_outputs/bb_20_bs_20_corr")

effectiveSize(bb_20_bs_20_corr$VCV)
heidel.diag(bb_20_bs_20_corr$VCV)

## setting up dataframe for correlation based on the model blups
bb_20_bs_20_coefs <- data_frame(Trait = attr(colMeans(bb_20_bs_20_corr$Sol), "names"),
                                Value = colMeans(bb_20_bs_20_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitBudSet", "traitcGDD")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bb_20_bs_20_coefs <- merge(x=bb_20_bs_20_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bb_20_bs_20_coefs$Region <- revalue(bb_20_bs_20_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bb_20_bs_20_coefs$Region <- factor(bb_20_bs_20_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bb_20_bs_20_coefs, aes(x = traitcGDD, y = traitBudSet, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Budbreak in 2020", y = "Budset in 2020") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.544", x = 0.2, y = -0.9)

cor.test(bb_20_bs_20_coefs$traitcGDD,bb_20_bs_20_coefs$traitBudSet) #0.5435233

## combined dataset-----------------------------------------------------------------------------------#
bb_20_bs_19_coefs$Dataset <- "2019"
bb_20_bs_20_coefs$Dataset <- "2020"
budbreak_budset_comb <- rbind(bb_20_bs_19_coefs,bb_20_bs_20_coefs)
# A data frame with labels for each facet
# f_labels2 <- data.frame(Dataset = c("a) Against previous year growth", "b) Against current year growth"), 
#                        label = c("Correlation between traits: 0.68", "Correlation between traits: -0.234"))

ggplot(budbreak_budset_comb, aes(x = traitBudSet, y = traitcGDD, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Budset in DOY", y = "Budbreak in cGDD") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  facet_wrap(. ~ Dataset) 

# within trait correlation-----------------------------------------------------------------------------#
######################################### Budset 2019-Budset 2020 ######################################
bs_20_bs_20_corr <- readRDS("./genetic_correlation_outputs/bs_20_bs_20_corr")

effectiveSize(bs_20_bs_20_corr$VCV)
heidel.diag(bs_20_bs_20_corr$VCV)

## setting up dataframe for correlation based on the model blups
bs_20_bs_20_coefs <- data_frame(Trait = attr(colMeans(bs_20_bs_20_corr$Sol), "names"),
                                Value = colMeans(bs_20_bs_20_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitBudSet2019", "traitBudSet2020")) %>%
  select(-Type) %>%
  spread(Trait, Value)


bs_20_bs_20_coefs <- merge(x=bs_20_bs_20_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
bs_20_bs_20_coefs$Region <- revalue(bs_20_bs_20_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
bs_20_bs_20_coefs$Region <- factor(bs_20_bs_20_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(bs_20_bs_20_coefs, aes(x = traitBudSet2020, y = traitBudSet2019, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Budset in 2020", y = "Budset in 2019") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.466", x = 0.3, y = -0.9)


cor.test(bs_20_bs_20_coefs$traitBudSet2020,bs_20_bs_20_coefs$traitBudSet2019) #0.4657773 

######################################### Growth 2019-Growth 2020 ######################################
gw_20_gw_20_corr <- readRDS("./genetic_correlation_outputs/gw_20_gw_20_corr")

effectiveSize(gw_20_gw_20_corr$VCV)
heidel.diag(gw_20_gw_20_corr$VCV)

## setting up dataframe for correlation based on the model blups
gw_20_gw_20_coefs <- data_frame(Trait = attr(colMeans(gw_20_gw_20_corr$Sol), "names"),
                                Value = colMeans(gw_20_gw_20_corr$Sol)) %>%
  separate(Trait, c("Trait","Type","Population"), sep = "\\.", fill = "right") %>%
  filter(Type == "Population") %>%
  filter(Trait %in% c("traitGrowth2019", "traitGrowth2020")) %>%
  select(-Type) %>%
  spread(Trait, Value)


gw_20_gw_20_coefs <- merge(x=gw_20_gw_20_coefs,
                           y=exome_meta[,c("Pop","Region")],
                           all.x=T,
                           by.x="Population",
                           by.y="Pop")

# rename regions
gw_20_gw_20_coefs$Region <- revalue(gw_20_gw_20_coefs$Region, c("C"="Core", "M"="Margin","E"="Edge")) 
# reorder regions
gw_20_gw_20_coefs$Region <- factor(gw_20_gw_20_coefs$Region,levels = c("Core", "Margin", "Edge"))

# plot
ggplot(gw_20_gw_20_coefs, aes(x = traitGrowth2019, y = traitGrowth2020, group = Population, color = Population)) +
  geom_point(aes(shape = Region), size=3, alpha = 0.7) + 
  labs(x = "Height Growth in 2019", y = "Height Growth in 2020") + 
  theme_bw() + theme(axis.text=element_text(size=14), 
                     axis.title=element_text(size=14,face="bold"),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 13)) +
  annotate("text", label = "Correlation between traits: 0.688", x = 0.8, y = -0.9)


cor.test(gw_20_gw_20_coefs$traitGrowth2019,gw_20_gw_20_coefs$traitGrowth2020) #0.6884412
