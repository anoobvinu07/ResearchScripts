# heritability and genetic correlation between traits

## notes
# gw - growth
# bb - budbreak
# bs - budset
# B - lmer family blups (means)
# P - plasticity value for the trait

# setwd("Z:/Anoob/MCMCglmm")

# packages
require(MCMCglmm)
require(lmerTest)
require(dplyr)
require(tidyr)

# data 
budset_19 <- read.csv("./trait_data/BudSet_2019.csv")
budset_20 <- read.csv("./trait_data/BudSet_2020.csv")
budbreak_20 <- read.csv("./trait_data/BudBreak_2020_cGDD.csv")
heightGrowth_19 <- read.csv("./trait_data/Growth_2019.csv")
heightGrowth_20 <- read.csv("./trait_data/Growth_2020.csv")

########################################################################################################
#-----------------------------------------------Priors-------------------------------------------------#
########################################################################################################
# set priors 
Prior <- list(R=list(V=1, n=0.002),           # R - prior on residual variance  
              G=list(G1=list(V=1, n=0.002),   # G prior for random variance # G1 = for first randomeffect, here its Family
                     G2=list(V=1, n=0.002)))  # G2 = for second random effect, here its Population

# parameter expansion  - trying a fix for budset and budbreak 2020 data  
extPrior1 <- list(R=list(V=1, n=1), 
                  G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=100),
                         G2=list(V=1,nu=1,alpha.mu=0,alpha.V=100)))



########################################################################################################
#--------------------------------MCMCglmm----Heritability----------------------------------------------#
########################################################################################################

############################################ Budset 2019 ###############################################


## Vermont
bs_19_mod_VT <- readRDS("./heritability_mcmcglmm_outputs/bs_19_mod_VT")
# heritability 
bs_19_H2_VT = (bs_19_mod_VT$VCV[,"Family"]+bs_19_mod_VT$VCV[,"Population"])/ (bs_19_mod_VT$VCV[,"Population"] + bs_19_mod_VT$VCV[,"Family"] + bs_19_mod_VT$VCV[,"Bed"] + bs_19_mod_VT$VCV[,"units"])
mean(bs_19_H2_VT)
bs_19_H2_VT <- (bs_19_mod_VT$VCV[,"Family"]+bs_19_mod_VT$VCV[,"Population"])/rowSums(bs_19_mod_VT[["VCV"]])
mean(bs_19_H2_VT)

HPDinterval(bs_19_H2_VT)

hist(bs_19_H2_VT)

plot(density(bs_19_H2_VT))

## Maryland
bs_19_mod_MD <- readRDS("./heritability_mcmcglmm_outputs/bs_19_mod_MD")
# heritability 
bs_19_H2_MD <- (bs_19_mod_MD$VCV[,"Family"]+bs_19_mod_MD$VCV[,"Population"])/rowSums(bs_19_mod_MD[["VCV"]])
mean(bs_19_H2_MD)

## North Carolina
bs_19_mod_NC <- readRDS("./heritability_mcmcglmm_outputs/bs_19_mod_NC")
# heritability 
bs_19_H2_NC <- (bs_19_mod_NC$VCV[,"Family"]+bs_19_mod_NC$VCV[,"Population"])/rowSums(bs_19_mod_NC[["VCV"]])
mean(bs_19_H2_NC)

############################################ Budset 2020 ###############################################

############################################ normal prior ##############################################
## Vermont
bs_20_mod_VT <- readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_VT")

# diagnostics
plot(bs_20_mod_VT)
effectiveSize(bs_20_mod_VT$VCV)
heidel.diag(bs_20_mod_VT$VCV)

# heritability 
bs_20_H2_VT <- (bs_20_mod_VT$VCV[,"Family"]+bs_20_mod_VT$VCV[,"Population"])/rowSums(bs_20_mod_VT[["VCV"]])
mean(bs_20_H2_VT)

## Maryland
bs_20_mod_MD <- readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_MD")

# diagnostics
plot(bs_20_mod_MD)
effectiveSize(bs_20_mod_MD$VCV)
heidel.diag(bs_20_mod_MD$VCV)

# heritability 
bs_20_H2_MD <- (bs_20_mod_MD$VCV[,"Family"]+bs_20_mod_MD$VCV[,"Population"])/rowSums(bs_20_mod_MD[["VCV"]])
mean(bs_20_H2_MD)

## North Carolina
bs_20_mod_NC <-readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_NC")

# diagnostics
plot(bs_20_mod_NC)
effectiveSize(bs_20_mod_NC$VCV)
heidel.diag(bs_20_mod_NC$VCV)

# heritability 
bs_20_H2_NC <- (bs_20_mod_NC$VCV[,"Family"]+bs_20_mod_NC$VCV[,"Population"])/rowSums(bs_20_mod_NC[["VCV"]])
mean(bs_20_H2_NC)

############################################ extended prior ############################################

## Vermont
bs_20_mod_VT_ext <- readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_VT_ext")

# diagnostics
plot(bs_20_mod_VT_ext)
effectiveSize(bs_20_mod_VT_ext$VCV)
heidel.diag(bs_20_mod_VT_ext$VCV)


# heritability 
bs_20_H2_VT_ext <- (bs_20_mod_VT_ext$VCV[,"Family"]+bs_20_mod_VT_ext$VCV[,"Population"])/rowSums(bs_20_mod_VT_ext[["VCV"]])
mean(bs_20_H2_VT_ext)

posterior.mode(bs_20_H2_VT_ext) #0.04652454 
HPDinterval(bs_20_H2_VT_ext) # 0.01004681 0.09119632
plot(bs_20_H2_VT_ext)
hist(bs_20_H2_VT_ext)


## Maryland
bs_20_mod_MD_ext <- readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_MD_ext")

# diagnostics
plot(bs_20_mod_MD_ext)
effectiveSize(bs_20_mod_MD_ext$VCV)
heidel.diag(bs_20_mod_MD_ext$VCV)

# heritability 
bs_20_H2_MD_ext <- (bs_20_mod_MD_ext$VCV[,"Family"]+bs_20_mod_MD_ext$VCV[,"Population"])/rowSums(bs_20_mod_MD_ext[["VCV"]])
mean(bs_20_H2_MD_ext)

posterior.mode(bs_20_H2_MD_ext) #0.1143763 
HPDinterval(bs_20_H2_MD_ext) # 0.06623493 0.1676811
plot(bs_20_H2_MD_ext)
hist(bs_20_H2_MD_ext)


## North Carolina
bs_20_mod_NC_ext <- readRDS("./heritability_mcmcglmm_outputs/bs_20_mod_NC_ext")

# diagnostics
plot(bs_20_mod_NC_ext)
effectiveSize(bs_20_mod_NC_ext$VCV)
heidel.diag(bs_20_mod_NC_ext$VCV)

# heritability 
bs_20_H2_NC_ext <- (bs_20_mod_NC_ext$VCV[,"Family"]+bs_20_mod_NC_ext$VCV[,"Population"])/rowSums(bs_20_mod_NC_ext[["VCV"]])
mean(bs_20_H2_NC_ext)

posterior.mode(bs_20_H2_NC_ext) #0.1222343
HPDinterval(bs_20_H2_NC_ext) # 0.06608991 0.1768585
plot(bs_20_H2_NC_ext)
hist(bs_20_H2_NC_ext)

########################################### Budbreak 2020 ##############################################
## Vermont
bb_20_mod_VT <- readRDS("./heritability_mcmcglmm_outputs/bb_20_mod_VT")
# heritability 
bb_20_H2_VT <- (bb_20_mod_VT$VCV[,"Family"]+bb_20_mod_VT$VCV[,"Population"])/rowSums(bb_20_mod_VT[["VCV"]])
mean(bb_20_H2_VT)


## Maryland
bb_20_mod_MD <- readRDS("./heritability_mcmcglmm_outputs/bb_20_mod_MD")
# heritability 
bb_20_H2_MD <- (bb_20_mod_MD$VCV[,"Family"]+bb_20_mod_MD$VCV[,"Population"])/rowSums(bb_20_mod_MD[["VCV"]])
mean(bb_20_H2_MD)

## North Carolina
bb_20_mod_NC <- readRDS("./heritability_mcmcglmm_outputs/bb_20_mod_NC")
# heritability 
bb_20_H2_NC <- (bb_20_mod_NC$VCV[,"Family"]+bb_20_mod_NC$VCV[,"Population"])/rowSums(bb_20_mod_NC[["VCV"]])
mean(bb_20_H2_NC)

############################################ Growth 2019 ###############################################

## Vermont
gw_19_mod_VT <- readRDS("./heritability_mcmcglmm_outputs/gw_19_mod_VT")
# heritability 
gw_19_H2_VT <- (gw_19_mod_VT$VCV[,"Family"]+gw_19_mod_VT$VCV[,"Population"])/rowSums(gw_19_mod_VT[["VCV"]])
mean(gw_19_H2_VT)

## Maryland
gw_19_mod_MD <- readRDS("./heritability_mcmcglmm_outputs/gw_19_mod_MD")
# heritability 
gw_19_H2_MD <- (gw_19_mod_MD$VCV[,"Family"]+gw_19_mod_MD$VCV[,"Population"])/rowSums(gw_19_mod_MD[["VCV"]])
mean(gw_19_H2_MD)

## North Carolina
gw_19_mod_NC <- readRDS("./heritability_mcmcglmm_outputs/gw_19_mod_NC")
# heritability 
gw_19_H2_NC <- (gw_19_mod_NC$VCV[,"Family"]+gw_19_mod_NC$VCV[,"Population"])/rowSums(gw_19_mod_NC[["VCV"]])
mean(gw_19_H2_NC)

############################################ Growth 2020 ###############################################

## Vermont
gw_20_mod_VT <- readRDS("./heritability_mcmcglmm_outputs/gw_20_mod_VT")
# heritability 
gw_20_H2_VT <- (gw_20_mod_VT$VCV[,"Family"]+gw_20_mod_VT$VCV[,"Population"])/rowSums(gw_20_mod_VT[["VCV"]])
mean(gw_20_H2_VT)

## Maryland
gw_20_mod_MD <- readRDS("./heritability_mcmcglmm_outputs/gw_20_mod_MD")
# heritability 
gw_20_H2_MD <- (gw_20_mod_MD$VCV[,"Family"]+gw_20_mod_MD$VCV[,"Population"])/rowSums(gw_20_mod_MD[["VCV"]])
mean(gw_20_H2_MD)

## North Carolina
gw_20_mod_NC <- readRDS("./heritability_mcmcglmm_outputs/gw_20_mod_NC")
# heritability 
gw_20_H2_NC <- (gw_20_mod_NC$VCV[,"Family"]+gw_20_mod_NC$VCV[,"Population"])/rowSums(gw_20_mod_NC[["VCV"]])
mean(gw_20_H2_NC)



########################################################################################################
#--------------------------------lmer--models----Heritability------------------------------------------#
########################################################################################################

############################################ budset 2019 ###############################################

## Vermont
bs_19_VT <- budset_19[budset_19$Garden=="Vermont",]
# lmer model
bs_19_lmer_VT <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_19_VT)
# heritability 
as.data.frame(VarCorr(bs_19_lmer_VT)) 

bs_19_lmer_H2_VT <- as.data.frame(VarCorr(bs_19_lmer_VT)) 
# or
bs_19_lmer_H2_VT <- as.data.frame(print(bs_19_lmer_H2_VT,comp=c("Variance","Std.Dev")))

bs_19_lmer_H2_VT <- (bs_19_lmer_H2_VT[bs_19_lmer_H2_VT$grp=="Family","vcov"]+bs_19_lmer_H2_VT[bs_19_lmer_H2_VT$grp=="Population","vcov"])/sum(bs_19_lmer_H2_VT$vcov)
bs_19_lmer_H2_VT

## Maryland
bs_19_MD <- budset_19[budset_19$Garden=="Maryland",]
# lmer model
bs_19_lmer_MD <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_19_MD)
# heritability 
as.data.frame(VarCorr(bs_19_lmer_MD)) 
bs_19_lmer_H2_MD <- as.data.frame(VarCorr(bs_19_lmer_MD)) 
bs_19_lmer_H2_MD <- (bs_19_lmer_H2_MD[bs_19_lmer_H2_MD$grp=="Family","vcov"]+bs_19_lmer_H2_MD[bs_19_lmer_H2_MD$grp=="Population","vcov"])/sum(bs_19_lmer_H2_MD$vcov)
bs_19_lmer_H2_MD

## North Carolina
bs_19_NC <- budset_19[budset_19$Garden=="North_Carolina",]
# lmer model
bs_19_lmer_NC <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_19_NC)
# heritability 
as.data.frame(VarCorr(bs_19_lmer_NC)) 
bs_19_lmer_H2_NC <- as.data.frame(VarCorr(bs_19_lmer_NC)) 
bs_19_lmer_H2_NC <- (bs_19_lmer_H2_NC[bs_19_lmer_H2_NC$grp=="Family","vcov"]+bs_19_lmer_H2_NC[bs_19_lmer_H2_NC$grp=="Population","vcov"])/sum(bs_19_lmer_H2_NC$vcov)
bs_19_lmer_H2_NC


############################################ budset 2020 ###############################################

## Vermont
bs_20_VT <- budset_20[budset_20$Garden=="Vermont",]
# lmer model
bs_20_lmer_VT <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_20_VT)
# heritability 
as.data.frame(VarCorr(bs_20_lmer_VT)) 

bs_20_lmer_H2_VT <- as.data.frame(VarCorr(bs_20_lmer_VT)) 
# or
bs_20_lmer_H2_VT <- as.data.frame(print(bs_20_lmer_H2_VT,comp=c("Variance","Std.Dev")))

bs_20_lmer_H2_VT <- (bs_20_lmer_H2_VT[bs_20_lmer_H2_VT$grp=="Family","vcov"]+bs_20_lmer_H2_VT[bs_20_lmer_H2_VT$grp=="Population","vcov"])/sum(bs_20_lmer_H2_VT$vcov)
bs_20_lmer_H2_VT

## Maryland
bs_20_MD <- budset_20[budset_20$Garden=="Maryland",]
# lmer model
bs_20_lmer_MD <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_20_MD)
# heritability 
as.data.frame(VarCorr(bs_20_lmer_MD)) 
bs_20_lmer_H2_MD <- as.data.frame(VarCorr(bs_20_lmer_MD)) 
bs_20_lmer_H2_MD <- (bs_20_lmer_H2_MD[bs_20_lmer_H2_MD$grp=="Family","vcov"]+bs_20_lmer_H2_MD[bs_20_lmer_H2_MD$grp=="Population","vcov"])/sum(bs_20_lmer_H2_MD$vcov)
bs_20_lmer_H2_MD

## North Carolina
bs_20_NC <- budset_20[budset_20$Garden=="North_Carolina",]
# lmer model
bs_20_lmer_NC <- lmer(BudSet ~ (1|Family) + (1|Population) + (1|Bed), data=bs_20_NC)
# heritability 
as.data.frame(VarCorr(bs_20_lmer_NC)) 
bs_20_lmer_H2_NC <- as.data.frame(VarCorr(bs_20_lmer_NC)) 
bs_20_lmer_H2_NC <- (bs_20_lmer_H2_NC[bs_20_lmer_H2_NC$grp=="Family","vcov"]+bs_20_lmer_H2_NC[bs_20_lmer_H2_NC$grp=="Population","vcov"])/sum(bs_20_lmer_H2_NC$vcov)
bs_20_lmer_H2_NC


############################################ budbreak 2020 ###############################################

## Vermont
bb_20_VT <- budbreak_20[budbreak_20$Garden=="Vermont",]
# lmer model
bb_20_lmer_VT <- lmer(cGDD ~ (1|Family) + (1|Population) + (1|Bed), data=bb_20_VT)
# heritability 
as.data.frame(VarCorr(bb_20_lmer_VT)) 

bb_20_lmer_H2_VT <- as.data.frame(VarCorr(bb_20_lmer_VT)) 
# or
bb_20_lmer_H2_VT <- as.data.frame(print(bb_20_lmer_H2_VT,comp=c("Variance","Std.Dev")))

bb_20_lmer_H2_VT <- (bb_20_lmer_H2_VT[bb_20_lmer_H2_VT$grp=="Family","vcov"]+bb_20_lmer_H2_VT[bb_20_lmer_H2_VT$grp=="Population","vcov"])/sum(bb_20_lmer_H2_VT$vcov)
bb_20_lmer_H2_VT

## Maryland
bb_20_MD <- budbreak_20[budbreak_20$Garden=="Maryland",]
# lmer model
bb_20_lmer_MD <- lmer(cGDD ~ (1|Family) + (1|Population) + (1|Bed), data=bb_20_MD)
# heritability 
as.data.frame(VarCorr(bb_20_lmer_MD)) 
bb_20_lmer_H2_MD <- as.data.frame(VarCorr(bb_20_lmer_MD)) 
bb_20_lmer_H2_MD <- (bb_20_lmer_H2_MD[bb_20_lmer_H2_MD$grp=="Family","vcov"]+bb_20_lmer_H2_MD[bb_20_lmer_H2_MD$grp=="Population","vcov"])/sum(bb_20_lmer_H2_MD$vcov)
bb_20_lmer_H2_MD

## North Carolina
bb_20_NC <- budbreak_20[budbreak_20$Garden=="North_Carolina",]
# lmer model
bb_20_lmer_NC <- lmer(cGDD ~ (1|Family) + (1|Population) + (1|Bed), data=bb_20_NC)
# heritability 
as.data.frame(VarCorr(bb_20_lmer_NC)) 
bb_20_lmer_H2_NC <- as.data.frame(VarCorr(bb_20_lmer_NC)) 
bb_20_lmer_H2_NC <- (bb_20_lmer_H2_NC[bb_20_lmer_H2_NC$grp=="Family","vcov"]+bb_20_lmer_H2_NC[bb_20_lmer_H2_NC$grp=="Population","vcov"])/sum(bb_20_lmer_H2_NC$vcov)
bb_20_lmer_H2_NC



############################################ Growth 2019 ###############################################

## Vermont
gw_19_VT <- heightGrowth_19[heightGrowth_19$Garden=="Vermont",]
# lmer model
gw_19_lmer_VT <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_19_VT)
# heritability 
as.data.frame(VarCorr(gw_19_lmer_VT)) 

gw_19_lmer_H2_VT <- as.data.frame(VarCorr(gw_19_lmer_VT)) 
# or
gw_19_lmer_H2_VT <- as.data.frame(print(gw_19_lmer_H2_VT,comp=c("Variance","Std.Dev")))

gw_19_lmer_H2_VT <- (gw_19_lmer_H2_VT[gw_19_lmer_H2_VT$grp=="Family","vcov"]+gw_19_lmer_H2_VT[gw_19_lmer_H2_VT$grp=="Population","vcov"])/sum(gw_19_lmer_H2_VT$vcov)
gw_19_lmer_H2_VT

## Maryland
gw_19_MD <- heightGrowth_19[heightGrowth_19$Garden=="Maryland",]
# lmer model
gw_19_lmer_MD <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_19_MD)
# heritability 
as.data.frame(VarCorr(gw_19_lmer_MD)) 
gw_19_lmer_H2_MD <- as.data.frame(VarCorr(gw_19_lmer_MD)) 
gw_19_lmer_H2_MD <- (gw_19_lmer_H2_MD[gw_19_lmer_H2_MD$grp=="Family","vcov"]+gw_19_lmer_H2_MD[gw_19_lmer_H2_MD$grp=="Population","vcov"])/sum(gw_19_lmer_H2_MD$vcov)
gw_19_lmer_H2_MD

## North Carolina
gw_19_NC <- heightGrowth_19[heightGrowth_19$Garden=="North_Carolina",]
# lmer model
gw_19_lmer_NC <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_19_NC)
# heritability 
as.data.frame(VarCorr(gw_19_lmer_NC)) 
gw_19_lmer_H2_NC <- as.data.frame(VarCorr(gw_19_lmer_NC)) 
gw_19_lmer_H2_NC <- (gw_19_lmer_H2_NC[gw_19_lmer_H2_NC$grp=="Family","vcov"]+gw_19_lmer_H2_NC[gw_19_lmer_H2_NC$grp=="Population","vcov"])/sum(gw_19_lmer_H2_NC$vcov)
gw_19_lmer_H2_NC


############################################ Growth 2020 ###############################################

## Vermont
gw_20_VT <- heightGrowth_20[heightGrowth_20$Garden=="Vermont",]
# lmer model
gw_20_lmer_VT <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_20_VT)
# heritability 
as.data.frame(VarCorr(gw_20_lmer_VT)) 

gw_20_lmer_H2_VT <- as.data.frame(VarCorr(gw_20_lmer_VT)) 
# or
gw_20_lmer_H2_VT <- as.data.frame(print(gw_20_lmer_H2_VT,comp=c("Variance","Std.Dev")))

gw_20_lmer_H2_VT <- (gw_20_lmer_H2_VT[gw_20_lmer_H2_VT$grp=="Family","vcov"]+gw_20_lmer_H2_VT[gw_20_lmer_H2_VT$grp=="Population","vcov"])/sum(gw_20_lmer_H2_VT$vcov)
gw_20_lmer_H2_VT

## Maryland
gw_20_MD <- heightGrowth_20[heightGrowth_20$Garden=="Maryland",]
# lmer model
gw_20_lmer_MD <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_20_MD)
# heritability 
as.data.frame(VarCorr(gw_20_lmer_MD)) 
gw_20_lmer_H2_MD <- as.data.frame(VarCorr(gw_20_lmer_MD)) 
gw_20_lmer_H2_MD <- (gw_20_lmer_H2_MD[gw_20_lmer_H2_MD$grp=="Family","vcov"]+gw_20_lmer_H2_MD[gw_20_lmer_H2_MD$grp=="Population","vcov"])/sum(gw_20_lmer_H2_MD$vcov)
gw_20_lmer_H2_MD

## North Carolina
gw_20_NC <- heightGrowth_20[heightGrowth_20$Garden=="North_Carolina",]
# lmer model
gw_20_lmer_NC <- lmer(Growth ~ (1|Family) + (1|Population) + (1|Bed), data=gw_20_NC)
# heritability 
as.data.frame(VarCorr(gw_20_lmer_NC)) 
gw_20_lmer_H2_NC <- as.data.frame(VarCorr(gw_20_lmer_NC)) 
gw_20_lmer_H2_NC <- (gw_20_lmer_H2_NC[gw_20_lmer_H2_NC$grp=="Family","vcov"]+gw_20_lmer_H2_NC[gw_20_lmer_H2_NC$grp=="Population","vcov"])/sum(gw_20_lmer_H2_NC$vcov)
gw_20_lmer_H2_NC







