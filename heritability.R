# heritability and genetic correlation between traits

## notes
# gw - growth
# bb - budbreak
# bs - budset
# B - lmer family blups (means)
# P - plasticity value for the trait

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
                     G2=list(V=1, n=0.002),
                     G3=list(V=1, n=0.002)))  # G2 = for second random effect, here its Population

# parameter expansion  - trying a fix for budset and budbreak 2020 data  
extPrior1 <- list(R=list(V=1, n=1), 
                  G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=100),
                         G2=list(V=1,nu=1,alpha.mu=0,alpha.V=100),
                         G3=list(V=1,nu=1,alpha.mu=0,alpha.V=100)))



########################################################################################################
#--------------------------------------------Heritability----------------------------------------------#
########################################################################################################

############################################ Budset 2019 ###############################################

# per site
budset_19_VT <- budset_19 %>% filter(Garden == "Vermont")
budset_19_MD <- budset_19 %>% filter(Garden == "Maryland")
budset_19_NC <- budset_19 %>% filter(Garden == "North_Carolina")

## Vermont
bs_19_mod_VT <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_19_VT,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_19_mod_VT,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_19_mod_VT")

## Maryland
bs_19_mod_MD <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_19_MD,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_19_mod_MD,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_19_mod_MD")

## North Carolina
bs_19_mod_NC <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_19_NC,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_19_mod_NC,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_19_mod_NC")

############################################ Budset 2020 ###############################################

############################################ normal prior ##############################################
# per site
budset_20_VT <- budset_20 %>% filter(Garden == "Vermont")
budset_20_MD <- budset_20 %>% filter(Garden == "Maryland")
budset_20_NC <- budset_20 %>% filter(Garden == "North_Carolina")

## Vermont
bs_20_mod_VT <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_VT,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_20_mod_VT,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_VT")

## Maryland
bs_20_mod_MD <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_MD,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_20_mod_MD,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_MD")

## North Carolina
bs_20_mod_NC <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_NC,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bs_20_mod_NC,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_NC")

############################################ extended prior ############################################

## Vermont
bs_20_mod_VT_ext <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_VT,
                          prior=extPrior1, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)


saveRDS(bs_20_mod_VT_ext,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_VT_ext")

## Maryland
bs_20_mod_MD_ext <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_MD,
                          prior=extPrior1, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)


saveRDS(bs_20_mod_MD_ext,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_MD_ext")

## North Carolina
bs_20_mod_NC_ext <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budset_20_NC,
                          prior=extPrior1, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)


saveRDS(bs_20_mod_NC_ext,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bs_20_mod_NC_ext")

########################################### Budbreak 2020 ##############################################

# per site
budbreak_20_VT <- budbreak_20 %>% filter(Garden == "Vermont")
budbreak_20_MD <- budbreak_20 %>% filter(Garden == "Maryland")
budbreak_20_NC <- budbreak_20 %>% filter(Garden == "North_Carolina")

## Vermont
bb_20_mod_VT <-  MCMCglmm(cGDD ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budbreak_20_VT,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bb_20_mod_VT,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bb_20_mod_VT")

## Maryland
bb_20_mod_MD <-  MCMCglmm(cGDD ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budbreak_20_MD,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bb_20_mod_MD,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bb_20_mod_MD")

## North Carolina
bb_20_mod_NC <-  MCMCglmm(cGDD ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=budbreak_20_NC,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(bb_20_mod_NC,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/bb_20_mod_NC")

############################################ Growth 2019 ###############################################

# per site
heightGrowth_19_VT <- heightGrowth_19 %>% filter(Garden == "Vermont")
heightGrowth_19_MD <- heightGrowth_19 %>% filter(Garden == "Maryland")
heightGrowth_19_NC <- heightGrowth_19 %>% filter(Garden == "North_Carolina")

## Vermont
gw_19_mod_VT <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_19_VT,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_19_mod_VT,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_19_mod_VT")

## Maryland
gw_19_mod_MD <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_19_MD,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_19_mod_MD,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_19_mod_MD")

## North Carolina
gw_19_mod_NC <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_19_NC,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_19_mod_NC,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_19_mod_NC")

############################################ Growth 2020 ###############################################

# per site
heightGrowth_20_VT <- heightGrowth_20 %>% filter(Garden == "Vermont")
heightGrowth_20_MD <- heightGrowth_20 %>% filter(Garden == "Maryland")
heightGrowth_20_NC <- heightGrowth_20 %>% filter(Garden == "North_Carolina")

## Vermont
gw_20_mod_VT <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_20_VT,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_20_mod_VT,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_20_mod_VT")

## Maryland
gw_20_mod_MD <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_20_MD,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_20_mod_MD,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_20_mod_MD")

## North Carolina
gw_20_mod_NC <-  MCMCglmm(Growth ~ 1, 
                          random= ~ Family + Population + Bed,
                          family="gaussian",
                          data=heightGrowth_20_NC,
                          prior=Prior, pr=TRUE, burnin=10000 , nitt=10000000 , thin=1000)
# prior=Prior, pr=TRUE, burnin=10000 , nitt=1000000 , thin=1000)

saveRDS(gw_20_mod_NC,"/home/Anoob/mydata/Anoob/MCMCglmm/heritability_mcmcglmm_outputs/gw_20_mod_NC")












