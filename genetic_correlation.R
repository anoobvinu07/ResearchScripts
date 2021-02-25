# genetic correlation between traits

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
# Prior <- list(R=list(V=1, n=0.002),           # R - prior on residual variance  
#               G=list(G1=list(V=1, n=0.002),   # G prior for random variance # G1 = for first randomeffect, here its Family
#                      G2=list(V=1, n=0.002)))  # G2 = for second random effect, here its Population
# 
# # parameter expansion  - trying a fix for budset and budbreak 2020 data for heritability estimation 
# extPrior1 <- list(R=list(V=1, n=1), 
#                   G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=100),
#                          G2=list(V=1,nu=1,alpha.mu=0,alpha.V=100)))

# prior for genetic correlation for just one random effect 
prior_corr <- list(R = list(V = diag(2), nu = 0.002),
                   G = list(G1 = list(V = diag(2), nu = 2,
                                      alpha.mu = rep(0,2),
                                      alpha.V = diag(25^2,2,2))))

# prior for genetic correlation for two random effect 
# prior_corr_2 <- list(R = list(V = diag(2), nu = 0.002),
#                    G = list(G1 = list(V = diag(2), nu = 2,
#                                       alpha.mu = rep(0,2),
#                                       alpha.V = diag(25^2,2,2)),
#                             G2 = list(V = diag(2), nu = 2,
#                                       alpha.mu = rep(0,2),
#                                       alpha.V = diag(25^2,2,2))))

# prior_corr_2 <- list(R = list(V = diag(2), nu = 0.002),
#                      G = list(G1 = list(V = diag(2), nu = 3,
#                                         alpha.mu = rep(0,2),
#                                         alpha.V = diag(25^2,2,2)),
#                               G2 = list(V = diag(2), nu = 3,
#                                         alpha.mu = rep(0,2),
#                                         alpha.V = diag(25^2,2,2))))




# Phenology vs growth----------------------------------------------------------------------------------#

######################################### Budset 2019-Growth 2019 ######################################

# against current year
bs_19_gw_19 <- merge(budset_19,heightGrowth_19[,c("plant_ID","Growth")])
bs_19_gw_19 <- bs_19_gw_19[!is.na(bs_19_gw_19$BudSet),]
bs_19_gw_19 <- bs_19_gw_19[!is.na(bs_19_gw_19$Growth),]

bs_19_gw_19 <- bs_19_gw_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_19_gw_19_corr <- MCMCglmm(cbind(scale(BudSet),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bs_19_gw_19,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=100)

saveRDS(bs_19_gw_19_corr,"./genetic_correlation_outputs/bs_19_gw_19_corr")





######################################### Budset 2020-Growth 2020 ######################################

# against current year
bs_20_gw_20 <- merge(budset_20,heightGrowth_20[,c("plant_ID","Growth")])
bs_20_gw_20 <- bs_20_gw_20[!is.na(bs_20_gw_20$BudSet),]
bs_20_gw_20 <- bs_20_gw_20[!is.na(bs_20_gw_20$Growth),]

bs_20_gw_20 <- bs_20_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_20_gw_20_corr <- MCMCglmm(cbind(scale(BudSet),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bs_20_gw_20,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=100)

saveRDS(bs_20_gw_20_corr,"./genetic_correlation_outputs/bs_20_gw_20_corr")


######################################## Budbreak 2020-Growth 2019 #####################################

# against previous year
bb_20_gw_19 <- merge(budbreak_20,heightGrowth_19[,c("plant_ID","Growth")])
bb_20_gw_19 <- bb_20_gw_19[!is.na(bb_20_gw_19$cGDD),]
bb_20_gw_19 <- bb_20_gw_19[!is.na(bb_20_gw_19$Growth),]

bb_20_gw_19 <- bb_20_gw_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_gw_19_corr <- MCMCglmm(cbind(scale(cGDD),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bb_20_gw_19,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=100)

saveRDS(bb_20_gw_19_corr,"./genetic_correlation_outputs/bb_20_gw_19_corr")






######################################## Budbreak 2020-Growth 2020 #####################################

# against current year
bb_20_gw_20 <- merge(budbreak_20,heightGrowth_20[,c("plant_ID","Growth")])
bb_20_gw_20 <- bb_20_gw_20[!is.na(bb_20_gw_20$cGDD),]
bb_20_gw_20 <- bb_20_gw_20[!is.na(bb_20_gw_20$Growth),]

bb_20_gw_20 <- bb_20_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_gw_20_corr <- MCMCglmm(cbind(scale(cGDD),scale(Growth)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_gw_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=100)

saveRDS(bb_20_gw_20_corr,"./genetic_correlation_outputs/bb_20_gw_20_corr")



# Phenology vs phenology-------------------------------------------------------------------------------#

######################################## Budbreak 2020-Budset 2019 #####################################

# against previous year
bb_20_bs_19 <- merge(budbreak_20,budset_19[,c("plant_ID","BudSet")])
bb_20_bs_19 <- bb_20_bs_19[!is.na(bb_20_bs_19$BudSet),]
bb_20_bs_19 <- bb_20_bs_19[!is.na(bb_20_bs_19$cGDD),]

bb_20_bs_19 <- bb_20_bs_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_bs_19_corr <- MCMCglmm(cbind(scale(cGDD),scale(BudSet)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_bs_19,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=100)

saveRDS(bb_20_bs_19_corr,"./genetic_correlation_outputs/bb_20_bs_19_corr")






######################################## Budbreak 2020-Budset 2020 #####################################

# against current year
bb_20_bs_20 <- merge(budbreak_20,budset_20[,c("plant_ID","BudSet")])
bb_20_bs_20 <- bb_20_bs_20[!is.na(bb_20_bs_20$BudSet),]
bb_20_bs_20 <- bb_20_bs_20[!is.na(bb_20_bs_20$cGDD),]

bb_20_bs_20 <- bb_20_bs_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_bs_20_corr <- MCMCglmm(cbind(scale(cGDD),scale(BudSet)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_bs_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=100)

saveRDS(bb_20_bs_20_corr,"./genetic_correlation_outputs/bb_20_bs_20_corr")

# within trait correlation-----------------------------------------------------------------------------#
######################################### Budset 2019-Budset 2020 ######################################
budset_2019 <- budset_19
colnames(budset_2019)[colnames(budset_2019)=="BudSet"] <- "BudSet2019"
budset_2020 <- budset_20
colnames(budset_2020)[colnames(budset_2020)=="BudSet"] <- "BudSet2020"

bs_20_bs_20 <- merge(budset_2019,budset_2020[,c("plant_ID","BudSet2020")])
bs_20_bs_20 <- bs_20_bs_20[!is.na(bs_20_bs_20$BudSet2020),]
bs_20_bs_20 <- bs_20_bs_20[!is.na(bs_20_bs_20$BudSet2019),]

bs_20_bs_20 <- bs_20_bs_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_20_bs_20_corr <- MCMCglmm(cbind(scale(BudSet2019),scale(BudSet2020)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bs_20_bs_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=100)

saveRDS(bs_20_bs_20_corr,"./genetic_correlation_outputs/bs_20_bs_20_corr")

######################################### Growth 2019-Growth 2020 ######################################
heightGrowth_2019 <- heightGrowth_19
colnames(heightGrowth_2019)[colnames(heightGrowth_2019)=="Growth"] <- "Growth2019"
heightGrowth_2020 <- heightGrowth_20
colnames(heightGrowth_2020)[colnames(heightGrowth_2020)=="Growth"] <- "Growth2020"

gw_20_gw_20 <- merge(heightGrowth_2019,heightGrowth_2020[,c("plant_ID","Growth2020")])
gw_20_gw_20 <- gw_20_gw_20[!is.na(gw_20_gw_20$Growth2020),]
gw_20_gw_20 <- gw_20_gw_20[!is.na(gw_20_gw_20$Growth2020),]

gw_20_gw_20 <- gw_20_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

gw_20_gw_20_corr <- MCMCglmm(cbind(scale(Growth2019),scale(Growth2020)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=gw_20_gw_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=100)

saveRDS(gw_20_gw_20_corr,"./genetic_correlation_outputs/gw_20_gw_20_corr")

