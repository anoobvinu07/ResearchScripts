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
budset_2019 <- read.csv("./trait_data/BudSet_2019.csv")
budset_2020 <- read.csv("./trait_data/BudSet_2020.csv")
budbreak_2020 <- read.csv("./trait_data/BudBreak_2020_cGDD.csv")
heightGrowth_2019 <- read.csv("./trait_data/Growth_2019.csv")
heightGrowth_2020 <- read.csv("./trait_data/Growth_2020.csv")
Plasticity <- read.csv("./trait_data/Plasticity.csv")




#########################################################
## prior
#########################################################
prior_corr <- list(R = list(V = diag(2), nu = 0.002),
                   G = list(G1 = list(V = diag(2), nu = 2,
                                      alpha.mu = rep(0,2),
                                      alpha.V = diag(25^2,2,2))))
# An uninformative prior for the correlations is an improper prior with V=diag(dim(V))???0 and nu=dim(V)+1.

#---------------------------------------Height growth-----------------------------------------------#

#########################################################
## Height growth of 2019 - Means vs Plasticity  
#########################################################

# blups for height growth
heightGrowth_2019_2 <- heightGrowth_2019
heightGrowth_2019_2 <- heightGrowth_2019_2%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2019_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2019_2)

HeightGrowth_2019_blup <- ranef(HeightGrowth_2019_mod)
HeightGrowth_2019_blup <- HeightGrowth_2019_blup$Family
HeightGrowth_2019_blup <- cbind(Family=rownames(HeightGrowth_2019_blup),HeightGrowth_2019_blup)
rownames(HeightGrowth_2019_blup) <- NULL
names(HeightGrowth_2019_blup)[2] <- "Height_Growth"

# plasticity based on height growth 
HeightGrowth_2019_blup <- merge(x=HeightGrowth_2019_blup,
                                y=Plasticity[,c("Family","growth2019_Plasticity")],
                                by="Family")
# rename plasticity column
colnames(HeightGrowth_2019_blup)[colnames(HeightGrowth_2019_blup)=="growth2019_Plasticity"] <- "Plasticity"

# remove NAs
HeightGrowth_2019_blup <- HeightGrowth_2019_blup[!is.na(HeightGrowth_2019_blup$Plasticity),]

# Population column creation
HeightGrowth_2019_blup$Population <- gsub("\\_.*", "", HeightGrowth_2019_blup$Family)

# Correlation between the two
gw_19_B_gw_19_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(Plasticity)) ~ trait - 1, 
                                     # fixed=Region,
                                     random=~us(trait):Population,
                                     rcov=~us(trait):units, # or idh
                                     family=c("gaussian","gaussian"),
                                     prior=prior_corr, 
                                     # pedigree=Ped, 
                                     data=HeightGrowth_2019_blup,
                                     pr=TRUE, # saves the posterior distribution of the individual random effects
                                     # analagous to the BLUP from the REML analysis, so we can visualise them later
                                     verbose=TRUE,
                                     nitt=10000000, 
                                     burnin=10000, 
                                     thin=1000)

# save the output
saveRDS(gw_19_B_gw_19_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_19_B_gw_19_P_corr")

#########################################################
## Height growth of 2020 - Means vs Plasticity  
#########################################################

# blups for height growth
heightGrowth_2020_2 <- heightGrowth_2020
heightGrowth_2020_2 <- heightGrowth_2020_2%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2020_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2020_2)

HeightGrowth_2020_blup <- ranef(HeightGrowth_2020_mod)
HeightGrowth_2020_blup <- HeightGrowth_2020_blup$Family
HeightGrowth_2020_blup <- cbind(Family=rownames(HeightGrowth_2020_blup),HeightGrowth_2020_blup)
rownames(HeightGrowth_2020_blup) <- NULL
names(HeightGrowth_2020_blup)[2] <- "Height_Growth"

# plasticity based on height growth 
HeightGrowth_2020_blup <- merge(x=HeightGrowth_2020_blup,
                                y=Plasticity[,c("Family","growth2020_Plasticity")],
                                by="Family")
# rename plasticity column
colnames(HeightGrowth_2020_blup)[colnames(HeightGrowth_2020_blup)=="growth2020_Plasticity"] <- "Plasticity"

# remove NAs
HeightGrowth_2020_blup <- HeightGrowth_2020_blup[!is.na(HeightGrowth_2020_blup$Plasticity),]

# Population column creation
HeightGrowth_2020_blup$Population <- gsub("\\_.*", "", HeightGrowth_2020_blup$Family)

# Correlation between the two
gw_20_B_gw_20_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(Plasticity)) ~ trait - 1, 
                                 # fixed=Region,
                                 random=~us(trait):Population,
                                 rcov=~us(trait):units, # or idh
                                 family=c("gaussian","gaussian"),
                                 prior=prior_corr, 
                                 # pedigree=Ped, 
                                 data=HeightGrowth_2020_blup,
                                 pr=TRUE, 
                                 verbose=TRUE,
                                 nitt=10000000, 
                                 burnin=10000, 
                                 thin=1000)

# save the output
saveRDS(gw_20_B_gw_20_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_20_B_gw_20_P_corr")

#----------------------------------Height growth vs budbreak-------------------------------------#

##########################################################
## Height growth of 2019 means vs budbreak 2020 Plasticity  
##########################################################

# blups for height growth
heightGrowth_2019_2 <- heightGrowth_2019
heightGrowth_2019_2 <- heightGrowth_2019_2%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2019_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2019_2)

HeightGrowth_2019_blup <- ranef(HeightGrowth_2019_mod)
HeightGrowth_2019_blup <- HeightGrowth_2019_blup$Family
HeightGrowth_2019_blup <- cbind(Family=rownames(HeightGrowth_2019_blup),HeightGrowth_2019_blup)
rownames(HeightGrowth_2019_blup) <- NULL
names(HeightGrowth_2019_blup)[2] <- "Height_Growth"

# plasticity based on budbreak 2020 
gw_19_B_bb_20_P <- merge(x=HeightGrowth_2019_blup,
                         y=Plasticity[,c("Family","BudBreak2020_Plasticity")],
                         by="Family")
# rename plasticity column if needed
colnames(gw_19_B_bb_20_P)[colnames(gw_19_B_bb_20_P)=="BudBreak2020_Plasticity"] <- "BudBreak2020_Plasticity"

# remove NAs
gw_19_B_bb_20_P <- gw_19_B_bb_20_P[!is.na(gw_19_B_bb_20_P$BudBreak2020_Plasticity),]

# Population column creation
gw_19_B_bb_20_P$Population <- gsub("\\_.*", "", gw_19_B_bb_20_P$Family)

# Correlation between the two
gw_19_B_bb_20_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(BudBreak2020_Plasticity)) ~ trait - 1, 
                                 # fixed=Region,
                                 random=~us(trait):Population,
                                 rcov=~us(trait):units, # or idh
                                 family=c("gaussian","gaussian"),
                                 prior=prior_corr, 
                                 # pedigree=Ped, 
                                 data=gw_19_B_bb_20_P,
                                 pr=TRUE, 
                                 verbose=TRUE,
                                 nitt=10000000, 
                                 burnin=10000, 
                                 thin=1000)
# save outputs
saveRDS(gw_19_B_bb_20_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_19_B_bb_20_P_corr")

##########################################################
## Height growth of 2020 means vs budbreak 2020 Plasticity
##########################################################

# blups for height growth
heightGrowth_2020_2 <- heightGrowth_2020
heightGrowth_2020_2 <- heightGrowth_2020_2 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2020_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2020_2)

HeightGrowth_2020_blup <- ranef(HeightGrowth_2020_mod)
HeightGrowth_2020_blup <- HeightGrowth_2020_blup$Family
HeightGrowth_2020_blup <- cbind(Family=rownames(HeightGrowth_2020_blup),HeightGrowth_2020_blup)
rownames(HeightGrowth_2020_blup) <- NULL
names(HeightGrowth_2020_blup)[2] <- "Height_Growth"

# plasticity based on budbreak 2020 
gw_20_B_bb_20_P <- merge(x=HeightGrowth_2020_blup,
                         y=Plasticity[,c("Family","BudBreak2020_Plasticity")],
                         by="Family")
# rename plasticity column if needed
colnames(gw_20_B_bb_20_P)[colnames(gw_20_B_bb_20_P)=="BudBreak2020_Plasticity"] <- "BudBreak2020_Plasticity"

# remove NAs
gw_20_B_bb_20_P <- gw_20_B_bb_20_P[!is.na(gw_20_B_bb_20_P$BudBreak2020_Plasticity),]

# Population column creation
gw_20_B_bb_20_P$Population <- gsub("\\_.*", "", gw_20_B_bb_20_P$Family)

# Correlation between the two
gw_20_B_bb_20_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(BudBreak2020_Plasticity)) ~ trait - 1, 
                                 # fixed=Region,
                                 random=~us(trait):Population,
                                 rcov=~us(trait):units, # or idh
                                 family=c("gaussian","gaussian"),
                                 prior=prior_corr, 
                                 # pedigree=Ped, 
                                 data=gw_20_B_bb_20_P,
                                 pr=TRUE, 
                                 verbose=TRUE,
                                 nitt=10000000, 
                                 burnin=10000, 
                                 thin=1000)

saveRDS(gw_20_B_bb_20_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_20_B_bb_20_P_corr")


#------------------------------------Height growth vs budset---------------------------------------#

##########################################################
## Height growth of 2019 means vs budset 2019 Plasticity  
##########################################################

# blups for height growth
heightGrowth_2019_2 <- heightGrowth_2019
heightGrowth_2019_2 <- heightGrowth_2019_2%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2019_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2019_2)

# HeightGrowth_2019_blup <- coef(HeightGrowth_2019_mod)

HeightGrowth_2019_blup <- ranef(HeightGrowth_2019_mod)

HeightGrowth_2019_blup <- HeightGrowth_2019_blup$Family
HeightGrowth_2019_blup <- cbind(Family=rownames(HeightGrowth_2019_blup),HeightGrowth_2019_blup)
rownames(HeightGrowth_2019_blup) <- NULL
names(HeightGrowth_2019_blup)[2] <- "Height_Growth"

# plasticity based on budset 2019 
gw_19_B_bs_19_P <- merge(x=HeightGrowth_2019_blup,
                         y=Plasticity[,c("Family","BudSet2019_Plasticity")],
                         by="Family")
# rename plasticity column if needed
# colnames(gw_19_B_bs_19_P)[colnames(gw_19_B_bs_19_P)=="BudSet2019_Plasticity"] <- "BudSet2019_Plasticity"

# remove NAs
gw_19_B_bs_19_P <- gw_19_B_bs_19_P[!is.na(gw_19_B_bs_19_P$BudSet2019_Plasticity),]

# Population column creation
gw_19_B_bs_19_P$Population <- gsub("\\_.*", "", gw_19_B_bs_19_P$Family)

# Correlation between the two
gw_19_B_bs_19_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(BudSet2019_Plasticity)) ~ trait - 1, 
                                 # fixed=Region,
                                 random=~us(trait):Population,
                                 rcov=~us(trait):units, # or idh
                                 family=c("gaussian","gaussian"),
                                 prior=prior_corr, 
                                 # pedigree=Ped, 
                                 data=gw_19_B_bs_19_P,
                                 pr=TRUE, 
                                 verbose=TRUE,
                                 nitt=10000000, 
                                 burnin=10000, 
                                 thin=1000)
# save outputs
saveRDS(gw_19_B_bs_19_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_19_B_bs_19_P_corr")

##########################################################
## Height growth of 2020 means vs budset 2020 Plasticity
##########################################################

# blups for height growth
heightGrowth_2020_2 <- heightGrowth_2020
heightGrowth_2020_2 <- heightGrowth_2020_2 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

HeightGrowth_2020_mod <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_2020_2)

HeightGrowth_2020_blup <- ranef(HeightGrowth_2020_mod)
HeightGrowth_2020_blup <- HeightGrowth_2020_blup$Family
HeightGrowth_2020_blup <- cbind(Family=rownames(HeightGrowth_2020_blup),HeightGrowth_2020_blup)
rownames(HeightGrowth_2020_blup) <- NULL
names(HeightGrowth_2020_blup)[2] <- "Height_Growth"

# plasticity based on budset 2020 
gw_20_B_bs_20_P <- merge(x=HeightGrowth_2020_blup,
                         y=Plasticity[,c("Family","BudSet2020_Plasticity")],
                         by="Family")
# rename plasticity column if needed
colnames(gw_20_B_bs_20_P)[colnames(gw_20_B_bs_20_P)=="BudSet2020_Plasticity"] <- "BudSet2020_Plasticity"

# remove NAs
gw_20_B_bs_20_P <- gw_20_B_bs_20_P[!is.na(gw_20_B_bs_20_P$BudSet2020_Plasticity),]

# Population column creation
gw_20_B_bs_20_P$Population <- gsub("\\_.*", "", gw_20_B_bs_20_P$Family)

# Correlation between the two
gw_20_B_bs_20_P_corr <- MCMCglmm(cbind(scale(Height_Growth),scale(BudSet2020_Plasticity)) ~ trait - 1, 
                                 # fixed=Region,
                                 random=~us(trait):Population,
                                 rcov=~us(trait):units, # or idh
                                 family=c("gaussian","gaussian"),
                                 prior=prior_corr, 
                                 # pedigree=Ped, 
                                 data=gw_20_B_bs_20_P,
                                 pr=TRUE, 
                                 verbose=TRUE,
                                 nitt=10000000, 
                                 burnin=10000, 
                                 thin=1000)

saveRDS(gw_20_B_bs_20_P_corr,"/home/Anoob/mydata/Anoob/MCMCglmm/plasticity_outputs/gw_20_B_bs_20_P_corr")