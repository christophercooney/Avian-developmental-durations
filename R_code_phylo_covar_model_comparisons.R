
## Phylogenetic covariance model comparisons

rm(list=ls())

library(ape)
library(phylolm)

setwd("~/Dropbox/Projects/Incubation/Comparative_analysis/")

# ---------- #

# get tree and data

phy <- read.tree("Outputs/MCC_tree.tre")
dat <- read.csv("Outputs/Comparative_dataset.csv"); rownames(dat) <- dat$binomial; head(dat)

# ---------- #

dat$dev_mod <- factor(dat$dev_mod, levels = c("precocial","semi-precocial","altricial"))
dat$parental_care_unibi <- factor(dat$parental_care_unibi, levels = c("bi","uni"))
dat$brood_parasite <- factor(dat$brood_parasite, levels = c("nonparasite","parasite"))
dat$nest_type <- factor(dat$nest_type, levels = c("cavity","closed","open","mixed"))
dat$habitat <- factor(dat$habitat, levels = c("none","low","med","high"))
dat$diet <- factor(dat$diet, levels = c("Omnivore","FruiNect","Invertebrate","PlantSeed","VertFishScav"))
dat$pelagic <- factor(dat$pelagic, levels = c("nonpelagic","pelagic"))
dat$nocturnal <- factor(dat$nocturnal, levels = c("diurnal","noctural"))
dat$migration <- factor(dat$migration, levels = c("sed","mig"))
dat$insularity <- factor(dat$insularity, levels = c("continental","insular"))

# ---------- #

# log_developmental_d

datx <- dat[,c("binomial","log_developmental_d", "log_adult_body_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.1.BM <- phylolm(log_developmental_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "BM"); summary(mx1.1.BM); mx1.1.BM[c("logLik","aic")]
mx1.1.OU <- phylolm(log_developmental_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.1.OU); mx1.1.OU[c("logLik","aic")]
mx1.1.LA <- phylolm(log_developmental_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "lambda"); summary(mx1.1.LA); mx1.1.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "log_egg_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.2.BM <- phylolm(log_developmental_d ~ log_egg_mass, data = datx, phy = phyx, model = "BM"); summary(mx1.2.BM); mx1.2.BM[c("logLik","aic")]
mx1.2.OU <- phylolm(log_developmental_d ~ log_egg_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.2.OU); mx1.2.OU[c("logLik","aic")]
mx1.2.LA <- phylolm(log_developmental_d ~ log_egg_mass, data = datx, phy = phyx, model = "lambda"); summary(mx1.2.LA); mx1.2.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "log_generation_length")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.3.BM <- phylolm(log_developmental_d ~ log_generation_length, data = datx, phy = phyx, model = "BM"); summary(mx1.3.BM); mx1.3.BM[c("logLik","aic")]
mx1.3.OU <- phylolm(log_developmental_d ~ log_generation_length, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.3.OU); mx1.3.OU[c("logLik","aic")]
mx1.3.LA <- phylolm(log_developmental_d ~ log_generation_length, data = datx, phy = phyx, model = "lambda"); summary(mx1.3.LA); mx1.3.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "log_clutch_size")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.4.BM <- phylolm(log_developmental_d ~ log_clutch_size, data = datx, phy = phyx, model = "BM"); summary(mx1.4.BM); mx1.4.BM[c("logLik","aic")]
mx1.4.OU <- phylolm(log_developmental_d ~ log_clutch_size, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.4.OU); mx1.4.OU[c("logLik","aic")]
mx1.4.LA <- phylolm(log_developmental_d ~ log_clutch_size, data = datx, phy = phyx, model = "lambda"); summary(mx1.4.LA); mx1.4.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "dev_mod")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.5.BM <- phylolm(log_developmental_d ~ dev_mod, data = datx, phy = phyx, model = "BM"); summary(mx1.5.BM); mx1.5.BM[c("logLik","aic")]
mx1.5.OU <- phylolm(log_developmental_d ~ dev_mod, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.5.OU); mx1.5.OU[c("logLik","aic")]
mx1.5.LA <- phylolm(log_developmental_d ~ dev_mod, data = datx, phy = phyx, model = "lambda"); summary(mx1.5.LA); mx1.5.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "parental_care_unibi")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.6.BM <- phylolm(log_developmental_d ~ parental_care_unibi, data = datx, phy = phyx, model = "BM"); summary(mx1.6.BM); mx1.6.BM[c("logLik","aic")]
mx1.6.OU <- phylolm(log_developmental_d ~ parental_care_unibi, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.6.OU); mx1.6.OU[c("logLik","aic")]
mx1.6.LA <- phylolm(log_developmental_d ~ parental_care_unibi, data = datx, phy = phyx, model = "lambda"); summary(mx1.6.LA); mx1.6.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "brood_parasite")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.7.BM <- phylolm(log_developmental_d ~ brood_parasite, data = datx, phy = phyx, model = "BM"); summary(mx1.7.BM); mx1.7.BM[c("logLik","aic")]
mx1.7.OU <- phylolm(log_developmental_d ~ brood_parasite, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.7.OU); mx1.7.OU[c("logLik","aic")]
mx1.7.LA <- phylolm(log_developmental_d ~ brood_parasite, data = datx, phy = phyx, model = "lambda"); summary(mx1.7.LA); mx1.7.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "nest_height_min")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.8.BM <- phylolm(log_developmental_d ~ nest_height_min, data = datx, phy = phyx, model = "BM"); summary(mx1.8.BM); mx1.8.BM[c("logLik","aic")]
mx1.8.OU <- phylolm(log_developmental_d ~ nest_height_min, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.8.OU); mx1.8.OU[c("logLik","aic")]
mx1.8.LA <- phylolm(log_developmental_d ~ nest_height_min, data = datx, phy = phyx, model = "lambda"); summary(mx1.8.LA); mx1.8.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "nest_type")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.9.BM <- phylolm(log_developmental_d ~ nest_type, data = datx, phy = phyx, model = "BM"); summary(mx1.9.BM); mx1.9.BM[c("logLik","aic")]
mx1.9.OU <- phylolm(log_developmental_d ~ nest_type, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.9.OU); mx1.9.OU[c("logLik","aic")]
mx1.9.LA <- phylolm(log_developmental_d ~ nest_type, data = datx, phy = phyx, model = "lambda"); summary(mx1.9.LA); mx1.9.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "habitat")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.10.BM <- phylolm(log_developmental_d ~ habitat, data = datx, phy = phyx, model = "BM"); summary(mx1.10.BM); mx1.10.BM[c("logLik","aic")]
mx1.10.OU <- phylolm(log_developmental_d ~ habitat, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.10.OU); mx1.10.OU[c("logLik","aic")]
mx1.10.LA <- phylolm(log_developmental_d ~ habitat, data = datx, phy = phyx, model = "lambda"); summary(mx1.10.LA); mx1.10.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "diet")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.11.BM <- phylolm(log_developmental_d ~ diet, data = datx, phy = phyx, model = "BM"); summary(mx1.11.BM); mx1.11.BM[c("logLik","aic")]
mx1.11.OU <- phylolm(log_developmental_d ~ diet, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.11.OU); mx1.11.OU[c("logLik","aic")]
mx1.11.LA <- phylolm(log_developmental_d ~ diet, data = datx, phy = phyx, model = "lambda"); summary(mx1.11.LA); mx1.11.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "pelagic")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.12.BM <- phylolm(log_developmental_d ~ pelagic, data = datx, phy = phyx, model = "BM"); summary(mx1.12.BM); mx1.12.BM[c("logLik","aic")]
mx1.12.OU <- phylolm(log_developmental_d ~ pelagic, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.12.OU); mx1.12.OU[c("logLik","aic")]
mx1.12.LA <- phylolm(log_developmental_d ~ pelagic, data = datx, phy = phyx, model = "lambda"); summary(mx1.12.LA); mx1.12.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "nocturnal")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.13.BM <- phylolm(log_developmental_d ~ nocturnal, data = datx, phy = phyx, model = "BM"); summary(mx1.13.BM); mx1.13.BM[c("logLik","aic")]
mx1.13.OU <- phylolm(log_developmental_d ~ nocturnal, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.13.OU); mx1.13.OU[c("logLik","aic")]
mx1.13.LA <- phylolm(log_developmental_d ~ nocturnal, data = datx, phy = phyx, model = "lambda"); summary(mx1.13.LA); mx1.13.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "migration")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.14.BM <- phylolm(log_developmental_d ~ migration, data = datx, phy = phyx, model = "BM"); summary(mx1.14.BM); mx1.14.BM[c("logLik","aic")]
mx1.14.OU <- phylolm(log_developmental_d ~ migration, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.14.OU); mx1.14.OU[c("logLik","aic")]
mx1.14.LA <- phylolm(log_developmental_d ~ migration, data = datx, phy = phyx, model = "lambda"); summary(mx1.14.LA); mx1.14.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "temperature")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.15.BM <- phylolm(log_developmental_d ~ temperature, data = datx, phy = phyx, model = "BM"); summary(mx1.15.BM); mx1.15.BM[c("logLik","aic")]
mx1.15.OU <- phylolm(log_developmental_d ~ temperature, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.15.OU); mx1.15.OU[c("logLik","aic")]
mx1.15.LA <- phylolm(log_developmental_d ~ temperature, data = datx, phy = phyx, model = "lambda"); summary(mx1.15.LA); mx1.15.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "precipitation")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.16.BM <- phylolm(log_developmental_d ~ precipitation, data = datx, phy = phyx, model = "BM"); summary(mx1.16.BM); mx1.16.BM[c("logLik","aic")]
mx1.16.OU <- phylolm(log_developmental_d ~ precipitation, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.16.OU); mx1.16.OU[c("logLik","aic")]
mx1.16.LA <- phylolm(log_developmental_d ~ precipitation, data = datx, phy = phyx, model = "lambda"); summary(mx1.16.LA); mx1.16.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.17.BM <- phylolm(log_developmental_d ~ latitude, data = datx, phy = phyx, model = "BM"); summary(mx1.17.BM); mx1.17.BM[c("logLik","aic")]
mx1.17.OU <- phylolm(log_developmental_d ~ latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.17.OU); mx1.17.OU[c("logLik","aic")]
mx1.17.LA <- phylolm(log_developmental_d ~ latitude, data = datx, phy = phyx, model = "lambda"); summary(mx1.17.LA); mx1.17.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.18.BM <- phylolm(log_developmental_d ~ insularity, data = datx, phy = phyx, model = "BM"); summary(mx1.18.BM); mx1.18.BM[c("logLik","aic")]
mx1.18.OU <- phylolm(log_developmental_d ~ insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.18.OU); mx1.18.OU[c("logLik","aic")]
mx1.18.LA <- phylolm(log_developmental_d ~ insularity, data = datx, phy = phyx, model = "lambda"); summary(mx1.18.LA); mx1.18.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "log_adult_body_mass","log_generation_length","log_clutch_size","parental_care_unibi","brood_parasite","nest_height_min","habitat","diet","pelagic","migration","temperature","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.19.BM <- phylolm(log_developmental_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + nest_height_min + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx1.19.BM); mx1.19.BM[c("logLik","aic")]
mx1.19.OU <- phylolm(log_developmental_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + nest_height_min + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.19.OU); mx1.19.OU[c("logLik","aic")]
mx1.19.LA <- phylolm(log_developmental_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + nest_height_min + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx1.19.LA); mx1.19.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_developmental_d", "log_egg_mass","log_generation_length","log_clutch_size","parental_care_unibi","brood_parasite","habitat","diet","pelagic","migration","temperature","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx1.20.BM <- phylolm(log_developmental_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx1.20.BM); mx1.20.BM[c("logLik","aic")]
mx1.20.OU <- phylolm(log_developmental_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx1.20.OU); mx1.20.OU[c("logLik","aic")]
mx1.20.LA <- phylolm(log_developmental_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + brood_parasite + habitat + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx1.20.LA); mx1.20.LA[c("logLik","aic")]

# prop_incubation_sqrt

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_adult_body_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.1.BM <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = datx, phy = phyx, model = "BM"); summary(mx2.1.BM); mx2.1.BM[c("logLik","aic")]
mx2.1.OU <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.1.OU); mx2.1.OU[c("logLik","aic")]
mx2.1.LA <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = datx, phy = phyx, model = "lambda"); summary(mx2.1.LA); mx2.1.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_egg_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.2.BM <- phylolm(prop_incubation_sqrt ~ log_egg_mass, data = datx, phy = phyx, model = "BM"); summary(mx2.2.BM); mx2.2.BM[c("logLik","aic")]
mx2.2.OU <- phylolm(prop_incubation_sqrt ~ log_egg_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.2.OU); mx2.2.OU[c("logLik","aic")]
mx2.2.LA <- phylolm(prop_incubation_sqrt ~ log_egg_mass, data = datx, phy = phyx, model = "lambda"); summary(mx2.2.LA); mx2.2.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_generation_length")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.3.BM <- phylolm(prop_incubation_sqrt ~ log_generation_length, data = datx, phy = phyx, model = "BM"); summary(mx2.3.BM); mx2.3.BM[c("logLik","aic")]
mx2.3.OU <- phylolm(prop_incubation_sqrt ~ log_generation_length, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.3.OU); mx2.3.OU[c("logLik","aic")]
mx2.3.LA <- phylolm(prop_incubation_sqrt ~ log_generation_length, data = datx, phy = phyx, model = "lambda"); summary(mx2.3.LA); mx2.3.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_clutch_size")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.4.BM <- phylolm(prop_incubation_sqrt ~ log_clutch_size, data = datx, phy = phyx, model = "BM"); summary(mx2.4.BM); mx2.4.BM[c("logLik","aic")]
mx2.4.OU <- phylolm(prop_incubation_sqrt ~ log_clutch_size, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.4.OU); mx2.4.OU[c("logLik","aic")]
mx2.4.LA <- phylolm(prop_incubation_sqrt ~ log_clutch_size, data = datx, phy = phyx, model = "lambda"); summary(mx2.4.LA); mx2.4.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "dev_mod")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.5.BM <- phylolm(prop_incubation_sqrt ~ dev_mod, data = datx, phy = phyx, model = "BM"); summary(mx2.5.BM); mx2.5.BM[c("logLik","aic")]
mx2.5.OU <- phylolm(prop_incubation_sqrt ~ dev_mod, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.5.OU); mx2.5.OU[c("logLik","aic")]
mx2.5.LA <- phylolm(prop_incubation_sqrt ~ dev_mod, data = datx, phy = phyx, model = "lambda"); summary(mx2.5.LA); mx2.5.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "parental_care_unibi")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.6.BM <- phylolm(prop_incubation_sqrt ~ parental_care_unibi, data = datx, phy = phyx, model = "BM"); summary(mx2.6.BM); mx2.6.BM[c("logLik","aic")]
mx2.6.OU <- phylolm(prop_incubation_sqrt ~ parental_care_unibi, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.6.OU); mx2.6.OU[c("logLik","aic")]
mx2.6.LA <- phylolm(prop_incubation_sqrt ~ parental_care_unibi, data = datx, phy = phyx, model = "lambda"); summary(mx2.6.LA); mx2.6.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "brood_parasite")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.7.BM <- phylolm(prop_incubation_sqrt ~ brood_parasite, data = datx, phy = phyx, model = "BM"); summary(mx2.7.BM); mx2.7.BM[c("logLik","aic")]
mx2.7.OU <- phylolm(prop_incubation_sqrt ~ brood_parasite, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.7.OU); mx2.7.OU[c("logLik","aic")]
mx2.7.LA <- phylolm(prop_incubation_sqrt ~ brood_parasite, data = datx, phy = phyx, model = "lambda"); summary(mx2.7.LA); mx2.7.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "nest_height_min")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.8.BM <- phylolm(prop_incubation_sqrt ~ nest_height_min, data = datx, phy = phyx, model = "BM"); summary(mx2.8.BM); mx2.8.BM[c("logLik","aic")]
mx2.8.OU <- phylolm(prop_incubation_sqrt ~ nest_height_min, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.8.OU); mx2.8.OU[c("logLik","aic")]
mx2.8.LA <- phylolm(prop_incubation_sqrt ~ nest_height_min, data = datx, phy = phyx, model = "lambda"); summary(mx2.8.LA); mx2.8.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "nest_type")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.9.BM <- phylolm(prop_incubation_sqrt ~ nest_type, data = datx, phy = phyx, model = "BM"); summary(mx2.9.BM); mx2.9.BM[c("logLik","aic")]
mx2.9.OU <- phylolm(prop_incubation_sqrt ~ nest_type, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.9.OU); mx2.9.OU[c("logLik","aic")]
mx2.9.LA <- phylolm(prop_incubation_sqrt ~ nest_type, data = datx, phy = phyx, model = "lambda"); summary(mx2.9.LA); mx2.9.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "habitat")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.10.BM <- phylolm(prop_incubation_sqrt ~ habitat, data = datx, phy = phyx, model = "BM"); summary(mx2.10.BM); mx2.10.BM[c("logLik","aic")]
mx2.10.OU <- phylolm(prop_incubation_sqrt ~ habitat, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.10.OU); mx2.10.OU[c("logLik","aic")]
mx2.10.LA <- phylolm(prop_incubation_sqrt ~ habitat, data = datx, phy = phyx, model = "lambda"); summary(mx2.10.LA); mx2.10.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "diet")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.11.BM <- phylolm(prop_incubation_sqrt ~ diet, data = datx, phy = phyx, model = "BM"); summary(mx2.11.BM); mx2.11.BM[c("logLik","aic")]
mx2.11.OU <- phylolm(prop_incubation_sqrt ~ diet, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.11.OU); mx2.11.OU[c("logLik","aic")]
mx2.11.LA <- phylolm(prop_incubation_sqrt ~ diet, data = datx, phy = phyx, model = "lambda"); summary(mx2.11.LA); mx2.11.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "pelagic")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.12.BM <- phylolm(prop_incubation_sqrt ~ pelagic, data = datx, phy = phyx, model = "BM"); summary(mx2.12.BM); mx2.12.BM[c("logLik","aic")]
mx2.12.OU <- phylolm(prop_incubation_sqrt ~ pelagic, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.12.OU); mx2.12.OU[c("logLik","aic")]
mx2.12.LA <- phylolm(prop_incubation_sqrt ~ pelagic, data = datx, phy = phyx, model = "lambda"); summary(mx2.12.LA); mx2.12.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "nocturnal")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.13.BM <- phylolm(prop_incubation_sqrt ~ nocturnal, data = datx, phy = phyx, model = "BM"); summary(mx2.13.BM); mx2.13.BM[c("logLik","aic")]
mx2.13.OU <- phylolm(prop_incubation_sqrt ~ nocturnal, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.13.OU); mx2.13.OU[c("logLik","aic")]
mx2.13.LA <- phylolm(prop_incubation_sqrt ~ nocturnal, data = datx, phy = phyx, model = "lambda"); summary(mx2.13.LA); mx2.13.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "migration")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.14.BM <- phylolm(prop_incubation_sqrt ~ migration, data = datx, phy = phyx, model = "BM"); summary(mx2.14.BM); mx2.14.BM[c("logLik","aic")]
mx2.14.OU <- phylolm(prop_incubation_sqrt ~ migration, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.14.OU); mx2.14.OU[c("logLik","aic")]
mx2.14.LA <- phylolm(prop_incubation_sqrt ~ migration, data = datx, phy = phyx, model = "lambda"); summary(mx2.14.LA); mx2.14.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "temperature")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.15.BM <- phylolm(prop_incubation_sqrt ~ temperature, data = datx, phy = phyx, model = "BM"); summary(mx2.15.BM); mx2.15.BM[c("logLik","aic")]
mx2.15.OU <- phylolm(prop_incubation_sqrt ~ temperature, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.15.OU); mx2.15.OU[c("logLik","aic")]
mx2.15.LA <- phylolm(prop_incubation_sqrt ~ temperature, data = datx, phy = phyx, model = "lambda"); summary(mx2.15.LA); mx2.15.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "precipitation")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.16.BM <- phylolm(prop_incubation_sqrt ~ precipitation, data = datx, phy = phyx, model = "BM"); summary(mx2.16.BM); mx2.16.BM[c("logLik","aic")]
mx2.16.OU <- phylolm(prop_incubation_sqrt ~ precipitation, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.16.OU); mx2.16.OU[c("logLik","aic")]
mx2.16.LA <- phylolm(prop_incubation_sqrt ~ precipitation, data = datx, phy = phyx, model = "lambda"); summary(mx2.16.LA); mx2.16.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.17.BM <- phylolm(prop_incubation_sqrt ~ latitude, data = datx, phy = phyx, model = "BM"); summary(mx2.17.BM); mx2.17.BM[c("logLik","aic")]
mx2.17.OU <- phylolm(prop_incubation_sqrt ~ latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.17.OU); mx2.17.OU[c("logLik","aic")]
mx2.17.LA <- phylolm(prop_incubation_sqrt ~ latitude, data = datx, phy = phyx, model = "lambda"); summary(mx2.17.LA); mx2.17.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.18.BM <- phylolm(prop_incubation_sqrt ~ insularity, data = datx, phy = phyx, model = "BM"); summary(mx2.18.BM); mx2.18.BM[c("logLik","aic")]
mx2.18.OU <- phylolm(prop_incubation_sqrt ~ insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.18.OU); mx2.18.OU[c("logLik","aic")]
mx2.18.LA <- phylolm(prop_incubation_sqrt ~ insularity, data = datx, phy = phyx, model = "lambda"); summary(mx2.18.LA); mx2.18.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_adult_body_mass","log_generation_length","log_clutch_size","parental_care_unibi","nest_height_min","diet","pelagic","nocturnal","migration","latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.19.BM <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "BM"); summary(mx2.19.BM); mx2.19.BM[c("logLik","aic")]
mx2.19.OU <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.19.OU); mx2.19.OU[c("logLik","aic")]
mx2.19.LA <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "lambda"); summary(mx2.19.LA); mx2.19.LA[c("logLik","aic")]

datx <- dat[,c("binomial","prop_incubation_sqrt", "log_egg_mass","log_generation_length","log_clutch_size","parental_care_unibi","diet","pelagic","nocturnal","migration","latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx2.20.BM <- phylolm(prop_incubation_sqrt ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "BM"); summary(mx2.20.BM); mx2.20.BM[c("logLik","aic")]
mx2.20.OU <- phylolm(prop_incubation_sqrt ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx2.20.OU); mx2.20.OU[c("logLik","aic")]
mx2.20.LA <- phylolm(prop_incubation_sqrt ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + nocturnal + migration + latitude, data = datx, phy = phyx, model = "lambda"); summary(mx2.20.LA); mx2.20.LA[c("logLik","aic")]

# log_incubation_d

datx <- dat[,c("binomial","log_incubation_d", "log_adult_body_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.1.BM <- phylolm(log_incubation_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "BM"); summary(mx3.1.BM); mx3.1.BM[c("logLik","aic")]
mx3.1.OU <- phylolm(log_incubation_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.1.OU); mx3.1.OU[c("logLik","aic")]
mx3.1.LA <- phylolm(log_incubation_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "lambda"); summary(mx3.1.LA); mx3.1.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "log_egg_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.2.BM <- phylolm(log_incubation_d ~ log_egg_mass, data = datx, phy = phyx, model = "BM"); summary(mx3.2.BM); mx3.2.BM[c("logLik","aic")]
mx3.2.OU <- phylolm(log_incubation_d ~ log_egg_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.2.OU); mx3.2.OU[c("logLik","aic")]
mx3.2.LA <- phylolm(log_incubation_d ~ log_egg_mass, data = datx, phy = phyx, model = "lambda"); summary(mx3.2.LA); mx3.2.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "log_generation_length")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.3.BM <- phylolm(log_incubation_d ~ log_generation_length, data = datx, phy = phyx, model = "BM"); summary(mx3.3.BM); mx3.3.BM[c("logLik","aic")]
mx3.3.OU <- phylolm(log_incubation_d ~ log_generation_length, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.3.OU); mx3.3.OU[c("logLik","aic")]
mx3.3.LA <- phylolm(log_incubation_d ~ log_generation_length, data = datx, phy = phyx, model = "lambda"); summary(mx3.3.LA); mx3.3.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "log_clutch_size")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.4.BM <- phylolm(log_incubation_d ~ log_clutch_size, data = datx, phy = phyx, model = "BM"); summary(mx3.4.BM); mx3.4.BM[c("logLik","aic")]
mx3.4.OU <- phylolm(log_incubation_d ~ log_clutch_size, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.4.OU); mx3.4.OU[c("logLik","aic")]
mx3.4.LA <- phylolm(log_incubation_d ~ log_clutch_size, data = datx, phy = phyx, model = "lambda"); summary(mx3.4.LA); mx3.4.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "dev_mod")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.5.BM <- phylolm(log_incubation_d ~ dev_mod, data = datx, phy = phyx, model = "BM"); summary(mx3.5.BM); mx3.5.BM[c("logLik","aic")]
mx3.5.OU <- phylolm(log_incubation_d ~ dev_mod, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.5.OU); mx3.5.OU[c("logLik","aic")]
mx3.5.LA <- phylolm(log_incubation_d ~ dev_mod, data = datx, phy = phyx, model = "lambda"); summary(mx3.5.LA); mx3.5.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "parental_care_unibi")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.6.BM <- phylolm(log_incubation_d ~ parental_care_unibi, data = datx, phy = phyx, model = "BM"); summary(mx3.6.BM); mx3.6.BM[c("logLik","aic")]
mx3.6.OU <- phylolm(log_incubation_d ~ parental_care_unibi, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.6.OU); mx3.6.OU[c("logLik","aic")]
mx3.6.LA <- phylolm(log_incubation_d ~ parental_care_unibi, data = datx, phy = phyx, model = "lambda"); summary(mx3.6.LA); mx3.6.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "brood_parasite")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.7.BM <- phylolm(log_incubation_d ~ brood_parasite, data = datx, phy = phyx, model = "BM"); summary(mx3.7.BM); mx3.7.BM[c("logLik","aic")]
mx3.7.OU <- phylolm(log_incubation_d ~ brood_parasite, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.7.OU); mx3.7.OU[c("logLik","aic")]
mx3.7.LA <- phylolm(log_incubation_d ~ brood_parasite, data = datx, phy = phyx, model = "lambda"); summary(mx3.7.LA); mx3.7.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "nest_height_min")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.8.BM <- phylolm(log_incubation_d ~ nest_height_min, data = datx, phy = phyx, model = "BM"); summary(mx3.8.BM); mx3.8.BM[c("logLik","aic")]
mx3.8.OU <- phylolm(log_incubation_d ~ nest_height_min, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.8.OU); mx3.8.OU[c("logLik","aic")]
mx3.8.LA <- phylolm(log_incubation_d ~ nest_height_min, data = datx, phy = phyx, model = "lambda"); summary(mx3.8.LA); mx3.8.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "nest_type")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.9.BM <- phylolm(log_incubation_d ~ nest_type, data = datx, phy = phyx, model = "BM"); summary(mx3.9.BM); mx3.9.BM[c("logLik","aic")]
mx3.9.OU <- phylolm(log_incubation_d ~ nest_type, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.9.OU); mx3.9.OU[c("logLik","aic")]
mx3.9.LA <- phylolm(log_incubation_d ~ nest_type, data = datx, phy = phyx, model = "lambda"); summary(mx3.9.LA); mx3.9.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "habitat")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.10.BM <- phylolm(log_incubation_d ~ habitat, data = datx, phy = phyx, model = "BM"); summary(mx3.10.BM); mx3.10.BM[c("logLik","aic")]
mx3.10.OU <- phylolm(log_incubation_d ~ habitat, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.10.OU); mx3.10.OU[c("logLik","aic")]
mx3.10.LA <- phylolm(log_incubation_d ~ habitat, data = datx, phy = phyx, model = "lambda"); summary(mx3.10.LA); mx3.10.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "diet")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.11.BM <- phylolm(log_incubation_d ~ diet, data = datx, phy = phyx, model = "BM"); summary(mx3.11.BM); mx3.11.BM[c("logLik","aic")]
mx3.11.OU <- phylolm(log_incubation_d ~ diet, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.11.OU); mx3.11.OU[c("logLik","aic")]
mx3.11.LA <- phylolm(log_incubation_d ~ diet, data = datx, phy = phyx, model = "lambda"); summary(mx3.11.LA); mx3.11.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "pelagic")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.12.BM <- phylolm(log_incubation_d ~ pelagic, data = datx, phy = phyx, model = "BM"); summary(mx3.12.BM); mx3.12.BM[c("logLik","aic")]
mx3.12.OU <- phylolm(log_incubation_d ~ pelagic, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.12.OU); mx3.12.OU[c("logLik","aic")]
mx3.12.LA <- phylolm(log_incubation_d ~ pelagic, data = datx, phy = phyx, model = "lambda"); summary(mx3.12.LA); mx3.12.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "nocturnal")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.13.BM <- phylolm(log_incubation_d ~ nocturnal, data = datx, phy = phyx, model = "BM"); summary(mx3.13.BM); mx3.13.BM[c("logLik","aic")]
mx3.13.OU <- phylolm(log_incubation_d ~ nocturnal, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.13.OU); mx3.13.OU[c("logLik","aic")]
mx3.13.LA <- phylolm(log_incubation_d ~ nocturnal, data = datx, phy = phyx, model = "lambda"); summary(mx3.13.LA); mx3.13.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "migration")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.14.BM <- phylolm(log_incubation_d ~ migration, data = datx, phy = phyx, model = "BM"); summary(mx3.14.BM); mx3.14.BM[c("logLik","aic")]
mx3.14.OU <- phylolm(log_incubation_d ~ migration, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.14.OU); mx3.14.OU[c("logLik","aic")]
mx3.14.LA <- phylolm(log_incubation_d ~ migration, data = datx, phy = phyx, model = "lambda"); summary(mx3.14.LA); mx3.14.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "temperature")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.15.BM <- phylolm(log_incubation_d ~ temperature, data = datx, phy = phyx, model = "BM"); summary(mx3.15.BM); mx3.15.BM[c("logLik","aic")]
mx3.15.OU <- phylolm(log_incubation_d ~ temperature, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.15.OU); mx3.15.OU[c("logLik","aic")]
mx3.15.LA <- phylolm(log_incubation_d ~ temperature, data = datx, phy = phyx, model = "lambda"); summary(mx3.15.LA); mx3.15.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "precipitation")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.16.BM <- phylolm(log_incubation_d ~ precipitation, data = datx, phy = phyx, model = "BM"); summary(mx3.16.BM); mx3.16.BM[c("logLik","aic")]
mx3.16.OU <- phylolm(log_incubation_d ~ precipitation, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.16.OU); mx3.16.OU[c("logLik","aic")]
mx3.16.LA <- phylolm(log_incubation_d ~ precipitation, data = datx, phy = phyx, model = "lambda"); summary(mx3.16.LA); mx3.16.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.17.BM <- phylolm(log_incubation_d ~ latitude, data = datx, phy = phyx, model = "BM"); summary(mx3.17.BM); mx3.17.BM[c("logLik","aic")]
mx3.17.OU <- phylolm(log_incubation_d ~ latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.17.OU); mx3.17.OU[c("logLik","aic")]
mx3.17.LA <- phylolm(log_incubation_d ~ latitude, data = datx, phy = phyx, model = "lambda"); summary(mx3.17.LA); mx3.17.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.18.BM <- phylolm(log_incubation_d ~ insularity, data = datx, phy = phyx, model = "BM"); summary(mx3.18.BM); mx3.18.BM[c("logLik","aic")]
mx3.18.OU <- phylolm(log_incubation_d ~ insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.18.OU); mx3.18.OU[c("logLik","aic")]
mx3.18.LA <- phylolm(log_incubation_d ~ insularity, data = datx, phy = phyx, model = "lambda"); summary(mx3.18.LA); mx3.18.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "log_adult_body_mass","log_generation_length","log_clutch_size","brood_parasite","nest_height_min","habitat","diet","pelagic","nocturnal","migration","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.19.BM <- phylolm(log_incubation_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + brood_parasite + nest_height_min + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx3.19.BM); mx3.19.BM[c("logLik","aic")]
mx3.19.OU <- phylolm(log_incubation_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + brood_parasite + nest_height_min + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.19.OU); mx3.19.OU[c("logLik","aic")]
mx3.19.LA <- phylolm(log_incubation_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + brood_parasite + nest_height_min + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx3.19.LA); mx3.19.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_incubation_d", "log_egg_mass","log_generation_length","log_clutch_size","brood_parasite","habitat","diet","pelagic","nocturnal","migration","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx3.20.BM <- phylolm(log_incubation_d ~ log_egg_mass + log_generation_length + log_clutch_size + brood_parasite + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx3.20.BM); mx3.20.BM[c("logLik","aic")]
mx3.20.OU <- phylolm(log_incubation_d ~ log_egg_mass + log_generation_length + log_clutch_size + brood_parasite + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx3.20.OU); mx3.20.OU[c("logLik","aic")]
mx3.20.LA <- phylolm(log_incubation_d ~ log_egg_mass + log_generation_length + log_clutch_size + brood_parasite + habitat + diet + pelagic + nocturnal + migration + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx3.20.LA); mx3.20.LA[c("logLik","aic")]

# log_fledging_d

datx <- dat[,c("binomial","log_fledging_d", "log_adult_body_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.1.BM <- phylolm(log_fledging_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "BM"); summary(mx4.1.BM); mx4.1.BM[c("logLik","aic")]
mx4.1.OU <- phylolm(log_fledging_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.1.OU); mx4.1.OU[c("logLik","aic")]
mx4.1.LA <- phylolm(log_fledging_d ~ log_adult_body_mass, data = datx, phy = phyx, model = "lambda"); summary(mx4.1.LA); mx4.1.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "log_egg_mass")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.2.BM <- phylolm(log_fledging_d ~ log_egg_mass, data = datx, phy = phyx, model = "BM"); summary(mx4.2.BM); mx4.2.BM[c("logLik","aic")]
mx4.2.OU <- phylolm(log_fledging_d ~ log_egg_mass, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.2.OU); mx4.2.OU[c("logLik","aic")]
mx4.2.LA <- phylolm(log_fledging_d ~ log_egg_mass, data = datx, phy = phyx, model = "lambda"); summary(mx4.2.LA); mx4.2.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "log_generation_length")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.3.BM <- phylolm(log_fledging_d ~ log_generation_length, data = datx, phy = phyx, model = "BM"); summary(mx4.3.BM); mx4.3.BM[c("logLik","aic")]
mx4.3.OU <- phylolm(log_fledging_d ~ log_generation_length, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.3.OU); mx4.3.OU[c("logLik","aic")]
mx4.3.LA <- phylolm(log_fledging_d ~ log_generation_length, data = datx, phy = phyx, model = "lambda"); summary(mx4.3.LA); mx4.3.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "log_clutch_size")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.4.BM <- phylolm(log_fledging_d ~ log_clutch_size, data = datx, phy = phyx, model = "BM"); summary(mx4.4.BM); mx4.4.BM[c("logLik","aic")]
mx4.4.OU <- phylolm(log_fledging_d ~ log_clutch_size, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.4.OU); mx4.4.OU[c("logLik","aic")]
mx4.4.LA <- phylolm(log_fledging_d ~ log_clutch_size, data = datx, phy = phyx, model = "lambda"); summary(mx4.4.LA); mx4.4.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "dev_mod")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.5.BM <- phylolm(log_fledging_d ~ dev_mod, data = datx, phy = phyx, model = "BM"); summary(mx4.5.BM); mx4.5.BM[c("logLik","aic")]
mx4.5.OU <- phylolm(log_fledging_d ~ dev_mod, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.5.OU); mx4.5.OU[c("logLik","aic")]
mx4.5.LA <- phylolm(log_fledging_d ~ dev_mod, data = datx, phy = phyx, model = "lambda"); summary(mx4.5.LA); mx4.5.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "parental_care_unibi")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.6.BM <- phylolm(log_fledging_d ~ parental_care_unibi, data = datx, phy = phyx, model = "BM"); summary(mx4.6.BM); mx4.6.BM[c("logLik","aic")]
mx4.6.OU <- phylolm(log_fledging_d ~ parental_care_unibi, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.6.OU); mx4.6.OU[c("logLik","aic")]
mx4.6.LA <- phylolm(log_fledging_d ~ parental_care_unibi, data = datx, phy = phyx, model = "lambda"); summary(mx4.6.LA); mx4.6.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "brood_parasite")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.7.BM <- phylolm(log_fledging_d ~ brood_parasite, data = datx, phy = phyx, model = "BM"); summary(mx4.7.BM); mx4.7.BM[c("logLik","aic")]
mx4.7.OU <- phylolm(log_fledging_d ~ brood_parasite, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.7.OU); mx4.7.OU[c("logLik","aic")]
mx4.7.LA <- phylolm(log_fledging_d ~ brood_parasite, data = datx, phy = phyx, model = "lambda"); summary(mx4.7.LA); mx4.7.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "nest_height_min")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.8.BM <- phylolm(log_fledging_d ~ nest_height_min, data = datx, phy = phyx, model = "BM"); summary(mx4.8.BM); mx4.8.BM[c("logLik","aic")]
mx4.8.OU <- phylolm(log_fledging_d ~ nest_height_min, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.8.OU); mx4.8.OU[c("logLik","aic")]
mx4.8.LA <- phylolm(log_fledging_d ~ nest_height_min, data = datx, phy = phyx, model = "lambda"); summary(mx4.8.LA); mx4.8.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "nest_type")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.9.BM <- phylolm(log_fledging_d ~ nest_type, data = datx, phy = phyx, model = "BM"); summary(mx4.9.BM); mx4.9.BM[c("logLik","aic")]
mx4.9.OU <- phylolm(log_fledging_d ~ nest_type, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.9.OU); mx4.9.OU[c("logLik","aic")]
mx4.9.LA <- phylolm(log_fledging_d ~ nest_type, data = datx, phy = phyx, model = "lambda"); summary(mx4.9.LA); mx4.9.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "habitat")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.10.BM <- phylolm(log_fledging_d ~ habitat, data = datx, phy = phyx, model = "BM"); summary(mx4.10.BM); mx4.10.BM[c("logLik","aic")]
mx4.10.OU <- phylolm(log_fledging_d ~ habitat, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.10.OU); mx4.10.OU[c("logLik","aic")]
mx4.10.LA <- phylolm(log_fledging_d ~ habitat, data = datx, phy = phyx, model = "lambda"); summary(mx4.10.LA); mx4.10.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "diet")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.11.BM <- phylolm(log_fledging_d ~ diet, data = datx, phy = phyx, model = "BM"); summary(mx4.11.BM); mx4.11.BM[c("logLik","aic")]
mx4.11.OU <- phylolm(log_fledging_d ~ diet, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.11.OU); mx4.11.OU[c("logLik","aic")]
mx4.11.LA <- phylolm(log_fledging_d ~ diet, data = datx, phy = phyx, model = "lambda"); summary(mx4.11.LA); mx4.11.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "pelagic")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.12.BM <- phylolm(log_fledging_d ~ pelagic, data = datx, phy = phyx, model = "BM"); summary(mx4.12.BM); mx4.12.BM[c("logLik","aic")]
mx4.12.OU <- phylolm(log_fledging_d ~ pelagic, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.12.OU); mx4.12.OU[c("logLik","aic")]
mx4.12.LA <- phylolm(log_fledging_d ~ pelagic, data = datx, phy = phyx, model = "lambda"); summary(mx4.12.LA); mx4.12.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "nocturnal")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.13.BM <- phylolm(log_fledging_d ~ nocturnal, data = datx, phy = phyx, model = "BM"); summary(mx4.13.BM); mx4.13.BM[c("logLik","aic")]
mx4.13.OU <- phylolm(log_fledging_d ~ nocturnal, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.13.OU); mx4.13.OU[c("logLik","aic")]
mx4.13.LA <- phylolm(log_fledging_d ~ nocturnal, data = datx, phy = phyx, model = "lambda"); summary(mx4.13.LA); mx4.13.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "migration")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.14.BM <- phylolm(log_fledging_d ~ migration, data = datx, phy = phyx, model = "BM"); summary(mx4.14.BM); mx4.14.BM[c("logLik","aic")]
mx4.14.OU <- phylolm(log_fledging_d ~ migration, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.14.OU); mx4.14.OU[c("logLik","aic")]
mx4.14.LA <- phylolm(log_fledging_d ~ migration, data = datx, phy = phyx, model = "lambda"); summary(mx4.14.LA); mx4.14.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "temperature")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.15.BM <- phylolm(log_fledging_d ~ temperature, data = datx, phy = phyx, model = "BM"); summary(mx4.15.BM); mx4.15.BM[c("logLik","aic")]
mx4.15.OU <- phylolm(log_fledging_d ~ temperature, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.15.OU); mx4.15.OU[c("logLik","aic")]
mx4.15.LA <- phylolm(log_fledging_d ~ temperature, data = datx, phy = phyx, model = "lambda"); summary(mx4.15.LA); mx4.15.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "precipitation")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.16.BM <- phylolm(log_fledging_d ~ precipitation, data = datx, phy = phyx, model = "BM"); summary(mx4.16.BM); mx4.16.BM[c("logLik","aic")]
mx4.16.OU <- phylolm(log_fledging_d ~ precipitation, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.16.OU); mx4.16.OU[c("logLik","aic")]
mx4.16.LA <- phylolm(log_fledging_d ~ precipitation, data = datx, phy = phyx, model = "lambda"); summary(mx4.16.LA); mx4.16.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "latitude")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.17.BM <- phylolm(log_fledging_d ~ latitude, data = datx, phy = phyx, model = "BM"); summary(mx4.17.BM); mx4.17.BM[c("logLik","aic")]
mx4.17.OU <- phylolm(log_fledging_d ~ latitude, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.17.OU); mx4.17.OU[c("logLik","aic")]
mx4.17.LA <- phylolm(log_fledging_d ~ latitude, data = datx, phy = phyx, model = "lambda"); summary(mx4.17.LA); mx4.17.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.18.BM <- phylolm(log_fledging_d ~ insularity, data = datx, phy = phyx, model = "BM"); summary(mx4.18.BM); mx4.18.BM[c("logLik","aic")]
mx4.18.OU <- phylolm(log_fledging_d ~ insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.18.OU); mx4.18.OU[c("logLik","aic")]
mx4.18.LA <- phylolm(log_fledging_d ~ insularity, data = datx, phy = phyx, model = "lambda"); summary(mx4.18.LA); mx4.18.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "log_adult_body_mass","log_generation_length","log_clutch_size","parental_care_unibi","nest_height_min","diet","pelagic","migration","temperature","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.19.BM <- phylolm(log_fledging_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx4.19.BM); mx4.19.BM[c("logLik","aic")]
mx4.19.OU <- phylolm(log_fledging_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.19.OU); mx4.19.OU[c("logLik","aic")]
mx4.19.LA <- phylolm(log_fledging_d ~ log_adult_body_mass + log_generation_length + log_clutch_size + parental_care_unibi + nest_height_min + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx4.19.LA); mx4.19.LA[c("logLik","aic")]

datx <- dat[,c("binomial","log_fledging_d", "log_egg_mass","log_generation_length","log_clutch_size","parental_care_unibi","diet","pelagic","migration","temperature","precipitation","latitude","insularity")]
datx <- datx[complete.cases(datx),]
phyx <- drop.tip(phy, setdiff(phy$tip.label, datx$binomial))
mx4.20.BM <- phylolm(log_fledging_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "BM"); summary(mx4.20.BM); mx4.20.BM[c("logLik","aic")]
mx4.20.OU <- phylolm(log_fledging_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "OUfixedRoot"); summary(mx4.20.OU); mx4.20.OU[c("logLik","aic")]
mx4.20.LA <- phylolm(log_fledging_d ~ log_egg_mass + log_generation_length + log_clutch_size + parental_care_unibi + diet + pelagic + migration + temperature + precipitation + latitude + insularity, data = datx, phy = phyx, model = "lambda"); summary(mx4.20.LA); mx4.20.LA[c("logLik","aic")]

##########################################
