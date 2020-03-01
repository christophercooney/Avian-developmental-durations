
## Allometry analyses

rm(list=ls())

library(ape)
library(caper)
library(phylolm)
library(rr2)

setwd("~/Dropbox/Projects/Incubation/Comparative_analysis/")

# ---------- #

# get tree and data

phy <- read.tree("Outputs/MCC_tree.tre")
dat <- read.csv("Outputs/Comparative_dataset.csv", stringsAsFactors = F); rownames(dat) <- dat$binomial; head(dat)

# ---------- #

orders.egg <- sort(table(dat$Taxon_subgroup[!is.na(dat$log_egg_mass)]))
orders.body <- sort(table(dat$Taxon_subgroup[!is.na(dat$log_adult_body_mass)]))
orders <- orders.body[orders.body >= 20]
background <- orders.body[!orders.body %in% orders]
dat$Taxon_subgroup[dat$Taxon_subgroup %in% names(background)] <- "Background"
orders <- c(orders, sum(background))
names(orders) <- c(names(orders)[-length(orders)], "Background")

# ---------- #


##### PHYLOGENETIC SIGNAL WITH CONFIDENCE INTERVALS VIA BOOTSTRAPPING

# Developmental duration

fx <- dat[,c("binomial","log_developmental_d")]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
#library(future); plan(multiprocess)
mx <- phylolm(log_developmental_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", boot = 1000); print(summary(mx))


# Incubation fraction

fx <- dat[,c("binomial","prop_incubation_sqrt")]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
#library(future); plan(multiprocess)
mx <- phylolm(prop_incubation_sqrt ~ 1, data = cx$data, phy = cx$phy, model = "lambda", boot = 1000); print(summary(mx))


# ---------- #


### DEVELOPMENTAL DURATION

# adult body mass

# single predictor (all data)
mx <- NA
fx <- dat[,c("binomial","log_developmental_d","log_adult_body_mass")]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
mx <- phylolm(log_developmental_d ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda"); #print(summary(mx))
mlpar <- data.frame(lambda=mx$optpar)
mx.r <- phylolm(log_developmental_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda); #print(summary(mx.r))
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(log_developmental_d ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda); #print(summary(mx.r.phy))
R2.lik(mx) # full model R2
R2.lik(mx.r) # phylogeny only model
R2.lik(mx.r.phy) # body size only model
R2.lik(mx, mx.r) # partial R2 (body mass)
R2.lik(mx, mx.r.phy) # partial R2 (phylogeny)

# different intercepts (all orders with 10 or more spp)
mx <- NA
fx <- dat[,c("binomial","log_developmental_d","log_adult_body_mass","Taxon_subgroup")]
fx <- fx[complete.cases(fx),]
fx <- fx[fx$Taxon_subgroup %in% names(orders),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
mx <- phylolm(log_developmental_d ~ log_adult_body_mass + Taxon_subgroup, data = cx$data, phy = cx$phy, model = "lambda"); #print(summary(mx))


# ---------- #


### INCUBATION FRACTION

# adult body mass

# single predictor (all data)
mx <- NA
fx <- dat[,c("binomial","prop_incubation_sqrt","log_adult_body_mass")]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
mx <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda"); #print(summary(mx))
mlpar <- data.frame(lambda=mx$optpar)
mx.r <- phylolm(prop_incubation_sqrt ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda); #print(summary(mx.r))
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda); #print(summary(mx.r.phy))
R2.lik(mx) # full model R2
R2.lik(mx.r) # phylogeny only model
R2.lik(mx.r.phy) # body size only model
R2.lik(mx, mx.r) # partial R2 (body mass)
R2.lik(mx, mx.r.phy) # partial R2 (phylogeny)

# different intercepts (all orders with 10 or more spp)
mx <- NA
fx <- dat[,c("binomial","prop_incubation_sqrt","log_adult_body_mass","Taxon_subgroup")]
fx <- fx[complete.cases(fx),]
fx <- fx[fx$Taxon_subgroup %in% names(orders),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
mx <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass + Taxon_subgroup, data = cx$data, phy = cx$phy, model = "lambda"); #print(summary(mx))


# ====================== #

