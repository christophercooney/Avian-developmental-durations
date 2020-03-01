
## Single and multipredictor models

rm(list=ls())

library(ape)
library(caper)
library(phylolm)
library(rr2)
library(car)

setwd("~/Dropbox/Projects/Incubation/Comparative_analysis/")

# ---------- #

# get tree and data

phy <- read.tree("Outputs/MCC_tree.tre")
dat <- read.csv("Outputs/Comparative_dataset.csv", stringsAsFactors = F); rownames(dat) <- dat$binomial; head(dat)

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

pred.vars <- c("log_adult_body_mass", "log_egg_mass", "log_generation_length", "log_clutch_size", "dev_mod",
               "parental_care_unibi", "brood_parasite",
               "nest_height_min", "nest_type",
               "habitat", "diet", "pelagic", "nocturnal", "migration",
               "temperature", "precipitation",
               "latitude", "insularity")

# ---------- #


### DEVELOPMENTAL DURATION

## Single predictor models
res.single <- as.list(rep(NA, length(pred.vars)))
names(res.single) <- pred.vars
out.single <- c()
for (i in 1:length(pred.vars)) {
  cat("\r", i, "of", length(pred.vars))
  fx <- dat[,c("binomial","log_developmental_d", pred.vars[i])]
  fx <- fx[complete.cases(fx),]
  px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
  cx <- comparative.data(px, fx, "binomial")
  mx <- phylolm(formula(paste0("log_developmental_d ~ ", pred.vars[i])), data = cx$data, phy = cx$phy, model = "lambda")
  mlpar <- data.frame(lambda=mx$optpar)
  mx.r <- phylolm(log_developmental_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  out.single <- rbind(out.single, c(pred.vars[i], nrow(cx$data), mx$p - mx.r$p, mx.r$aic - mx$aic, R2.lik(mx, mx.r)))
  res.single[[i]] <- list(mx, mx.r)
}
out.single <- data.frame(out.single, stringsAsFactors=F); colnames(out.single) <- c("predictor","n","df","dAIC","R2.lik")

## Multipredictor models (with body mass)
single <- out.single[out.single$predictor != "log_egg_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_developmental_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_developmental_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_developmental_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_developmental_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_developmental_d ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)

## Multipredictor models (with egg mass)
single <- out.single[out.single$predictor != "log_adult_body_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_developmental_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_developmental_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_developmental_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_developmental_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_developmental_d ~ log_egg_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)


# --------- #


### INCUBATION FRACTION

## Single predictor models
res.single <- as.list(rep(NA, length(pred.vars)))
names(res.single) <- pred.vars
out.single <- c()
for (i in 1:length(pred.vars)) {
  cat("\r", i, "of", length(pred.vars))
  fx <- dat[,c("binomial","prop_incubation_sqrt", pred.vars[i])]
  fx <- fx[complete.cases(fx),]
  px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
  cx <- comparative.data(px, fx, "binomial")
  mx <- phylolm(formula(paste0("prop_incubation_sqrt ~ ", pred.vars[i])), data = cx$data, phy = cx$phy, model = "lambda")
  mlpar <- data.frame(lambda=mx$optpar)
  mx.r <- phylolm(prop_incubation_sqrt ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  out.single <- rbind(out.single, c(pred.vars[i], nrow(cx$data), mx$p - mx.r$p, mx.r$aic - mx$aic, R2.lik(mx, mx.r)))
  res.single[[i]] <- list(mx, mx.r)
}
out.single <- data.frame(out.single, stringsAsFactors=F); colnames(out.single) <- c("predictor","n","df","dAIC","R2.lik")

## Multipredictor models (with body mass)
single <- out.single[out.single$predictor != "log_egg_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","prop_incubation_sqrt", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("prop_incubation_sqrt ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("prop_incubation_sqrt ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(prop_incubation_sqrt ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(prop_incubation_sqrt ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)

## Multipredictor models (with egg mass)
single <- out.single[out.single$predictor != "log_adult_body_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","prop_incubation_sqrt", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("prop_incubation_sqrt ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("prop_incubation_sqrt ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(prop_incubation_sqrt ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(prop_incubation_sqrt ~ log_egg_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)


# --------- #


### INCUBATION DURATION

## Single predictor models
res.single <- as.list(rep(NA, length(pred.vars)))
names(res.single) <- pred.vars
out.single <- c()
for (i in 1:length(pred.vars)) {
  cat("\r", i, "of", length(pred.vars))
  fx <- dat[,c("binomial","log_incubation_d", pred.vars[i])]
  fx <- fx[complete.cases(fx),]
  px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
  cx <- comparative.data(px, fx, "binomial")
  mx <- phylolm(formula(paste0("log_incubation_d ~ ", pred.vars[i])), data = cx$data, phy = cx$phy, model = "lambda")
  mlpar <- data.frame(lambda=mx$optpar)
  mx.r <- phylolm(log_incubation_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  out.single <- rbind(out.single, c(pred.vars[i], nrow(cx$data), mx$p - mx.r$p, mx.r$aic - mx$aic, R2.lik(mx, mx.r)))
  res.single[[i]] <- list(mx, mx.r)
}
out.single <- data.frame(out.single, stringsAsFactors=F); colnames(out.single) <- c("predictor","n","df","dAIC","R2.lik")

## Multipredictor models (with body mass)
single <- out.single[out.single$predictor != "log_egg_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_incubation_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_incubation_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_incubation_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_incubation_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_incubation_d ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)

## Multipredictor models (with egg mass)
single <- out.single[out.single$predictor != "log_adult_body_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_incubation_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_incubation_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_incubation_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_incubation_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_incubation_d ~ log_egg_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)


# --------- #


### FLEDGING DURATION

## Single predictor models
res.single <- as.list(rep(NA, length(pred.vars)))
names(res.single) <- pred.vars
out.single <- c()
for (i in 1:length(pred.vars)) {
  cat("\r", i, "of", length(pred.vars))
  fx <- dat[,c("binomial","log_fledging_d", pred.vars[i])]
  fx <- fx[complete.cases(fx),]
  px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
  cx <- comparative.data(px, fx, "binomial")
  mx <- phylolm(formula(paste0("log_fledging_d ~ ", pred.vars[i])), data = cx$data, phy = cx$phy, model = "lambda")
  mlpar <- data.frame(lambda=mx$optpar)
  mx.r <- phylolm(log_fledging_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  out.single <- rbind(out.single, c(pred.vars[i], nrow(cx$data), mx$p - mx.r$p, mx.r$aic - mx$aic, R2.lik(mx, mx.r)))
  res.single[[i]] <- list(mx, mx.r)
}
out.single <- data.frame(out.single, stringsAsFactors=F); colnames(out.single) <- c("predictor","n","df","dAIC","R2.lik")

## Multipredictor models (with body mass)
single <- out.single[out.single$predictor != "log_egg_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_fledging_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_fledging_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_fledging_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_fledging_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_fledging_d ~ log_adult_body_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)

## Multipredictor models (with egg mass)
single <- out.single[out.single$predictor != "log_adult_body_mass",]
multi.pred.vars <- single$predictor[as.numeric(single$dAIC) > 2]
fx <- dat[,c("binomial","log_fledging_d", multi.pred.vars)]
fx <- fx[complete.cases(fx),]
px <- drop.tip(phy, setdiff(phy$tip.label, fx$binomial))
cx <- comparative.data(px, fx, "binomial")
form <- formula(paste0("log_fledging_d ~ ", paste(multi.pred.vars, collapse=" + "))); vif(lm(form, cx$data))
mx <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda"); summary(mx)
mlpar <- data.frame(lambda=mx$optpar)
pred.out <- data.frame(predictor=multi.pred.vars, n=NA, df=NA, dAIC=NA, R2.lik=NA)
for (i in 1:length(multi.pred.vars)) {
  form.r <- paste0("log_fledging_d ~ ", paste(multi.pred.vars[-i], collapse=" + "))
  mx.r <- phylolm(form.r, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda)
  pred.out$n[i] <- nrow(mx$X)
  pred.out$df[i] <- mx$p - mx.r$p
  pred.out$dAIC[i] <- mx.r$aic - mx$aic
  pred.out$R2.lik[i] <- R2.lik(mx, mx.r)
}
# r2 values
mx.r <- phylolm(log_fledging_d ~ 1, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # intercept-only model
R2.lik(mx) # full model r2
R2.lik(mx, mx.r) # partial r2 (all predictors)
# 'Ecology' effect
mx.r.ecol <- phylolm(log_fledging_d ~ log_egg_mass, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.ecol) # partial r2 (ecology)
# Phylogeny effect
mlpar <- data.frame(lambda=0)
mx.r.phy <- phylolm(form, data = cx$data, phy = cx$phy, model = "lambda", starting.value = mlpar, lower.bound = mlpar$lambda, upper.bound = mlpar$lambda) # predictors only, no phylo
R2.lik(mx, mx.r.phy) # partial r2 (phylo)


# ---------------------------------------- #

