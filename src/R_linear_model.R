install.packages("lmerTest")
install.packages("tidyverse")
install.packages("emmeans")
install.packages("xtable")
install.packages("car")
install.packages("Hmisc")
install.packages("rstanarm")
library(lme4)
library(pbkrtest)
library(tidyverse)
#library(xtable)
library(emmeans)
library(car)
library(Hmisc)
library(lmerTest)
library(rstanarm)
options(contrasts = c("contr.helmert", "contr.poly"))
si_data <- as.data.frame(read.csv("data/real_data/biomass_si_data.csv"))
si_filtered  <- si_data[si_data$isDamaged != "Undamaged",]
si_filtered$Genotype <- as.factor(si_filtered$Genotype)
si_filtered$Induction <- as.factor(si_filtered$Induction)
si_filtered$Species <- as.factor(si_filtered$Species)
si_mod <- lmerTest::lmer(log(Si_ppm) ~ Induction * Species + mass_g + (1|Species/Genotype), si_filtered)
si_mod <- lme4::lmer(log(Si_ppm) ~ Induction * Species + mass_g + (1|Genotype), si_filtered)
anova(si_mod, type = 2)
summary(aov(Si_ppm/10000 ~ Induction * Species + mass_g, si_filtered))
summary(si_mod)
lmerTest::anova(si_mod)
Anova(si_mod, type =2, test.statistic = "F")
si_mod0 <- lmer(log(Si_ppm) ~  mass_g + (1|Genotype), si_filtered)
si_mod1 <- lmer(log(Si_ppm) ~ Induction + mass_g + (1|Genotype), si_filtered)
si_mod2 <- lmer(log(Si_ppm) ~ Species + mass_g + (1|Genotype), si_filtered)
si_mod3 <- lmer(log(Si_ppm) ~ Induction + Species + mass_g + (1|Genotype), si_filtered)
si_mod4 <- lmer(log(Si_ppm) ~ Induction * Species + mass_g + (1|Genotype), si_filtered)
si_mod5 <- lmer(log(Si_ppm) ~ Induction:Species + mass_g + (1|Genotype), si_filtered)
fit.formula <- log(Si_ppm) ~ Induction+Species + mass_g
X <- model.matrix(fit.formula, si_filtered)
pb01 <- PBmodcomp(si_mod1, si_mod0, nsim = 10000, cl = 4)
pb02 <- PBmodcomp(si_mod2, si_mod0, nsim = 10000, cl = 4)
pb03 <- PBmodcomp(si_mod3, si_mod0, nsim = 10000, cl = 4)
pb32 <- PBmodcomp(si_mod3, si_mod2, nsim=10000, cl=4)
pb34 <- PBmodcomp(si_mod4, si_mod3, nsim = 10000, cl = 4)
pbinonly <- PBmodcomp(si_mod5, si_mod2, nsim= 1000, cl = 4)


emgrid1  <- emmeans(si_mod, list(pairwise ~ Species), adjust = "tukey", type = "response")
emgrid2 <- emmeans(si_mod,  ~ Induction * Species, adjust = "mvt", type = "response")
my_emm <- pairs(emmeans(si_mod, ~ Induction * Species))
foo <- pwpm(my_emm, by = NULL, adjust = "bonferroni")
xtable::xtable(foo)
emmeans(si_mod, list(pairwise ~ Species|Induction), adjust = "tukey", type = "response")
ggplot(si_filtered, aes(x = Induction, y = Si_ppm)) +
geom_col() +
theme_classic(base_size = 22)

insect_data <- si_data[si_data$Induction == "Insect",]
insect_data$isDamaged <- as.factor(insect_data$isDamaged)
insect_data$Si_percent <- insect_data$Si_ppm/10000
logit_mod <- glmer(Si_percent ~ Si_percent + Species + (1|Genotype), insect_data, family = "binomial")
summary(logit_mod)
anova(logit_mod)



si_mod0 <- lmer(log(Si_ppm) ~  mass_g + (1|Genotype), si_filtered)
si_mod1 <- lmer(log(Si_ppm) ~ Induction + mass_g + (1|Genotype), si_filtered)
si_mod2 <- lmer(log(Si_ppm) ~ Species + mass_g + (1|Genotype), si_filtered)
si_mod3 <- lmer(log(Si_ppm) ~ Induction + Species + mass_g + (1|Genotype), si_filtered)
si_mod4 <- lmer(log(Si_ppm) ~ Induction * Species + mass_g + (1|Genotype), si_filtered)
si_mod5 <- lmer(log(Si_ppm) ~ Induction:Species + mass_g + (1|Genotype), si_filtered)
fit.formula <- log(Si_ppm) ~ Induction+Species + mass_g
X <- model.matrix(fit.formula, si_filtered)
pb01 <- PBmodcomp(si_mod1, si_mod0, nsim = 10000, cl = 4)
pb02 <- PBmodcomp(si_mod2, si_mod0, nsim = 10000, cl = 4)
pb03 <- PBmodcomp(si_mod3, si_mod0, nsim = 10000, cl = 4)
pb34 <- PBmodcomp(si_mod4, si_mod3, nsim = 10000, cl = 4)
pbinonly <- PBmodcomp(si_mod5, si_mod2, nsim= 1000, cl = 4)

absorbancedf <- as.data.frame(read.csv("data/real_data/si_absorbance_data.csv"))
abdf_filtered  <- absorbancedf[absorbancedf$isDamaged != "Undamaged",]
abdf_filtered$Genotype <- as.factor(abdf_filtered$Genotype)
abdf_filtered$Species <- as.factor(abdf_filtered$Species)
abdf_filtered$Induction <- as.factor(abdf_filtered$Induction)
ab_mod <- lme4::lmer(mcabsorbance ~ Induction * Species + mass_g + (1|Genotype), abdf_filtered)
summary(ab_mod)
Anova(ab_mod)
summary(lm(mcabsorbance ~ mass_g, absorbancedf))
absorbance_mod <- lmer(scale(mcabsorbance) ~ scale(Si_ppm) + scale(mass_g) + Species+ Induction + (1|Genotype) ,REML = FALSE, abdf_filtered)
control_absorbance_only <- absorbancedf[absorbancedf$Induction == "Control",]
summary(lmerTest::lmer(scale(log(mcabsorbance)) ~ scale(log(Si_ppm)) + Species + mass_g + (1|Genotype), control_absorbance_only))
summary(lm(scale(mcabsorbance) ~ scale(Si_ppm), control_absorbance_only))
absi_mod0 <- lme4::lmer(log(mcabsorbance) ~  mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
absi_mod1 <- lme4::lmer(log(mcabsorbance) ~ Induction + mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
absi_mod2 <- lme4::lmer(log(mcabsorbance) ~ Species + mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
absi_mod3 <- lme4::lmer(log(mcabsorbance) ~ Induction + Species + mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
absi_mod4 <- lme4::lmer(log(mcabsorbance) ~ Induction * Species + mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
absi_mod5 <- lme4::lmer(log(mcabsorbance) ~ Induction:Species + mass_g + (1|Genotype), abdf_filtered, REML = FALSE)
fit.formula <- log(Si_ppm) ~ Induction+Species + mass_g
X <- model.matrix(fit.formula, si_filtered)
abpb01 <- PBmodcomp(si_mod1, si_mod0, nsim = 10000, cl = 4)
abpb02 <- PBmodcomp(si_mod2, si_mod0, nsim = 10000, cl = 4)
abpb03 <- PBmodcomp(si_mod3, si_mod0, nsim = 10000, cl = 4)
abpb34 <- PBmodcomp(si_mod4, si_mod3, nsim = 10000, cl = 4)
pbinonly <- PBmodcomp(si_mod5, si_mod2, nsim = 1000, cl = 4)



si_mod0 <- lme4::lmer(log(Si_ppm) ~  mass_g + (1|Species/Genotype), si_filtered)
si_mod1 <- lme4::lmer(log(Si_ppm) ~ Induction + mass_g + (1|Species/Genotype), si_filtered)
si_mod2 <- lme4::lmer(log(Si_ppm) ~ Species + mass_g + (1|Species/Genotype), si_filtered)
si_mod3 <- lme4::lmer(log(Si_ppm) ~ Induction + Species + mass_g + (1|Species/Genotype), si_filtered)
si_mod4 <- lme4::lmer(log(Si_ppm) ~ Induction * Species + mass_g + (1|Species/Genotype), si_filtered)
si_mod5 <- lme4::lmer(log(Si_ppm) ~ Induction:Species + mass_g + (1|Species/Genotype), si_filtered)


stan_mod <- stan_lmer((log(Si_ppm) ~ Induction * Species + mass_g + (1|Genotype)), si_filtered)
lmerTestFullMod <- lmerTest::lmer((log(Si_ppm) ~ Induction * Species + mass_g + (1|Genotype/Species)), si_filtered)
summary(lmerTestFullMod)
anova(lmerTestFullMod, type = 2)
with(si_filtered, table(Genotype, Species, Induction))
with(si_filtered, table(Genotype, Induction, Species))

test_dat <- data.frame(
    "response" = c(rnorm(288)), 
    "Induction" = c(rep(c(rep(c(rep("1", 2), rep("2",2), rep("3",2)),12)),4)), 
    "Species" = c(rep(c(rep("A",18), rep("B", 18), rep("C",18), rep("D",18)),4)), 
    "GenotypeID" = c(rep(c(rep(c(rep("1", 6), rep("2",6), rep("3", 6)),4)),4)),
    "biomass" = rnorm(288))
test_dat$UniqueGeno <- paste(test_dat$Species,test_dat$GenotypeID, sep="")
test_dat$Induction <- as.factor(test_dat$Induction)
test_dat$Species <- as.factor(test_dat$Species)
test_dat$UniqueGeno <- as.factor(test_dat$Species)
ex_mod <- lmer(response ~ Induction * Species + biomass + (1|UniqueGeno), test_dat)
summary(ex_mod)
anova(ex_mod, type = 2)
