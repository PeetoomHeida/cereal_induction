install.packages("lmerTest")
install.packages("tidyverse")
install.packages("emmeans")
install.packages("xtable")
install.packages("car")
install.packages("Hmisc")
library(lme4)
library(tidyverse)
#library(xtable)
library(emmeans)
library(car)
library(Hmisc)
si_data <- as.data.frame(read.csv("data/real_data/si_absorbance_data.csv"))
si_filtered  <- si_data[si_data$isDamaged != "Undamaged",]
si_mod <- lme4::lmer(Si_ppm/10000 ~ Induction * Species +mass_g + (1|Genotype), si_filtered)
summary(si_mod)
anova(si_mod)
Anova(si_mod, type ="III", test.statistic = "F")
emgrid1  <- emmeans(si_mod, list(pairwise ~ Species), adjust = "tukey", type = "response")
emgrid2 <- emmeans(si_mod,  ~ Induction * Species, adjust = "mvt", type = "response")
my_emm <- pairs(emmeans(si_mod, ~ Induction * Species))
foo <- pwpm(my_emm, by = NULL, adjust = "bonferroni")
xtable::xtable(foo)
emmeans(si_mod, list(pairwise ~ Induction), adjust = "tukey", type = "response")
ggplot(si_filtered, aes(x = Induction, y = Si_ppm)) +
geom_col() +
theme_classic(base_size = 22)

insect_data <- si_data[si_data$Induction == "Insect",]
insect_data$isDamaged <- as.factor(insect_data$isDamaged)
insect_data$Si_percent <- insect_data$Si_ppm/10000
logit_mod <- glmer(Si_percent ~ Si_percent + Species + (1|Genotype), insect_data, family = "binomial")
summary(logit_mod)
anova(logit_mod)

pigsint.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
pigs.sav <- as.list(ref_grid(pigsint.lm))
pigs.anew <- as.emmGrid(pigs.sav)
pigsint.emm <- pairs(emmeans(pigsint.lm, ~ percent | source))
test(pigsint.emm, by = NULL, adjust = "bonferroni")
typeof(pigsint.emm)
xtable::xtable(pigsint.emm, type = "response")


absorbancedf <- as.data.frame(read.csv("data/real_data/si_absorbance_data.csv"))
abdf_filtered  <- absorbancedf[absorbancedf$isDamaged != "Undamaged",]
ab_mod <- lme4::lmer(mcabsorbance ~ Induction * Species + mass_g + (1|Genotype), abdf_filtered)
summary(ab_mod)
Anova(ab_mod)
ggplot(abdf_filtered, aes(x = Species, y = mcabsorbance, group=Induction)) +
    #geom_col(position = "dodge") +
    stat_summary(fun.y = "mean", geom = "col", width = .9, fill = "gray69", position = "dodge") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width = 0.90), fun.args = list(mult = 1.96))
