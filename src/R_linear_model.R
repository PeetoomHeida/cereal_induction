install.packages("lmerTest")
install.packages("tidyverse")
install.packages("emmeans")
install.packages("xtable")
install.packages("car")
library(lme4)
library(tidyverse)
#library(xtable)
library(emmeans)
library(car)
si_data <- as.data.frame(read.csv("data/real_data/biomass_si_data.csv"))
si_filtered  <- si_data[si_data$isDamaged != "Undamaged",]
my_mod <- lme4::lmer(Si_ppm ~ Induction * Species +mass_g + (1|Genotype), si_filtered)
summary(my_mod)
anova(my_mod)
Anova(my_mod, type ="III", test.statistic = "F")
emgrid1  <- emmeans(my_mod, list(pairwise ~ Species), adjust = "tukey", type = "response")
emgrid2 <- emmeans(my_mod,  ~ Induction * Species, adjust = "mvt", type = "response")
my_emm <- pairs(emmeans(my_mod, ~ Induction * Species))
foo <- pwpm(my_emm, by = NULL, adjust = "bonferroni")
xtable::xtable(foo)
emmeans(my_mod, list(pairwise ~ Induction), adjust = "tukey", type = "response")
ggplot(si_filtered, aes(x = Induction, y = Si_ppm)) +
geom_col() +
theme_classic(base_size = 22)

insect_data <- si_data[si_data$Induction == "Insect",]
insect_data$isDamaged <- as.factor(insect_data$isDamaged)
insect_data$Si_percent <- insect_data$Si_ppm/10000
logit_mod <- glmer(isDamaged ~ Si_percent + Species + (1|Genotype), insect_data, family = "binomial")
summary(logit_mod)
anova(logit_mod)

pigsint.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
pigs.sav <- as.list(ref_grid(pigsint.lm))
pigs.anew <- as.emmGrid(pigs.sav)
pigsint.emm <- pairs(emmeans(pigsint.lm, ~ percent | source))
test(pigsint.emm, by = NULL, adjust = "bonferroni")
typeof(pigsint.emm)
xtable::xtable(pigsint.emm, type = "response")

