# Results

# Author: Amadeus
# Version: 2021-11-25

# Libraries
library(tidyverse)
library(rstatix)
library(corrplot)
library(Rfit)
library(mediation)
library(WRS2)
library(sfsmisc)
library(psych)
library(MBESS)
library(robustbase)

# Sources
dat <- read_csv("processed_data.csv")

# Parameters

# ============================================================================

# Code

## reliability

dat_minik <- dat %>%
  dplyr::select(minik1:minik20)

omega(dat_minik)
ci.reliability(data = dat_minik, 
               type = "omega", 
               conf.level = 0.95, 
               interval.type = "bca", B = 1000) # 0.88, 95%CI 0.8576-0.9073

dat_FTP <- dat %>%
  dplyr::select(FTP1:FTP10)

omega(dat_FTP)
ci.reliability(data = dat_FTP, 
               type = "omega", 
               conf.level = 0.95, 
               interval.type = "bca", B = 1000) # 0.88, 0.8549-0.8989

dat_SES <- dat %>%
  dplyr::select(SES1:SES3)

ci.reliability(data = dat_SES, 
               type = "omega", 
               conf.level = 0.95, 
               interval.type = "bca", B = 1000) # 0.9148 95%CI 0.88637-0.93338

## Table 1 and distributions

### Age: Distribution

ggplot(dat) +
  geom_histogram(aes(x = age, y = ..density..), 
                 col = "#D53E4F", fill = "#D53E4F", alpha = .5,
                 binwidth = 1) +
  geom_density(aes(x = age, y = ..density..), 
               col = "#D53E4F", size = 1.5) +
  theme_bw()

dat %>%
  group_by(age_group) %>%
  summarise(across(age, list(mean = mean, sd = sd, range = range)))

### subjective health

dat %>%
  group_by(age_group) %>%
  summarise(across(contains("health"), list(mean = mean, sd = sd)))

## Vertical and Horizontal Trade-offs

dat %>%
  mutate(logk = log(k)) %>%
  cor_test(minik, logk)

### horizontal life history tradeoffs

rlm <- lmrob(logk ~ minik + age + age_condition + age:age_condition, 
             data = dat %>%
               filter(age_group != "middle") %>%
               mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
               mutate(logk = log(k)))

summary(rlm)


## Mechanisms underlying horizontal trade-offs

### age bias: antecedents
rlm <-  lmrob(age_bias_phys ~ minik + age + health_phys, 
               data = dat %>%
                 mutate(age_bias_phys = age_phys - age,
                        logk = log(k)) %>%
                filter(age_group == "young"))

summary(rlm)

rlm <-  lmrob(health_phys ~ age, 
              data = dat %>%
                mutate(age_bias_phys = age_phys - age,
                       logk = log(k)) %>%
                filter(age_group == "old"))

summary(rlm)

dat %>%
  mutate(age_bias_phys = age_phys - age,
         logk = log(k)) %>%
  group_by(age_group) %>%
  cor_test(age, health_phys , method = "kendall")

### age bias and health condition
dat %>%
  mutate(age_bias_phys = age_phys - age,
         logk = log(k)) %>%
  cor_test(age_phys, age_bias_phys, method = "kendall")

### age bias and interaction
dat %>%
  mutate(age_bias_phys = age_phys - age,
              logk = log(k)) %>%
  cor_test(age_bias_phys, age)
  

rlm2 <-  lmrob(logk ~ minik + age + age_bias_phys + age_condition + age_condition:age_bias_phys + age:age_condition, 
               data = dat %>%
                 mutate(age_bias_phys = age_phys - age,
                        logk = log(k)) %>%
                 filter(age_group != "middle", !is.na(logk)) %>%
                 mutate(age_condition = ifelse(age_group == "young", 0, 1)))

summary(rlm2)

### mediation effects of age bias scores
robu_mediation <- with(dat %>%
                         drop_na() %>%
                         mutate(logk = log(k),
                                age_bias_phys = age_phys - age) %>%
                         filter(age_group == "young"), 
                       ZYmediate(health_phys, logk, age_bias_phys, nboot = 1999))
robu_mediation
robu_mediation$a.est
robu_mediation$b.est
robu_mediation$ab.est
robu_mediation

dat %>%
  group_by(platform, age_group) %>%
  summarise(across(minik, list(mean = mean, sd = sd)))
