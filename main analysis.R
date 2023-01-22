# Results

# Author: Junsong Lu
# Version: 2021-11-24

# Libraries

required_pack <- c('tidyverse','rstatix','corrplot', 'Rfit', 'mediation', 
                   'WRS2', 'sfsmisc', 'psych', 'MBESS', 'robustbase')

for(p in required_pack){
  if(!require(p, character.only = TRUE)){
    install.packages(p)
  }else{
    library(p, character.only = TRUE)
  }
} 

# Sources
dat <- read_csv("data.csv")

dat %>%
  group_by(age_group) %>%
  summarise(across(age, list(mean = mean, sd = sd,
                             rang1 = range)))

# Parameters

# ============================================================================

# Code

## Reliability

dat_minik <- dat %>%
  dplyr::select(minik1:minik20)

ci.reliability(data = dat_minik, 
               type = "omega", 
               conf.level = 0.95, 
               interval.type = "bca", B = 1000) # 0.88, 95%CI 0.8576-0.9073

dat_FTP <- dat %>%
  dplyr::select(FTP1:FTP10)

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

### Delay discounting 
summary(dat$k)

ggplot(dat %>%
         filter(!is.na(k))) +
  geom_histogram(aes(x = k), 
                 col = "#D53E4F", fill = "#D53E4F", alpha = .5,
                 binwidth = 0.01) +
  geom_density(aes(x = k, y = ..density..), 
               col = "#D53E4F", size = 1.5) +
  theme_bw()

ggplot(dat %>%
         rowwise() %>%
         mutate(ll_count = sum(c_across(starts_with("DD")) == 2),
                ss_count = sum(c_across(starts_with("DD")) == 1))) +
  geom_histogram(aes(x = ss_count), 
                 col = "#D53E4F", fill = "#D53E4F", alpha = .5,
                 binwidth = 1) +
  geom_density(aes(x = ss_count, y = ..density..), 
               col = "#D53E4F", size = 1.5) +
  theme_bw()

### Subjective health

dat %>%
  group_by(age_group) %>%
  summarise(across(contains("health"), list(mean = mean, sd = sd)))

## Delay Discounting as Vertical and Horizontal Life History Trade-offs

dat %>%
  cor_test(minik, SES, method = "kendall")

dat %>%
  cor_test(minik, k, method = "kendall")

dat %>%
  group_by(age_group) %>%
  cor_test(k, age, method = "kendall")

dat %>%
  group_by(age_group) %>%
  cor_test(k, age_phys, method = "kendall")

dat %>%
  group_by(age_group) %>%
  cor_test(k, age_psyc, method = "kendall")

## Differential Effects of Vertical and Horizontal Life History Trade-offs on Delay Discounting

rlm <- lmrob(logk ~ minik + age + age_condition + age:age_condition, 
             data = dat %>%
               filter(age_group != "middle") %>%
               mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
               mutate(logk = log(k)))

summary(rlm)


## Mechanisms underlying vertical trade-offs
set.seed(111)
robu_mediation <- with(dat %>%
                         drop_na() %>%
                         mutate(logk = log(k)), ZYmediate(SES, logk,
                                                          minik, nboot = 1999))
robu_mediation
robu_mediation$a.est
robu_mediation$b.est
robu_mediation$ab.est
robu_mediation

set.seed(111)
robu_mediation2 <- with(dat %>%
                         drop_na() %>%
                         mutate(logk = log(k)), ZYmediate(minik, logk,
                                                          beta, nboot = 1999))
robu_mediation2
robu_mediation2$a.est
robu_mediation2$b.est
robu_mediation2$ab.est
robu_mediation2

dat %>%
  cor_test(alpha, minik, method = "kendall")

dat %>%
  cor_test(beta, minik, method = "kendall")

dat %>%
  cor_test(FTP, minik, method = "kendall")

dat %>%
  cor_test(FTP, k, method = "kendall")

## Mechanisms underlying horizontal trade-offs

dat %>%
  cor_test(k, FTP, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(k, alpha, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(k, beta, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(FTP, beta, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(FTP, alpha, method = "kendall")

dat %>%
  cor_test(FTP, minik, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(minik, alpha, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 1.323433) %>%
  cor_test(minik, beta, method = "kendall")

rlm2 <-  lmrob(logk ~ minik + age + age_bias_phys + age_condition + age_condition:age_bias_phys + age:age_condition, 
               data = dat %>%
                 mutate(age_bias_phys = age_phys - age,
                        logk = log(k)) %>%
                 filter(age_group != "middle", !is.na(logk)) %>%
                 mutate(age_condition = ifelse(age_group == "young", 0, 1)))

summary(rlm2)

rlm3 <-  lmrob(logk ~ minik + age + age_bias_psyc + age_condition + age_condition:age_bias_psyc + age:age_condition, 
               data = dat %>%
                 mutate(age_bias_psyc = age_psyc - age,
                        logk = log(k)) %>%
                 filter(age_group != "middle", !is.na(logk)) %>%
                 mutate(age_condition = ifelse(age_group == "young", 0, 1)))

summary(rlm3)


dat %>%
  mutate(age_bias_phys = age_phys - age) %>%
  group_by(age_group) %>%
  cor_test(health_phys, age_phys, method = "kendall")

dat %>%
  mutate(age_bias_phys = age_phys - age) %>%
  group_by(age_group) %>%
  cor_test(age_bias_phys, k , method = "kendall")

set.seed(111)
robu_mediation3 <- with(dat %>%
                         drop_na() %>%
                          filter(age_group == "old") %>%
                         mutate(logk = log(k),
                                age_bias_phys = age_phys - age), ZYmediate(health_phys, logk,
                                                                           age_phys, nboot = 1999))
robu_mediation3
robu_mediation3$a.est
robu_mediation3$b.est
robu_mediation3$ab.est
robu_mediation3

set.seed(1111)
robu_mediation4 <- with(dat %>%
                          drop_na() %>%
                          filter(age_group == "middle") %>%
                          mutate(logk = log(k),
                                 age_bias_phys = age_phys - age), ZYmediate(health_phys, logk,
                                                                            age_phys, nboot = 1999))
robu_mediation4
robu_mediation4$a.est
robu_mediation4$b.est
robu_mediation4$ab.est
robu_mediation4

dat %>%
  drop_na() %>%
  group_by(age_group) %>%
  mutate(logk = log(k),
         age_bias_phys = age_phys - age) %>%
  cor_test(age_bias_phys, k, method = "kendall")

### gender issues

rlm4 <- lmrob(logk ~ minik + age + age:gender, 
               data = dat %>%
                 mutate(logk = log(k)) %>%
                 filter(age_group == "old", !is.na(logk)) %>%
                 mutate(gender = ifelse(gender == 1, 0, 1)) %>%
                 mutate(gender = factor(gender)))

summary(rlm4)

rlm5 <- lmrob(logk ~ minik + age, 
              data = dat %>%
                mutate(logk = log(k)) %>%
                filter(age_group == "old", !is.na(logk)) %>%
                mutate(gender = ifelse(gender == 1, 0, 1)) %>%
                filter(gender == 1))

summary(rlm5) # female

rlm6 <- lmrob(logk ~ minik + age, 
               data = dat %>%
                 mutate(logk = log(k)) %>%
                 filter(age_group == "old", !is.na(logk)) %>%
                 mutate(gender = ifelse(gender == 1, 0, 1)) %>%
                 filter(gender == 0))

summary(rlm6) # male

dat %>%
  filter(age_group == "old") %>%
  pivot_longer(cols = c(age, age_psyc),
               names_to = "age_type",
               values_to = "age_value") %>%
  t_test(age_value ~ age_type, paired = TRUE)

dat %>%
  filter(age_group == "old") %>%
  pivot_longer(cols = c(age, age_psyc),
               names_to = "age_type",
               values_to = "age_value") %>%
  group_by(age_type) %>%
  summarise(across(age_value, list(mean = mean,
                                   sd = sd)))

dat %>%
  filter(age_group == "old") %>%
  group_by(gender) %>%
  cor_test(k, age, method = "kendall")

dat %>%
  filter(age_group == "young") %>%
  mutate(logk = log(k)) %>%
  anova_test(logk ~ gender + minik)

# reply
dat <- dat %>%
  mutate(logk = log(k))

shapiro.test(dat$logk)
shapiro.test(dat$age)

shapiro.test(dat$minik)
shapiro.test(dat$FTP)
shapiro.test(dat$SES)
shapiro.test(dat$health_phys)
shapiro.test(dat$health_psyc)

# two-lines analysis
png("twoline2.png", units = "in",
    width = 7.5, height = 5.5, res = 700)
twolines(logk ~ age + minik + SES, data = dat %>%
                 filter(!is.na(k)) %>%
                 mutate(logk = log(k),
                        age_sq = age*age) %>%
                 as.data.frame)
dev.off()
summary(tl)

lm1 <- rlm(logk ~ age + age_sq + minik + SES, data = dat %>%
                  filter(!is.na(k)) %>%
                  mutate(logk = log(k),
                         age_sq = age*age) %>%
                  as.data.frame)
summary(lm1)
f.robftest(lm1, var = "age_sq")
