# Data Analysis

# Author: 不知年
# Version: 2021-11-16

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

## Age: Distribution

ggplot(dat) +
  geom_histogram(aes(x = age, y = ..density..), 
                 col = "#D53E4F", fill = "#D53E4F", alpha = .5,
                 binwidth = 1) +
  geom_density(aes(x = age, y = ..density..), 
               col = "#D53E4F", size = 1.5) +
  theme_bw()

ggplot(dat) +
  geom_histogram(aes(log(k)), fill = "#FF9900") +
  theme_bw()

ggplot(dat) +
  geom_histogram(aes(beta), fill = "#FF9900") +
  theme_bw()

summary(dat$age[dat$age_group == "old"])
sd(dat$age)

dat %>%
  group_by(age_group, gender) %>%
  summarise(count = n())

dat %>%
  group_by(age_group) %>%
  cor_test(FTP, age, method = "kendall")

dat %>%
  group_by(age_group) %>%
  mutate(logk = log(k)) %>%
  filter(!is.na(logk)) %>%
  summarise(across(logk, list(mean = mean, sd = sd)))

dat %>%
  group_by(age_group) %>%
  summarise(count = n(), range_c = c("lower", "upper"), range = range(age))

ggplot(dat) +
  geom_point(aes(age, k), 
             position = position_jitter(width = .3, height = .08),
             fill = "#FF9900") +
  theme_bw()

ggplot(dat) +
  geom_point(aes(age, alpha), 
             position = position_jitter(width = .3, height = .08),
             fill = "#FF9900") +
  theme_bw()

dat %>%
  cor_test(k, age, method = "kendall")


rlm1 <- lqs(alpha ~ age, data = dat %>%
              filter(alpha != 5, age <= 34), 
            method = "lms")

rlm1

rlm1 <- rlm(alpha ~ age + I(age^2), data = dat %>%
              filter(alpha != 5), 
            method = "M", maxit = 50)

summary(rlm1)

lm1 <- lm(alpha ~ age, data = dat %>%
              filter(alpha != 5) %>%
            mutate(age2 = (age-mean(age))*(age-mean(age))))

summary(lm1)

# correlation matrix

cor_matrix <- cor(dat %>%
      drop_na() %>%
      filter(age <= 35) %>%
      dplyr::select(k, age, age_phys, age_psyc, minik,
             SES, FTP, alpha, beta, health_phys,
             health_psyc), method = c("kendall"))

corrplot(cor_matrix, order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

dat %>%
  cor_test(age_phys, FTP, method = "kendall")

lm1 <- dat %>%
  filter(age <= 50) %>%
  mutate(age_12 = ifelse(age <= 35, 0, 1)) %>%
  rlm(k ~ age_12 + alpha + age_12:alpha, data = ., method = "MM")

summary(lm1)

# anova

t1way(formula = logk ~ age_group, data = dat %>%
        mutate(logk = log(k)),
      tr = 0.2)

t1way(formula = minik ~ age_group, data = dat,
      tr = 0.2)

dat %>%
  filter(age_group == "young") %>%
  cor_test(k, age, method = "kendall")


fit <- rfit(k ~ age, data = dat %>%
              filter(age_group == "young"))

summary(fit)

# Results

## age effect - U shape

dat %>%
  group_by(gender, age_group) %>%
  summarise(count = n())

dat %>%
  filter(age_group == "young") %>%
  mutate(logk = log(k)) %>%
  cor_test(logk, age, method = "pearson")

dat %>%
  filter(age_group == "young") %>%
  cor_test(k, age, method = "kendall")

dat %>%
  filter(age_group == "old", gender == 1) %>%
  mutate(logk = log(k)) %>%
  cor_test(logk, age_psyc, method = "pearson")

dat %>%
  filter(age_group == "old") %>%
  mutate(logk = log(k)) %>%
  cor_test(logk, age, method = "pearson")

dat %>%
  filter(age_group == "middle") %>%
  mutate(logk = log(k)) %>%
  cor_test(logk, age, method = "kendall")

dat %>%
  filter(age_group == "old") %>%
  cor_test(k, age_psyc, method = "kendall")

dat %>%
  filter(!is.na(k)) %>%
  group_by(age_group) %>%
  summarise(k_median = mest(k),
            k_mad = mad(k))

anova_test(data = dat,
           k ~ age_group)

t1way(k ~ age_group, data = dat,
      tr = 0.05)

lincon(k ~ age_group, data = dat, tr = 0.05)

## age_phys

dat %>%
  filter(age_group == "young") %>%
  cor_test(k, age_phys, method = "kendall")

dat %>%
  filter(age_group == "old") %>%
  cor_test(k, age_phys, method = "kendall")

## age_psyc

dat %>%
  filter(age_group == "young") %>%
  cor_test(k, age_phys, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(beta, age_psyc, method = "kendall")


## minik

anova_test(data = dat,
           minik ~ age_group)

dat %>%
  group_by(age_group) %>%
  summarise(across(c(minik, FTP), list(mean = mean, sd = sd)))

dat %>%
  cor_test(minik, FTP, method = "kendall")

dat %>%
  filter(age_group == "young") %>%
  cor_test(minik, k, method = "kendall")

dat %>%
  filter(age_group == "middle") %>%
  cor_test(minik, k, method = "kendall")

dat %>%
  filter(age_group == "old") %>%
  cor_test(minik, k, method = "kendall")

dat %>%
  mutate(logk = log(k)) %>%
  cor_test(minik, SES, method = "kendall")

dat %>%
  cor_test(FTP, k, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(minik, alpha, method = "kendall") # no association

dat %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik, data = .) %>%
  summary()

dat %>%
  filter(age_group == "young") %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik + age + SES, data = .) %>%
  summary()

dat %>%
  filter(age_group == "old") %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik + age + SES, data = .) %>%
  summary()

dat %>%
  filter(age_group != "middle") %>%
  mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik + age + age_condition + age:age_condition, data = .) %>%
  summary()

dat %>%
  mutate(logk = log(k)) %>%
  rlm(logk ~ minik + age + age_group + age:age_group, data = .) %>%
  summary() # not significant

rlm <- rlm(logk ~ minik + age + age_condition + age:age_condition, 
           data = dat %>%
                              filter(age_group != "middle") %>%
                              mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
                              mutate(logk = log(k)), 
                            psi = psi.hampel, method = "M")

summary(rlm)

f.robftest(rlm, var = -1)
f.robftest(rlm, var = "minik")
f.robftest(rlm, var = "age")
f.robftest(rlm, var = "age_condition")
f.robftest(rlm, var = "age:age_condition")

rlm <- lmrob(logk ~ minik + age + age_condition + age:age_condition, 
           data = dat %>%
             filter(age_group != "middle") %>%
             mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
             mutate(logk = log(k)))

summary(rlm)

dat %>%
  filter(age_group != "middle") %>%
  mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
  mutate(logk = log(k)) %>%
  filter(!is.na(logk))

rlm_young <- lmrob(logk ~ age, 
             data = dat %>%
               filter(age_group == "young") %>%
               mutate(logk = log(k)) %>%
               filter(!is.na(logk)))

summary(rlm_young)

rlm_middle <- lmrob(logk ~ age, 
             data = dat %>%
               filter(age_group == "middle") %>%
               mutate(logk = log(k)))

summary(rlm_middle)

rlm_old <- lmrob(logk ~ age, 
                    data = dat %>%
                      filter(age_group == "old") %>%
                      mutate(logk = log(k)))

summary(rlm_old)

# P1 <- ggplot(dat %>%
#                filter(age_group != "middle") %>%
#                mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
#                mutate(logk = log(k)) %>%
#                filter(age_condition == 0), 
#              aes(x = age, y = logk)) +
#   geom_point(aes(color = minik),
#     position = position_jitter(width = .8, height = 1),
#              alpha = 0.8,
#              shape = 16,
#              size = 3) +
#   scale_colour_gradient(low = "white", high = "#2c7fb8") +
#   stat_smooth(method = lm,
#               se = FALSE,
#               size = 1.25,
#               color = "#2c7fb8") +
#   theme_bw() +  
#   labs(x = "Chronological Age",
#        y = "Delay Discounting Rate (log k)",
#        color = "Mini-K") +
#   theme(axis.title.x = element_text(size = 10, face = "bold"),
#         axis.title.y = element_text(size = 10, face = "bold"),
#         legend.title = element_text(size = 9, face = "bold"),
#         legend.text = element_text(size = 8.5),
#         legend.position = "top",
#         legend.box = "horizontal")
# 
# P1

curve_young <- function(x) -0.15359*x - 0.64163
curve_middle <- function(x) -0.05758*x - 1.33934
curve_old <- function(x) 0.16703*x - 13.87336

P1 <- ggplot(dat %>%
               filter(age_group == "young") %>%
               mutate(logk = log(k)) %>%
               filter(!is.na(logk)), 
             aes(x = age, y = logk)) +
  geom_point(aes(color = minik),
             position = position_jitter(width = .8, height = 0.5),
             alpha = 0.8,
             shape = 16,
             size = 3) +
  scale_colour_gradient(low = "white", high = "#2c7fb8") +
  stat_function(fun = curve_young, color = "#2c7fb8", size = 1.25) +
  theme_bw() +  
  labs(x = "Chronological Age",
       y = "Delay Discounting Rate (log k)",
       color = "Mini-K") +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8.5),
        legend.position = "top",
        legend.box = "horizontal")

P1

P2 <- ggplot(dat %>%
               filter(age_group == "old") %>%
               mutate(logk = log(k)) %>%
               filter(!is.na(logk)), 
             aes(x = age, y = logk)) +
  geom_point(aes(color = minik),
    position = position_jitter(width = .8, height = 1),
             alpha = 0.8,
             shape = 16,
             size = 3) +
  stat_function(fun = curve_old, color = "#dd1c77", size = 1.25) +
  scale_colour_gradient(low = "white", high = "#dd1c77") +
  theme_bw() +  
  labs(x = "Chronological Age",
       y = "Delay Discounting Rate (log k)",
       color = "Mini-K") +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8.5),
        legend.position = "top",
        legend.box = "horizontal")

P2

cowplot::plot_grid(P1, P2, nrow = 1)

P3 <- ggplot(dat %>%
               filter(age_group == "middle") %>%
               mutate(logk = log(k)) %>%
               filter(!is.na(logk)), 
             aes(x = age, y = logk)) +
  geom_point(aes(color = minik),
             position = position_jitter(width = .8, height = 1),
             alpha = 0.8,
             shape = 16,
             size = 3) +
  stat_function(fun = curve_middle, color = "#03808a", size = 1.25) +
  scale_colour_gradient(low = "white", high = "#03808a") +
  theme_bw() +  
  labs(x = "Chronological Age",
       y = "Delay Discounting Rate (log k)",
       color = "Mini-K") +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8.5),
        legend.position = "top",
        legend.box = "horizontal")

P3

cowplot::plot_grid(P1, P3, P2, nrow = 1)

png("age_effect_m_est.png", units = "in",
    width = 8, height = 5, res = 700)

cowplot::plot_grid(P1, P3, P2, nrow = 1)

dev.off()

dat_res_minik <- dat %>%
  drop_na() %>%
  select(ID, k, minik, age_group)

k_on_lhs <- dat_res_minik %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik, data = .)

res_after_lhs <- k_on_lhs$residuals %>% unname() %>% as_tibble()

dat_res_minik <- bind_cols(dat_res_minik, res_after_lhs)

dat_res_minik %>%
  group_by(age_group) %>%
  summarise(mean = mean(value),
            sd = sd(value))

dat_res_minik %>%
  t1way(value ~ age_group, data = .,
        tr = 0.2)
  
ancova_model <- aov(logk ~ minik + age_group, data = dat %>%
                      mutate(logk = log(k)))

Anova(ancova_model, type="III")

dat %>%
  mutate(logk = log(k)) %>%
  anova_test(data = .,
             logk ~ minik + age_group) %>%
  get_anova_table()

pwc <- dat %>%
  mutate(logk = log(k)) %>% 
  emmeans_test(
    logk ~ age_group, covariate = minik,
    p.adjust.method = "bonferroni"
  )

get_emmeans(pwc)

dat %>%
  mutate(logk = log(k)) %>%
  drop_na() %>%
  group_by(age_group) %>%
  summarise(across(logk, list(mean, sd)))

## alpha beta

dat %>%
  group_by(age_group) %>%
  summarise(across(c(alpha, beta), list(m = mest, mad = mad)))

dat %>%
  filter(age_group == "young", alpha != 5, beta <= 1.323433) %>%
  cor_test(alpha, age, method = "kendall")

dat %>%
  filter(age_group == "young", alpha != 5, beta <= 1.323433) %>%
  cor_test(beta, age, method = "kendall")


dat %>%
  filter(age_group == "middle", alpha != 5, beta <= 1.323433) %>%
  cor_test(alpha, age, method = "kendall")

dat %>%
  filter(age_group == "middle", alpha != 5, beta <= 1.323433) %>%
  cor_test(beta, age, method = "kendall")


dat %>%
  filter(age_group == "old", alpha != 5, beta <= 1.323433) %>%
  cor_test(alpha, age, method = "kendall")

dat %>%
  filter(age_group == "old", alpha != 5, beta <= 1.323433) %>%
  cor_test(beta, age, method = "kendall")

## alpha\beta and k
dat %>%
  mutate(logk = log(k)) %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(beta, minik, method = "kendall") # 不显著

dat %>%
  filter(alpha != 5) %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ alpha + age, data = .) %>% # alpha 预测k
  summary()

dat %>%
  mutate(logk = log(k)) %>%
  filter(alpha != 5, beta <= 2, age_group == "young") %>%
  lmrob(logk ~ alpha, data = .) %>% # alpha 预测k
  summary()

dat %>%
  mutate(logk = log(k)) %>%
  filter(alpha != 5, beta <= 2, age_group == "middle") %>%
  lm(logk ~ age, data = .) %>% # alpha 预测k
  summary()

## t test 

yuen(beta ~ age_group, data = dat %>%
       filter(age_group != "old"))

## FTP
dat %>%
  anova_test(FTP ~ age_group)

dat %>%
  group_by(age_group) %>%
  summarise(across(c(FTP), list(m = mest, mad = mad)))

dat %>%
  group_by(age_group) %>%
  cor_test(FTP, age, method = "kendall")

dat %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(FTP, alpha, method = "kendall") #不显著

dat %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(FTP, beta, method = "kendall")

dat %>%
  cor_test(FTP, age_phys, method = "kendall")

dat %>%
  cor_test(FTP, minik, method = "pearson")

dat %>%
  group_by(age_group) %>%
  cor_test(FTP, age, method = "kendall")

dat %>%
  lm(FTP ~ age + minik, data = .) %>% #都显著
  summary()

dat %>%
  filter(age_group == "young") %>%
  cor_test(FTP, age, method = "kendall")

dat %>%
  filter(age_group == "middle") %>%
  cor_test(FTP, age, method = "kendall")

dat %>%
  filter(age_group == "old") %>%
  cor_test(FTP, age, method = "kendall")

## mediation

dat %>%
  cor_test(k, SES, method = "kendall")

dat %>%
  cor_test(minik, k, method = "kendall")


fit.mx <- lm(minik ~ SES, data = dat %>%
               drop_na() %>%
               mutate(logk = log(k)))
fit.yxm <- lm(logk ~ minik + SES, data = dat %>%
                drop_na() %>%
                mutate(logk = log(k)))
set.seed(123)
fitmed <- mediation::mediate(fit.mx, fit.yxm, treat = "SES",
                             mediator = "minik", sims = 999, boot = TRUE, boot.ci.type = "bca")
summary(fitmed)

robu_mediation <- with(dat %>%
       drop_na() %>%
       mutate(logk = log(k)), ZYmediate(SES, logk,
                        minik, nboot = 1999))
robu_mediation
robu_mediation$a.est
robu_mediation$b.est
robu_mediation$ab.est
robu_mediation

robu_mediation2 <- with(dat %>%
                          mutate(logk = log(k)) %>%
                          filter(age_group == "young"), ZYmediate(age, logk, FTP, nboot = 1999))
robu_mediation2
robu_mediation2$a.est
robu_mediation2$b.est
robu_mediation2$ab.est
robu_mediation2

dat %>%
  mutate(across(c(SES, minik), scale)) %>%
  lm(SES ~ minik, data = .) %>%
  summary()

## health

dat %>%
  cor_test(k, health_phys, method = "kendall")

dat %>%
  cor_test(age, health_phys, method = "kendall")

dat %>%
  mutate(logk = log(k)) %>%
  filter(age_group == "old") %>%
  lm(logk ~ age + health_phys + age:health_phys, data = .) %>%
  summary()

dat %>%
  cor_test(minik, SES, method = "kendall")

## provide support for Sozou

dat %>%
  filter(age_group != "middle") %>%
  mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
  mutate(logk = log(k)) %>%
  lm(logk ~ minik + age_phys + age_condition + age_phys:age_condition, data = .) %>%
  summary()

dat %>%
  filter(age_group != "middle") %>%
  mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
  mutate(logk = log(k)) %>%
  pivot_longer(cols = c(age, age_phys),
               names_to = "age_type",
               values_to = "age_value") %>%
  mutate(across(c(age_value, minik), scale, scale = F)) %>%
  rlm(logk ~ minik + age_value + age_condition + age_type +
       age_value:age_condition + age_value:age_type + age_condition:age_type + 
       age_value:age_condition:age_type, 
      data = .,
      method = "M") %>%
  summary()

dat %>%
  mutate(logk = log(k)) %>%
  pivot_longer(cols = c(age, age_phys),
               names_to = "age_type",
               values_to = "age_value") %>%
  mutate(across(c(age_value, minik), scale)) %>%
  lm(logk ~ minik + age_value + age_group + age_type +
        age_value:age_group + age_value:age_type + age_group:age_type + 
        age_value:age_group:age_type, 
      data = .) %>%
  summary()

compKendall <- function(data, indices){
  df <- data[indices, ] 
  ks <- sapply(2:3, function(x) cor(df[ , 1], df[ , x], method = "kendall"))
  kPairs <- combn(ks, 2)
  
  return(kPairs[1, ] - kPairs[2, ])
}

df <- dat %>%
  filter(age_group == "young", !is.na(k)) %>%
  dplyr::select(k, age, age_phys)

kDiff <- boot(df, statistic = compKendall, R = 1999)

broom::tidy(kDiff, conf.int = TRUE, conf.method = "perc")

# compKendall <- function(data, indices){
#   
#   d1 <- data %>% filter(gender == 1, age_group == "old")
#   df1 <- d1[indices, ]
#   d2 <- data %>% filter(gender == 2, age_group == "old")
#   df2 <- d2[indices, ]
#   
#   cor1 <- df1 %>% cor_test(age, k, method = "kendall")
#   cor2 <- df2 %>% cor_test(age, k, method = "kendall")
#   cor_d <- cor1$cor - cor2$cor
#   return(cor_d)
# }
#   
# 
# kDiff <- boot(dat, statistic = compKendall, R = 999)
# 
# broom::tidy(kDiff, conf.int = TRUE, conf.method = "bca")

## EDA1

lm_eda2 <- dat %>%
  mutate(logk = log(k)) %>%
  filter(alpha != 5, beta <= 2, age_group == "young") %>%
  rlm(logk ~ alpha, data = ., psi = psi.hampel, method = "M")

lm_eda2

f.robftest(lm_eda2, var = "alpha")

rlm <- rlm(logk ~ minik + age + age_condition + age:age_condition, 
           data = dat %>%
             filter(age_group != "middle") %>%
             mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
             mutate(logk = log(k)), 
           psi = psi.hampel, method = "M")

summary(rlm)

## EDA 2

dat %>%
  filter(alpha != 5, beta <= 2) %>%
  cor_test(minik, alpha, method = "kendall")

## EDA3

dat %>%
  cor_test(age, age_phys, method = "kendall")

dat %>%
  mutate(age_bias_phys = age_phys - age) %>%
  group_by(age_group) %>%
  cor_test(k, age_bias_phys, method = "kendall")

dat %>%
  mutate(age_bias_phys = age_phys - age,
         logk = log(k)) %>%
  filter(age_group == "old") %>%
  lm(logk ~ minik + age + health_phys, data = .) %>%
  summary()

dat %>%
  mutate(age_bias_phys = age_phys - age) %>%
  filter(age_group == "young") %>%
  cor_test(k, age_bias_phys, method = "kendall")

dat %>%
  mutate(age_bias_phys = age_phys - age) %>%
  group_by(age_group) %>%
  summarise(across(c(age_bias_phys, age, age_phys), list(mean = mean,
                                             sd = sd)))

rlm <- rlm(logk ~ minik + age_bias_phys + age_condition + age_bias_phys:age_condition, 
           data = dat %>%
             filter(age_group != "middle") %>%
             mutate(age_condition = ifelse(age_group == "young", 0, 1)) %>%
             mutate(logk = log(k), age_bias_phys = age_phys - age), 
           psi = psi.hampel, method = "M")

summary(rlm)

f.robftest(rlm, var = "minik")
f.robftest(rlm, var = "age")
f.robftest(rlm, var = "age_condition")
f.robftest(rlm, var = "age_bias_phys:age_condition")


dat %>%
  mutate(age_bias_phys = age_phys - age,
         logk = log(k)) %>%
  filter(age_group == "old") %>%
  rlm(logk ~ minik + age + age_bias_phys, data = .) %>%
  summary()

dat %>%
  mutate(age_bias_phys = age_phys - age,
         age_bias_psyc = age_psyc - age,
         logk = log(k)) %>%
  filter(age_group == "old") %>%
  cor_test(age_bias_phys, minik, method = "kendall")

rlm2 <-  lm(age_bias_phys ~ minik + age + health_phys, 
             data = dat %>% mutate(age_bias_phys = age_phys - age))

summary(rlm2)


rlm2 <-  lmrob(logk ~ minik + age + age_bias_phys + age_condition + age_condition:age_bias_phys + age:age_condition, 
             data = dat %>%
               mutate(age_bias_phys = age_phys - age,
                      logk = log(k)) %>%
               filter(age_group != "middle") %>%
               mutate(age_condition = ifelse(age_group == "young", 0, 1)))

summary(rlm2)

f.robftest(rlm2, var = "age_bias_phys:age_condition")

dat %>%
  mutate(age_bias_phys = age_phys - age,
         logk = log(k)) %>%
  filter(age_group != "middle") %>%
  lm(health_phys ~ age + minik, data = .) %>%
  summary()

dat %>%
  group_by(age_group) %>%
  summarise(across(c(minik, health_phys, FTP), list(mean = mean, sd = sd)))

dat %>%
  anova_test(health_phys ~ age_group + minik)

summary(dat$minik)

minik_pts <- seq(4.7, 5.75, length.out = 5)

ancova(health_phys ~ age_group + minik, fr1 = 0.3, fr2 = 0.3, pts = minik_pts,
       data = dat %>% filter(age_group != "middle") %>% mutate(age_group = factor(age_group)))

ancova(logk ~ age_group + minik, fr1 = 0.3, fr2 = 0.3, pts = minik_pts,
       data = dat %>% filter(age_group != "old") %>% mutate(age_group = factor(age_group),
                                                               logk = log(k)))

ancova(FTP ~ age_group + minik, fr1 = 0.3, fr2 = 0.3, pts = minik_pts,
       data = dat %>% filter(age_group != "old") %>% mutate(age_group = factor(age_group),
                                                            logk = log(k)))

aov(logk ~ age_group + minik, dat %>% 
      filter(age_group != "young") %>% 
      mutate(age_group = factor(age_group), logk = log(k))) %>%
  summary()


dat %>%
  group_by(age_group) %>%
  summarise(across(SES, list(mean = mean, sd = sd)))
