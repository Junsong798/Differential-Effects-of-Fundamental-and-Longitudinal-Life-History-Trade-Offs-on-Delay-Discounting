# Data Cleaning

# Author: Amadeus
# Version: 2021-11-15

# Libraries
library(drc)
library(WRS2)
library(nlme)
library(aomisc)
library(rstatix)
library(tidyverse)


# Sources
dat <- read_csv("total.csv")

# Parameters

# ============================================================================

# Code

## scale

dat <- dat %>%
  mutate(across(Item19:Item21, ~ recode(.,
                                        '1' = 7,
                                        '2' = 6,
                                        '3' = 5,
                                        '4' = 4,
                                        '5' = 3,
                                        '6' = 2,
                                        '7' = 1))) %>%
  mutate(FTP = rowMeans(select(., Item12:Item21)),
         minik = rowMeans(select(., Item42:Item61)),
         SES = rowMeans(select(., Item62:Item64))) %>%
  mutate(ID = 1:242) %>%
  rowwise() %>%
  mutate(ss_count = sum(c_across(Item5:Item11) == 1))

## transformation of subjective time

dat <- dat %>%
  mutate(across(Item22:Item31, ~ ifelse(. <= 100, .,
                                         (. - 42648)/(87551 - 42648)*100)))

## rename variables

names(dat) <- c("gender", "age", 
                "age_phys", "age_psyc",
                paste0("DD", 1:7),
                paste0("FTP", 1:10),
                paste0("ST", 1:10),
                paste0("DDT", 1:10),
                paste0("minik", 1:20),
                paste0("SES", 1:3),
                "health_phys",
                "health_psyc",
                "platform",
                "FTP", "minik", "SES", "ID", "ss_count")


## calculate k values

dat_dd <- dat %>%
  select(ID, DD1:DD7)

k_item <- function(V, A, D){
  k <- (A - V)/(V*D)
  return(k)
}

D <- c(31, 4, 939, 183, 527, 25, 184)
A <- c(3480, 8190, 2311, 4620, 4200, 5160, 2760) # LL
V <- c(840, 2730, 1680, 3570, 3990, 3150, 1260) # SS


k_items7 <- k_item(V, A, D)

k_items7_order <- k_items7[order(k_items7)]

dat_dd_order <- dat_dd %>%
  select(ID, c(order(k_items7) + 1)) %>%
  as_tibble()

### MLE

# nll <- function(k, mu) {
#   
#   prob_SV1 <- vector("numeric", length = 7)
#   f <- vector("numeric", length = 7)
#   llvalue <- vector("numeric", length = 7)
#   
#   for (i in 1:7) {
#     
#     SV1 <- V[i]
#     SV2 <- A[i]/(1 + k*D[i])
#     
#     prob_SV1[i] <- 1/(1 + exp(-(SV1 - SV2)*mu))
#     
#     if(data[i] == 1){
#       f[i] <- prob_SV1[i]
#     }else{
#       f[i] <- 1 - prob_SV1[i]
#     }
#     
#   }
#   
#   -sum(log(f))
# }

nll2 <- function(pars, data) {
  
  k <- pars[1]
  mu <- pars[2]
  
  prob_SV1 <- vector("numeric", length = 7)
  f <- vector("numeric", length = 7)
  llvalue <- vector("numeric", length = 7)
  
  for (i in 1:7) {
    
    SV1 <- V[i]
    SV2 <- A[i]/(1 + k*D[i])
    
    prob_SV1[i] <- 1/(1 + exp(-(SV1 - SV2)*mu))
    
    if(data[i] == 1){
      f[i] <- prob_SV1[i]
    }else{
      f[i] <- 1 - prob_SV1[i]
    }
    
  }
  
  -sum(log(f))
}

# mle_res <- mle(nll, start = list(k = 0.1, mu = 0.02), method = "CG")
# mle_res
# 
# mle_res2 <- optim(par = c(k = 0.1, mu = 0.02), 
#                   fn = nll2, data = data,
#                   method = "CG",
#                   control = list(parscale = c(k = 0.1, mu = 0.02)))
# mle_res2$par


dat_dd_mle <- dat_dd %>%
  pivot_longer(cols = starts_with("DD"),
               names_to = "NO",
               values_to = "choice") %>%
  select(-NO) %>%
  nest(data = choice)

est_para <- function(data){
  
  data_unlist <- unlist(data) %>% unname()
  
  mle_res <- optim(par = c(k = 0.05, mu = 0.005), 
                   fn = nll2, data = data_unlist,
                   lower = c(0, 0),
                   upper = c(0.8, 1),
                   method = "L-BFGS-B",
                   control = list(parscale = c(k = 0.01, mu = 0.005)))
  
  return(list(k = unname(mle_res$par[1]),
              mu = unname(mle_res$par[2])))
  
}

edat_dd_mle <- dat_dd_mle %>%
  mutate(para = map(data, possibly(est_para, list(k = NA,
                                                  mu = NA)))) %>%
  unnest(para) %>%
  mutate(para_name = rep(c("k", "mu"), times = 242)) %>%
  pivot_wider(names_from = para_name,
              values_from = para)

dat <- left_join(dat, dat_dd_mle, by = "ID") %>%
  mutate(across(k:mu, unlist))


# subjective time estimates

# dat %>%
#   drop_na() %>%
#   filter(age < 35) %>%
#   cor_test(age, k, method = "kendall")
# 
# dat %>%
#   drop_na() %>%
#   cor_test(minik, k, method = "kendall")


ggplot(dat) +
  geom_histogram(aes(age)) +
  theme_bw()

dat <- dat %>%
  mutate(age_group = case_when(
    age <= 35 ~ 0,
    age > 35 & age <= 50 ~ 1,
    age > 50 ~ 2
  )) %>%
  mutate(age_group = factor(age_group,
                            levels = c(0, 1, 2),
                            labels = c("young",
                                       "middle",
                                       "old")))

# anova_test(dat,
#            ST2 ~ age_group,
#            type = 3)
# 
# mean(dat$ST1)
# 
# dat %>%
#   group_by(age_group) %>%
#   summarise(mean_sv = mean(ST9),
#             sd = sd(ST9))
mean(dat$ST1)

dat <- dat %>%
  mutate(across(ST1:ST10, ~ ./19.28291*3, .names = "{.col}_sub"))


OT <- c(3, 6, 9, 12, 15, 21, 24, 30, 36, 60)

dat_st <- dat %>%
  dplyr::select(ends_with("_sub"), ID) %>%
  pivot_longer(cols = starts_with("ST"),
               names_to = "NO",
               values_to = "sub_time") %>%
  mutate(obj_time = rep(OT, 242)) %>%
  dplyr::select(-NO) %>%
  nest(data = c(obj_time, sub_time))

est_time <- function(data){
  
  dat <- data %>% as.data.frame()
  
  model <- drm(sub_time ~ obj_time, fct = DRC.powerCurve(),
               type = "continuous",
               robust = "median",
               lowerl = c(0, 0),
               upperl = c(5, 5),
               data = dat)
  
  return(list(alpha = model$parmMat[1],
              beta = model$parmMat[2]))
  
}

dat_st <- dat_st %>%
  mutate(para_list = map(data, est_time)) %>%
  unnest(para_list) %>%
  mutate(para_name = rep(c("alpha", "beta"), times = 242)) %>%
  pivot_wider(names_from = para_name,
              values_from = para_list) %>%
  dplyr::select(-data)

dat <- dat %>%
  left_join(., dat_st, by = "ID") %>%
  mutate(across(alpha:beta, unlist))

# dat %>%
#   group_by(age_group) %>%
#   filter(!(alpha == 5)) %>%
#   summarise(across(c(alpha, beta), list(mean = mean, sd = sd)))
# 
# dat %>%
#   filter(age <= 40, alpha != 5) %>%
#   cor_test(alpha, age, method = "kendall")
# 
# dat %>%
#   lm(k ~ alpha + age + I(age^2), data = .) %>%
#   summary()

dat %>%
  group_by(age_group) %>%
  filter(alpha != 5) %>%
  summarise(across(c(alpha, beta), list(M = mest, sd = mad)))

dat_mest <- dat %>%
  group_by(age_group) %>%
  filter(alpha != 5) %>%
  summarise(across(ST1_sub:ST10_sub, list(M = median))) 

mest_long <- dat_mest %>%
  pivot_longer(cols = contains("sub_M"),
               names_to = "sub_delay",
               values_to = "value") %>%
  mutate(obj_delay = rep(OT, 3))

# obj_matrix <- matrix(rep(OT, times = 3), nrow = 3, byrow = TRUE) 
# colnames(obj_matrix) <- names(dat_mest)[2:11]
# 
# dat_obj_time <- obj_matrix %>%
#   as_tibble() %>%
#   mutate(age_group = c("young", "middle", "old")) %>%
#   select(age_group, everything())
# 
# dat_mest <- bind_rows(dat_mest, dat_obj_time)

ggplot(mest_long, aes(x = obj_delay, y = value, fill = age_group)) +
  geom_line() + 
  geom_point(size = 4.5, shape = 21) +
  theme_bw() +
  labs(x = "Objective Delay (Month)", 
       y = "Subjective Delay (Month)", 
       fill = "Age Group") +
  scale_fill_brewer()

F1 <- ggplot(mest_long %>%
               mutate(age_group = factor(age_group,
                                         levels = c("young",
                                                    "middle",
                                                    "old"))), 
             aes(x = obj_delay, y = value, fill = age_group)) +
  geom_line() + 
  geom_point(size = 2.75, shape = 21) +
  theme_bw() +
  labs(x = "Calendar Time (Month)", 
       y = "Subjective Time (Month)", 
       fill = "Age Group") +
  ylim(0, 16) +
  scale_fill_manual(values = c("#BDD8DA", "#FF9900", "#AD235E")) +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8.5))

obj_to_sub <- function(obj_delay, age_group){
  
  if(age_group == "young"){
    sub_delay <- 1.33*(obj_delay^0.593)
  }else if(age_group == "middle"){
    sub_delay <- 1.08*(obj_delay^0.614)
  }else{
    sub_delay <- 1.34*(obj_delay^0.578)
  }
  return(sub_delay)
}

dat_time_curve <- tibble(age_group = rep(c("young", "middle", "old"), each = 10),
                         obj_delay = rep(seq(1, 60, length.out = 10), times = 3),
                         sub_delay = map2_dbl(obj_delay, age_group, obj_to_sub))

F2 <- ggplot(dat_time_curve %>%
               mutate(age_group = factor(age_group,
                                         levels = c("young",
                                                    "middle",
                                                    "old"))), 
             aes(x = obj_delay, y = sub_delay, fill = age_group)) +
  geom_line() + 
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "Calendar Time (Month)", 
       y = "Subjective Time (Month)",
       fill = "Age Group") +
  ylim(0, 16) +
  scale_fill_manual(values = c("#BDD8DA", "#FF9900", "#AD235E")) +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8.5))

cowplot::plot_grid(F1, F2, nrow = 1)


# png("time_est.png", units = "in",
#     width = 8, height = 4, res = 700)
# 
# cowplot::plot_grid(F1, F2, nrow = 1)
# 
# dev.off()

write_csv(dat %>% select(-data), "processed_data.csv")



