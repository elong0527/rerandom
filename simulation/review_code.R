# ------------------------------------------------------- #
#    Curve of u(alpha, L), l(alpha, L), p_l, p_u over L   #
# ------------------------------------------------------- #
library(dplyr)
library(ggplot2)
tbl <- tibble::tribble(~L, ~m_lower, ~m_upper, ~p_lower, ~p_upper,
                       1000, 0, 6, 0.000000, 0.006000,
                       2000, 0, 6, 0.000000, 0.003000,
                       3000, 0, 7, 0.000000, 0.002333,
                       4000, 0, 7, 0.000000, 0.001750,
                       5000, 0, 7, 0.000000, 0.001400,
                       10000, 0, 8, 0.000000, 0.000800,
                       50000, 1, 15, 0.000020, 0.000300,
                       100000, 4, 22, 0.000040, 0.000220,
                       500000, 31, 76, 0.000062, 0.000152,
                       1000000, 70, 138, 0.000070, 0.000138,
                       2000000, 151, 258, 0.000076, 0.000129,
                       3000000, 234, 376, 0.000078, 0.000125,
                       4000000, 318, 492, 0.000080, 0.000123,
                       5000000, 403, 608, 0.000081, 0.000122,
                       6636000, 543, 796, 0.000082, 0.000120)

# l, u vs L
setEPS()
postscript("~/Rtemp/test_temp/L_vs_u_l.eps")
ggplot(data = tbl %>%
         select(L, m_lower, m_upper) %>%
         tidyr::pivot_longer(cols = 2:3, names_to = "Label"),
       aes(x = log(L), y = log(value), group = Label, color = Label)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
  ) +
  scale_color_discrete(labels = c(expression(l(alpha, L)), expression(u(alpha, L)))) +
  xlab(expression(log(L))) +
  ylab("Logarithm of the upper and lower bounds of the L")
dev.off()

# p_l, p_u vs L
setEPS()
postscript("~/Rtemp/test_temp/L_vs_pu_pl.eps")
ggplot(data = tbl %>%
         select(L, p_lower, p_upper) %>%
         tidyr::pivot_longer(cols = 2:3, names_to = "Label"),
       aes(x = log(L), y = value * 1000, group = Label, color = Label)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 20),
  ) +
  scale_color_discrete(labels = c(expression(p[l]), expression(p[u]))) +
  xlab(expression(log(L))) +
  ylab("Values of the lower/upper bound of the observed p-value * 1000")
dev.off()

# ------------------------------------------------------- #
#               Signal to noise of l_i                    #
# ------------------------------------------------------- #
library(rerandom)
library(dplyr)
library(survival)
library(gt)
set.seed(123)

fit_rerandom <- function(ana){
  
  fit_binary <- glm(binary ~ treatment + ecog + tmb, family = "binomial", data = ana)
  fit_continous <- lm(continous ~ treatment + ecog + tmb, data = ana)
  fit_survival <- survdiff(Surv(surv_time, surv_status) ~ treatment + strata(ecog) + strata(tmb), data=ana)
  
  binary <- list(est = summary(fit_binary)$coeff[2, 1], 
                 sd = summary(fit_binary)$coeff[2, 2], 
                 z = summary(fit_binary)$coeff[2, 3], 
                 wald_p = summary(fit_binary)$coeff[2, 4])
  
  continous <- list(est = summary(fit_continous)$coeff[2, 1], 
                    sd = summary(fit_continous)$coeff[2, 2], 
                    z = summary(fit_continous)$coeff[2, 3], 
                    wald_p = summary(fit_continous)$coeff[2, 4])
  
  
  survival <- list(est = NULL,
                   sd = NULL, 
                   z = sqrt(fit_survival$chisq), 
                   wald_p = 1 - pchisq(fit_survival$chisq, df = 1) )
  
  dplyr::bind_rows(binary = binary, 
                   continous = continous, 
                   survival = survival, .id = "method")
  
}

# set the number of subjects
n <- 100

# set the number of repetitions
n_rerandom <- 1000

# define 3 stratum
stratum <- define_stratum(
  site = rep(1, 50),
  ecog = c("0" = 1, "1" = 1),
  tmb = c("<=6" = 1, ">6 and <= 12" = 1, ">12" = 1) )

# define 2 arms
treatment <- c("drug", "placebo")

# randomize subjects into 2 arms
df <- n %>% 
  simu_stratum(stratum = stratum) %>%
  simu_treatment_minimization(treatment = treatment, 
                              prob = define_prob(0.9, treatment), 
                              ratio = c(1, 1),
                              imbalance_fun = imbalance_fun_range)

ans <- NULL
for (i in 1:1e4) {
  set.seed(i + 1)
  
  beta <- list(
    site = rnorm(n = 50, 0.2),
    ecog = c(0.3, 0.5), 
    tmb = c(-0.3, 0, 0.3),
    treatment = c(0.2, 0.8))
  
  l <- beta$site[as.numeric(df$site)] + 
    beta$ecog[as.numeric(df$ecog)] + 
    beta$tmb[as.numeric(df$tmb)] + 
    beta$treatment[as.numeric(df$treatment)]
  
  if(i %in% c(9, 49, 99, 499)){
    postscript(file.path(paste0("./l_hist_seed", i + 1, ".eps")), width = 5, height = 5, horizontal = FALSE)
    ggplot(data.frame(l = l), aes(x = l)) +
      geom_histogram(fill = "grey", color = "black", bins = 8) +
      labs(x = "l", y = "Percentage") +
      theme(axis.text.x = element_text(size = 60),
            axis.text.y = element_text(size = 60)) +
      theme_bw()
    dev.off()
  }
  
  R2 <- (lm(l ~ beta$site[as.numeric(df$site)]) %>% summary())$r.squared
  p <- exp(l) / (1 + exp(l))
  
  ans_new <- data.frame(mean_l = mean(l),
                        l_ge1 = sum(abs(l) > 1),
                        l_ls02 = sum(abs(l) < 0.5),
                        R2 = R2)
  ans <- union_all(ans, ans_new)
  print(i)
}

apply(ans, 2, mean)

# figure: histogram of percentage l_1 > 1
postscript("./l_ge_1_hist.eps", width = 5, height = 5, horizontal = FALSE)
ggplot(data.frame(l_ge1 = ans$l_ge1 / length(l)), aes(x = l_ge1)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 fill = "grey", color = "black", bins = 15) +
  labs(x = "Pergentage of elements of l greater than 1", y = "Percentage") +
  theme(axis.text.x = element_text(size = 80),
        axis.text.y = element_text(size = 80)) +
  theme_bw()
dev.off()

# figure: histogram of R2
postscript("./r2_hist.eps", width = 5, height = 5, horizontal = FALSE)
ggplot(data.frame(R2 = ans$R2), aes(x = R2)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 fill = "grey", color = "black", bins = 15) +
  labs(x = expression(R^2), y = "Percentage") +
  theme(axis.text.x = element_text(size = 80),
        axis.text.y = element_text(size = 80)) +
  theme_bw()
dev.off()

# figure: histogram of percentage l_i < 0.5
postscript("./l_ls_05_hist.eps", width = 5, height = 5, horizontal = FALSE)
ggplot(data.frame(l_ls02 = ans$l_ls02), aes(x = l_ls02)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 fill = "grey", color = "black", bins = 15) +
  labs(x = "Pergentage of elements of l smaller than 0.5", y = "Percentage") +
  theme(axis.text.x = element_text(size = 80),
        axis.text.y = element_text(size = 80)) +
  theme_bw()
dev.off()

# ------------------------------------------------------- #
#        running time of 10 repetition and                #
#   calculate the average time spend for each repletion   #
# ------------------------------------------------------- #
library(rerandom)
library(dplyr)
library(survival)
library(tictoc)

# simulation setup 
n <- 200
n_rerandom <- 10

tic()
stratum <- define_stratum(
  site = rep(1, 50),
  ecog = c("0" = 1, "1" = 1),
  tmb = c("<=6" = 1, ">6 and <= 12" = 1, ">12" = 1) )

treatment <- c("drug", "placebo")

df <- n %>% 
  simu_stratum(stratum = stratum) %>%
  simu_treatment_minimization(treatment = treatment, 
                              prob = define_prob(0.9, treatment), 
                              ratio = c(1, 1),
                              imbalance_fun = imbalance_fun_range)

beta <- list(
  site = rnorm(n = 50, 0.2),
  ecog = c(0.3, 0.5), 
  tmb = c(-0.3, 0, 0.3),
  treatment = c(0.2, 0.8))

l <- beta$site[as.numeric(df$site)] + 
  beta$ecog[as.numeric(df$ecog)] + 
  beta$tmb[as.numeric(df$tmb)] + 
  beta$treatment[as.numeric(df$treatment)]

p <- exp(l) / (1 + exp(l)) 

c_min <- 0
c_max <- 5

time <- exp(l + log(- log(runif(n)))) # PH model based on transformation model
cen <- runif(n, min = c_min, max = c_max) # Uniform censoring  

df$binary <- rbinom(n = n, size = 1, prob = p)
df$continous <- l + rnorm(n)
df$surv_time <- pmin(time, cen)
df$surv_status <- time < cen

# analysis of single data
fit_rerandom <- function(ana){
  
  fit_binary <- glm(binary ~ treatment + ecog + tmb, family = "binomial", data = ana)
  fit_continous <- lm(continous ~ treatment + ecog + tmb, data = ana)
  fit_survival <- survdiff(Surv(surv_time, surv_status) ~ treatment + strata(ecog) + strata(tmb), data=ana)
  
  binary <- list(est = summary(fit_binary)$coeff[2, 1], 
                 sd = summary(fit_binary)$coeff[2, 2], 
                 z = summary(fit_binary)$coeff[2, 3], 
                 wald_p = summary(fit_binary)$coeff[2, 4])
  
  continous <- list(est = summary(fit_continous)$coeff[2, 1], 
                    sd = summary(fit_continous)$coeff[2, 2], 
                    z = summary(fit_continous)$coeff[2, 3], 
                    wald_p = summary(fit_continous)$coeff[2, 4])
  
  
  survival <- list(est = NULL,
                   sd = NULL, 
                   z = sqrt(fit_survival$chisq), 
                   wald_p = 1 - pchisq(fit_survival$chisq, df = 1) )
  
  dplyr::bind_rows(binary = binary, 
                   continous = continous, 
                   survival = survival, .id = "method")
}

# observed data
res <- fit_rerandom(df)

# re-randomization 
res_rerandom <- replicate(n_rerandom, {
  
  df_rerandom_trt <- df[, c("site", "ecog", "tmb")] %>%
    simu_treatment_minimization(treatment = treatment, 
                                prob = define_prob(0.9, treatment), 
                                ratio = c(1, 1),
                                imbalance_fun = imbalance_fun_range)
  df_rerandom <- df
  df_rerandom$treatment <- df_rerandom_trt$treatment
  
  tmp <- fit_rerandom(df_rerandom)
  tmp$res_z <- res$z
  tmp$res_p <- res$wald_p
  tmp
}, simplify = FALSE)

res_rerandom <- dplyr::bind_rows(res_rerandom, .id = "id")

res1 <- res_rerandom %>% 
  group_by(method) %>% 
  summarise(sd_200 = sd(z[1:200]), 
            full_p = mean(abs(z) > abs(res_z)) ) 

res2 <- res %>% 
  left_join(res1) %>% 
  mutate(parat_p = 2*(1 - pt(q = z / sd_200, df = 200 - 1)), 
         paran_p = 2*(1 - pnorm(z, sd = sd_200)) )
toc()
