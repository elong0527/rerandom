# To install the SimuExample library from the source. 
# 1. Open the .Rproj file using Rstudio
# 2. Run the code devtools::install() in the Rstudio console

# Save all functions in the R folder with proper documentation using usethis::use_r() 
# Save all required dataset in the data folder using usethis::use_data()
library(rerandom)
library(dplyr)
library(survival)
# Get task id for each simulation job (not run in Rstudio Serve)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Set up Simulation Enviroment

# task_id <- 1  # used for debugging purpose
set.seed(task_id)

#---------- Simulation Setup --------------------#

n <- 200
n_rerandom <- 1e+6

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
  treatment = c(0.2, 0.8)
)

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

# Analysis of single data

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

res2

save.image(file = paste0(task_id, ".Rdata"))


