#' Simulate design matrix based on stratum information
#' 
#' @param n a numeric value of sample size
#' @param stratum a list of stratum inforamtion created by `define_stratum`
#' 
#' @examples 
#' simu_stratum(n = 10)
#' simu_stratum(n = 10, 
#'              stratum = define_stratum(sex = c("M" = 1, "F" = 1)))
#'  
#' @export
simu_stratum <- function(n, stratum = NULL){
  
  if(is.null(stratum)) return(data.frame(stratum = factor(rep(1, n))))
  
  res <- lapply(stratum, function(x) sample(x$level, size = n, replace = TRUE, prob = x$weight))
  dplyr::bind_cols(res) # TODO: For production, replace to a base R implementation.
  
}

#' Simulate treatment group using block randomization
#' 
#' @param df a design matrix contain all stratum.
#' @param treatment a character vector of unique treatment group name. 
#' @param ratio a integer vector of randomization ratio.
#' @param blocksize a numeric value of number of subjects in each block.
#' 
#' @examples 
#' library(magrittr)
#' simu_stratum(n = 10, 
#'              stratum = define_stratum(sex = c("M" = 1, "F" = 1))) %>% 
#' simu_treatment_block(treatment = c("case", "control"), 
#'                      ratio = c(3,2),
#'                      block = 2)
#'  
#' @export
simu_treatment_block <- function(df, 
                                 treatment, 
                                 ratio = rep(1L, length(treatment)), 
                                 blocksize = length(treatment)*2){
  
  treatment <- factor(treatment, levels = unique(treatment))
  element <- rep(treatment, ratio)
  
  n <- nrow(df)
  
  x <- apply(df, 1, paste, collapse = "+") 
  x <- factor(x)
  x_loc <- as.numeric(x)
  n_loc <- table(x_loc)
  x_level <- levels(x)
  
  x_trt <- lapply(x_level, function(trt){
    trt <- replicate( max(n_loc) %/% blocksize + 1 , 
                      sample(rep(element, length.out = blocksize)), 
                      simplify = FALSE)
    unlist(trt)
  })
  
  trt <- rep(NA, n)
  for(i in seq_along(x_level)){
    trt[x_loc == i] <- x_trt[[i]][1:n_loc[i]]
  }
  
  df$treatment <- factor(trt, labels = levels(treatment))
  
  df
}


#' Simulate treatment group using minimization randomization
#' 
#' @inheritParams simu_treatment_block
#' @param prob a numeric vector of minimization treatment assignment probability. 
#' The length is equal to the length of `treatment`. 
#' @param weight_df a numeric vector of weight to sum imbalance value of different stratum. 
#' The length is equal to the number of columns in `df`.
#' @param imbalance_fun a function to calculate imbalance. e.g. `imbalance_fun_range`.
#' 
#' @examples 
#' library(magrittr)
#' simu_stratum(n = 10, 
#'              stratum = define_stratum(sex = c("M" = 1, "F" = 1))) %>% 
#' simu_treatment_minimization(treatment = c("case", "control"), 
#'                             ratio = c(3,2))
#' 
#' @export
simu_treatment_minimization <- function(df, 
                                        treatment, 
                                        ratio, 
                                        prob = define_prob(0.9, treatment),
                                        weight_df = rep(1, ncol(df)), 
                                        imbalance_fun = imbalance_fun_range){
  
  df0 <- lapply(df, as.numeric)
  
  treatment <- factor(treatment, levels = unique(treatment))
  n_trt <- length(treatment)
  
  # Create a list of table object
  tbl <- list()
  for(j in 1:ncol(df)){
    # Create a table object with all 0
    tbl[[j]] <- table(treatment, df[1:n_trt, ][[j]])
    tbl[[j]][tbl[[j]] != 0] <- 0  
  }
  
  # Assign treatment group
  trt <- c()
  for(i in 1:nrow(df)){
    
    
    # Calculate imbalance 
    imb <- list()
    for(j in 1:ncol(df)){
      imb[[j]] <- imbalance_x(tbl[[j]], df0[[j]][i], ratio, imbalance_fun)
    }
    
    imb <- do.call(rbind, imb)
    imb <- apply(imb, 2, function(x) sum(x * weight_df))  
    
    # Calculate probability to treatment
    p <- prob[rank(imb, ties.method = "min")] 
    p <- p / sum(p)
    
    # Assign treatment group
    trt[i] <- sample(treatment, size = 1, prob = p)
    
    # Update table object
    i_index <- as.numeric(trt[i])
    for(j in 1:ncol(df)){
      tbl[[j]][i_index, df0[[j]][i]] <- tbl[[j]][i_index, df0[[j]][i]] + 1
    }
  }
  
  df$treatment <- factor(trt, labels = levels(treatment))
  
  df
}