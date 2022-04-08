#' Define stratum 
#' 
#' @param ... unique level of each variable 
#' 
#' @examples 
#' define_stratum(site = rep(1,10), sex = c("male" = 1, "female" = 1))
#' 
#' @export
define_stratum <- function(...){
  
  # TODO: input checking 
  # 1. each argument contain unique values 
  res <- list(...)
  res <- lapply(res, function(x){
    
    x_level <- names(x)
    if(is.null(x_level)) x_level <- 1:length(x)
    
    if(any(duplicated(x_level))){
      stop("There is duplicated value: \n", paste(x_level, collapse = ", "))
    }
    
    x <- x / sum(x) # standardize weight
    
    list(level = factor(x_level, levels = unique(x_level)), 
         weight = x)
  })
  
  res
  
}

#' Define treatment group
#' 
#' @param treatment a vector of treatment group and its randomization ratio
#' 
#' @examples 
#' define_treatment(c("case" = 1, "control" = 1))
#' 
#' @export
define_treatment <- function(treatment){
  
  level <- names(treatment)
  list(level = factor(level, levels = unique(level)), 
       weight = treatment / sum(treatment))    
  
}

#' Define probability of randomization 
#' 
#' @param p a value of probability for first treatment
#' @param treatment a vector of treatment group and its randomization ratio
#' 
#' @examples 
#' define_prob(0.9, c("group1", "group2", "group3")) 
#' 
#' @export
define_prob <- function(p, treatment){
  n <- length(treatment) - 1
  q <- (1 - p) / n
  c(p, rep(q, n))
}
