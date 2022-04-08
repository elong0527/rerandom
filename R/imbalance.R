#' Imbalance function using range
#' 
#' @param x a numeric value 
#' 
#' @examples 
#' imbalance_fun_range(1:10)
#' 
#' @export
imbalance_fun_range <- function(x){
  max(x) - min(x)
}

#' Imbalance function of a variable 
#' 
#' @inheritParams simu_treatment_minimization
#' @param tbl a `table` object. 
#' @param j a integer of new category 
imbalance_x <- function(tbl, 
                        j, 
                        ratio,
                        imbalance_fun){
  
  x0 <- tbl[, j] / ratio
  x1 <- (tbl[, j] + 1) / ratio
  
  res <- c()
  for(i in 1:nrow(tbl)){
    x <- x0
    x[i] <- x1[i]
    res[i] <- imbalance_fun_range(x)
  }
  names(res) <- rownames(tbl)
  
  res
  
}




