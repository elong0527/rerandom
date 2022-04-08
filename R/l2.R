#' L following Andrews & Buchinsky (2000)
#' 
#' @param a one sided alpha
#' 
#' @examples 
#' 
#' # Match results from Andrews &Buchinsky (2000) at page 37. 
#' l2(0.025, tau = 0.95) 
#' (267/10)^2
l2 <- function(a, pdb = 10, tau = 0.99, w = NULL){
  # pdb = 10
  # tau = 0.99
  
  z_tau = qnorm(1 - (1 - tau)/2 )
  
  if(is.null(w)){
    w = a * (1 - a) / (qnorm(1 - a) * dnorm(qnorm(1 - a))) ^2
  }
  
  n <- (10000 * z_tau^2 * w / pdb^2)
  
  as.integer(n) 
}

#' Step 1 in  Andrews &Buchinsky (2000)
#' 
#' @param a one sided alpha
#' 
#' @examples 
#' 
#' # Match results from Andrews &Buchinsky (2000) at page 37. 
#' l2_frac(0.025, tau = 0.95) 
#' (267/10)^2
l2_frac <- function(a, pdb = 10, tau = 0.99, w = NULL){
  # pdb = 10
  # tau = 0.99
  
  a_frac <- MASS:::.rat(a)$rat
  
  z_tau = qnorm(1 - (1 - tau)/2 )
  
  if(is.null(w)){
    w = a * (1 - a) / (qnorm(1 - a) * dnorm(qnorm(1 - a))) ^2
  }
  
  n <- (10000 * z_tau^2 * w / pdb^2 / a_frac[2])
  
  # as.integer(n) * a_frac[2] - 1
  ceiling(n) * a_frac[2] - 1
}

#' Step 2 and 3 to update L in Andrews & Buchinsky (2000)
#' 
#' @param a one sided alpha level
#' @param t_stat statistics approximately following standard normal distribution 
#' @param pdb value of pdb in Andrews & Buchinsky (2000)
#' @param tau value of tau in Andrews & Buchinsky (2000)
#' 
#' @examples 
#' b <- l2_frac(a = 0.01, pdb = 10, tau = 0.99)
#' l2_w(a = 0.01, t_stat = rnorm(b), pdb = 10, tau = 0.99)
#' 
#' @export
l2_w <- function(a, t_stat, pdb = 10, tau = 0.99){
  
  b <- l2_frac(a, pdb = pdb, tau = tau)
  t_stat <- sort(t_stat)
  
  nu <- as.integer((b + 1) * (1 - a)) 
  
  c_a <- (1.5 * qnorm(1 - a/2)^2 * dnorm(qnorm(1 - a))^2 / (2 * qnorm(1 - a/2)^2 + 1))^(1/3)
  m <- ceiling(c_a * b^(2/3))
  
  
  q <- t_stat[nu] 
  g <- 2 * m / b / (t_stat[nu + m] - t_stat[nu - m])
  
  w <- (a * (1 - a)) / q^2 / g^2
  
  l2_frac(a, pdb, tau, w)
}
