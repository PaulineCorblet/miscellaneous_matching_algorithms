
library('numDeriv')

################ DATA SET UP ##############

rm(list=ls())

seed = 777
nbX  = 2
nbY  = 2

n = rep(1,nbX)
m = rep(1,nbY)

set.seed(seed)

A_n = kronecker(matrix(1,1,nbY),diag(nbX))
A_m = kronecker(diag(nbY),matrix(1,1,nbX))

alpha_xy = runif(nbX*nbY)
gamma_xy = -runif(nbX*nbY)

alpha_x0 = rep(0, nbX)
gamma_0y = rep(0, nbY)

############# FUNCTIONS SET UP #############

# constraint_n <- function(pi, alpha, w, n){
#   return(c(cbind(A_n, diag(nbX),matrix(0,nbX,nbY))%*%pi-n))
# }
# 
# constraint_m <- function(pi, gamma, w, m){
#   return(c(cbind(A_m, matrix(0,nbY,nbX),diag(nbY))%*%pi-m))
# }

constraint_nm <- function(pi, utility, w, nm){
  return(c(cbind(rbind(A_n,A_m),diag(nbX+nbY))%*%pi - nm))
}


g <- function(x){x}
h <- function(x){x}

returns_firms <- function(pi, utility, w, nm){
  
  pi_xy = pi[1:(nbX*nbY)]
  pi_x0 = pi[(nbX*nbY+1):(nbX*nbY+nbX)]
  pi_0y = pi[(nbX*nbY+nbX+1):(nbX*nbY+nbX+nbY)]
  
  obj = unlist(lapply(utility-w,g)) 
  return(-(c(pi_xy,pi_x0)%*%obj)[1])
}

returns_workers <- function(pi, utility, w, nm){
  
  pi_xy = pi[1:(nbX*nbY)]
  pi_x0 = pi[(nbX*nbY+1):(nbX*nbY+nbX)]
  pi_0y = pi[(nbX*nbY+nbX+1):(nbX*nbY+nbX+nbY)]
  
  obj = unlist(lapply(utility+w,h)) 
  return(-(c(pi_xy,pi_0y)%*%obj)[1])
}

returns_firms_gradW <- function(pi, utility, w, nm){
  grad(returns_firms, w, pi = pi, utility = utility, nm=nm)
}

returns_workers_gradW <- function(pi, utility, w, nm){
  grad(returns_workers, w, pi = pi, utility = utility, nm=nm)
}


######## SOLVING FOR OPTIMAL WAGE ########

w0_xy = rep(0, nbX*nbY)
w0_x0 = rep(0, nbX)
w0_0y = rep(0, nbY)

res_P = solnp(pars = x0,
              fun = returns_firms,
              eqfun = constraint_nm,
              eqB = rep(0.0, nbX+nbY),
              LB = rep(0.0, nbX*nbY+nbX+nbY),
              utility = c(alpha_xy, alpha_x0),
              w =  c(w0_xy,w0_x0),
              nm = c(n,m),
              control = control)

pi_P = res_P$pars
G_P = -res_P$values[length(res_P$values)]

res_D = solnp(pars = x0,
              fun = returns_workers,
              eqfun = constraint_nm,
              eqB = rep(0.0, nbX+nbY),
              LB = rep(0.0, nbX*nbY+nbX+nbY),
              utility = c(gamma_xy, gamma_0y),
              w =  c(w0_xy,w0_0y),
              nm = c(n,m),
              control = control)

pi_D = res_D$pars
G_D = -res_D$values[length(res_D$values)]

w_prev = c(w0_xy,w0_x0,w0_0y)
pi_P_prev = pi_P

w = w_prev + (pi_P - pi_D) + (pi_P_prev-pi_P)
