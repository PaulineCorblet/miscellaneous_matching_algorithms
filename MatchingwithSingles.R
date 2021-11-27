#################################################
###############  Pauline Corblet  ###############
#################################################

################  September 2019  ###############
#################################################

library('Matrix')
library('tictoc')
library('Rsolnp')
library('numDeriv')

################# DATA SET UP ###################

rm(list=ls())

seed = 777
nbX  = 40
nbY  = 45
nb_tot = nbX*nbY+nbX+nbY

n = rep(1,nbX)
m = rep(1,nbY)

set.seed(seed)

# Two characteristics per type, generated with unform distribution and uncorrelated
charX = 2
charY = 2
Xvals = matrix(runif(charX*nbX), nrow = nbX)
Yvals = matrix(runif(charY*nbY), nrow = nbY)

# Generate affinity matrix
A_xy = 10* matrix(runif(charX*charY), nrow = charX, ncol = charY)
A_x0 = c(runif(charX))
A_0y = c(runif(charY))

A = c(A_xy, A_x0, A_0y)

# generate corresponding surplus
Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
Phi_x0 = Xvals %*% A_x0
Phi_0y = Yvals %*% A_0y

Phi= c(Phi_xy, Phi_x0,Phi_0y)

A_n = kronecker(matrix(1,1,nbY),sparseMatrix(1:nbX,1:nbX))
A_m = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,1,nbX))
  
########## FUNCTIONS FOR OPTIMIZATION #########
  
objective_func_reg <- function(pi, params, Phi){
    
    n = params[1:nbX]
    m = params[(nbX+1):(nbX+nbY)]
    
    pi_xy = pi[1:(nbX*nbY)]
    pi_x0 = pi[((nbX*nbY)+1):((nbX*nbY)+nbX)]
    pi_0y = pi[((nbX*nbY)+nbX+1):((nbX*nbY)+nbX+nbY)]
    
    entropy = 2*sum(pi_xy*log(pi_xy))+sum(pi_x0*log(pi_x0))+sum(pi_0y*log(pi_0y))-sum(n*log(n))-sum(m*log(m))
    return(-(t(Phi)%*%pi - entropy)[1])
}

constraint_func <- function(pi, params, Phi){
  
  n = params[1:nbX]
  m = params[(nbX+1):(nbX+nbY)]
  
  return(as.matrix(cbind(rbind(A_n,A_m),diag(nbX+nbY))%*%pi - c(n,m))[,])
}

IPFP <- function(Phi, params, tolIPFP = 1e-8, maxiterIPFP = 1e4){
  
  n = params[1:nbX]
  m = params[(nbX+1):(nbX+nbY)]
  
  Phi_xy = matrix(Phi[1:(nbX*nbY)], nrow = nbX)
  Phi_x0 = Phi[((nbX*nbY)+1):((nbX*nbY)+nbX)]
  Phi_0y = Phi[((nbX*nbY)+nbX+1):((nbX*nbY)+nbX+nbY)]
  
  contIPFP = TRUE
  iterIPFP = 0

  B = rep(1,nbY)
  K_xy = exp(Phi_xy/2)
  K_x0 = exp(Phi_x0)
  K_0y = exp(Phi_0y)
  
  while(contIPFP){
    
    iterIPFP = iterIPFP+1
    #print(paste0("IPFP iteration number ", iterIPFP))
    
    Delta_x = c(K_xy%*%B)^2+4*n*c(K_x0)
    A = (-c(K_xy%*%B)+sqrt(Delta_x))/(2*c(K_x0))
    
    Delta_y = c(A%*%K_xy)^2+4*m*c(K_0y)
    B = (-c(A%*%K_xy)+sqrt(Delta_y))/(2*c(K_0y))
    
    polynom_x = A*c(K_xy%*%B)+A^2*c(K_x0)-n
    polynom_y = B*c(A%*%K_xy)+B^2*c(K_0y)-m
    
    if(max(abs(polynom_x))<tolIPFP & max(abs(polynom_y))<tolIPFP){
      #print("IPFP converged")
      contIPFP = FALSE
    }
    
    if(iterIPFP>maxiterIPFP){
      print("IPFP maximum number of iterations reached")
      contIPFP = FALSE
    }
  }
  
  u = -2*log(A)
  v = -2*log(B)
  
  return(list("u"=u,"v"=v))
  
}

pi_A <- function(A, params, tolIPFP = 1e-8, maxiterIPFP = 1e4){

  A_xy = matrix(A[1:(charX*charY)], nrow = charX)
  A_x0 = A[((charX*charY)+1):((charX*charY)+charX)]
  A_0y = A[((charX*charY)+charX+1):((charX*charY)+charX+charY)]
  
  Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
  Phi_x0 = Xvals %*% A_x0
  Phi_0y = Yvals %*% A_0y
  
  Phi= c(Phi_xy, Phi_x0,Phi_0y)
  res = IPFP(Phi, params)
  
  pi_xy = exp((Phi_xy-matrix(rep(res$u, nbY),nrow=nbX)-t(matrix(rep(res$v, nbX),nrow=nbY)))/2)
  pi_x0 = exp(Phi_x0-res$u)
  pi_0y = exp(Phi_0y-res$v)
  
  return(c(pi_xy, pi_x0, pi_0y))
  
}

########### GENERATE BENCHMARK MATCHING ############

tic()
pi_hat=pi_A(c(A_xy, A_x0, A_0y),c(n,m))
toc()

pi_hat_xy = pi_hat[1:(nbX*nbY)]
pi_hat_x0 = pi_hat[((nbX*nbY)+1):((nbX*nbY)+nbX)]
pi_hat_0y = pi_hat[((nbX*nbY)+nbX+1):((nbX*nbY)+nbX+nbY)]

constraint_func(pi_hat, c(n,m),Phi)
objective_func_reg(pi_hat, c(n,m),Phi)

############### ESTIMATE AFFINITY MATRIX ############

## 1.By GMM/ Maximum Likelihood

Sigma = kronecker(Yvals,Xvals)

estA_xy = matrix(1, nrow = charX, ncol = charY)
estA_x0 = rep(1, charX)
estA_0y = rep(1, charY)

estA = c(estA_xy,estA_x0,estA_0y)

contGD = TRUE
iterGD = 0
maxiterGD = 1e6
tolGD = 1e-6

gamma = 1e-1

while(contGD){
  
  iterGD = iterGD+1
  print(paste0("Iteration number ", iterGD))

  estPhi_xy =  Xvals %*% estA_xy %*% t(Yvals)
  estPhi_x0 = Xvals %*% estA_x0
  estPhi_0y = Yvals %*% estA_0y
  
  estPhi= c(estPhi_xy, estPhi_x0, estPhi_0y)
  
  tic()
  estres = IPFP(estPhi, c(n,m))
  toc()
  
  estpi_xy = exp((estPhi_xy-matrix(rep(estres$u, nbY),nrow=nbX)-t(matrix(rep(estres$v, nbX),nrow=nbY)))/2)
  estpi_x0 = exp(estPhi_x0-estres$u)
  estpi_0y = exp(estPhi_0y-estres$v)
  
  if(max(abs(constraint_func(c(estpi_xy,estpi_x0,estpi_0y), c(n,m),estPhi)))>1e-4){
    print("Shrodinger equations not satisfied")
    contGD = FALSE
  }

  next_estA_xy = estA_xy - gamma*matrix(c(estpi_xy-pi_hat_xy)%*%Sigma, nrow = charX)
  next_estA_x0 = estA_x0 - gamma*c(estpi_x0-pi_hat_x0)%*%Xvals
  next_estA_0y = estA_0y - gamma*c(estpi_0y-pi_hat_0y)%*%Yvals
  
  next_estA = c(next_estA_xy, next_estA_x0, next_estA_0y)
  
  if(iterGD>maxiterGD){
    print('Max number of iterations reached')
    contGD = FALSE
  }
  
  print(paste0("Max error is ",max(abs((next_estA-estA)/next_estA))))
  
  if (max(abs((next_estA-estA)/next_estA))<tolGD){
    print('Gradient descent converged')
    contGD = FALSE
  }
  
  estA_prev = estA
  estA = next_estA
  
  estA_xy = next_estA_xy
  estA_x0 = c(next_estA_x0)
  estA_0y = c(next_estA_0y)

}

## 1.By OLS as in search

outcome = 2*c(log(pi_hat_xy/(t(matrix(1,nbY,1)%*%pi_hat_x0)*(matrix(1,nbX,1)%*%pi_hat_0y))^0.5))

lm(outcome ~ Sigma - 1)


