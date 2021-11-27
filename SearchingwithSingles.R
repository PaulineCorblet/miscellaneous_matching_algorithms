#####################################################
#################  Pauline Corblet  #################
#####################################################

##################  September 2019  #################
#####################################################

library('Matrix')
library('tictoc')
library('numDeriv')
library('crayon')

################### DATA SET UP #####################

rm(list=ls())

seed = 777
nbX  = 40
nbY  = 45

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

# Population densities
lf = rep(1/(nbX+nbY), nbX)
lm = rep(1/(nbX+nbY), nbY)

# Useful tools for algebra
A_n = kronecker(matrix(1,1,nbY),sparseMatrix(1:nbX,1:nbX))
A_m = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,1,nbX))


######### FUNCTION FOR FINDING STEADY STATE #########

meeting_rate <- function(U,V){
  return(U*V)
}

steady_state <- function(Phi_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  
  lf = params[1:nbX]
  lm = params[(nbX+1):(nbX+nbY)]
  
  cont = TRUE
  iter = 0
  
  U = nbX/(nbX+nbY)
  V = nbY/(nbX+nbY)
  H = 0
  
  # Densities for singles and couples
  u = rep(U/nbX, nbX)
  v = rep(V/nbY, nbY)
  h = matrix(0, nrow = nbX, ncol = nbY)
  
  # Meeting rate
  lambda = meeting_rate(U,V)/(U*V)
  
  while(cont){
    
    iter = iter+1
    
    # Flows of marriages and divorces
    alpha_xy = matrix(1/(1+exp(-c(Phi_xy))), nrow = nbX, ncol = nbY) 
    
    marriages = lambda*matrix(rep(u, nbY), nrow = nbX)*t(matrix(rep(v, nbX), nrow = nbY))*alpha_xy
    divorces = delta*h*(1-alpha_xy)
    
    # Update densities accordingly 
    h = h+marriages-divorces
    u = u-apply(marriages, 1, sum)+apply(divorces, 1, sum)
    v = v-apply(marriages, 2, sum)+apply(divorces, 2, sum)
    
    U = sum(u)
    V = sum(v)
    H = sum(h)
    
    # Meeting rate
    lambda = meeting_rate(U,V)/(U*V)
    
    if(iter > maxiterSS){
      cat(magenta('Stop: max number of iterations reached\n'))
      cont = FALSE
    }
    
    max_diff = max(abs(marriages-divorces))
    
    if(max_diff<tolSS){
      cat(green('Stop: steady state reached\n'))
      cont = FALSE
    }
  }
  
  cat(blue('Max difference between marriages and divorces is: ',max_diff,'\n'))
  return(list("u"=u,"v"=v,"h"=h))
  
}

steady_state_couples <- function(A_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
  return(steady_state(Phi_xy, params, delta, tolSS, maxiterSS)$h)
}

steady_state_single_women <- function(A_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
  return(steady_state(Phi_xy, params, delta, tolSS, maxiterSS)$u)
}

steady_state_single_men <- function(A_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
  return(steady_state(Phi_xy, params, delta, tolSS, maxiterSS)$v)
}

steady_state_singles <- function(A_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  Phi_xy =  Xvals %*% A_xy %*% t(Yvals)
  return(c(steady_state(Phi_xy, params, delta, tolSS, maxiterSS)$u,steady_state(Phi_xy, params, delta, tolSS, maxiterSS)$v))
}

next_period_flows <- function(h, u, v, Phi_xy, params, delta, tolSS = 1e-9, maxiterSS = 1e6){
  
  lf = params[1:nbX]
  lm = params[(nbX+1):(nbX+nbY)]
  
  cont = TRUE
  iter = 0
  
  U = sum(u)
  V = sum(v)
  H = sum(h)
  
  # Meeting rate
  lambda = meeting_rate(U,V)/(U*V)
    
  # Flows of marriages and divorces
  alpha_xy = matrix(1/(1+exp(-c(Phi_xy))), nrow = nbX, ncol = nbY) 
    
  marriages = lambda*matrix(rep(u, nbY), nrow = nbX)*t(matrix(rep(v, nbX), nrow = nbY))*alpha_xy
  divorces = delta*h*(1-alpha_xy)
    
  return(list('divorces' = divorces, 'marriages' = marriages))
  
}

################## MODEL SIMULATION #################

mu = (nbX*nbY)/(nbX+nbY)^2
delta = 1

matching = steady_state(Phi_xy, c(lf,lm), delta)

h_hat = matching$h
u_hat = matching$u
v_hat = matching$v

hist(u_hat)
hist(v_hat)
hist(h_hat)

library(plotly)
plot_ly(x=Xvals[,1], y = Xvals[,2], z = u_hat)


test_u <- function(vec){
  x=vec[1]
  y=vec[2]
  return(x+y+x^2+y^2+2*x*y)
}

z <- sapply(Xvals, FUN = test_u)

############### ESTIMATE AFFINITY MATRIX ############

## 1.a By OLS over one period

outcome = c(log(h_hat/(t(matrix(1,nbY,1)%*%u_hat)*(matrix(1,nbX,1)%*%v_hat))))
regressors = kronecker(Yvals,Xvals)

res1 = lm(outcome ~ regressors - 1)

max_diff_OLS = max(abs(c(A_xy)-res$coefficients))
cat(blue('Max difference between OLS coefficients and true values is ', max_diff_OLS, '\n'))

## 1.b By OLS over two period

next_flows = next_period_flows(h_hat, u_hat, v_hat, Phi_xy,c(lf,lm), delta)

divorces = next_flows$divorces
marriages = next_flows$marriages

MR = c(marriages/(t(matrix(1,nbY,1)%*%u_hat)*(matrix(1,nbX,1)%*%v_hat)))
DR = c(divorces/h_hat)

res2_rates = lm(rep(1, nbX*nbY) ~ MR + DR - 1 )

outcome = c(log(h_hat/(t(matrix(1,nbY,1)%*%u_hat)*(matrix(1,nbX,1)%*%v_hat))))
regressors = kronecker(Yvals,Xvals)

cons = log(res2_rates$coefficients[2]/res2_rates$coefficients[1])
res2 = lm((outcome - cons) ~ regressors - 1)

max_diff_OLS = max(abs(c(A_xy)-res2$coefficients))
cat(blue('Max difference between OLS coefficients and true values is ', max_diff_OLS, '\n'))

## 2. By maximum likelihood

Sigma = kronecker(Yvals,Xvals)

estA_xy = matrix(1, nrow = charX, ncol = charY)

contGD = TRUE
iterGD = 0
maxiterGD = 1e6
tolGD = 1e-4

gamma = 10

while(contGD){
  
  iterGD = iterGD+1
  print(paste0("Iteration number ", iterGD))
  
  estPhi_xy =  Xvals %*% estA_xy %*% t(Yvals)
  estmatching = steady_state(estPhi_xy, c(lf,lm), delta)
  
  estpi_xy = c(estmatching$h, estmatching$u, estmatching$v)
  
  jac_singles = jacobian(func=steady_state_singles, x=estA_xy, method = 'simple', params= c(lf,lm), delta = delta)
  jac_u = jac_singles[1:nbX,]
  jac_v = jac_singles[(nbX+1):(nbX+nbY),]
  
  grad_ll = (c(h_hat)%*%Sigma 
    + c(t(jac_u/estmatching$u)%*%apply(h_hat, 1, sum))
    + c(t(jac_u/estmatching$u)%*%u_hat)
    + c(t(jac_v/estmatching$v)%*%apply(h_hat, 2, sum))
    + c(t(jac_v/estmatching$v)%*%v_hat)
  )
  
  next_estA_xy = estA_xy + gamma*matrix(grad_ll, nrow = charX)
  
  if(iterGD>maxiterGD){
    print('Max number of iterations reached')
    contGD = FALSE
  }
  
  print(paste0("Max error is ",max(abs((next_estA_xy-estA_xy)/next_estA_xy))))
  
  cat(cyan('Affinity matrix is: \n'))
  cat(cyan(next_estA_xy, '\n'))
  
  if (max(abs((next_estA_xy-estA_xy)/next_estA_xy))<tolGD){
    print('Gradient descent converged')
    contGD = FALSE
  }
  
  estA_xy = next_estA_xy
  
}

## 1.By GMM/ Maximum Likelihood for matching model

Sigma = kronecker(Yvals,Xvals)

estA_xy = matrix(1, nrow = charX, ncol = charY)
estA_x0 = rep(1, charX)
estA_0y = rep(1, charY)

estA = c(estA_xy,estA_x0,estA_0y)

contGD = TRUE
iterGD = 0
maxiterGD = 1e6
tolGD = 1e-4

gamma = 1e-1

while(contGD){
  
  iterGD = iterGD+1
  print(paste0("Iteration number ", iterGD))
  
  estPhi_xy =  Xvals %*% estA_xy %*% t(Yvals)
  estPhi_x0 = Xvals %*% estA_x0
  estPhi_0y = Yvals %*% estA_0y
  
  estPhi= c(estPhi_xy, estPhi_x0, estPhi_0y)
  
  tic()
  estres = IPFP(estPhi, c(n,m), maxiterIPFP = 1e6)
  toc()
  
  estpi_xy = exp((estPhi_xy-matrix(rep(estres$u, nbY),nrow=nbX)-t(matrix(rep(estres$v, nbX),nrow=nbY)))/2)
  estpi_x0 = exp(estPhi_x0-estres$u)
  estpi_0y = exp(estPhi_0y-estres$v)
  
  if(max(abs(constraint_func(c(estpi_xy,estpi_x0,estpi_0y), c(n,m),estPhi)))>1e-4){
    print("Shrodinger equations not satisfied")
    contGD = FALSE
  }
  
  next_estA_xy = estA_xy - gamma*matrix(c(estpi_xy-c(h_hat))%*%Sigma, nrow = charX)
  next_estA_x0 = estA_x0 - gamma*c(estpi_x0-c(u_hat))%*%Xvals
  next_estA_0y = estA_0y - gamma*c(estpi_0y-c(v_hat))%*%Yvals
  
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
  
  cat(cyan('Affinity matrix is: \n'))
  cat(cyan(next_estA_xy, '\n'))
  
  estA_prev = estA
  estA = next_estA
  
  estA_xy = next_estA_xy
  estA_x0 = c(next_estA_x0)
  estA_0y = c(next_estA_0y)
  
}


