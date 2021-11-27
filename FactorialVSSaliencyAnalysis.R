###############################################################
################ Pauline Corblet - Sciences Po ################
###############################################################
#
###############################################################
# EXPLORING THE LINKS BETWEEN SALIENCY AND FACTORIAL ANALYSIS #
###############################################################

require('Matrix')
require('gurobi')
require('geigen')
library(FactoMineR)

path_dta = "C:/Users/pauli/Documents/01Work/00PhDProjects/03CEREQData/CEREQ/dta"
path_output = "C:/Users/pauli/Documents/01Work/00PhDProjects/07FactorialAnalysis/output/"
path_notes = "C:/Users/pauli/Documents/01Work/00PhDProjects/07FactorialAnalysis/notes/"

######################## EXAMPLE DATA #########################

#Td = matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,1,1), nrow = 10)
#Tc = matrix(c(3,2,1,0,1,3), nrow = 3)

#X = matrix(c(1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1), nrow = 3)
#Y = matrix(c(1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,1), nrow = 2)

nlow = 10
nmid = 3
nhigh = 2

X = t(cbind(matrix(rep(c(1,0,0), nlow), nrow = 3), matrix(rep(c(0,1,0), nmid), nrow = 3), matrix(rep(c(0,0,1), nhigh), nrow = 3)))

mlow = 5
mmid = 9
mhigh = 1

Y = t(cbind(matrix(rep(c(1,0,0), mlow), nrow = 3), matrix(rep(c(0,1,0), mmid), nrow = 3), matrix(rep(c(0,0,1), mhigh), nrow = 3)))

###################### SET UP MATCHING ########################

Xvals = X
Yvals = Y

colnames(Xvals) = c("low", "mid", "high")
colnames(Yvals) = c("low", "mid", "high")


# sdX = apply(Xvals,2,sd)
# sdY = apply(Yvals,2,sd)
# mX = apply(Xvals,2,mean)
# mY = apply(Yvals,2,mean)
# Xvals = t( (t(Xvals)-mX) / sdX)
# Yvals = t( (t(Yvals)-mY) / sdY)

nobs = dim(Xvals)[1]

#Aff = matrix(c(1,2,3,2,4,6,3,6,9), nrow = 3)

Aff = matrix(sample(-5:5, 9, replace = TRUE), nrow = 3)

Phi = Xvals %*% Aff %*% t(Yvals)

p = rep(1,nobs)
q = rep(1,nobs)

N = dim(Phi)[1]
M = dim(Phi)[2]

c = c(Phi)

A1 = kronecker(matrix(1,1,M),sparseMatrix(1:N,1:N))
A2 = kronecker(sparseMatrix(1:M,1:M),matrix(1,1,N))
A = rbind2(A1,A2)

d = c(p,q) 

result   = gurobi ( list(A=A,obj=c,modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0) ) 
pi = matrix(result$x,nrow=N)

match = which(pi==1, arr.ind=TRUE)
match = match[order(match[,1]),]

Xm = X
Ym = as.matrix(Y[match[,2],])

pim = diag(nobs)

phim =  (Xm %*% Aff %*% t(Ym))*pim

red_phim = t(Xm)%*%phim%*%Ym

########## DEDUCE CONTIGENCY AND DISJONCTIVE TABLES ###########

Tc = t(X)%*%pi%*%Y 
# same as:
t(Xm)%*%pim%*%Ym

n = dim(X)[1]

pi_n = pi/n
T_tilde = t(X)%*%pi_n%*%Y

Td = cbind(X,pi%*%Y)

# Another way to retrive Aff - Only where TC is not zero
other_Aff = red_phim/Tc
#other_Aff[is.nan(other_Aff)]=0

# Or Tc
other_Tc = red_phim/Aff
#other_Tc[is.nan(other_Tc)]=0

######################## COMPARE METHODS ######################

Z <- Tc/sum(Tc)

c <- colSums(Z)
r <- rowSums(Z)

diagr = diag((r)^(-1))
diagc = diag((c)^(-1))

c2 <- colSums(Z)/sum(Z)
r2 <- rowSums(Z)/sum(Z)

diagr2 = diag((r2)^(-0.5))
diagc2 = diag((c2)^(-0.5))

# Should all be the same
V = diagr%*%as.matrix(Z-r%*%t(c))%*%diagc
V2 =  t(t(Z/r)/c) - 1
V3 = diagr2%*%as.matrix(Z-r2%*%t(c2))%*%diagc2
V4 = t(t(V)*sqrt(c))*sqrt(r)

resSVD = svd.triplet(V,  row.w = r, col.w = c)
resSVD2 = svd.triplet(V2,  row.w = r, col.w = c)
resSVD3 = svd(V3)
resSVD4 = svd(V4)

rownames(Tc) <- c("blue", "green", "red")
colnames(Tc) <- c("square", "circle", "triangle")

resCA <- CA(Tc)

# png(paste0(path_notes,'/ex_T_svd.png'))
# plot.CA(resCA, title = "")
# dev.off()

D.CA = diag(resCA$eig[,1])

V.CA = t(t(resCA$col$coord)*(1/resCA$eig[,1]))
U.CA = t(t(resCA$row$coord)*(1/resCA$eig[,1]))

Aff.svd = svd(Aff, nu = nrow(Aff), nv = ncol(Aff))
D.match = diag(Aff.svd$d)
U.match = Aff.svd$u
V.match = Aff.svd$v

# t(matrix(Xvals[1,]))%*%Aff%*%matrix(Yvals[1,])
# # Should be the same as:
# t(matrix(Xvals[1,]))%*%U%*%D%*%t(V)%*%matrix(Yvals[1,])
# # and  
# Xvals %*% Aff %*% t(Yvals)
# #the same as:
# Xvals %*% U.match%*%D.match%*%t(V.match) %*% t(Yvals)

#source('C:/Users/pauli/Documents/01Work/00PhDProjects/07FactorialAnalysis/affinity.R')

s = 0.2
l = 0.005
Aff_est = affinity_old(Xm, Ym, s, l)

Aff_est_norelax = affinity_old(Xm, Ym, s, 0)

# Matrix rank has to be al least 2 for the rest of the code to work!
rankMatrix(Aff_est)

Aff_est.svd = svd(Aff_est, nu = nrow(Aff_est), nv = ncol(Aff_est))
D_est.match = diag(Aff_est.svd$d)
U_est.match = Aff_est.svd$u
V_est.match = Aff_est.svd$v

eig = (diag(D_est.match)^2)[1:3]

contrib_dim = eig/sum(eig)

coord_Aff_est <- list()
coord_Aff_est$col = t(t(V_est.match)*sqrt(eig))
coord_Aff_est$row = t(t(U_est.match)*sqrt(eig))

# range_x = range(c(resCA$col$coord[,1], resCA$row$coord[,1]))
# range_y = range(c(resCA$col$coord[,2], resCA$row$coord[,2]))

range_x = c(-2,1)
range_y = c(-2,2)

# png(paste0(path_notes,'/ex_Ahat_svd.png'))
# 
# plot(x=coord_Aff_est$col[,1], y=coord_Aff_est$col[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=17, col="red")
# abline(v=0,lty=2)
# abline(h=0,lty=2)
# text(x=coord_Aff_est$col[,1], y=coord_Aff_est$col[,2],labels=c("square", "circle", "triangle"), cex = 0.8, pos=c(3,3,3), col = "red")
# 
# par(new=TRUE)
# plot(x=coord_Aff_est$row[,1], y=coord_Aff_est$row[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=16, col="blue")
# text(x=coord_Aff_est$row[,1], y=coord_Aff_est$row[,2],labels=c("blue", "green", "red"), cex = 0.8, pos=c(1,1,1), col = "blue")
# legend('topleft',legend=c(parse(text=sprintf('sigma == %s',s)),parse(text=sprintf('lambda == %s',l))),bty='n');
# 
# dev.off()

resCA$row$coord[,1]/coord_Aff_est$row[,1]
resCA$row$coord[,2]/coord_Aff_est$row[,2]
resCA$col$coord[,1]/coord_Aff_est$col[,1]
resCA$col$coord[,2]/coord_Aff_est$col[,2]

# Calcul des distances entre les points

match_dist = matrix(NA, 3, 3)
CA_dist = matrix(NA, 3,3)
for (x in 1:3){
  for (y in 1:3){
    match_dist[x,y] = dist(rbind(coord_Aff_est$row[x,],coord_Aff_est$col[y,]))
    CA_dist[x,y] = dist(rbind(resCA$row$coord[x,],resCA$col$coord[y,]))
  }
}

prop_dist = CA_dist/match_dist

#prop_dist_ref = CA_dist/match_dist

dist_origin = cbind(sqrt(rowSums(coord_Aff_est$row^2)), sqrt(rowSums(coord_Aff_est$col^2)))

#dist_origin_ref = cbind(sqrt(rowSums(coord_Aff_est$row^2)), sqrt(rowSums(coord_Aff_est$col^2)))

dist_origin_ref/dist_origin

phim_est =  (Xm %*% Aff_est %*% t(Ym))*pim
red_phim_est = t(Xm)%*%phim_est%*%Ym

red_phim_est/Aff_est

# Compute pi^A
phis = kronecker(t(Ym),t(Xm))  
dX=dim(Xm)[2]
dY=dim(Ym)[2]
n=dim(Xm)[1]

p=rep(1/n,n)
q=rep(1/n,n)
IX=rep(1,n)
tIY=matrix(rep(1,n),nrow=1)
f = p %*% tIY
g = IX %*% t(q)
pihat = diag(n)/n
v=rep(0,n)

Phi = Xm %*% matrix(Aff_est,nrow=dX) %*% t(Ym)
contIpfp = TRUE
iterIpfp = 0
tolIpfp = 1E-6
maxiterIpfp = 1E3
while(contIpfp)
{
iterIpfp = iterIpfp+1
u = s*log(apply(g * exp( ( Phi - IX %*% t(v) ) / s ),1,sum))
vnext = s*log(apply(f * exp( ( Phi - u %*% tIY ) / s ),2,sum))
error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / s ),1,sum)-1))
if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
v=vnext
}

pi_est = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / s )

Tc_est = t(Xm)%*%pi_est%*%Ym

Z_est <- Tc_est/sum(Tc_est)

c2 <- colSums(Z_est)/sum(Z_est)
r2 <- rowSums(Z)/sum(Z)

diagr2 = diag((r2)^(-0.5))
diagc2 = diag((c2)^(-0.5))

V_est = diagr2%*%as.matrix(Z-r2%*%t(c2))%*%diagc2

Phi_norelax = Xm %*% matrix(Aff_est_norelax,nrow=dX) %*% t(Ym)
v=rep(0,n)
contIpfp = TRUE
iterIpfp = 0
tolIpfp = 1E-6
maxiterIpfp = 1E3
while(contIpfp)
{
  iterIpfp = iterIpfp+1
  u = s*log(apply(g * exp( ( Phi - IX %*% t(v) ) / s ),1,sum))
  vnext = s*log(apply(f * exp( ( Phi_norelax - u %*% tIY ) / s ),2,sum))
  error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / s ),1,sum)-1))
  if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp) | iterIpfp == 138) {contIpfp=FALSE}
  v=vnext
}

pi_est_norelax = f * g * exp( ( Phi_norelax - IX %*% t(v) - u %*% tIY ) / s )

Tc_est_norelaw = t(Xm)%*%pi_est_norelax%*%Ym


resCA_est = CA(Tc_est)

total_phim = sum(exp(Phi)/s)

Tc_est2 = Tc_est/total_phim
resCA_est2 = CA(Tc_est2)

# phim_est =  (Xm %*% Aff_est %*% t(Ym))*pim
phim_est =  (Xm %*% Aff_est %*% t(Ym))
red_phim_est = t(Xm)%*%phim_est%*%Ym

# Matrix rank has to be al least 2 for the rest of the code to work!
rankMatrix(red_phim_est)

red_phim_est.svd = svd(red_phim_est, nu = nrow(red_phim_est), nv = ncol(red_phim_est))
D_phim_est.match = diag(red_phim_est.svd$d)
U_phim_est.match = red_phim_est.svd$u
V_phim_est.match = red_phim_est.svd$v

eig_phim = (diag(D_phim_est.match)^2)

coord_phim_est <- list()
coord_phim_est$col = t(t(V_phim_est.match)*sqrt(eig_phim))
coord_phim_est$row = t(t(U_phim_est.match)*sqrt(eig_phim))

range_x = c(-20,20)
range_y = c(-20,20)

plot(x=coord_phim_est$col[,1], y=coord_phim_est$col[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=17, col="red")
abline(v=0,lty=2)
abline(h=0,lty=2)
text(x=coord_phim_est$col[,1], y=coord_phim_est$col[,2],labels=c("low", "mid", "high"), cex = 0.8, pos=c(3,3,3), col = "red")

par(new=TRUE)
plot(x=coord_phim_est$row[,1], y=coord_phim_est$row[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=17, col="blue")
text(x=coord_phim_est$row[,1], y=coord_phim_est$row[,2],labels=c("low", "mid", "high"), cex = 0.8, pos=c(1,1,1), col = "blue")
legend('topleft',legend=c(parse(text=sprintf('sigma == %s',s)),parse(text=sprintf('lambda == %s',l))),bty='n');




Aff_est_exp = exp(Aff_est)
rankMatrix(Aff_est_exp)

Aff_est_exp.svd = svd(Aff_est_exp, nu = nrow(Aff_est_exp), nv = ncol(Aff_est_exp))
D_est_exp.match = diag(Aff_est_exp.svd$d)
U_est_exp.match = Aff_est_exp.svd$u
V_est_exp.match = Aff_est_exp.svd$v

eig_est_exp = (diag(D_est_exp.match)^2)

coord_est_exp <- list()
coord_est_exp$col = t(t(V_est_exp.match)*sqrt(eig_est_exp))
coord_est_exp$row = t(t(U_est_exp.match)*sqrt(eig_est_exp))


range_x = c(-3,1)
range_y = c(-2,2)

plot(x=coord_est_exp$col[,1], y=coord_est_exp$col[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=17, col="red")
abline(v=0,lty=2)
abline(h=0,lty=2)
text(x=coord_est_exp$col[,1], y=coord_est_exp$col[,2],labels=c("low", "mid", "high"), cex = 0.8, pos=c(2,2,2), col = "red")

par(new=TRUE)
plot(x=coord_est_exp$row[,1], y=coord_est_exp$row[,2], xlim=range_x, ylim=range_y,xlab="Dim 1", ylab="Dim 2", main="",pch=17, col="blue")
text(x=coord_est_exp$row[,1], y=coord_est_exp$row[,2],labels=c("low", "mid", "high"), cex = 0.8, pos=c(4,4,4), col = "blue")
legend('topleft',legend=c(parse(text=sprintf('sigma == %s',s)),parse(text=sprintf('lambda == %s',l))),bty='n');


# Not working because of negative coefficients in Aff_est
resCA_aff = CA(Aff_est)

resCA_affexp = CA(Aff_est_exp)

# Taking only log of Tc yields -inf because of 0s
log_Tc_est = log(Tc_est)
# Completely different
resCA_Tclog = CA(log_Tc_est)

# Studying exp(A)

Z3 = Aff_est_exp/sum(Aff_est_exp)

c3 <- colSums(Z3)/sum(Z3)
r3 <- rowSums(Z3)/sum(Z3)

diagr3 = diag((r2)^(-0.5))
diagc3 = diag((c2)^(-0.5))

# Should all be the same
V5 = diagr3%*%as.matrix(Z3-r3%*%t(c3))%*%diagc3

rankMatrix(V5)
test_link = lm(c(V5) ~ c(Aff_est))

Aff_est_expred = Aff_est_exp[1:2,]

k = 2
Aff_est_red = Aff_est[1:2,]
Xm_red = cbind(Xm[,1], Xm[,2]+k*Xm[,3])

lm(Aff_est[3,] ~ Aff_est[1,] + Aff_est[2,])

res_gsvd = gsvd(Aff_est,Tc_est)

D1 = diag(res_gsvd$alpha)
D2 = diag(res_gsvd$beta)

res_gsvd$U%*%D1%*%res_gsvd$A%*%t(res_gsvd$Q)
res_gsvd$V%*%D2%*%res_gsvd$A%*%t(res_gsvd$Q)
