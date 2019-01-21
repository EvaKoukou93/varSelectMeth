# source('/home/koukouli/Desktop/Dissertation/variables_selection_varying_coef_models/Simulation Study/simulation_study.R')

library(grpreg); library(splines)

# simulate data
n = 100; n_i = 30; p = 23
missing  = F; uneven = F; err = T

ngrid = 101
Tgrid <- seq(from = 0, to = 30, length = ngrid)

'
  These are the coefficient functions according to which the data was created.
  All the rest coefficient functions are zero.
  beta_0 <- 15+20*sin(pi*t/60)
  beta_1 <- 2-3*cos(pi*(t-25)/15)
  beta_2 <- 6-0.2*t
  beta_3 <- -4 +((20-t)^3)/2000
'
beta <- matrix(0, nrow = p+1, ncol = 101)
beta[1,] <- 15+20*sin(pi*Tgrid/60)
beta[2,] <- 2-3*cos(pi*(Tgrid-25)/15)
beta[3,] <- 6-0.2*Tgrid
beta[4,] <- -4 +((20-Tgrid)^3)/2000

par(mfrow = c(2,2))
plot(Tgrid, beta[1,], type = 'l'); plot(Tgrid, beta[2,], type = 'l') 
plot(Tgrid, beta[3,], type = 'l'); plot(Tgrid, beta[4,], type = 'l')

plot(Tgrid, beta[5,], type = 'l'); plot(Tgrid, beta[6,], type = 'l') 
plot(Tgrid, beta[7,], type = 'l'); plot(Tgrid, beta[8,], type = 'l')
par(mfrow = c(1,1))


## generate data from 100 subjects, with 30 measurements for each subject and 10 covariates 
## of which only 3 are associated with the response 
data <- simul_dat(n, n_i, p,  missing  = F, uneven = F, err = T)


y <- data$responses
X <-data$design_matrix
time <-data$t
  
## plot data
plot(time, y)                              
rXsum = apply(X, 1, sum)
cXsum = apply(X, 2, sum)
plot(rXsum)
plot(cXsum)

## construct the Basis approximation matrix
## Create the basis approximation functions matrix
B0 = bs(time, df =  df ,intercept = int, degree = deg)
matplot(time, B0, type='l')
## generate the design matrix
modsplines = B0
for (i in 2:ncol(X)){
  Bi = diag(X[,i])%*%B0
  modsplines = cbind(modsplines, Bi)
}

rmsum = apply(modsplines, 1, sum)
cmsum = apply(modsplines, 2, sum)
cmsum[1:30]
par(mfrow=c(1,2))
plot(rmsum)
plot(cmsum)

par(mfrow=c(1,1))
ind = 1:10
#ind = 11:20
matplot(time, modsplines[,ind], type='l')

## groups specification - for each predictor define one group
l <- dim(modsplines)[2]
group <- rep(0:p, each = df)

## 5-fold cross validation to choose the optimal penalty parameter
Xmod = modsplines[,-1]
Xgroup = group[-1]
lambda <- cv.grpreg(Xmod, y, group = Xgroup, nfolds = 5)
lambda_opt <- lambda$lambda.min

fit <- grpreg(Xmod,y, Xgroup, penalty = "grLasso", family = "gaussian", 
              lambda = lambda_opt)
fit$group
fit$group.multiplier  # missing group 0: need to redefine group from 1

b <- as.matrix(fit$beta)
rownames(b) <- as.character(group)
print(b)

Beta <- matrix(0, nrow = p+1, ncol = ngrid)
BB=bs(Tgrid, df = df, intercept = int, degree = deg)
for (i in 1:(p+1)){
  ind = (i*df-(df-1)):(i*df)
  #print(ind)
  Beta[i,] <- BB %*%b[ind]
}
rownames(Beta) <- 1:(p+1)
colnames(Beta) <- Tgrid

par(mfrow=c(2,3))
for (i in 1:6){
  ymax = max(Beta[i,], beta[i,])
  ymin = min(Beta[i,], beta[i,])
  plot(Tgrid, beta[i,], type='l', ylim=c(1.1*ymin, 1.1*ymax), ylab=paste("beta", i))
  lines(Tgrid, Beta[i,], col=2)
} 

##########################################################################################
## check bias in intercept
ymean = mean(y)
y0 = y - ymean
## 5-fold cross validation to choose the optimal penalty parameter
Xmod = modsplines[,-1]
Xgroup = group[-1]
lambda1 <- cv.grpreg(Xmod, y0, group = Xgroup, nfolds = 5)
par(mfrow=c(1,1))
plot(lambda1)
lambda1_opt <- lambda1$lambda.min

fit1 <- grpreg(Xmod,y0, Xgroup, penalty = "grLasso", family = "gaussian", 
               lambda = lambda1_opt)
b1 = as.matrix(fit1$beta)

Beta1 <- matrix(0, nrow = p+1, ncol = ngrid)
BB=bs(Tgrid, df = df, intercept = int, degree = deg)
for (i in 1:(p+1)){
  ind = (i*df-(df-1)):(i*df)
  #print(ind)
  Beta1[i,] <- BB %*%b1[ind]
}

par(mfrow=c(2,3))
const = rep(0,p+1)
const[1]=ymean
for (i in 1:6){
ymax = max(Beta1[i,]+const[i], beta[i,])
ymin = min(Beta1[i,]+const[i], beta[i,])
plot(Tgrid, beta[i,], type='l', ylim=c(1.1*ymin, 1.1*ymax))
lines(Tgrid, Beta1[i,]+const[i], col=2)
}

par(mfrow=c(1,2))
plot(time, y)
lines(Tgrid, Beta[1,], col=2)
lines(Tgrid, Beta1[1,]+ymean, col=3)

plot(Tgrid, beta[1,], type='l', ylim=range(y))
lines(Tgrid, Beta[1,], col=2)
lines(Tgrid, Beta1[1,]+ymean, col=3)

###############################################################################################
## group weighting
#group <- rep(1:(p+1), each = df)
Xmod = modsplines[,-1]
Xgroup = group[-1]
cxsum = apply(Xmod, 2, sum)
ngr = length(unique(Xgroup))
gweight = rep(0, ngr)
for (i in 1:ngr){
  gind = which(Xgroup==(i-1))
  gweight[i] = sum(abs(cxsum[gind]))
}
gw = sqrt(gweight)

fit2 <- grpreg(Xmod,y, Xgroup, penalty = "grLasso", family = "gaussian", group.multiplier = gw[-1])

lambda2 <- cv.grpreg(Xmod, y, group = Xgroup, nfolds = 5, penalty="grLasso", group.multipler=gw)
plot(lambda2)
lambda_opt2 <- lambda2$lambda.min

fit2 <- grpreg(Xmod,y, Xgroup, penalty = "grLasso", family = "gaussian", group.multiplier = gw[-1],
              lambda = lambda_opt)

b2 <- as.matrix(fit2$beta)
rownames(b2) <- as.character(group)
print(b2)

Beta2 <- matrix(0, nrow = p+1, ncol = ngrid)
BB=bs(Tgrid, df = df, intercept = int, degree = deg)
for (i in 1:(p+1)){
  ind = (i*df-(df-1)):(i*df)
  #print(ind)
  Beta2[i,] <- BB %*%b2[ind]
}

par(mfrow=c(2,3))
for (i in 1:6){
  ymax = max(Beta2[i,], beta[i,])
  ymin = min(Beta2[i,], beta[i,])
  plot(Tgrid, beta[i,], type='l', ylim=c(1.1*ymin, 1.1*ymax), ylab=paste("beta", i))
  lines(Tgrid, Beta2[i,], col=2)
} 

