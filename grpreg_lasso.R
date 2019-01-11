# source('/home/koukouli/Desktop/Dissertation/variables_selection_varying_coef_models/Simulation Study/simulation_study.R')

library(grpreg); library(splines)

# simulate data
n = 100; n_i = 30; p = 23
missing  = F; uneven = F; err = T

Tgrid <- seq(from = 0, to = 30, length = 101)

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
par(mfrow = c(1,1))


## generate data from 100 subjects, with 30 measurements for each subject and 10 covariates 
## of which only 3 are associated with the response 
data <- simul_dat(n, n_i, p,  missing  = F, uneven = F, err = T)


y <- data$responses
X <-data$design_matrix
time <-data$t
  
## plot data
plot(time, y)                              

## construct the Basis approximation matrix
modsplines <- c()
int <- FALSE
df = 10 
for (i in 1:ncol(X)){
  modsplines <- cbind(modsplines, bs(X[,i], df =  df ,intercept = int, degree = 3))
  #print(i)
}


## groups specification - for each predictor define one group
l <- dim(modsplines)[2]
group <- rep(0:p, each = df)

## 5-fold cross validation to choose the optimal penalty parameter
lambda <- cv.grpreg(modsplines, y, group = group, nfolds = 5)
lambda_opt <- lambda$lambda.min

fit <- grpreg(modsplines,y, group, penalty = "grLasso", family = "gaussian", 
              lambda = lambda_opt)

b <- as.matrix(fit$beta)
rownames(b) <- as.character(c(0,group))
print(b)

Beta <- matrix(0, nrow = p+1, ncol = 101)
for (i in 1:(p+1)){
  Beta[i,] <- bs(Tgrid, df = df, intercept = int, degree = 2) %*%
    b[(i*df-(df-1)):(i*df)]
}
rownames(Beta) <- 0:23
colnames(Beta) <- Tgrid


