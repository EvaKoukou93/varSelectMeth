simul_dat <- function(n = 1L, n_i = 1L, p = 3L, missing = TRUE, prob = 0.5, uneven = TRUE, err = TRUE){
  
  require(mgcv); require(Matrix); require(dplyr)
  
  # n        : stands for the total number of subjects in the study
  # n_i      : stands for the total number of measurements per subject
  # p        : stands for the number of predictors in the study
  # missing  : this argument stands for longitudinal data generation with missing data
  # prob     : stands for the probability of each skipped slot for a subject
  
  # time-points generation according to the number of measuremetns 
  # specified by the user:
  
  
  # create time points for each subject
  time_points <- data.frame(time = rep(1:n_i, n),
                            subj = rep(1:n, each = n_i))
  if (missing == TRUE){
    # add missing values: each scheduled time had 60% probability of being skipped
    for (i in 1:n){
      skip <- rbinom(n_i, 1, prob)
      time_points$time[(i*n_i-(n_i-1)):(i*n_i)] <- ifelse(skip == 0, NA,
                                                          time_points$time[(i*n_i-(n_i-1)):(i*n_i)])
    }
    # remove NA's
    time_points <- time_points[which(complete.cases(time_points) == TRUE),]
  }
  if (uneven == TRUE){
    time_points$time <- time_points$time + runif(length(time_points$time),
                                                 min = -0.5, max = 0.5)
  }
  t <- time_points$time
  
  # generate coefficient functions
  b <- matrix(ncol = p+1, nrow = length(t))
  rownames(b) <- t 
  
  b[,1] <- 15+20*sin(pi*t/60)
  b[,2] <- 2-3*cos(pi*(t-25)/15)
  b[,3] <- 6-0.2*t
  b[,4] <- -4 +((20-t)^3)/2000
  if (p > 3L){
    for (i in 5:(p+1)){b[,i] = 0}  
  }
  b <- t(b)
  
  # generation of covariate processes
  X <- list()
  lt <- length(t)
  X[[1]] <- rep(1, lt)
  X[[2]] <- runif(lt, min = t/10, max = 2+t/10)
  X[[3]] <- rnorm(lt, mean = 0, sd = sqrt((1+X[[2]])/(2+X[[2]])))
  X[[4]] <- rbinom(lt, size = 1, prob = 0.6)
  
  t_u <- unique(t)
  t_l <- length(t_u)
  if (p >= 3L){
    for (i in 4:(p+1)){
      # sigma_1 = 1
      sigma_1 <- 1
      cov_1 <- matrix(nrow = t_l, ncol = t_l)
      for (j in 1:t_l){
        for (k in 1:t_l){
          cov_1[j,k] <- (sigma_1^2)*exp(-abs(t_u[j]-t_u[k]))
        }
      }
      X[[i]] <- rmvn(n, mu = numeric(t_l), V = cov_1)
      colnames(X[[i]]) <- t_u
      rownames(X[[i]]) <- 1:n
      newX <- c()
      
      X[[i]] <- as.data.frame(X[[i]])
      names(X[[i]]) <- t_u
      X[[i]]$id <- 1:n
      long_X <- reshape(data = X[[i]], direction = 'long', varying = list(1:t_l), v.names = 'x',
                        timevar = 'time', times = t_u)
      names(long_X)[1] <- 'subj'
      newX <- merge(time_points, long_X, by = c('subj', 'time'))
      newX <- newX[order(newX$subj),]
      X[[i]] <- newX$x
    }  
  }
  
  X_all <- X[[1]]
  for (i in 2:(p+1)){
    X_all <- cbind(X_all,X[[i]])
  }
  
  # responses generation
  # error process generation
  Z <- rmvn(1, mu = numeric(t_l), V = cov_1)
  # sigma_2 = 2
  u <- rmvn(1, mu = numeric(t_l), V = diag(x=2^2,
                                           ncol = t_l, nrow = t_l))
  if (err == T){error <- Z+u}else{error <- 0*Z}
  responses <- diag(X_all %*% b)+error
  
  # return dataset
  dfr <- list(id = time_points$subj, responses = responses, 
              design_matrix = X_all, t = time_points$time)
  return(dfr)
}
