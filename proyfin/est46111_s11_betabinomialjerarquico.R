##############################################################
#
# Rat tumor problem
#
#
rm(list=ls())
 
data <- read.table("~/Downloads/est46111_betabinomial_data.txt", header=T)

y <- data$y
n <- data$N
J <- length(y)

##############################################################

log.prior <- function(alpha,beta) {
  {-2.5}*log(alpha + beta)
}

draw.thetas <- function(alpha,beta) {
  return(rbeta(J,alpha+y,beta+n-y))
}

draw.alpha <- function(alpha,beta,theta,prop.sd) {
  alpha.star <- rnorm(1,alpha,prop.sd)
  if (alpha.star<0) { alpha.star <- 0 }
  num <- J*(lgamma(alpha.star+beta) - lgamma(alpha.star)) +
    alpha.star*sum(log(theta)) + log.prior(alpha.star,beta)
  den <- J*(lgamma(alpha+beta)      - lgamma(alpha)) +
    alpha     *sum(log(theta)) + log.prior(alpha,beta)
# print(c(alpha,alpha.star,num,den))
  acc <- ifelse(log(runif(1))<=num - den,1,0)
  alpha.acc <<- alpha.acc + acc
  return(ifelse(acc,alpha.star,alpha))
}

draw.beta <- function(alpha,beta,theta,prop.sd) {
  beta.star <- rnorm(1,beta,prop.sd)
  if (beta.star<0) { beta.star <- 0 }
  num <- J*(lgamma(alpha+beta.star) - lgamma(beta.star)) +
    beta.star*sum(log(1-theta)) + log.prior(alpha,beta.star)
  den <- J*(lgamma(alpha+beta)      - lgamma(beta)) +
    beta     *sum(log(1-theta)) + log.prior(alpha,beta)
# print(c(beta,beta.star,num,den))
  acc <- ifelse(log(runif(1))<=num - den,1,0)
  beta.acc <<- beta.acc + acc

  return(ifelse(acc,beta.star,beta))
}

################################################################  

# B <- 0
# M <- 1000

run.chain <- function(a.start=1,b.start=1,B=0,M=1000) {
  
  MM <- B + M
  
  alpha <- matrix(NA,MM)
  beta <- alpha
  theta <- matrix(NA,nrow=MM,ncol=J)
  
                                        # Metropolis tuning parameters
  alpha.prop.sd <- 0.5
  beta.prop.sd <- 3
  
                                        # Initial values for the chain
  alpha[1] <- a.start
  beta[1] <- b.start
  theta[1,] <- draw.thetas(alpha[1],beta[1]) # or theta[1,] <- (y+.5)/(n+.5)
  
                                        # Monitor acceptance frequency
  alpha.acc <<- 0
  beta.acc <<- 0
  
                                        # MCMC simulation
  for (m in 2:MM) {
    alpha[m] <- draw.alpha(alpha[m-1],beta[m-1],theta[m-1,],alpha.prop.sd)
    beta[m] <- draw.beta(alpha[m],beta[m-1],theta[m-1,],beta.prop.sd)
    theta[m,] <- draw.thetas(alpha[m],beta[m])
  }
  
  good <- (B+1):MM
  
  return(list(alpha=alpha[good],beta=beta[good],theta=theta[good,],
         alpha.rate=alpha.acc/MM,beta.rate=beta.acc/MM))
  
}

test <- run.chain(M=1000)

alpha.mcmc <- test$alpha
beta.mcmc <- test$beta
theta.mcmc <- test$theta
alpha.rate <- test$alpha.rate
beta.rate <- test$beta.rate

par(mfrow=c(2,2))
plot(alpha.mcmc,type="l")
plot(beta.mcmc,type="l")
acf(alpha.mcmc,1000)
acf(beta.mcmc,1000)

print(round(c(alpha.rate,beta.rate),2))

#######################################################################

test <- run.chain(M=10000)

alpha.mcmc <- test$alpha
beta.mcmc <- test$beta
theta.mcmc <- test$theta
alpha.rate <- test$alpha.rate
beta.rate <- test$beta.rate

par(mfrow=c(2,2))
plot(alpha.mcmc,type="l")
plot(beta.mcmc,type="l")
acf(alpha.mcmc,1000)
acf(beta.mcmc,1000)

print(round(c(alpha.rate,beta.rate),2))

#######################################################################

R.hat <- function(phi) {

  M <- dim(phi)[1]
  R <- dim(phi)[2]

  phi.dot <- apply(phi,2,mean)
  phi.dotdot <- mean(phi)

# print(round(c(pd=phi.dot,pdd=phi.dotdot),2))
# scan()

  B <- (M/(R-1))*sum((phi.dot - phi.dotdot)^2)

  s2 <- (sweep(phi,2,phi.dot,"-"))^2

  W <- sum(s2)/(R*(M-1))

  varplus <- (M-1)*W/M + B/M
  varminus <- W

# print(round(c(B=B,W=W,vp=varplus,vm=varminus),2))
# scan()

  return(sqrt(varplus/varminus))
  
}

R <- 3
M <- 10000
B <- 0

alpha <- array(NA,c(M,R))
beta <- array(NA,c(M,R))
theta <- array(NA,c(M,J,R))

for (r in 1:R) {

  alpha.start <- 0
  while(alpha.start<=0) { alpha.start <- rt(1,2) + 2 }
  beta.start <- 0
  while(beta.start<=0) { beta.start <- rt(1,2) + 14 }

  output <- run.chain(alpha.start,beta.start,B=B,M=M)
  alpha[,r] <- output$alpha
  beta[,r] <- output$beta
  theta[,,r] <- output$theta
}

# R.hat(alpha)
#
# R.hat(beta)

R.th <- rep(NA,J)

for (j in 1:J) {
  R.th[j] <- R.hat(theta[,j,])
}

round(c(R.alpha=R.hat(alpha),R.beta=R.hat(beta),R.theta=R.th),4)


#######################################################################

# MC Standard Error...

round(sqrt(var(alpha[,1])/M),2)
# [1] 0.01
round(sqrt(var(alpha[,1])*(2*sum(acf(alpha[,1],M,plot=F)$acf) - 1)/M),2)
# [1] 0

# whoah!

# Note: Long-lag correlations in the ACF function are severely biased.
# you can see this by considering something like this:

bozo <- rnorm(1000)
acf(bozo,1000,plot=F)$acf[1:6]
# [1]  1.00000000 -0.01441198  0.01562296 -0.03066300  0.02334996 -0.01119617
sum(acf(bozo,1000,plot=F)$acf)  # should be 1!
# [1]  0.5

# A potentially better approach is to deliberately use only the
# first k observations and the first k lags, in both the raw
# variance and the 

for (k in 1000*(1:5)) {
  print(round(
     sqrt(var(alpha[(1:k),1])*(2*sum(acf(alpha[,1],k,plot=F)$acf) - 1)/k),2))
}


#######################################################################

# Thinning ...

thin <- .01

subsample <- round(1/thin*(1:(M*thin)))
SM <- length(subsample)

thinned <- alpha[subsample,1]

par(mfrow=c(2,2))

plot(alpha[,1],type="l")
plot(thinned,type="l")
acf(alpha[,1],1000)
acf(thinned,1000)

#######################################################################

# Batching  / batch means...

thin <- .01

subsample <- round(1/thin*(1:(M*thin)))
SM <- length(subsample)

batched <- rep(NA,SM-1)

for (i in 1:(SM-1)) {
  batched[i] <- mean(alpha[(subsample[i]+1):subsample[i+1],1])
}

par(mfrow=c(2,2))

plot(alpha[,1],type="l")
plot(batched,type="l")
acf(alpha[,1],1000)
acf(batched,1000)

#######################################################################






