# simulate data
set.seed(1)
size = 3
n <- rep(1,nrow(frax))#sample(1:50, size, replace = TRUE)
a <- 1
b <- 15
p = rbeta(size, a, b)
y = frax$fraudRisk
data <- data.frame(y = y, n = n)


# log prior evaluation function
log.prior <- function(alpha, beta) {
  -2.5 * log(alpha + beta)
}
# sample proposal value
getstar = function(current, lambda) {
  return(current * exp(lambda * (runif(1) -
                                   0.5)))
}
# conditional sample of p
draw.p <- function(alpha, beta) {
  return(rbeta(size, alpha + y, beta +
                 n - y))
}

draw.alpha <- function(alpha, beta, p, lambda) {
  alpha.star <- getstar(alpha, lambda)
  num <- size * (lgamma(alpha.star + beta) -
                   lgamma(alpha.star)) + alpha.star *
    sum(log(p))
  num <- num + log.prior(alpha.star, beta)
  num <- num + log(alpha.star)
  den <- size * (lgamma(alpha + beta) -
                   lgamma(alpha)) + alpha * sum(log(p))
  den <- den + log.prior(alpha, beta)
  den <- den + log(alpha)
  acc <- ifelse((log(runif(1)) <= num -
                   den) && (alpha.star > 0), 1, 0)
  return(c(acc, ifelse(acc, alpha.star,
                       alpha)))
}



draw.beta <- function(alpha, beta, p, lambda) {
  beta.star <- getstar(beta, lambda)
  num <- size * (lgamma(alpha + beta.star) -
                   lgamma(beta.star)) + beta.star *
    sum(log(1 - p))
  num <- num + log.prior(alpha, beta.star)
  num <- num + log(beta.star)
  den <- size * (lgamma(alpha + beta) -
                   lgamma(beta)) + beta * sum(log(1 -
                                                    p))
  den <- den + log.prior(alpha, beta)
  den <- den + log(beta)
  acc <- ifelse((log(runif(1)) <= num -
                   den) && (beta.star > 0), 1, 0)
  return(c(acc, ifelse(acc, beta.star,
                       beta)))
}

Nitr <- 200
a.draw <- b.draw <- matrix(NA, Nitr, 2)
ps <- matrix(NA, nrow = Nitr, ncol = size)
# Metropolis tuning parameters
lambda.alpha <- 0.7
lambda.beta <- 0.7
# Initial values for the chain
a.draw[1, 2] <- 1
b.draw[1, 2] <- 1
ps[1, ] <- draw.p(a.draw[1, 2], b.draw[1,
                                       2])
# MCMC simulation
for (m in 2:Nitr) {
  a.draw[m, ] <- draw.alpha(a.draw[m -
                                     1, 2], b.draw[m - 1, 2], ps[m - 1,
                                                                 ], lambda.alpha)
  b.draw[m, ] <- draw.beta(a.draw[m, 2],
                           b.draw[m - 1, 2], ps[m - 1, ], lambda.beta)
  ps[m, ] <- draw.p(a.draw[m, 2], b.draw[m,
                                         2])
}
# thinning
draws <- seq(10, Nitr, by = 5)
a.draw <- a.draw[draws, ]
b.draw <- b.draw[draws, ]
ps <- ps[draws, ]


par(mfcol = c(2, 3))
plot(a.draw[, 2], type = "l", main = "Trace plot of alpha")
abline(h = a, col = "red")
plot(b.draw[, 2], type = "l", main = "Trace plot of beta")
abline(h = b, col = "red")
plot(density(a.draw[, 2]), main = "Posterior distribution of alpha")
abline(v = a, col = "red")
plot(density(b.draw[, 2]), main = "Posterior distribution of beta")
abline(v = b, col = "red")
plot(p, y/n, xlim = c(0, 0.4), ylim = c(0,
                                        0.4), xlab = "Truth", ylab = "MLE", main = "Prevalence: MLE")
abline(c(0, 1))
plot(p, apply(ps, 2, median), xlim = c(0,
                                       0.4), ylim = c(0, 0.4), xlab = "Truth",
     ylab = "Posterior median", main = "Prevalence: Bayes")
abline(c(0, 1))
print(round(c(alpha.rate = mean(a.draw[,
                                       1]), beta.rate = mean(b.draw[, 1])),
            2))




###################################################################################


log_fc_alpha = function(theta, alpha, beta) {
  if (alpha<0) return(-Inf)
  n = length(theta)
  (alpha-1)*sum(log(theta))-n*lbeta(alpha,beta)-5/2*(alpha+beta)
}
log_fc_beta = function(theta, alpha, beta) {
  if (beta<0) return(-Inf)
  n = length(theta)
  (beta-1)*sum(log(1-theta))-n*lbeta(alpha,beta)-5/2*(alpha+beta)
}



mcmc = function(n_sims, dat, inits, tune) {
  n_groups = nrow(dat)
  alpha = inits$alpha
  beta = inits$beta
  # Recording structure
  theta_keep = matrix(NA, nrow=n_sims, ncol=n_groups)
  alpha_keep = rep(alpha, n_sims)
  beta_keep = rep(beta , n_sims)
  for (i in 1:n_sims) {
    # Sample thetas
    theta = with(dat, rbeta(length(y), alpha+y, beta+n-y))
    # Sample alpha
    alpha_prop = rnorm(1, alpha, tune$alpha)
    logr = log_fc_alpha(theta, alpha_prop, beta)-log_fc_alpha(theta, alpha, beta)
    alpha = ifelse(log(runif(1))<logr, alpha_prop, alpha)
    # Sample beta
    beta_prop = rnorm(1, beta, tune$beta)
    logr = log_fc_beta(theta, alpha, beta_prop)-log_fc_beta(theta, alpha, beta)
    beta = ifelse(log(runif(1))<logr, beta_prop, beta)
    # Record parameter values
    theta_keep[i,] = theta
    alpha_keep[i] = alpha
    beta_keep[ i] = beta
  }
  return(data.frame(iteration=1:n_sims,
                    parameter=rep(c("alpha","beta",paste("theta[",1:n_groups,"]",sep="")),each=n_sims),
                    value=c(alpha_keep,beta_keep,as.numeric(theta_keep))))
}

dat=data.frame(y=y, n=n)
inits = list(alpha=1, beta=1)
# Run the MCMC
r = mcmc(200, dat=dat, inits=inits, tune=list(alpha=1,beta=1))



'model{
  for(t in 1:N){
    y[t] ~ dbern(theta[seg[t]])
  }
  for(j in 1:nCoins){
    theta[j] ~ dbeta(a, b)
  }
  a <- mu * kappa
  b <- (1 - mu) * kappa
  # A_mu = 2, B_mu = 2
  mu ~ dbeta(.1, .1)
  kappa ~ dgamma(1, 0.1)
}
'

#y<-frax$fraudRisk
m<-3
N<-100000
a<-matrix(0L, N, m)
b<-matrix(0L, N, m)
th<-matrix(0L, N, m)
wy<-matrix(0L, N, m)

for(j in 1:m){
  ky<-frax[frax$Segmento==j,8]
  unos<-sum(ky)
  dif<-length(ky)-unos  
  #a<-rep(0,N)
  #b<-rep(0,N)
  #th<-rep(0,N)
for (t in 1:N){
  k<-rgamma(1,.01,.01)
  mu<-rbeta(1,.1,.1)
  a[t,j]<-mu*k
  b[t,j]<-(1-mu)*k
  th[t,j]<-rbeta(1,a[t,j]+unos,b[t,j]+dif)
  wy[t,j]<-rbernoulli(1,th[t,j])
  
}
  
}
























