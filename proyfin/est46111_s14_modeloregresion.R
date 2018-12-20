rm(list=ls())

set.seed(1)

if(!require('LearnBayes')){install.packages("LearnBayes")}
library("LearnBayes")

data(birdextinct)
attach(birdextinct)

summary(birdextinct)
dim(birdextinct)

plot(birdextinct[,c("time","nesting","size","status")])

# Figura uno
logtime <- log(time)
plot(nesting,logtime)

# Filtro de mediciones
out <- (logtime > 3)
text(nesting[out], logtime[out], label=species[out], pos = 2)	

# Figura dos
plot(jitter(size),logtime,xaxp=c(0,1,1))

# Figura tres
plot(jitter(status),logtime,xaxp=c(0,1,1))

# Aprendizaje de referencia
# logtime = b_0 = b_1 nesting + b_2 size + b_3 status
fit <- lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
summary(fit)

# ----------------------
# Aprendizaje bayesiano
# ----------------------

#
# Prior difusa
# \pi^{o}(beta,lambda) 

# Simulamos datos de la distribucion posterior
# 5 mil datos
beta.sample <- blinreg(fit$y,fit$x,5000)

par(mfrow=c(2,2))
hist(beta.sample$beta[,2],main="NESTING",
 xlab=expression(beta[1]))
hist(beta.sample$beta[,3],main="SIZE",
 xlab=expression(beta[2]))
hist(beta.sample$beta[,4],main="STATUS",
 xlab=expression(beta[3]))
hist((beta.sample$sigma)^(-1),main="PRECISION",
 xlab=expression(lambda))

# Resumen de la posterior
apply(beta.sample$beta,2,quantile,c(.05,.5,.95))

quantile((beta.sample$sigma)^(-1),c(.05,.5,.95))


# Estimacion de los tiempos de extincion
cov1 <- c(1,4,0,0)
cov2 <- c(1,4,1,0)
cov3 <- c(1,4,0,1)
cov4 <- c(1,4,1,1)
X1 <- rbind(cov1,cov2,cov3,cov4)
mean.draws <- blinregexpected(X1,beta.sample)

c.labels <- c("A","B","C","D")
par(mfrow=c(2,2))
j <- 1
for (j in 1:4)
   hist(mean.draws[,j],
      main <- paste("Covariates",c.labels[j]),xlab="log TIME")
 
######## Predicting extinction times

 cov1 <- c(1,4,0,0)
 cov2 <- c(1,4,1,0)
 cov3 <- c(1,4,0,1)
 cov4 <- c(1,4,1,1)
 X1 <- rbind(cov1,cov2,cov3,cov4)
 pred.draws <- blinregpred(X1,beta.sample)

 c.labels <- c("A","B","C","D")
 par(mfrow=c(2,2))
 for (j in 1:4)
   hist(pred.draws[,j],
      main=paste("Covariates",c.labels[j]),xlab="log TIME")

# Validacion del modelo via prediccion
pred.draws <- blinregpred(fit$x,beta.sample)
pred.sum <- apply(pred.draws,2,quantile,c(.05,.95))
par(mfrow=c(1,1))
 ind <- 1:length(logtime)
 matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,
 xlab="INDEX",ylab="log TIME")
points(ind,logtime,pch=19)
out=(logtime>pred.sum[2,])
text(ind[out], logtime[out], label=species[out], pos = 4)

# Validacion del modelo via probabilidades
prob.out <- bayesresiduals(fit,beta.sample,2)
 par(mfrow=c(1,1))
 plot(nesting,prob.out)
 out = (prob.out > 0.35)
 text(nesting[out], prob.out[out], label=species[out], pos = 4)	

