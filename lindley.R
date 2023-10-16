library(lamW)

mu=.7
q=.8
d<-function(y){
  ((1-mu)^2/(mu*(1-y)^3))*exp((-y*(1-mu))/(mu*(1-y)))
}
integrate(d,0,q)

u=1-(1-(1-mu)*q/(q-1))*exp(
  -(1-mu)*q/(mu*(1-q))
)

theta<-1/mu-1

(1+theta+lambertWm1((1+theta)*(u-1)*exp(-(1+theta))))/
  (1+lambertWm1((1+theta)*(u-1)*exp(-(1+theta))))

rep<- 1/mu
q<-(rep+lambertWm1((rep)*(u-1)*exp(-(rep))))/
  (1+lambertWm1((rep)*(u-1)*exp(-(rep))))
q

#########################################################################

# density function
dUL<-function(y, mu=mu){
  ((1-mu)^2/(mu*(1-y)^3))*exp((-y*(1-mu))/(mu*(1-y)))
}

# quantile function
qUL<-function(u,mu){
  rep<- 1/mu
  q<-(rep+lambertWm1((rep)*(u-1)*exp(-(rep))))/
    (1+lambertWm1((rep)*(u-1)*exp(-(rep))))
}

# inversion method for random generation
rUL<-function(n,mu)
{
  u<- runif(n)
  y<- qUL(u=u,mu=mu)
}

# cumulative distribution function
pUL<-function(q, mu=mu, lower.tail = TRUE, log.p = FALSE)
{
  cdf1<- 1-(1-(1-mu)*q/(q-1))*exp(
    -(1-mu)*q/(mu*(1-q))
  )
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<-1-cdf1
  if(log.p==FALSE) cdf<-cdf else cdf<-log(cdf)
  cdf
}

##########################################################################
#log-likelihood of the unit Lindley

t = sum(x/1-x)

UL<-expression(2*n*log(theta)-n*log(1+theta)-theta*t)

m1UL<-D(UG,"mu")
s1UL<-D(UG,"sigma")
ms2UL<-D(m1UG,"sigma")

UG<-function (mu.link = "logit", sigma.link = "identity")
{
  mstats <- checklink("mu.link", "UG", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UG", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("UG", "Unit-Gamma"),
                 parameters = list(mu = TRUE, sigma = TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 dldm = function(y, mu, sigma) {
                   dldm <- eval(m1UG)
                   dldm
                 },
                 d2ldm2 = function(y,mu, sigma) {
                   dldm <- eval(m1UG)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma) {
                   dldd <- eval(s1UG)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   dldd <- eval(s1UG)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   dldm <- eval(m1UG)
                   dldd <- eval(s1UG)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dUG(y=y, mu=mu, sigma=sigma)),
                 rqres = expression(
                   rqres(pfun = "pUG", type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 sigma.initial = expression(sigma<- rep(4, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}
#------------------------------------------------------------------------------------------
# density function
dUG<-function(y, mu = 0.7, sigma = 2.1, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  fy1 <- 1/y*dgamma(-log(y),sigma,mu^(1/sigma)/(1-mu^(1/sigma)))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
# integrate(dUG,0,1) # checking the pdf
#------------------------------------------------------------------------------------------
# cumulative distribution function
pUG<-function(q, mu = 0.7, sigma = 2.1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1<- 1-pgamma(-log(q),sigma,mu^(1/sigma)/(1-mu^(1/sigma)))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
# pUG(.5)
# integrate(dUG,0,.5) # checking the cdf with the pdf
#------------------------------------------------------------------------------------------
# quantile function
qUG<-function(u,mu,sigma)
{
  q<- exp(-qgamma(1-u,sigma,mu^(1/sigma)/(1-mu^(1/sigma))))
  q
}
# u=pUG(.5)
# qUG(u,mu=.7,sigma=2.1) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for random generation
rUG<-function(n,mu,sigma)
{
  u<- runif(n)
  y<- qUG(u,mu =mu, sigma =sigma)
  y
}
# Checking the results
library(gamlss)

set.seed(10)
n<-1000
# Case 1: without regressors
mu_true<-.7
sigma_true<-7
mu_result<-sigma_result<-c()
for (i in 1:100) {
  y<-rUG(n,mu_true,sigma_true)
  fit1<-gamlss(y~1, family="UG", trace = F)
  logit_link<-make.link("logit")
  mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)
  sigma_result[i]<-fit1$sigma.coefficients
}
result1<- matrix(c(mu_true, mean(mu_result),
                   sigma_true, mean(sigma_result)),2,2)
colnames(result1)<-c("mu","sigma")
rownames(result1)<-c("true value","mean")
print(round(result1,2))

##########################################################################

UL<-expression(log(
  sigma*log(2)/(y*(-log(mu)))*(log(y)/log(mu))^(sigma-1)*2^(-(log(y)/log(mu))^(sigma))
)
)
m1UW<-D(UW,"mu")
s1UW<-D(UW,"sigma")
ms2UW<-D(m1UW,"sigma")


# density function
dUW<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  fy1 <- sigma*log(2)/(y)*(-log(mu))^(-1)*(log(y)/log(mu))^(sigma-1)*
    2^(-(log(y)/log(mu))^sigma)
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}

pUW<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1<-  .5^((log(q)/log(mu))^sigma)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}