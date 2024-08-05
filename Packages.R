rm(list=ls())
# Packages
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")
}
library("tidyverse")

### VIF
if("car" %in% rownames(installed.packages()) == FALSE) {install.packages("car")
}
library(car)

### Shrinkage regression
if("glmnet" %in% rownames(installed.packages()) == FALSE) {install.packages("glmnet")
}
library(glmnet)


### Dunnet
if("DescTools" %in% rownames(installed.packages()) == FALSE) {install.packages("DescTools")
}
library(DescTools)

### lme4
#oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
if("lme4" %in% rownames(installed.packages()) == FALSE) {install.packages("lme4")
}
library(lme4)
library(lmerTest)
#options(oo)

# Regression Analaysis



## 2. Calculation
n <- 100
x <- rnorm(n, mean=1)
y <- 1 + 0.5 * x + rnorm(n)
fit1 <- lm(y ~ x)
summary(fit1)

### beta hat
beta0_hat <- fit1$coefficients[1]
beta1_hat <- fit1$coefficients[2]

### sigmahat
sigmasqhat <- sum(fit1$residuals^2)/fit1$df.residual
sqrt(sigmasqhat)
ssx <- sum(scale(x, center=T, scale=F)^2)
se_beta1hat <- sqrt(sigmasqhat/ssx)
se_beta1hat
se_beta0hat <- sqrt(sigmasqhat*(1/n + mean(x)^2/ssx))
se_beta0hat

### upper and lower bound
alpha <- 0.05
beta0_hat  - se_beta0hat * qt(1-alpha/2, df=fit1$df.residual); beta0_hat  + se_beta0hat * qt(1-alpha/2, df=fit1$df.residual)
beta1_hat  - se_beta1hat * qt(1-alpha/2, df=fit1$df.residual); beta1_hat  + se_beta1hat * qt(1-alpha/2, df=fit1$df.residual)

### p_value beta (H_0: beta1=0, two sided test)
2*pt((beta1_hat -0)/se_beta1hat, df=fit1$df.residual, lower.tail = F)



### hat matrix
X <- matrix(c(rep(1,5), 8,4,0,-4,-8), ncol=2) 
Y <- c(7.8,9.0,10.2,11.0,11.7)
H <- X%*%solve(t(X)%*%X)%*%t(X)
H%*%Y
lm(Y~X)$fitted.values

### Working-Hotelling Covariance
p <- 2
W <- sqrt(2*qf(1-alpha, df1=p, df2=n-p))

### Joint confidence intervals
beta0_hat  - se_beta0hat * W; beta0_hat  + se_beta0hat * W
beta1_hat  - se_beta1hat * W; beta1_hat  + se_beta1hat * W

### Bonferroni correction 
alpha2 <- alpha/p
beta0_hat - se_beta0hat * qt(1-alpha2/2, df=fit1$df.residual); beta0_hat  + se_beta0hat * qt(1-alpha2/2, df=fit1$df.residual)
beta1_hat - se_beta1hat * qt(1-alpha2/2, df=fit1$df.residual); beta1_hat  + se_beta1hat * qt(1-alpha2/2, df=fit1$df.residual)

### Mean response
hat_y2 <- predict(fit1, data.frame(x=20))
se_haty2 <- sqrt(sigmasqhat*(1/n + (2-mean(x))^2/ssx))
se_haty2

hat_y2 - se_haty2*qt(1-alpha/2, df=fit1$df.residual); hat_y2 + se_haty2*qt(1-alpha/2, df=fit1$df.residual)
hat_y2 - se_haty2*W; hat_y2 + se_haty2*W

### Diagnosis
#### Plots
par(mfrow=c(1,2))
plot(fit1, which = c(1,2))

#### DFFITS
dffit1 <- dffits(fit1)

dffit1[which.max(dffit1)]
dfbeta1 <- dfbetas(fit1)
dfbeta1


#### Box-cox
library(MASS)
bc <- boxcox(Y~X)
lambda <- bc$x[which.max(bc$y)]
lambda

Y_star <- (Y^lambda-1)/lambda

### 3. Example 1
data("mtcars")
df_ex1 <- mtcars %>% select(qsec, disp, cyl, wt, mpg)
df_ex1$cyl <- as.factor(df_ex1$cyl)
attach(df_ex1)
grid <- seq(from=min(df_ex1$disp), to=max(disp), l=1e2) 
pred_mean <- mtcars %>% select(wt, mpg) %>% apply(2, FUN=mean) 
newdata <- cbind(matrix(c(grid, rep(pred_mean, each=1e2)),
                        ncol=1+length(pred_mean))) %>% as.data.frame()
colnames(newdata) = c("disp", "wt", "mpg")

fit2 <- lm(qsec ~ disp + mpg + wt, data=df_ex1)
                                                                                                                                                                              
pred_line <- predict(fit2, as.data.frame(newdata))
par(mar = c(4,4,1,1))
plot(qsec ~ disp, pch=16, data=df_ex1, xlab = "disp",
     ylab="qsec")

#### the curve of the predicted line is added
lines(grid, pred_line, lty=1, lwd=2, col=2)

#### the 95% Working-Hotelling Confidence band
alpha <- 0.05
n <- nrow(mtcars)
p <- 3
MSE1_2 <- sum(fit2$residuals^2)/fit2$df.residual
summary(fit2)

X <- model.matrix(~ disp + mpg + wt, data=df_ex1)
X_grid <- model.matrix(~ disp + mpg + wt, data=newdata)
H_grid <- X_grid%*%solve(t(X)%*%X)%*%t(X_grid)
se_pred <- sqrt(MSE1_2*diag(H_grid))
W2 <- sqrt(2*qf(1-alpha, df1=p, df2=n-p))
se_pred
H_grid 

pred_line_W_upper <- pred_line + se_pred*W2
pred_line_W_lower <- pred_line - se_pred*W2
lines(grid, pred_line_W_upper , lty=2, lwd=1, col=2) 
lines(grid, pred_line_W_lower, lty=2, lwd=1, col=2)
legend("topright", legend = c("Fitted Values", "the 95% Working-Hoteling Confidence Band"), 
       col=2, lty = 1:2, lwd = 1)

### Grid expansion matrix
fit3 <- lm(qsec ~ disp + mpg + wt + cyl, data=df_ex1)
grid <- seq(from=min(disp), to=max(disp), l=1e2) 
nlevel_cyl <- length(levels(df_ex1$cyl))
newdata2 <- cbind(expand.grid(disp=grid, cyl=levels(df_ex1$cyl)), matrix(rep(pred_mean, each=1e2*nlevel_cyl), ncol=length(pred_mean), byrow=F))
colnames(newdata2) = c("disp", "cyl", "wt", "mpg")
pred_mat <- matrix(predict(fit3, newdata2), ncol=nlevel_cyl)

plot(qsec ~ disp, pch=16, col=cyl, data=df_ex1, xlab = "disp",
     ylab="qsec")
legend("topright", pch=16, col=1:nlevel_cyl,
       legend=levels(df_ex1$cyl), title="cyl")
matlines(grid, pred_mat, type = "l", lty=1, lwd=2, col=1:nlevel_cyl, xlab = "disp", ylab="qsec")


X <- model.matrix(qsec~., data=df_ex1) 
X_grid <- model.matrix(~., data=newdata2) 
H_grid <- X_grid%*%solve(t(X)%*%X)%*%t(X_grid) 
df1 <- ncol(X)
df2 <- fit3$df.residual
MSE1_3 <- sum(fit3$residuals^2)/df2
se_pred <- sqrt(MSE1_3*diag(H_grid)) 
alpha <- 0.05
W2 <- sqrt(2 * qf(1-alpha, df1, df2))

pred_line_W_upper <- pred_mat + matrix(se_pred*W2, ncol=nlevel_cyl, byrow=F)
pred_line_W_lower <- pred_mat - matrix(se_pred*W2, ncol=nlevel_cyl, byrow=F)

for(j in 1:nlevel_cyl){
  polygon(x=c(grid, rev(grid)), y=c(pred_line_W_lower[,j], rev(pred_line_W_upper[,j])), col=adjustcolor(j, 0.3))
}

### Diagnosis
vif(fit3)

### Variable selection
back3 <- step(fit3, direction="backward") ## AIC
stepwise3 <- step(fit3, direction="both")

### Ridge and Lasso
coefficients(fit3) ### Ordinary LS

ridge_cvfit <- cv.glmnet(X[,-1], df_ex1$qsec, nfolds = 5, alpha = 0) ## Ridge pen
lambda <- ridge_cvfit$lambda.min; cv_error <- min(ridge_cvfit$cvsd)
lambda
bestridge_fit <-glmnet(X[,-1], df_ex1$qsec, nfolds = 5, alpha = 0, lambda=lambda) 
coef.glmnet(bestridge_fit)

lasso_cvfit <- cv.glmnet(X[,-1], df_ex1$qsec, nfolds = 5, alpha = 1) ## LASSO pen
lambda <- lasso_cvfit$lambda.min; cv_error <- min(lasso_cvfit$cvsd)
lambda
bestlasso_fit <-glmnet(X[,-1], df_ex1$qsec, nfolds = 5, alpha = 1, lambda=lambda) 
coef.glmnet(bestlasso_fit)


# Experimental Design

## 1. Sampling distributions

ybars <- sapply(1:200, function(rep){
  mean(rnorm(n=50, mean=2, sd=3))
})
hist(ybars, probability = T, breaks="fd")
lines(seq(0, 4, l=1e2), dnorm(seq(0,4, l=1e2), mean=2, sd=3/sqrt(50)))


## 2. Testing difference in means
data(iris)
df1 <- iris$Sepal.Length
t.test(df1[iris$Species=="setosa"], df1[iris$Species=="versicolor"], var.equal=F)


### Checking assumptions

Fstat <- var(df1[iris$Species=="setosa"]) / var(df1[iris$Species=="versicolor"])
n1 = length(df1[iris$Species=="setosa"])
n2 = length(df1[iris$Species=="versicolor"])

### p-value
2*min(pf(Fstat, df1=n1-1, df2=n2-1), pf(1/Fstat, df1=n1-1, df2=n2-1))
pf(Fstat, df1=n1-1, df2=n2-1)


## 3. ANOVA
data("mtcars")

fit2 <- lm(qsec~as.factor(cyl), data=mtcars)
summary(fit2)

### Checking assumptions
qqnorm(resid(fit2))
qqline(resid(fit2))
plot(fit2, which=2)

plot(fitted.values(fit2), resid(fit2))
plot(fit2, which=1)

bartlett.test(qsec~as.factor(cyl), data=mtcars)
library(car)
leveneTest(qsec~as.factor(cyl), data=mtcars)

### Mean vs sds plot
sds <- aggregate(qsec~as.factor(cyl), data=mtcars, FUN=sd)
means <- aggregate(qsec~as.factor(cyl), data=mtcars, FUN=mean)
par(mar=c(4,4,1,1))
plot(log(means$qsec), log(sds$qsec), xlab="log-mean", ylab="log-SD")

### Variance-stabilizing Transformation

### alpha : slope, lambda = 1 - alpha : power transformation


### Post Hoc
library(DescTools)
PostHocTest(aov(qsec~as.factor(cyl), data=mtcars), method = "hsd")
TukeyHSD(aov(qsec~as.factor(cyl), data=mtcars), which=1, conf.level=1-alpha)

### ANOVA test
library(car)
anova(fit2, fit3)


### Refit with transformed data
fit3 <- lm(sqrt(qsec)~as.factor(cyl), data=mtcars)
summary(fit3)

par(mfrow=c(1,2))
plot(fit3, which=c(1,2))


### Non-parametric inference
kruskal.test(qsec~as.factor(cyl), data=mtcars)

### Example 2: Determining Sample Size
df_ex2 <- data.frame(times = c(65, 81, 57, 66, 82, 82, 67, 59, 75, 70, 
                            64, 71, 83, 59, 65, 56, 69, 74, 82, 79), 
                  type=factor(rep(1:2, each=10)))
fit3_2 <- aov(times~type, data=df_ex2)
fit3_2_sum <- summary(fit3_2)
MSE3 <- unlist(fit3_2_sum)[6]
n <- 10
alpha <- 0.05
F_alpha <- qf(1-alpha, df1=1, df2=18)

ncp_3 <- ((1)^2+(-1)^2) / (MSE3/n) ## Non-centrality parameter with an error variance estimate
power3 <- pf(F_alpha, df1=1, df2=18, ncp=ncp_3, lower.tail = F) 

power3_2 <- c()
for(n in 1801:1850){
  F_alpha2 <- qf(1-alpha, df1=1, df2=2*(n-1))
  ncp_3_2 <- ((0.5)^2+(-0.5)^2) / (MSE3/n) ## Non-centrality parameter with an error variance estimate
  power3_2 <- c(power3_2, pf(F_alpha2, df1=1, df2=2*(n-1), ncp=ncp_3_2, lower.tail = F))
}

power.anova.test(groups=2, between.var=((0.5)^2+(-0.5)^2), within.var = MSE3, sig.level = 0.05, power=0.9)

which(power3_2 > 0.9)

## 4. Random effect models

### MME
### Fixed model
#### MS_A / MS_E

### Mixed model if block effect is random
#### MS_A / MS_AB

#### hat sigma^2_block = 1/(an) * (MS_{B} - MS_{AB})
#### hat sigma^2_{A X B} = 1/n * (MS_{AB} - MS_{E})

### REML
library(lme4)
fit4 <- lmer(Yield~as.factor(Pressure)+as.factor(Temperature)+(1 | Day))
ranova()

## 5. Factorial Design
fit4 <- lm(Yield~A*B*C*D) ## All orthogonal design 
summary(fit4)

## Basic -> Interaction plot -> Interaction test -> Multiplicative model fit
## Multiplicative model fit -> Reduce term (Especially Single Replicate Experiment)

## 6. Nested Design

lmer(Y~Treatment+(1|Batch:Treatment))

aov(Y~Treatment+Batch:Treatment)


# Others

df1[df1$a %in% c(1,3), ]