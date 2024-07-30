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

sds <- aggregate(qsec~as.factor(cyl), data=mtcars, FUN=sd)
means <- aggregate(qsec~as.factor(cyl), data=mtcars, FUN=mean)
par(mar=c(4,4,1,1))
plot(log(means$qsec), log(sds$qsec), xlab="log-mean", ylab="log-SD")


### Refit with transformed data
fit3 <- lm(sqrt(qsec)~as.factor(cyl), data=mtcars)
summary(fit3)

par(mfrow=c(1,2))
plot(fit3, which=c(1,2))

### Non-parametric inference
kruskal.test(qsec~as.factor(cyl), data=mtcars)

### Determining Sample Size
df2 <- data.frame(times = c(65, 81, 57, 66, 82, 82, 67, 59, 75, 70, 
                              64, 71, 83, 59, 65, 56, 69, 74, 82, 79), 
                    type=factor(rep(1:2, each=10)))
fit4 <- aov(times~type, data=df2)
fit4_sum <- summary(fit4)
MSE4 <- unlist(fit4_sum)[6]
n <- 10
alpha <- 0.05
F_alpha <- qf(1-alpha, df1=1, df2=18)

ncp_4 <- ((1)^2+(-1)^2) / (MSE4/n) ## Non-centrality parameter with an error variance estimate
power4 <- pf(F_alpha, df1=1, df2=18, ncp=ncp_4, lower.tail = F) 

power4_2 <- c()
for(n in 1801:1850){
  F_alpha2 <- qf(1-alpha, df1=1, df2=2*(n-1))
  ncp_4_2 <- ((0.5)^2+(-0.5)^2) / (MSE4/n) ## Non-centrality parameter with an error variance estimate
  power4_2 <- c(power4_2, pf(F_alpha2, df1=1, df2=2*(n-1), ncp=ncp_4_2, lower.tail = F))
}

which(power4_2 > 0.9)

