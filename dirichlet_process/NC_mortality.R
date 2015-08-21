dr <- read.csv("mortality.csv")
dr <- na.omit(dr)
View(dr)
dr$Year <- dr$Year - 1999
dr$Hospital <- 1e5 * dr$Total_Hospital /dr$Population
dr$Density <- dr$Population_Density/1e3
dr$Income <- dr$Median_Household_Income/1e4

drNC <- subset(dr, State == "NC")

nc.fips <- split(drNC, drNC$FIPS_Code)
hosp <- sapply(nc.fips, function(d) return(d$Hospital[1]))
dens <- sapply(nc.fips, function(d) return(d$Density[1]))
incm <- sapply(nc.fips, function(d) return(d$Income[1]))

dr.train <- subset(drNC, Year < 13) ## 2011 and before
dr.test <- subset(drNC, Year == 13) ## 2012

dr.cty <- split(dr.train, dr.train$FIPS_Code)
length(dr.cty)

lms.coef <- t(sapply(dr.cty, function(dd) return(coef(lm(Crude_Rate ~ Year, data = dd)))))

require(maps)
data(county.fips)
county.fips$state <- sapply(strsplit(as.character(county.fips$polyname), ","), function(ll) return(ll[1]))
colors = heat.colors(4)

cty.mapper <- function(x, fips, colors, st = 'north carolina', ...){
  m <- length(colors)
  q.x <- quantile(x, 0:m/m)
  x.bins <- as.numeric(cut(x, q.x, include.lowest = TRUE)) 
  ct.fips <- county.fips$fips
  if(nchar(st) > 0) ct.fips <- subset(county.fips, state == st)$fips
  match.county <- match(ct.fips, fips)
  map('county', st, col = colors[x.bins[match.county]], fill = TRUE, lty = 0, ...)
  q.x <- signif(q.x, 3)
  x.txt <- paste(q.x[-length(q.x)], rep(" - ",length(q.x)-1), q.x[-1], sep = "")
  legend("bottomright", x.txt, horiz = FALSE, fill = colors, cex = 0.7, bty = "n", ncol = 1)
}

par(mfrow = c(2,2))
xx <- lms.coef[,1:2]
for(i in 1:2){
  hist(xx[,i], 20, freq = FALSE, col = "gray", border = "white", main = "", ann = FALSE)
  g <- seq(min(xx[,i]), max(xx[,i]), len = 501)
  lines(g, dnorm(g, mean(xx[,i]), sd(xx[,i])))  
  title(main = dimnames(xx)[[2]][i])
  cty.mapper(xx[,i], dimnames(xx)[[1]], colors, myborder = 0.1)
  title(main = dimnames(xx)[[2]][i])
}

par(mfrow = c(1,1), mar = 5:2 + .1, bg = "cyan")
plot(dens, incm, col = gray(pnorm(scale(lms.coef[,1]))), pch = 19, cex = 0.8, bty = "n")


## lm full
lmf <- lm(Crude_Rate ~ Year*(Hospital+Density+Income), data = dr.train)


### nlme > lme()

require(nlme)
fit.lmm <- lme(fixed = Crude_Rate ~ Hospital+Density+Income, random = ~Year|FIPS_Code, data = dr.train, method = "ML")
lmm.coef <- coef(fit.lmm)

par(mfrow = c(2,2), bg = "white")
xx <- lmm.coef[,c(1,5)]
for(i in c(1,2)){
  hist(xx[,i], 20, freq = FALSE, col = "gray", border = "white", main = "", ann = FALSE)
  g <- seq(min(xx[,i]), max(xx[,i]), len = 501)
  lines(g, dnorm(g, mean(xx[,i]), sd(xx[,i])))  
  title(main = dimnames(xx)[[2]][i])
  cty.mapper(xx[,i], dimnames(xx)[[1]], colors, myborder = 0.1)
  title(main = dimnames(xx)[[2]][i])
}

### DPpackage > DPMlmm()

require(DPpackage)

m0 <- colMeans(lmm.coef)[c(1,5)]
S0 <- as.matrix(cov(lmm.coef))[c(1,5), c(1,5)]
prior <- list(alpha = 1e-2, tau1 = 0.01, tau2 = 0.01, nu0 = 4.01, tinv = S0,
              nub = 4.01, tbinv = S0, mb = m0, Sb = S0, beta0 = colMeans(lmm.coef)[2:4], Sbeta0 = diag(100,3))

mcmc <- list(nburn = 200, nsave = 2e3, nskip = 10, ndisplay = 200)
fit1 <- DPMlmm(fixed = Crude_Rate ~ Hospital+Density+Income, random = ~Year|FIPS_Code, prior=prior, mcmc=mcmc, state=NULL, status=TRUE, data = dr.train)
theta.samp <- fit1$save.state$randsave
dpm.coef <- t(matrix(colMeans(theta.samp), nrow = 2))
dimnames(dpm.coef)[[2]] <- c("(Intercept)", "Year")
dimnames(dpm.coef)[[1]] <- c(dimnames(lmm.coef)[[1]], "Pred")
betaF <- colMeans(fit1$save.state$thetasave[,3:5])  

par(mfrow = c(2,2), bg = "white")
xx <- dpm.coef[,c(1,2)]
for(i in c(1,2)){
  plot(density(xx[,i]), col = "royalblue", ann = FALSE, bty = "n", main = "")
  g <- seq(min(xx[,i]), max(xx[,i]), len = 501)
  lines(g, dnorm(g, mean(xx[,i]), sd(xx[,i])))  
  title(main = dimnames(xx)[[2]][i])
  cty.mapper(xx[,i], dimnames(xx)[[1]], colors, myborder = 0.1)
  title(main = dimnames(xx)[[2]][i])
}


## point prediction for 2012

lms.2012 <- c(lms.coef %*% c(1, 13))
lmf.2012 <- predict(lmf, newdata = dr.test)
lmm.2012 <- predict(fit.lmm, newdata = dr.test)
match.st <- match(dr.test$FIPS_Code, dimnames(dpm.coef)[[1]])
dpm.2012 <- sapply(1:nrow(dr.test), function(j) return(sum(as.numeric(c(1, dr.test[j,c("Year","Hospital","Density","Income")]))*c(dpm.coef[match.st[j],],betaF))))

plot(dr.test$Crude_Rate, lms.2012, pch = 19, cex = 0.4, col = "royalblue", bty = "n", ann = FALSE)
abline(0,1,col = 2)
plot(dr.test$Crude_Rate, lmf.2012, pch = 19, cex = 0.4, col = "royalblue", bty = "n", ann = FALSE)
abline(0,1,col = 2)
plot(dr.test$Crude_Rate, lmm.2012, pch = 19, cex = 0.4, col = "royalblue", bty = "n", ann = FALSE)
abline(0,1,col = 2)
plot(dr.test$Crude_Rate, dpm.2012, pch = 19, cex = 0.4, col = "royalblue", bty = "n", ann = FALSE)
abline(0,1,col = 2)

rmse <- function(x, y) return(sqrt(mean((x-y)^2)))
rmse(dr.test$Crude_Rate, lms.2012)
rmse(dr.test$Crude_Rate, lmf.2012)
rmse(dr.test$Crude_Rate, lmm.2012)
rmse(dr.test$Crude_Rate, dpm.2012)
