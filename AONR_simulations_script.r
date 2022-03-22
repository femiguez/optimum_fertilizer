## This scripts contains the code used in a publication
## "How can we estimate optimum fertilizer rates with accuracy and precision?"

library(ggplot2)
library(nlraa)
library(car)
library(patchwork)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
options(warn = -1)

#### Chunk: quadp3xs-true ----
## Set up nitrogen rates
nrates <- seq(0, 300, by = 50)

set.seed(123) ## Set the seed for reproducibility

## Here I'm assuming that both the yield and nrate are in kg/ha
yield <- replicate(5, quadp3xs(nrates, 6000, 50, 200) + rnorm(length(nrates), 0, 1000))

dat <- data.frame(nrate = nrates, yield = c(yield))

fm.LP <- nls(yield ~ SSlinp(nrate, a, b, xs), data = dat)
fm.QP <- nls(yield ~ SSquadp3xs(nrate, a, b, xs), data = dat)

new.nrates <- data.frame(nrate = seq(0, 300))
prds.LP <- predict_nls(fm.LP, newdata = new.nrates)
prds.QP <- predict_nls(fm.QP, newdata = new.nrates)

## Bootstrap confidence intervals
fm.LP.bt <- quiet(boot_nls(fm.LP))
fm.LP.bt.ci <- suppressWarnings(confint(fm.LP.bt)[3,])
fm.QP.bt <- quiet(boot_nls(fm.QP))
fm.QP.bt.ci <- suppressWarnings(confint(fm.QP.bt)[3,])

#### Chunk: fig-quadp3xs-true ----
gp1 <- ggplot() + 
  geom_point(data = dat, aes(nrate, yield)) + 
  geom_point(aes(x = 200, y = 5500), color = "black", shape = 2) + 
  geom_text(aes(x = 200, y = 5800, label = "True AONR")) + 
  geom_point(aes(x = coef(fm.LP)[3], y = 5000), color = "#F8766D") + 
  geom_errorbarh(aes(y = 5000, xmin = fm.LP.bt.ci[1], xmax = fm.LP.bt.ci[2]), color = "#F8766D") + 
  geom_point(aes(x = coef(fm.QP)[3], y = 5250), color = "#00BFC4") +
  geom_errorbarh(aes(y = 5250, xmin = fm.QP.bt.ci[1], xmax = fm.QP.bt.ci[2]), color = "#00BFC4") + 
  geom_line(aes(y = prds.LP, x = new.nrates$nrate, color = "linear plateau")) +
  geom_line(aes(y = prds.QP, x = new.nrates$nrate, color = "quadratic plateau")) +
  xlab("Nitrogen rate (kg/ha)") + 
  ylab("Yield (kg/ha)") + 
  xlim(c(0,400)) + 
  ggtitle("a) Quadratic-plateau is the true model") + 
  guides(color = guide_legend(title = element_blank())) + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.45)) 

#### Chunk: linp-true ----
## Set up nitrogen rates
nrates <- seq(0, 300, by = 50)

set.seed(123) ## Set the seed for reproducibility

## Here I'm assuming that both the yield and nrate are in kg/ha
yield <- replicate(5, linp(nrates, 6000, 50, 200) + rnorm(length(nrates), 0, 1000))

dat <- data.frame(nrate = nrates, yield = c(yield))

fm2.LP <- nls(yield ~ SSlinp(nrate, a, b, xs), data = dat)
fm2.QP <- nls(yield ~ SSquadp3xs(nrate, a, b, xs), data = dat)

new.nrates <- data.frame(nrate = seq(0, 300))
prds2.LP <- predict_nls(fm2.LP, newdata = new.nrates)
prds2.QP <- predict_nls(fm2.QP, newdata = new.nrates)

## Bootstrap confidence intervals
fm2.LP.bt <- quiet(boot_nls(fm2.LP))
fm2.LP.bt.ci <- suppressWarnings(confint(fm2.LP.bt)[3,])
fm2.QP.bt <- quiet(boot_nls(fm2.QP))
fm2.QP.bt.ci <- suppressWarnings(confint(fm2.QP.bt)[3,])

#### Chunk: fig-linp-true ----
gp2 <- ggplot() + 
  geom_point(data = dat, aes(nrate, yield)) + 
  geom_point(aes(x = 200, y = 5500), color = "black", shape = 2) + 
  geom_text(aes(x = 200, y = 6000, label = "True AONR")) + 
  geom_point(aes(x = coef(fm2.LP)[3], y = 5000), color = "#F8766D") + 
  geom_errorbarh(aes(y = 5000, xmin = fm2.LP.bt.ci[1], xmax = fm2.LP.bt.ci[2]), color = "#F8766D") + 
  geom_point(aes(x = coef(fm2.QP)[3], y = 5250), color = "#00BFC4") +
  geom_errorbarh(aes(y = 5250, xmin = fm2.QP.bt.ci[1], xmax = fm2.QP.bt.ci[2]), color = "#00BFC4") + 
  geom_line(aes(y = prds2.LP, x = new.nrates$nrate, color = "linear plateau")) +
  geom_line(aes(y = prds2.QP, x = new.nrates$nrate, color = "quadratic plateau")) +
  xlab("Nitrogen rate (kg/ha)") + 
  ylab("Yield (kg/ha)") + 
  ggtitle("b) Linear-plateau is the true model") + 
  guides(color = guide_legend(title = element_blank())) + 
  xlim(c(0, 400)) + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.45)) 
gp1 / gp2

#### Chunk: sims-quadp3xs ----
## These simulations assume that the correct model is the QP3 model
n.nrates <- c(4, 6, 8, 10)
n.reps <- 3:10
n.sim <- 500

res <- data.frame(n.reps = NA, n.nrates = NA, i = NA, aonr = NA, aonr.se = NA, 
                  chosen.model = NA)

set.seed(123) ## Set the seed for reproducibility
k <- 1
for(i in seq_along(n.reps)){
  for(j in seq_along(n.nrates)){
    ## cat("N reps:", i, "N nrates:", j, "\n")
    for(q in 1:n.sim){
      nrates <- seq(0, 300, length.out = n.nrates[j])
      ## Here I'm assuming that both the yield and nrate are in kg/ha
      yield <- replicate(n.reps[i], quadp3xs(nrates, 6000, 50, 200) + rnorm(length(nrates), 0, 1000))
      dat <- data.frame(nrate = nrates, yield = c(yield))
      fm.LP <- try(nls(yield ~ SSlinp(nrate, a, b, xs), data = dat), silent = TRUE)
      fm.QP <- try(nls(yield ~ SSquadp3xs(nrate, a, b, xs), data = dat), silent = TRUE)
      
      res[k, "n.reps"] <- n.reps[i]
      res[k, "n.nrates"] <- n.nrates[j]
      res[k, "i"] <- q
      
      if(inherits(fm.LP, "try-error") || inherits(fm.QP, "try-error")){
        k <- k + 1
        next
      }
      ## Pick results from the chosen model
      if(AIC(fm.QP) < AIC(fm.LP)){
        best.model <- fm.QP
        res[k, "chosen.model"] <- 1 
      }else{
        best.model <- fm.LP
        res[k, "chosen.model"] <- 0 
      }
      res[k, "aonr"] <- coef(best.model)[3]  
      res[k, "aonr.se"] <- summary(best.model)$coefficients[3,2]  
      k <- k + 1
    }
  }
}

#### Chunk: (blank) ----
res.na <- na.omit(res)
pc.ag <- aggregate(chosen.model ~ n.reps + n.nrates, data = res, FUN = function(x) sum(x, na.rm = TRUE))
total.ag <- aggregate(chosen.model ~ n.reps + n.nrates, data = res, FUN = function(x) length(na.omit(x)))
pc.ag$total <- total.ag$chosen.model
pc.ag$percent.correct <- pc.ag$chosen.model / pc.ag$total * 100
avg.percent.correct.QP <- round(mean(pc.ag$percent.correct))

#### Chunk: fig-sims-quadp3xs ----
res2 <- res.na
res2$n.nrates.f <- as.factor(res2$n.nrates)
ggplot(data = res2, aes(n.reps, aonr.se, color = n.nrates.f)) + 
  geom_smooth(se = FALSE, formula = y ~ s(x, k = 4, bs = "cs"), method = "gam") + 
  ##geom_point(alpha = 1/3) + 
  ## ylim(c(0, 30)) + 
  xlab("Number of replications") + 
  ylab("Standard error of AONR (kg/ha)") + 
  guides(color = guide_legend(title = "Number N rates")) + 
  theme_bw() + 
  theme(legend.position = c(0.85, 0.8)) 

#### Chunk: fig-sims-quadp3xs-percent-correct ----
pc.ag$n.nrates.f <- as.factor(pc.ag$n.nrates)
ggplot(data = pc.ag, aes(n.reps, percent.correct, color = n.nrates.f)) + 
  geom_smooth(se = FALSE, formula = y ~ s(x, k = 4, bs = "cs"), method = "gam") + 
  ##geom_point(alpha = 1/3) + 
  ## ylim(c(0, 30)) + 
  xlab("Number of replications") + 
  ylab("Choosing correct model (%)") + 
  ggtitle("QP is the true model") + 
  guides(color = guide_legend(title = "Number N rates")) + 
  theme_bw() + 
  theme(legend.position = c(0.85, 0.5)) 

#### Chunk: quadp4 ----
set.seed(1234)

niter <- 5e3

bias.LP <- numeric(niter)
bias.QP <- numeric(niter)
bias.best <- numeric(niter)
bias.avg <- numeric(niter)
best.model <- numeric(niter)

dat.comb <- NULL

true.aonr <- 200

start <- Sys.time()
for(i in 1:niter){
  
  highest.nrate <- 350
  nrates <- seq(0, highest.nrate, length.out = 8)
  a.int <- runif(1, min = 4000, max = 10000)
  b.slp <- runif(1, min = 5, max = 70)
  c.qud <- runif(1, min = -0.05, max = -0.01) 
  yld <- replicate(5, quadp(nrates, a.int, b.slp, c.qud, true.aonr) + rnorm(length(nrates), 0, 1000))
  
  dat <- data.frame(nrate = nrates, yield = c(yld))
  
  fit.LP <- try(nls(yield ~ SSlinp(nrate, a, b, xs), data = dat), silent = TRUE)
  fit.QP <- try(nls(yield ~ SSquadp3xs(nrate, a, b, xs), data = dat), silent = TRUE)
  
  if(inherits(fit.LP, "try-error") || inherits(fit.QP, "try-error") || coef(fit.QP)[3] > 400){
    bias.LP[i] <- NA
    bias.QP[i] <- NA
    bias.avg[i] <- NA
    bias.best[i] <- NA
    next
  } 
  
  bias.LP[i] <- coef(fit.LP)[3] - true.aonr
  bias.QP[i] <- coef(fit.QP)[3] - true.aonr
  
  ic.weights <- IC_tab(fit.LP, fit.QP, sort = FALSE)
  
  bias.avg[i] <- c(coef(fit.LP)[3], coef(fit.QP)[3]) %*% ic.weights$weight - true.aonr
  
  if(AIC(fit.LP) < AIC(fit.QP)){
    bias.best[i] <- coef(fit.LP)[3] - true.aonr
    best.model[i] <- 1
  }else{
    bias.best[i] <- coef(fit.QP)[3] - true.aonr
    best.model[i] <- 2
  }
  
  dat.comb <- rbind(dat.comb, data.frame(i = i, dat))
}
## beepr::beep()
end <- Sys.time() ## Apparently this only takes about 2 minutes

bias.comb <- data.frame(model = rep(c("LP", "QP", "avg", "best"), each = niter), 
                        bias = c(bias.LP, bias.QP, bias.avg, bias.best))

mean.bias.comb <- aggregate(bias ~ model, data = bias.comb, FUN = mean, na.rm = TRUE)
median.bias.comb <- aggregate(bias ~ model, data = bias.comb, FUN = median, na.rm = TRUE)
sd.bias.comb <- aggregate(bias ~ model, data = bias.comb, FUN = sd, na.rm = TRUE)

#### Chunk: model-averaging-bias ----
ggplot(data = bias.comb, aes(x = bias, color = model)) + 
  geom_density() + 
  geom_vline(xintercept = 0) + 
  xlab("Bias (kg/ha)") + 
  theme_bw()