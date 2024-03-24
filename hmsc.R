library(Hmsc)

rm(list = ls())

# import data
det <- read.csv('detection.csv',row.names = 1)
env <- read.csv('env_new.csv',row.names = 1)
tr <- read.csv('trait.csv',row.names = 1)
tr$logm <- log(tr$mass)
tr$type <- as.factor(tr$type)

y <- as.matrix(det)
xraw <- env[,3:9]
mychart.Correlation(xraw,method = 'pearson',histogram=TRUE, pch=20, col = 'grey35')
xdata <- env[,c(3:7,9)]
colnames(xdata) <- c('s1','s2','s3','s4','s5','s6')
trdata <- tr[,c(1,3)]

Design <- data.frame(PA = as.factor(env$Reserve))
rL <- HmscRandomLevel(units = Design$PA)

# fit the model
x_fm = ~s1+s2+s3+s4+s5+s6
tr_fm = ~type+logm

mod <- Hmsc(Y = y, 
            XData = xdata, XFormula = x_fm,
            TrData = trdata, TrFormula = tr_fm,
            studyDesign = Design, ranLevels = list(PA = rL), 
            distr='probit')

# MCMC
nchains = 4
thin = 10
samples = 10000
transient = 50000
verbose = T

set.seed(704)
fit <- sampleMcmc(mod, thin = thin, samples = samples, transient = transient,
           nChains = nchains, verbose = verbose)

# explanatory power
preds <- computePredictedValues(fit)
eva <- evaluateModelFit(fit,preds)

# check convergence
post <- convertToCodaObject(fit)
#ess.beta & psrf.beta
es.beta <- effectiveSize(post$Beta)
ge.beta <- gelman.diag(post$Beta,multivariate=FALSE)$psrf
es.gamma <- effectiveSize(post$Gamma)
ge.gamma <- gelman.diag(post$Gamma,multivariate=FALSE)$psrf
#es.rho <- effectiveSize(post$Rho)
#ge.rho <- gelman.diag(post$Rho,multivariate=FALSE)$psrf
es.V <- effectiveSize(post$V)
ge.V <- gelman.diag(post$V,multivariate=FALSE)$psrf
hist(es.beta)
hist(es.gamma)
hist(es.V)
mean(ge.beta)
mean(ge.gamma)
mean(ge.V)

# predictive power
part <- createPartition(fit, nfolds = 3, column = 'PA')
preds1 <- computePredictedValues(fit, expected=FALSE, partition = part)
eva1 <- evaluateModelFit(fit, preds1)

# contribution of variance
var_part <- computeVariancePartitioning(fit)

# plot
post_beta <- getPostEstimate(fit, 'Beta')
plotBeta(fit, post=post_beta, param="Support", spNamesNumbers=c(TRUE,FALSE))

plotGamma(fit, post=post_beta, param="Support", trNamesNumbers=c(TRUE,TRUE))


library(corrplot)
OmegaCor <- computeAssociations(fit)
sptLevel <- 0.90
for (r in 1:fit$nr){
  plotOrder <- corrMatOrder(OmegaCor[[r]]$mean, order = "AOE")
  toPlot <- ((OmegaCor[[r]]$support > sptLevel) +
              (OmegaCor[[r]]$support < (1-sptLevel)) > 0) * OmegaCor[[r]]$mean
  par(xpd = T)
  colnames(toPlot) = rownames(toPlot) = gsub("_", " ", x = colnames(toPlot))
  corrplot(toPlot[plotOrder, plotOrder], method = "color",
           col = colorRampPalette(c("blue","white","red"))(200),
           title = "", type = "lower", tl.col = "black", tl.cex = .7, mar = c(0,0,6,0))
}



corrplot(OmegaCor[[r]]$mean, method = "color",order = 'FPC',
         col = colorRampPalette(c("blue","white","red"))(200),
         title = "", type = "lower", tl.col = "black", tl.cex = .7, mar = c(0,0,6,0))



Gradient <- constructGradient(fit, focalVariable = "s4")
predY_s4 <- predict(fit, XData = Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)

par(mfrow=c(1,3))
plotGradient(fit, Gradient, pred=predY_s4, measure="S", las=1,
             showData = TRUE, main='Species richness')
plotGradient(fit, Gradient, pred=predY_s4, measure="Y", index=1, las=1,
             showData = TRUE, main='Focal species occurrence')
plotGradient(fit, Gradient, pred=predY_s4, measure="T", index=1, las=1,
             showData = TRUE, main='Mean trait value')

