require(MuMIn)
library(unmarked)
library(AICcmodavg)
library(lattice)
library(ROCR)
library(arm)
library(MASS)
library(Matrix)
library(lme4)
library(car)
library(psych)
library(PerformanceAnalytics)
library(pROC)

rm(list = ls())

##section1----correlation test####

#correlation test
data_dhole <- read.csv('dhole.csv')
data_mdeer <- read.csv('lesser mouse deer.csv')
data_serow <- read.csv('Chinese serow.csv')
data_wboar <- read.csv('wild boar.csv')
data_muntj <- read.csv('red muntjac.csv')
data_sambar <- read.csv('sambar.csv')

data_dhole$label <- 1
for (i in 1:238) {
  if (is.na(data_dhole[i,2]) &&is.na(data_dhole[i,3]) &&is.na(data_dhole[i,4]) &&
      is.na(data_dhole[i,5]) &&is.na(data_dhole[i,6]) &&is.na(data_dhole[i,7]) &&
      is.na(data_dhole[i,8]) &&is.na(data_dhole[i,9]) &&is.na(data_dhole[i,10]) &&
      is.na(data_dhole[i,11]) &&is.na(data_dhole[i,12]) &&is.na(data_dhole[i,13]) &&
      is.na(data_dhole[i,14]) &&is.na(data_dhole[i,15]) &&is.na(data_dhole[i,16]) &&
      is.na(data_dhole[i,17]) &&is.na(data_dhole[i,18]) &&is.na(data_dhole[i,19]) &&
      is.na(data_dhole[i,20]) &&is.na(data_dhole[i,21])) {
    data_dhole$label[i] <- 0
  }
}

env <- read.csv('env_factors.csv')
env$label <- data_dhole$label

#for (i in 1:238) {
#  if (is.na(env[i,8])) {
#    env$label[i] <- 0
#  }
#}

#data_dhole$label <- env$label
data_mdeer$label <- data_dhole$label
data_muntj$label <- data_dhole$label
data_sambar$label <- data_dhole$label
data_serow$label <- data_dhole$label
data_wboar$label <- data_dhole$label

label <- env$label
new_env <- env[label == 1,]
new_dhole <- data_dhole[label == 1,2:21]
new_mdeer <- data_mdeer[label == 1,2:21]
new_muntj <- data_muntj[label == 1,2:21]
new_sambar <- data_sambar[label == 1,2:21]
new_serow <- data_serow[label == 1,2:21]
new_wboar <- data_wboar[label == 1,2:21]

data_cortest <- new_env[,c(5:18,21)]
#mychart.Correlation <- edit(chart.Correlation)
mychart.Correlation(data_cortest, method = 'pearson',histogram=TRUE, pch=20, col = 'grey35')

env_data <- new_env[,c(5,6,12,15,18,21)]
mychart.Correlation(env_data, method = 'pearson',histogram=TRUE, pch=20, col = 'grey35')

##section2----import and preprocess raw data####

#import raw data
y <- list(dhole = as.matrix(new_dhole), mdeer = as.matrix(new_mdeer),
          muntj = as.matrix(new_muntj), sambar = as.matrix(new_sambar),
          serow = as.matrix(new_serow), wboar = as.matrix(new_wboar))
xraw <- env_data
obsCovs <- list(dhole = as.matrix(new_dhole))

y0 <- list(mdeer = as.matrix(new_mdeer),
           muntj = as.matrix(new_muntj), sambar = as.matrix(new_sambar),
           serow = as.matrix(new_serow), wboar = as.matrix(new_wboar))

colnames(xraw) <- c('s1','s2','s3','s4','s5','s6')

#standardize raw data
xraw$s1 <- (xraw$s1-min(xraw$s1))/(max(xraw$s1)-min(xraw$s1))
xraw$s2 <- (xraw$s2-min(xraw$s2))/(max(xraw$s2)-min(xraw$s2))
xraw$s3 <- (xraw$s3-min(xraw$s3))/(max(xraw$s3)-min(xraw$s3))
xraw$s4 <- (xraw$s4-min(xraw$s4))/(max(xraw$s4)-min(xraw$s4))
xraw$s5 <- (xraw$s5-min(xraw$s5))/(max(xraw$s5)-min(xraw$s5))
xraw$s6 <- (xraw$s6-min(xraw$s6))/(max(xraw$s6)-min(xraw$s6))

summary(xraw)
##section3----select det_var####

#select Detection variable with global site variables
testD <- unmarkedFrameOccuMulti(y=y0, siteCovs=xraw, obsCovs=obsCovs)
plot(testD)
summary(testD)

summary(testD@fDesign)

det_fm <- rep('~dhole',5)
det_fm0 <- c('~1','~1','~1','~1','~1')
occu_fm1 <- c('~.','~.','~.','~.','~.')


occu1 <- occuMulti(det_fm, occu_fm1, data=testD, maxOrder = 1, control=list(maxit=10000))
null <- occuMulti(det_fm0, occu_fm1, data=testD, maxOrder = 1, control=list(maxit=10000))
modList <- fitList(D1 = occu1, Null = null)
modSel(modList,nullmod = 'Null')
#null is the best fit

set.seed(123)
mod_penalty <- optimizePenalty(null)
summary(mod_penalty)

occu_fm <- c('~.','~.','~.','~.','~.')
nx <- dim(xraw)[2]
for (i in 1:nx) {
#  i <- 1
  index <- combn(1:nx,i)
  list.out <- list()
  model.names <- character(dim(index)[2])
  
  for (j in 1:dim(index)[2]) {
    UMF.t <- unmarkedFrameOccuMulti(y=y0, siteCovs=data.frame(xraw[,index[,j]]), 
                                    obsCovs=obsCovs)
#    list.out[j] <- optimizePenalty(occuMulti(det_fm0, occu_fm, UMF.t, maxOrder = 1, control=list(maxit=10000)))
    list.out[j] <- occuMulti(det_fm0, occu_fm, UMF.t, maxOrder = 1, 
                             control=list(maxit=10000))
    model.names[j] <- paste(paste('s', index[,j], sep=''), collapse=",")
    print(paste(j,'/',dim(index)[2]))
  }
  
  
  if (i==1) {aic.table <- aictab(list.out, modnames=model.names)}
  else {aic.table <- rbind(aic.table, aictab(list.out, modnames=model.names))}
}
aic.table2 <- aic.table[order(aic.table$AICc),]
aic.table2 <- data.frame(aic.table2)
names(aic.table2)[1] <- 'model'
head(aic.table2)
aic.table3 <- aic.table2[,c('model','K','LL','AICc')]
aic.table.out <- transform(aic.table3, delta.AICc=AICc-min(AICc))
aic.table.out[c(1:20),]



predUMF <- unmarkedFrameOccuMulti(y=y0, siteCovs=xraw, obsCovs=obsCovs)
occupred_fm1 <- rep('~s1+s2+s4',5)
occupred_fm2 <- rep('~s1+s2+s4+s5',5)
occupred_fm3 <- rep('~s1+s4',5)
occupred_fm4 <- rep('~s1+s4+s5',5)

pred1 <- occuMulti(det_fm0, occupred_fm1, data=predUMF, maxOrder = 1)
pred2 <- occuMulti(det_fm0, occupred_fm2, data=predUMF, maxOrder = 1)
pred3 <- occuMulti(det_fm0, occupred_fm3, data=predUMF, maxOrder = 1)
pred4 <- occuMulti(det_fm0, occupred_fm4, data=predUMF, maxOrder = 1)

set.seed(123)
penalty1 <- optimizePenalty(pred1)
penalty2 <- optimizePenalty(pred2)
penalty3 <- optimizePenalty(pred3)
penalty4 <- optimizePenalty(pred4)

modList <- fitList(pred1, pred2, pred3, pred4,
                   penalty1, penalty2, penalty3, penalty4)
coef(modList)


psi1 <- predict(pred1, type = 'state', species = 'mdeer')
psi1$Predicted

modavg(list(pred1,pred2,pred3,pred4), 
       parm = c('[mdeer] (Intercept)'), 
       parm.type = 'psi', exclude = list('s1','s2','s4','s5'))

modavg(list(penalty1,penalty2,penalty3,penalty4), 
       parm = c('[mdeer] (Intercept)'), 
       parm.type = 'psi', exclude = list('s1','s2','s4','s5'))


####for mdeer####
x <- xraw
predUMF <- unmarkedFrameOccuMulti(y=y0, siteCovs=x, obsCovs=obsCovs)
par(mfrow=c(2,3))
#Hold everything constant except for elev
newDats1 <- data.frame(s1 = seq(min(x$s1), max(x$s1), length=100), 
                       s2 = mean(x$s2), s3 = mean(x$s3), s4 = mean(x$s4), 
                       s5 = mean(x$s5), s6 = mean(x$s6))
occu_fm <- rep('~s1+s2+s4',5)
pred<-occuMulti(det_fm0, occu_fm, data=predUMF, maxOrder = 1)
Epreds1 <- as.data.frame(predict(pred, newdata=newDats1, type="state", species = 'mdeer'))
with(Epreds1, {
  s1 <- newDats1$s1*(max(env_data$elev)-min(env_data$elev))+min(env_data$elev)
  plot(s1, Predicted, xlab="Elevation", ylab="Occurrence probability_mouse deer", 
       pch=16, cex.lab=1, xlim=c(800,1800), ylim=c(0,1),type="l", lwd=2)
  lines(s1, Predicted-SE, lty=3)
  lines(s1, Predicted+SE, lty=3)
})

######
pdata <- read.csv('pred_data.csv')
head(pdata)

praw <- pdata[,5:8]
praw$road <- (praw$road-min(env_data$dist_road))/(max(env_data$dist_road)-min(env_data$dist_road))
praw$elev <- (praw$elev-min(env_data$elev))/(max(env_data$elev)-min(env_data$elev))
praw$slope <- (praw$slope-min(env_data$slope))/(max(env_data$slope)-min(env_data$slope))
praw$map <- (praw$map-min(env_data$map))/(max(env_data$map)-min(env_data$map))
summary(praw)

predUMF <- unmarkedFrameOccuMulti(y=y0, siteCovs=x, obsCovs=obsCovs)
set.seed(123)
occu_fm_pred <- rep('~s1+s2+s4+s5',5)
Ppred <- occuMulti(det_fm0, occu_fm_pred, data=predUMF, maxOrder = 1)
Ppred0 <- optimizePenalty(Ppred)

praw$s1 <- praw$elev
praw$s2 <- praw$slope
praw$s4 <- praw$road
praw$s5 <- praw$map

Poccu_mdeer <- predict(Ppred0, newdata=praw, type="state", species = 'mdeer')
Poccu_muntj <- predict(Ppred0, newdata=praw, type="state", species = 'muntj')
Poccu_sambar <- predict(Ppred0, newdata=praw, type="state", species = 'sambar')
Poccu_serow <- predict(Ppred0, newdata=praw, type="state", species = 'serow')
Poccu_wboar <- predict(Ppred0, newdata=praw, type="state", species = 'wboar')

write.csv(cbind.data.frame(Poccu_mdeer$Predicted,Poccu_muntj$Predicted,Poccu_sambar$Predicted,
            Poccu_serow$Predicted,Poccu_wboar$Predicted), file = 'pred_result.csv')

######
testS <- unmarkedFrameOccuMulti(y=y, siteCovs=xraw, obsCovs=obsCovs)

s1normal<-occu(~1~s1, data=civettestS,control=list(maxit=10000))
s1log<-occu(~1~log(s1+1), data=civettestS,control=list(maxit=10000))
s1sqrt<-occu(~1~sqrt(s1), data=civettestS,control=list(maxit=10000))
modList <- fitList(S1normal = s1normal, S1log=s1log, S1sqrt = s1sqrt)
modSel(modList)

civettestD@siteCovs$s1 <- sqrt(civettestD@siteCovs$s1)

summary(civettestD)