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

##section1----correlation test####

#correlation test
data_survey <- read.csv('civet.csv')
s_survey <- data_survey[,c('s1','s2','s3','s4','s5','s6','s7','s8','s9')]

#mychart.Correlation <- edit(chart.Correlation)
chart.Correlation(s_survey, method = 'spearman',histogram=TRUE, pch=20, col = 'grey35')

##section2----import and preprocess raw data####

#import raw data
modeldata <- read.csv('civet.csv')

y <- modeldata[,c('ob1','ob2','ob3', 'ob4', 'ob5', 'ob6','ob7','ob8', 'ob9', 'ob10')]
xraw <- modeldata[,c('s1','s2','s3','s4','s5','s6','s7','s8','s9')]
obsCovs <- list(wind1=modeldata[,c('cs1', 'cs2', 'cs3', 'cs4', 'cs5', 
                                   'cs6', 'cs7', 'cs8', 'cs9', 'cs10')])

#as.factor
xraw$s7<-as.factor(xraw$s7)

#standardize raw data
xraw$s1 <- (xraw$s1-min(xraw$s1))/(max(xraw$s1)-min(xraw$s1))
xraw$s2 <- (xraw$s2-min(xraw$s2))/(max(xraw$s2)-min(xraw$s2))
xraw$s3 <- (xraw$s3-min(xraw$s3))/(max(xraw$s3)-min(xraw$s3))
xraw$s4 <- (xraw$s4-min(xraw$s4))/(max(xraw$s4)-min(xraw$s4))
xraw$s5 <- (xraw$s5-min(xraw$s5))/(max(xraw$s5)-min(xraw$s5))
xraw$s6 <- (xraw$s6-min(xraw$s6))/(max(xraw$s6)-min(xraw$s6))

##section3----select det_var####

#select Detection variable with global site variables
civettestD <- unmarkedFrameOccu(y=y, siteCovs=xraw, obsCovs=obsCovs)
summary(civettestD)

occu1 <- occu(~wind1~., data=civettestD,control=list(maxit=1000))
null <- occu(~1~., data=civettestD, control=list(maxit=1000))
modList <- fitList(D1 = occu1, Null = null)
modSel(modList,nullmod = 'Null')
#none is the best fit

civettestS <- unmarkedFrameOccu(y=y, siteCovs=civettestD@siteCovs, obsCovs=obsCovs)

s1normal<-occu(~1~s1, data=civettestS,control=list(maxit=10000))
s1log<-occu(~1~log(s1+1), data=civettestS,control=list(maxit=10000))
s1sqrt<-occu(~1~sqrt(s1), data=civettestS,control=list(maxit=10000))
modList <- fitList(S1normal = s1normal, S1log=s1log, S1sqrt = s1sqrt)
modSel(modList)

s2normal<-occu(~1~s2, data=civettestS,control=list(maxit=10000))
s2log<-occu(~1~log(s2+1), data=civettestS,control=list(maxit=10000))
s2sqrt<-occu(~1~sqrt(s2), data=civettestS,control=list(maxit=10000))
modList <- fitList(S2normal = s2normal, S2log=s2log, S2sqrt = s2sqrt)
modSel(modList)

s3normal<-occu(~1~s3, data=civettestS,control=list(maxit=10000))
s3log<-occu(~1~log(s3+1), data=civettestS,control=list(maxit=10000))
s3sqrt<-occu(~1~sqrt(s3), data=civettestS,control=list(maxit=10000))
modList <- fitList(S3normal = s3normal, S3log=s3log, S3sqrt = s3sqrt)
modSel(modList)

s4normal<-occu(~1~s4, data=civettestS,control=list(maxit=10000))
s4log<-occu(~1~log(s4+1), data=civettestS,control=list(maxit=10000))
s4sqrt<-occu(~1~sqrt(s4), data=civettestS,control=list(maxit=10000))
modList <- fitList(S4normal = s4normal, S4log=s4log, S4sqrt = s4sqrt)
modSel(modList)

s5normal<-occu(~1~s5, data=civettestS,control=list(maxit=10000))
s5log<-occu(~1~log(s5+1), data=civettestS,control=list(maxit=10000))
s5sqrt<-occu(~1~sqrt(s5), data=civettestS,control=list(maxit=10000))
modList <- fitList(S5normal = s5normal, S5log=s5log, S5sqrt = s5sqrt)
modSel(modList)

s6normal<-occu(~1~s6, data=civettestS,control=list(maxit=10000))
s6log<-occu(~1~log(s6+1), data=civettestS,control=list(maxit=10000))
s6sqrt<-occu(~1~sqrt(s6), data=civettestS,control=list(maxit=10000))
modList <- fitList(S6normal = s6normal, S6log=s6log, S6sqrt = s6sqrt)
modSel(modList)
#sqrt-transform:s1

civettestD@siteCovs$s1 <- sqrt(civettestD@siteCovs$s1)
civettestD@siteCovs$s2 <- sqrt(civettestD@siteCovs$s2)
civettestD@siteCovs$s6 <- sqrt(civettestD@siteCovs$s6)

summary(civettestD)

##section5----build up occupancy model####

#Run the occupancy model with every site variable combinations
y <- civettestD@y
x <- civettestD@siteCovs
obsCovs <- list(wind1=modeldata[,c('cs1', 'cs2', 'cs3', 'cs4', 'cs5', 
                                   'cs6', 'cs7', 'cs8', 'cs9', 'cs10')])


nx <- dim(x)[2]
for (i in 1:nx) {
  index <- combn(1:nx,i)
  
  list.out <- list()
  model.names <- character(dim(index)[2])
  
  for (j in 1:dim(index)[2]) {
    civetUMF.t <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x[,index[,j]]), 
                                    obsCovs=obsCovs)
    list.out[j] <- occu(~1~., civetUMF.t, se=F, control=list(maxit=1000))
    model.names[j] <- paste(paste('s', index[,j], sep=''), collapse=",")
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

#head(aic.table.out)


##section6----model selection####
testdata <- x

pred1 <- occu(~1~s2+s3+s6+s7+s8+s9, data=predUMF)
pred2 <- occu(~1~s2+s3+s5+s6+s7+s8+s9, data=predUMF)
pred3 <- occu(~1~s2+s6+s7+s8+s9, data=predUMF)
pred4 <- occu(~1~s1+s2+s3+s4+s5+s6+s7+s8+s9, data=predUMF)
pred5 <- occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
pred6 <- occu(~1~s2+s5+s6+s7+s8+s9, data=predUMF)

avg <- model.avg(pred1,pred2,pred3,pred4,pred5,pred6)
summary(avg)

#predict the response curve

plot(predUMF,panel=4)

pred <- occu(~1~s2+s3+s6+s7+s8+s9, data=predUMF)
pred
plot(pred)


#expect site psi&p

site_psi <- round(predict(pred, type="state"), 2)  # Expected psi at each site

site_p <- round(predict(pred, type="det"), 2)  # Expected p at each site

##ROC curve####
testdata$label <- data_civet$label

summary(testdata)
head(testdata)


tests7 <- model.matrix(~s7,data=testdata)
testdata$s70.5 <- tests7[,2]
testdata$s71 <- tests7[,3]

testdata$pred<-( 38.36 + -5.28*testdata$s1 + 39.58*testdata$s2 + -27.91*testdata$s3 +
                   -2.54*testdata$s4 + -6.26*testdata$s5 + 44.67*testdata$s6 + -20.82*testdata$s70.5 + 
                   -47.01*testdata$s71 + -40.95*testdata$s8 + -44.69*testdata$s9 )

testdata$Psi<-invlogit(testdata$pred)

round(testdata$Psi,4)

perf <- performance(prediction(testdata$Psi, testdata$label), 'tpr', 'fpr')
pred <- prediction(testdata$Psi, testdata$label)
plot(perf)

perf_cost = performance(pred, "cost")
perf_err = performance(pred, "err")
perf_tpr = performance(pred, "tpr")
perf_sn_sp = performance(pred, "sens", "spec")

plot(perf_cost)

roc = performance(pred,"tpr","fpr")
plot(roc, colorize = T, lwd = 2,xlab='False Positive Rate',ylab='True Positive Rate')
abline(a = 0, b = 1)
text(0.8,0.1,'AUC = 0.902',cex = 1.5)

auc = performance(pred, measure = "auc")
print(auc@y.values)

#
cost.perf = performance(pred, "cost")
pred@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]


opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc, pred))


##species-environment relationship#####

predUMF <- unmarkedFrameOccu(y=y, siteCovs=x, obsCovs=obsCovs)

par(mfrow=c(3,3))
#Hold everything constant except for RAI_prey
newDats1 <- data.frame(s1 = seq(min(x$s1), max(x$s1), length=100), 
                       s2 = mean(x$s2), s3 = mean(x$s3), s4 = mean(x$s4), 
                       s5 = mean(x$s5), s6 = mean(x$s6), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats1$s7 <- as.factor(newDats1$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds1 <- as.data.frame(predict(pred, newdata=newDats1, type="state"))
with(Epreds1, {
  s1 <- newDats1$s1^2*(max(modeldata$s1)-min(modeldata$s1))+min(modeldata$s1)
  plot(s1, Predicted, xlab="Slope", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(4,55), ylim=c(0,1),type="l", lwd=2)
  lines(s1, Predicted-SE, lty=3)
  lines(s1, Predicted+SE, lty=3)
})

newDats2 <- data.frame(s2 = seq(min(x$s2), max(x$s2), length=100), 
                       s1 = mean(x$s1), s3 = mean(x$s3), s4 = mean(x$s4), 
                       s5 = mean(x$s5), s6 = mean(x$s6), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats2$s7 <- as.factor(newDats2$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds2 <- as.data.frame(predict(pred, newdata=newDats2, type="state"))
with(Epreds2, {
  s2 <- newDats2$s2^2*(max(modeldata$s2)-min(modeldata$s2))+min(modeldata$s2)
  plot(s2, Predicted, xlab="Prey abundance", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(0,18),ylim=c(0,1),type="l", lwd=2)
  lines(s2, Predicted-SE, lty=3)
  lines(s2, Predicted+SE, lty=3)
})

newDats3 <- data.frame(s3 = seq(min(x$s3), max(x$s3), length=100), 
                       s1 = mean(x$s1), s2 = mean(x$s2), s4 = mean(x$s4), 
                       s5 = mean(x$s5), s6 = mean(x$s6), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats3$s7 <- as.factor(newDats3$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds3 <- as.data.frame(predict(pred, newdata=newDats3, type="state"))
with(Epreds3, {
  s3 <- newDats3$s3*(max(modeldata$s3)-min(modeldata$s3))+min(modeldata$s3)
  plot(s3, Predicted, xlab="Distance to river", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(0,700),ylim=c(0,1),type="l", lwd=2)
  lines(s3, Predicted-SE, lty=3)
  lines(s3, Predicted+SE, lty=3)
})

newDats4 <- data.frame(s4 = seq(min(x$s4), max(x$s4), length=100), 
                       s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3),
                       s5 = mean(x$s5), s6 = mean(x$s6), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats4$s7 <- as.factor(newDats4$s7)
pred<-occu(~1~s1+s2+s3+s4+s5+s6+s7+s8+s9, data=predUMF)
Epreds4 <- as.data.frame(predict(pred, newdata=newDats4, type="state"))
with(Epreds4, {
  s4 <- newDats4$s4*(max(modeldata$s4)-min(modeldata$s4))+min(modeldata$s4)
  plot(s4, Predicted, xlab="Distance to road", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(0,3000),ylim=c(0,1),type="l", lwd=2)
  lines(s4, Predicted-SE, lty=3)
  lines(s4, Predicted+SE, lty=3)
})

newDats5 <- data.frame(s5 = seq(min(x$s5), max(x$s5), length=100), 
                       s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3),
                       s4 = mean(x$s4), s6 = mean(x$s6), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats5$s7 <- as.factor(newDats5$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds5 <- as.data.frame(predict(pred, newdata=newDats5, type="state"))
with(Epreds5, {
  s5 <- newDats5$s5*(max(modeldata$s5)-min(modeldata$s5))+min(modeldata$s5)
  plot(s5, Predicted, xlab="Dog abundance", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(0,1.2),ylim=c(0,1),type="l", lwd=2)
  lines(s5, Predicted-SE, lty=3)
  lines(s5, Predicted+SE, lty=3)
})

newDats6 <- data.frame(s6 = seq(min(x$s6), max(x$s6), length=100), 
                       s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3), s4 = mean(x$s4), 
                       s5 = mean(x$s5), s7 = 0.5, s8 = mean(x$s8), s9 = 0)
newDats6$s7 <- as.factor(newDats6$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds6 <- as.data.frame(predict(pred, newdata=newDats6, type="state"))
with(Epreds6, {
  s6 <- newDats6$s6^2*(max(modeldata$s6)-min(modeldata$s6))+min(modeldata$s6)
  plot(s6, Predicted, xlab="Human abundance", ylab="Occupancy probability", 
       pch=16, cex.lab=1, xlim=c(0,3), ylim=c(0,1),type="l", lwd=2)
  lines(s6, Predicted-SE, lty=3)
  lines(s6, Predicted+SE, lty=3)
})

newDats7 <- data.frame(s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3), 
                       s4 = mean(x$s4), s5 = mean(x$s5), s6 = mean(x$s6), 
                       s7 = c(0,0.5,1), s8 = mean(x$s8), s9 = 0)
newDats7$s7 <- as.factor(newDats7$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds7 <- as.data.frame(predict(pred, newdata=newDats7, type="state"))
with(Epreds7, {
  bp <- barplot(Predicted, ylim=c(0, 1.1), names.arg=newDats7$s7, 
                xlab="Forest structure",ylab="Occupancy probability", cex.lab=1, space=0.5,
                axis.lty=1)
  box()
  arrows(bp, Predicted, bp, Predicted+SE, code=2, angle=90, length=0.05)
})

newDats8 <- data.frame(s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3), 
                       s4 = mean(x$s4), s5 = mean(x$s5), s6 = mean(x$s6), 
                       s7 = 0.5, s8 = c(0,0.5,1), s9 = 0)
newDats8$s7 <- as.factor(newDats8$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds8 <- as.data.frame(predict(pred, newdata=newDats8, type="state"))
with(Epreds8, {
  bp <- barplot(Predicted, ylim=c(0, 1.1), names.arg=newDats8$s8, 
                xlab="Height of tree",ylab="Occupancy probability", cex.lab=1, space=0.5,
                axis.lty=1)
  box()
  arrows(bp, Predicted, bp, Predicted+SE, code=2, angle=90, length=0.05)
})

newDats9 <- data.frame(s1 = mean(x$s1), s2 = mean(x$s2), s3 = mean(x$s3),
                       s4 = mean(x$s4), s5 = mean(x$s5), s6 = mean(x$s6), 
                       s7 = 0.5, s8 = mean(x$s8), s9 = c(0,1))
newDats9$s7 <- as.factor(newDats9$s7)
pred<-occu(~1~s1+s2+s3+s5+s6+s7+s8+s9, data=predUMF)
Epreds9 <- as.data.frame(predict(pred, newdata=newDats9, type="state"))
with(Epreds9, {
  bp <- barplot(Predicted, ylim=c(0, 1.1), names.arg=newDats9$s9, 
                xlab="Bamboo presence",ylab="Occupancy probability", cex.lab=1, space=0.5,
                axis.lty=1)
  box()
  arrows(bp, Predicted, bp, Predicted+SE, code=2, angle=90, length=0.05)
})

