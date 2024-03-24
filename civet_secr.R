library(secr)
library(rgdal)

#rm(list = ls())

# import data
## traps

CVtraps1 <- read.traps(data = read.csv('CVtraps1.csv', row.names = 1),
                       detector = 'count')
CVtraps2 <- read.traps(data = read.csv('CVtraps2.csv', row.names = 1),
                       detector = 'count')
CVtraps3 <- read.traps(data = read.csv('CVtraps3.csv', row.names = 1),
                       detector = 'count')
CVtraps4 <- read.traps(data = read.csv('CVtraps4.csv', row.names = 1),
                       detector = 'count')


## capture history
CVcaptures2107 <- read.csv('CVcaps1_l3.csv')
CVdata1 <- make.capthist(captures = CVcaptures2107, traps = CVtraps1, fmt = 'trapID') 

CVcaptures2201 <- read.csv('CVcaps2_l3.csv')
CVdata2 <- make.capthist(captures = CVcaptures2201, traps = CVtraps2, fmt = 'trapID') 

CVcaptures2207 <- read.csv('CVcaps3_l3.csv')
CVdata3 <- make.capthist(captures = CVcaptures2207, traps = CVtraps3, fmt = 'trapID') 

CVcaptures2301 <- read.csv('CVcaps4_l3.csv')
CVdata4 <- make.capthist(captures = CVcaptures2301, traps = CVtraps4, fmt = 'trapID') 


## detection of 4 seasons
#CVtraps <- read.traps(data = read.csv('CVsites_new.csv', row.names = 1),
#                      detector = 'count')
#CVcaptures <- read.csv('CVcaptures4s_l2.csv')
#CVdata <- make.capthist(captures = CVcaptures, traps = CVtraps, fmt = 'trapID') 
#summary(CVdata)

# closure.test(CVdata)
# p > 0.05, closure can be assumed
# closure.test(CVdata1)
# closure.test(CVdata2)
# closure.test(CVdata3)
# closure.test(CVdata4)


# initial model
#CVdist <- RPSV(CVdata, CC = TRUE)
CVdist1 <- RPSV(CVdata1, CC = TRUE)
CVdist2 <- RPSV(CVdata2, CC = TRUE)
CVdist3 <- RPSV(CVdata3, CC = TRUE)
CVdist4 <- RPSV(CVdata4, CC = TRUE)
CV_RPSV <- c(CVdist1, CVdist2, CVdist3, CVdist4)

CVdist01 <- MMDM(CVdata1)
CVdist02 <- MMDM(CVdata2)
CVdist03 <- MMDM(CVdata3)
CVdist04 <- MMDM(CVdata4)
CV_MMDM <- c(CVdist01, CVdist02, CVdist03, CVdist04)

CVmask0 <- make.mask(CVtraps1, buffer = 4*max(CV_RPSV), type ='trapbuffer')

#---
CVcounts_hn01 <- secr.fit(CVdata1, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)
coefficients(CVcounts_hn01)
CVpred01 <- predict(CVcounts_hn01)
CVnumb01 <- region.N(CVcounts_hn01)
#---
CVcounts_hn02 <- secr.fit(CVdata2, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)
coefficients(CVcounts_hn02)
CVpred02 <- predict(CVcounts_hn02)
CVnumb02 <- region.N(CVcounts_hn02)
#---
CVcounts_hn03 <- secr.fit(CVdata3, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)
coefficients(CVcounts_hn03)
CVpred03 <- predict(CVcounts_hn03)
CVnumb03 <- region.N(CVcounts_hn03)
#---
CVcounts_hn04 <- secr.fit(CVdata4, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)
coefficients(CVcounts_hn04)
CVpred04 <- predict(CVcounts_hn04)
CVnumb04 <- region.N(CVcounts_hn04)


# corrected null model
#CVbuffer <- suggest.buffer(CVdata, detectfn = 'HN',
#                            detectpar = list(
#                              g0 = CVpred0$estimate[2],
#                              sigma = CVpred0$estimate[3]))
#
#CVmask <-  make.mask(traps(CVdata), buffer = CVbuffer, type ='trapbuffer')
#CVcounts_hn <- secr.fit(CVdata, detectfn = 'HN', mask = CVmask)
#
#coefficients(CVcounts_hn)
#CVpred <- predict(CVcounts_hn)
#CVnumb <- region.N(CVcounts_hn)


# import habitats
hbt1 <- readOGR('.', 'CV_mask') # buffer = 0.5*MMDM
CVmask1 <- make.mask(traps(CVdata1), type='polybuffer', 
                     poly = hbt1, poly.habitat = T)
plot(CVmask1)



# for capthist 202107
CVcounts_hn1 <- secr.fit(CVdata1, detectfn = 'HN', mask = CVmask1, start = autoini(CVdata1,CVmask1))
coefficients(CVcounts_hn1)
CVpred1 <- predict(CVcounts_hn1)
CVnumb1 <- region.N(CVcounts_hn1)

CVsur1 <- fx.total(CVcounts_hn1)
plot(CVsur1, covariate = 'D.sum', scale = 1000)


# for capthist 202201
CVcounts_hn2 <- secr.fit(CVdata2, detectfn = 'HN', mask = CVmask1, start = autoini(CVdata2,CVmask1))
coefficients(CVcounts_hn2)
CVpred2 <- predict(CVcounts_hn2)
CVnumb2 <- region.N(CVcounts_hn2)

CVsur2 <- fx.total(CVcounts_hn2)
plot(CVsur2, covariate = 'D.sum', scale = 1000)


# for capthist 202207
CVcounts_hn3 <- secr.fit(CVdata3, detectfn = 'HN', mask = CVmask1, start = autoini(CVdata3,CVmask1))
coefficients(CVcounts_hn3)
CVpred3 <- predict(CVcounts_hn3)
CVnumb3 <- region.N(CVcounts_hn3)

CVsur3 <- fx.total(CVcounts_hn3)
plot(CVsur3, covariate = 'D.sum', scale = 1000)


# for capthist 202301
CVcounts_hn4 <- secr.fit(CVdata4, detectfn = 'HN', mask = CVmask1, start = autoini(CVdata4,CVmask1))
coefficients(CVcounts_hn4)
CVpred4 <- predict(CVcounts_hn4)
CVnumb4 <- region.N(CVcounts_hn4)

CVsur4 <- fx.total(CVcounts_hn4)
plot(CVsur4, covariate = 'D.sum', scale = 1000)

##
coef(CVcounts_hn1)
coef(CVcounts_hn2)
coef(CVcounts_hn3)
coef(CVcounts_hn4)

## D-beta
m1 <- coef(CVcounts_hn1)[1,1]
s1 <- coef(CVcounts_hn1)[1,2]
m2 <- coef(CVcounts_hn2)[1,1]
s2 <- coef(CVcounts_hn2)[1,2]
m3 <- coef(CVcounts_hn3)[1,1]
s3 <- coef(CVcounts_hn3)[1,2]
m4 <- coef(CVcounts_hn4)[1,1]
s4 <- coef(CVcounts_hn4)[1,2]

m_D <- c(m1,m2,m3,m4)
s_D <- c(s1,s2,s3,s4)

## g0-beta
m1 <- coef(CVcounts_hn1)[2,1]
s1 <- coef(CVcounts_hn1)[2,2]
m2 <- coef(CVcounts_hn2)[2,1]
s2 <- coef(CVcounts_hn2)[2,2]
m3 <- coef(CVcounts_hn3)[2,1]
s3 <- coef(CVcounts_hn3)[2,2]
m4 <- coef(CVcounts_hn4)[2,1]
s4 <- coef(CVcounts_hn4)[2,2]

m_g <- c(m1,m2,m3,m4)
s_g <- c(s1,s2,s3,s4)

## sigma-beta
m1 <- coef(CVcounts_hn1)[3,1]
s1 <- coef(CVcounts_hn1)[3,2]
m2 <- coef(CVcounts_hn2)[3,1]
s2 <- coef(CVcounts_hn2)[3,2]
m3 <- coef(CVcounts_hn3)[3,1]
s3 <- coef(CVcounts_hn3)[3,2]
m4 <- coef(CVcounts_hn4)[3,1]
s4 <- coef(CVcounts_hn4)[3,2]

m_s <- c(m1,m2,m3,m4)
s_s <- c(s1,s2,s3,s4)

##
m <- mean(m_s)
s <- sqrt(mean(s_s^2))

m <- mean(m_D)
s <- sqrt(mean(s_D^2))

# for log-link
m_org = exp(m)
#s_org = sqrt(exp(s^2)-1)*exp(m)
lcl_org = exp(m-1.96*s)
ucl_org = exp(m+1.96*s)

##
m <- mean(m_g)
s <- sqrt(mean(s_g^2))

# for logit-link
m_org = plogis(m)
#s_org = ?
lcl_org = plogis(m-1.96*s)
ucl_org = plogis(m+1.96*s)


# Cohen's d
d12 = (m1-m2)/sqrt((s1^2+s2^2)/2) # 0.53
d13 = (m1-m3)/sqrt((s1^2+s3^2)/2) # -0.46
d14 = (m1-m4)/sqrt((s1^2+s4^2)/2) # 0.89
d23 = (m2-m3)/sqrt((s2^2+s3^2)/2) # -0.96
d24 = (m2-m4)/sqrt((s2^2+s4^2)/2) # 0.38
d34 = (m3-m4)/sqrt((s3^2+s4^2)/2) # 1.30

