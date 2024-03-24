library(dplyr)
library(rstan)

rm(list = ls())
#### observation data ########################################
setwd('test/Data/All_widedata')
files <- list.files(pattern = '.csv')
files <- files[-7]

sp_data <- list()
for (i in 1:length(files)) {
  sp <- read.csv(files[i])
  sp <- as.matrix(sp[-1])
  sp_data[[i]] <- sp
}

segment <- (ncol(cam_hour)-1)%/%5
det <- list()
df <- data.frame()
colname <- vector()

for (i in 1:length(sp_data)) {
  print(i)
  for (k in 1:segment) {
    for (j in 1:nrow(cam_hour)) {
      df[j,k] <- sum(sp_data[[i]][j,(5*k-4):(5*k)])
      if (!is.na(df[j,k])) {
        if (df[j,k] >= 1) {
          df[j,k] <- 1
        }
      }
    }
    colname[k] <- paste('seg',k,sep = '')
  }
  colnames(df) <- colname
  det[[i]] <- df
}

Y_ob <- array(unlist(det), dim = c(j,k,i))

#### data augmentation ########################################
un_ob = 20
Y_zero <- Y_ob[,,1]
Y_zero[which(Y_zero > 0)] <- 0
Y_aug <- array(0, dim = c(j,k,i+un_ob))
Y_aug[,,1:i] <- Y_ob
Y_aug[,,(i+1):(i+un_ob)] <- rep(Y_zero, un_ob)

site <- cam_hour$station
seg <- colname
sp <- c(1:(i+un_ob))
dimnames(Y_aug) <- list(site, seg, sp)

#### detection covariates ########################################
cam_hour <- read.csv('datawide_camhour.csv')
setwd('C:/Users/William/Desktop/Rproject/rstan/test')
env <- read.csv('Data/MSOM_sitecov.csv')
park <- read.csv('Data/park_level_cov.csv')

cam_hour5 <- Y_zero
for (k in 1:segment) {
  cam_hour5[,k] <- rowSums(cam_hour[,(5*k-3):(5*k+1)])/5
}
hour <- (cam_hour5 - mean(cam_hour5, na.rm = T))/sd(cam_hour5, na.rm = T)

angle <- Y_zero
for (j in 1:nrow(hour)) {
  if (env$Cam_angle[j] == 1) {
    angle[j,] <- angle[j,] + 1
  }
}

#### site level covariates ########################################
env_std <- read.csv('Data/RN_model_final/STDED.DEC.sitecov.csv', header = T, row.names = 1)

#### PA level covariates ########################################
pasize <- env_std %>% select(park.ind,pasize) %>% distinct() %>% arrange(park.ind)
punish_vill <- env_std %>% select(park.ind,punishment) %>% distinct() %>% arrange(park.ind)
reach_vill <- env_std %>% select(park.ind,outreach) %>% distinct() %>% arrange(park.ind)
punish_park <- env_std %>%  select(park.ind,park_punishment) %>% distinct() %>% arrange(park.ind)
reach_park <- env_std %>%  select(park.ind,park_outreach) %>% distinct() %>% arrange(park.ind)


#### data preparation ###########################################################
new_hour <- hour
new_angle <- angle
new_Y <- Y_aug
new_hour[is.na(new_hour)] <- 9999
new_angle[is.na(new_angle)] <- 9999
new_Y[is.na(new_Y)] <- 9999

occu_data1 <- list(n1=dim(Y_ob)[3], n0=un_ob, n = dim(Y_aug)[3], J=dim(Y_aug)[1], K=rowSums(!is.na(Y_zero)), 
                 park=env_std$park.ind, PA=max(env_std$park.ind), 
                 elev=env_std$elevation, pop=env_std$population, dist=env_std$distance,
                 pasize=pasize[,2], punish=punish_vill[,2], reach=reach_vill[,2],
                 hour=new_hour, angle=new_angle, Y=new_Y)


occu_data2 <- list(n1=dim(Y_ob)[3], n0=un_ob, n = dim(Y_aug)[3], J=dim(Y_aug)[1], K=rowSums(!is.na(Y_zero)),
                  park=env_std$park.ind, PA=max(env_std$park.ind),  
                  elev=env_std$elevation, pop=env_std$population, dist=env_std$distance,
                  pasize=pasize[,2], punish=punish_park[,2], reach=reach_park[,2],
                  hour=new_hour, angle=new_angle, Y=new_Y)



#####################################################################
model <- stan_model(file = 'test/MSOM_code.stan')

fit <- sampling(model, data = occu_data1, chains = 3, iter = 10000, warmup = 5000)

