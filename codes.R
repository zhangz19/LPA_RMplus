
rm(list=ls())
library('MplusAutomation')
library(MASS)

M <- 5000   #total replications 
# simulation method: sampling from reduced population 334 cases with 
# full observation to match the missing pattern

#--------------------- data processing
dat0 <- read.table('train8.txt',sep='\t',na='.')
tmp <- dat0[,-c(1,ncol(dat0))]   # exclude ID info (1st col) and gender info (last col)
mis <- numeric(); for(i in 1:nrow(tmp)) mis <- c(mis, sum(is.na(tmp[i,])))
print(table(mis))  # summary missing cases
dat <- na.omit(dat0)  # for training data set, keep cases with full observations
n <- nrow(dat); p <- ncol(dat)
dat[,p] <- factor(dat[,p], levels=c('M','F'))

#--------------------- simulation of missing counts according to X
dat1 <- read.table('test8.txt',sep='\t',na='.')
X <- dat1[,2:ncol(dat1)]  # with no ID
n1 <- nrow(X); p1 <- ncol(X)
mis1 <- numeric(); for(i in 1:n1) mis1 <- c(mis1, sum(is.na(X[i,])))
print(table(mis1))  # summary missing cases

#--------------------- define util functions
foo <- function(x) sum(is.na(x))
# missin values simulation function: exactly the same missing pattern
simuMis <- function(X0, X, n1){  
  tmp <- X[sample(1:nrow(X)),]  # random permutation
  for(i in 1:n1) X0[i,which(is.na(tmp[i,]))] <- NA
  return(X0)
}

#--------------------- run simulations

start.time <- Sys.time(); 
MissRate <- numeric(0)
mat1 <- mat2 <- matrix(numeric(M*4), ncol=4)
err1 <- err2 <- numeric(M)

set.seed(200)

for(repl in 1:M){
  cat(repl,' ')
  X0 <- dat[sample(1:nrow(dat),size=nrow(X)),2:ncol(dat)]
  sex <- X0[,ncol(X0)]; X0 <- X0[,-ncol(X0)]

  #--------------------- for complete observation case
  write.table(X0, file='dat.txt',quote=F,na='.',sep='\t',row.names=F,col.names=F)
  runModels('./')  # run model, call "simu.inp" under the same directory
  out <- read.table('datOutput.txt', sep='') #collect output results from Mplus running
  if(nrow(out)!=nrow(X0)) stop('nrow(out)!=nrow(X0)')
  prediction <- out[,ncol(out)]
  # now need to identify the group
  m1 <- colMeans(X0[g1 <- which(prediction==1),],na.rm=T)
  m2 <- colMeans(X0[g2 <- which(prediction==2),],na.rm=T)
  if(mean(m1>m2, na.rm=T)>0.5){ #the empirical rule: if half of the meansures dominate over group2
     prediction[g1] <- "M"; prediction[g2] <- "F" }
  else{prediction[g1] <- "F"; prediction[g2] <- "M"}
  prediction <- factor(prediction, levels=c("M","F"))
  tab <- table(sex, prediction)
  N <- sum(tab); nc <- sum(diag(tab)); K <- 2
  mat1[repl,] <- as.vector(tab)
  err1[repl] <- (N-nc)/N
  
  #--------------------- for simulated missing observation case
  X0 <- simuMis(X0, X, n1)
  write.table(X0, file='dat.txt',quote=F,na='.',sep='\t',row.names=F,col.names=F)
  runModels('./')  # run model
  out <- read.table('datOutput.txt', sep='')
  if(nrow(out)!=nrow(X0)) stop('nrow(out)!=nrow(X0)')
  prediction <- out[,ncol(out)]
  # now need to identify the group
  m1 <- colMeans(X0[g1 <- which(prediction==1),],na.rm=T)
  m2 <- colMeans(X0[g2 <- which(prediction==2),],na.rm=T)
  if(mean(m1>m2, na.rm=T)>0.5){ #the empirical rule: if half of the meansures dominate over group2
     prediction[g1] <- "M"; prediction[g2] <- "F" }
  else{prediction[g1] <- "F"; prediction[g2] <- "M"}
  prediction <- factor(prediction, levels=c("M","F"))

  tab <- table(sex, prediction)
  N <- sum(tab); nc <- sum(diag(tab)); K <- 2
  mat2[repl,] <- as.vector(tab)
  err2[repl] <- (N-nc)/N
  
  MissRate <- rbind(MissRate, apply(X0,2,foo)/nrow(X0))
}

# report running time
print(times <- Sys.time() - start.time )

#--------------------- save the results
 save(times,M,mat1,err1,mat2,err2,file='out_120203_simu9_m2')

# comparison of full vs missing
q1 <- as.character(round(quantile(1-err1, c(.025, .5, .975)),3))
q2 <- as.character(round(quantile(1-err2, c(.025, .5, .975)),3))
q12 <- as.character(round(quantile(err2-err1, c(.025, .5, .975)),3))

par(mfrow=c(1,3))
truehist(1-err1, prob=F, xlab='accuracy (full)',main='histogram', sub='5000 replicates'); legend('topleft',c('2.5%,50%,97.5%:',q1), inset=-.01,bty='n')
truehist(1-err2, prob=F, xlab='accuracy (missing)',main='histogram', sub='5000 replicates'); legend('topleft',c('2.5%,50%,97.5%:',q2), inset=-.01,bty='n')
truehist(err2-err1, prob=F, xlab='difference of accuracy (full - missing)',main='histogram', sub='5000 replicates'); legend('topleft',c('2.5%,50%,97.5%:',q12), inset=-.01,bty='n')

## not run

