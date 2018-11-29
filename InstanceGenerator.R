#========================================================================================
# File:     Random sample generation for Correlation induction with Java
#           & benchmark using IC from mc2d package
#           Multivariate case
# Created:  01/08/2016
# Updated:  29/11/2018
# Author:   Oscar Guaje (oo.guaje10@uniandes.edu.co)
#           http://copa.uniandes.edu.co
#=======================================================================================


library(mc2d)
library(dplyr)
library(MASS)
library(DiscreteWeibull)

#setwd("C:/Users/oo.guaje10/Dropbox/01 - Tesis/Multivariate Java/data")

rm(list = ls())

generador5 <- function(obs, discretas, continuas, semilla, target5){
  
  for(i in 1:5){
    
    set.seed(12+(semilla+i))
    D1 <- sample(0:99, obs, replace = T)
    
    set.seed(12*(semilla+i))
    D2 <- rpois(n = obs, lambda = 5)
    
    set.seed(1000+(semilla+i))
    D3 <- rbinom(n = obs, size = 100, prob = 0.3)
    
    set.seed(1000*(semilla+i))
    D4 <- rhyper(nn = obs, m = 20, n = 30, k = 10)
    
    set.seed(8000+(semilla+i))
    D5 <- rdweibull(n = obs, q = 0.4, beta = 0.7)
    
    set.seed(8000*(semilla+i))
    C1 <- rnorm(obs, mean = 5, sd = 3)
    
    set.seed(52300+(semilla+i))
    C2 <- rexp(obs, rate = 4)
    
    set.seed(52300*(semilla+i))
    C3 <- runif(obs)
    
    set.seed(510+(semilla+i))
    C4 <- rlnorm(n = obs, meanlog = log(5.14^2/sqrt(2.76^2 + 5.14^2)), sdlog = sqrt(log(((2.76^2)/(5.14^2))+1)))
    
    set.seed(510*(semilla+i))
    C5 <- rgamma(n = obs, shape = 8, scale = 5)
    
    ddiscretas <- cbind(D1,D2,D3,D4,D5)
    dcontinuas <- cbind(C1,C2,C3,C4,C5)
    
    if(continuas == 5){
      datos <- dcontinuas
    }else if(discretas == 5){
      datos <- ddiscretas
    }else{
      datos <- as.matrix(cbind(ddiscretas[,1:discretas], dcontinuas[,1:continuas])) 
    }
    
    IC <- cornode(datos, target = target5)
    IC_out=cor(IC,method="pearson")
    obj <- abs(target5-IC_out)
    obj[upper.tri(obj, diag = FALSE)]
    sum(obj[upper.tri(obj, diag = FALSE)])
    
    
    filename <- paste('5vars', obs, sep='_')
    filename <- paste(filename,discretas,sep='_')
    filename <- paste(filename,continuas,sep='_')
    filename <- paste(filename,i,sep='_')
    dir.create(filename)
    filename <- paste(filename,'.txt',sep='')
    
    paste('K:',ncol(IC),sep='') %>% write(filename,append=T)
    paste('N:',nrow(IC),sep='') %>% write(filename,append=T)
    paste('Gap:',sum(obj[upper.tri(obj, diag = FALSE)]),sep='') %>% write(filename,append=T)
    write('Rho: ', filename, append=T)
    write.table(target5, file = filename,append = T, quote=F,col.names=F,row.names=F)
    write('X:', filename, append=T)
    write.table(IC,filename,sep=',',quote=F,col.names=F,row.names=F,append=T)  
  }
}

generador10 <- function(obs, discretas, continuas, semilla, target10){
  
  for(i in 1:5){
    
    set.seed(12+(semilla+i))
    D1 <- sample(0:99, obs, replace = T)
    
    set.seed(12*(semilla+i))
    D2 <- rpois(n = obs, lambda = 5)
    
    set.seed(1000+(semilla+i))
    D3 <- rbinom(n = obs, size = 100, prob = 0.3)
    
    set.seed(1000*(semilla+i))
    D4 <- rhyper(nn = obs, m = 20, n = 30, k = 10)
    
    set.seed(8000+(semilla+i))
    D5 <- rdweibull(n = obs, q = 0.4, beta = 0.7)
    
    
    set.seed(60+(semilla+i))
    D6 <- sample(0:99, obs, replace = T)
    
    set.seed(60*(semilla+i))
    D7 <- rpois(n = obs, lambda = 8)
    
    set.seed(2011+(semilla+i))
    D8 <- rbinom(n = obs, size = 100, prob = 0.3)
    
    set.seed(2011*(semilla+i))
    D9 <- rhyper(nn = obs, m = 70, n = 30, k = 30)
    
    set.seed(7369+(semilla+i))
    D10 <- rdweibull(n = obs, q = 0.5, beta = 0.9)
    
    
    set.seed(8000*(semilla+i))
    C1 <- rnorm(obs, mean = 5, sd = 3)
    
    set.seed(52300+(semilla+i))
    C2 <- rexp(obs, rate = 4)
    
    set.seed(52300*(semilla+i))
    C3 <- runif(obs)
    
    set.seed(510+(semilla+i))
    C4 <- rlnorm(n = obs, meanlog = log(5.14^2/sqrt(2.76^2 + 5.14^2)), sdlog = sqrt(log(((2.76^2)/(5.14^2))+1)))
    
    set.seed(510*(semilla+i))
    C5 <- rgamma(n = obs, shape = 8, scale = 5)
    
    
    set.seed(7369*(semilla+i))
    C6 <- rnorm(obs, mean = 8, sd = 16)
    
    set.seed(100000+(semilla+i))
    C7 <- rexp(obs, rate = 10)
    
    set.seed(100000*(semilla+i))
    C8 <- runif(obs)
    
    set.seed(666+(semilla+i))
    C9 <- rlnorm(n = obs, meanlog = log(9.74^2/sqrt(6.58^2 + 9.74^2)), sdlog = sqrt(log(((6.58^2)/(9.74^2))+1)))
    
    set.seed(666*(semilla+i))
    C10 <- rgamma(n = obs, shape = 10, scale = 7)
    
    
    ddiscretas <- cbind(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10)
    dcontinuas <- cbind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10)
    
    if(continuas == 10){
      datos <- dcontinuas
    }else if(discretas == 10){
      datos <- ddiscretas
    }else{
      datos <- as.matrix(cbind(ddiscretas[,1:discretas], dcontinuas[,1:continuas])) 
    }
    
    IC <- cornode(datos, target = target10)
    IC_out=cor(IC,method="pearson")
    obj <- abs(target10-IC_out)
    obj[upper.tri(obj, diag = FALSE)]
    sum(obj[upper.tri(obj, diag = FALSE)])
    
    
    filename <- paste('10vars', obs, sep='_')
    filename <- paste(filename,discretas,sep='_')
    filename <- paste(filename,continuas,sep='_')
    filename <- paste(filename,i,sep='_')
    dir.create(filename)
    filename <- paste(filename,'.txt',sep='')
    
    paste('K:',ncol(IC),sep='') %>% write(filename,append=T)
    paste('N:',nrow(IC),sep='') %>% write(filename,append=T)
    paste('Gap:',sum(obj[upper.tri(obj, diag = FALSE)]),sep='') %>% write(filename,append=T)
    write('Rho: ', filename, append=T)
    write.table(target10, file = filename,append = T, quote=F,col.names=F,row.names=F)
    write('X:', filename, append=T)
    write.table(IC,filename,sep=',',quote=F,col.names=F,row.names=F,append=T)  
  }
}

# Para 5 variables:

target5 <- matrix(c(1,0.8,0.45,0.8,1,0.45,0.45,0.45,1),nrow=3)
target5 <- cbind(target5,c(0.2,0.5,0.1))
target5 <- rbind(target5,c(0.2,0.5,0.1,1))

target5 <- cbind(target5,c(-0.3,-0.1,-0.3,-0.5))
target5 <- rbind(target5,c(-0.3,-0.1,-0.3,-0.5,1))


# Para 10 variables

target10 <- matrix(c(1, 	0.8, 	0.45, 	0.2, 	-0.3, 	0.7, 	0.75, 	0.5, 	-0.3, 	-0.1, 
                     0.8, 	1, 	0.45, 	0.5, 	-0.1, 	0.55, 	0.4, 	0.6, 	-0.4, 	-0.2, 
                     0.45, 	0.45, 	1, 	0.1, 	-0.3, 	0.6, 	0.3, 	0.4, 	0, 	-0.3, 
                     0.2, 	0.5, 	0.1, 	1, 	-0.5, 	0.35, 	0.1, 	0.1, 	-0.5, 	-0.45, 
                     -0.3, 	-0.1, 	-0.3, 	-0.5, 	1, 	-0.4, 	-0.3, 	-0.1, 	0.2, 	0.2, 
                     0.7, 	0.55, 	0.6, 	0.35, 	-0.4, 	1, 	0.8, 	0.45, 	-0.2, 	-0.45, 
                     0.75, 	0.4, 	0.3, 	0.1, 	-0.3, 	0.8, 	1, 	0.45, 	-0.1, 	-0.3, 
                     0.5, 	0.6, 	0.4, 	0.1, 	-0.1, 	0.45, 	0.45, 	1, 	0.1, 	-0.2, 
                     -0.3, 	-0.4, 	0, 	-0.5, 	0.2, 	-0.2, 	-0.1, 	0.1, 	1, 	0.4, 
                     -0.1, 	-0.2, 	-0.3, 	-0.45, 	0.2, 	-0.45, 	-0.3, 	-0.2, 	0.4, 	1
),nrow=10)


observaciones2 <- c(500, 1000, 2000, 3000)


for(observaciones in observaciones2){
  semilla = 28
  generador5(observaciones, 5, 0 , semilla, target5)
  generador5(observaciones, 4, 1 , semilla, target5)
  generador5(observaciones, 3, 2 , semilla, target5)
  generador5(observaciones, 2, 3 , semilla, target5)
  generador5(observaciones, 1, 4 , semilla, target5)
  generador5(observaciones, 0, 5 , semilla, target5)
  
  generador10(observaciones, 10, 0 , semilla, target10)
  generador10(observaciones, 8, 2 , semilla, target10)
  generador10(observaciones, 6, 4 , semilla, target10)
  generador10(observaciones, 4, 6 , semilla, target10)
  generador10(observaciones, 2, 8 , semilla, target10)
  generador10(observaciones, 0, 10 , semilla, target10)
}


