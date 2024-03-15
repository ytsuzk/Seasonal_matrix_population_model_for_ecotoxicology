##############################################################################
# R script for the Three-stage model 1
# 
# Author: Yoichi Tsuzuki
# Affiliation: National Institute for Environmental Studies, Japan (NIES)
# Contact: tsuzuki.yoichi@nies.go.jp, 
#          yoichi.tsuzuki.95@gmail.com
#
##############################################################################

rm(list=ls())

# ============================================================================
#  1. calculating population statistics under no chemical exposure
# ============================================================================
library(expm)
library(popbio)

# demographic rates
s1 <- 0.95  # survival of stage 1
s2 <- 0.95  # survival of stage 2
s3 <- 0.95  # survival of stage 3
f3 <- 300   # fecundity of stage 3

# result storage for aseasonal model
lambda.null <- array(dim = c(2,6))
elast.null <- array(dim = c(2,3,3))
elast.4rates.null <- array(dim = c(2,4))

# result storage for seasonal model
mpm.A <- array(dim = c(2,10,10,10,6,6,3,3,60))
lambda <- array(dim = c(2,10,10,10,6,6,60))
sensit <- array(dim = c(2,10,10,10,6,6,3,3,60))
elast <- array(dim = c(2,10,10,10,6,6,3,3,60))
elast.4rates <- array(dim = c(2,10,10,10,6,6,4,60))

for(g in 1:2){  # variations in growth probability
  g1 <- c(0.1, 0.9)[g]  # growth of stage 1
  g2 <- c(0.1, 0.9)[g]  # growth of stage 2
  
  # per-phase matrix during the breeding season
  mpm.r <- matrix(c(s1*(1-g1), 0,         f3,
                    s1*g1,     s2*(1-g2), 0,
                    0,         s2*g2,     s3), 3, byrow = T)
  
  # population statistics of the aseasonal model
    # annaul population growth rate
    lambda.null[g,] <- Re(eigen(mpm.r)$values[1]) ^ (1:6 * 10)
    # elasticities of per-phase population growth rate to matrix elements
    elast.null[g,,] <- elasticity(mpm.r) 
    # elasticities of per-phase population growth rate to demographic rates
    elast.4rates.null[g,] <- c(sum(elast.null[g,,1:2]), 
                              sum(elast.null[g,,] * 
                                    matrix(c(-s1*g1, 0,      0,
                                              s1*g1,  -s2*g2, 0,
                                              0,      s2*g2,  0), 
                                            3, byrow = T) / 
                                    mpm.r, na.rm = T),
                              elast.null[g,1,3],
                              elast.null[g,3,3])
  
  for(wj in 1:10){  # variations in seasonal coefficient of juvenile survival
    w1 <- 0.1 * wj 
    w2 <- 0.1 * wj 
    
    for(wa in 1:10){  # variations in seasonal coefficient of adult survival
      w3 <- 0.1 * wa 
      
      for(wg1 in 1:10){  # variations in seasonal coefficient of growth
        wg <- 0.1 * wg1
        
        # per-phase matrix during the non-breeding season
        mpm.w <- matrix(c(s1*w1*(1-g1*wg), 0,               0,
                          s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                          0,               s2*w2*g2*wg,     s3*w3),
                        3, byrow = T)
        
        for(n in 1:6){  # variations in the number of phases
          n.season <- 10 * n 
          
          for(p in 1:6){  # variations in the proportion of the breeding season
            p.repseason <- 0.1 * p 
            
            # duration of breeding
            n.replength <- n.season * p.repseason 
            # duration of non-breeding
            n.nonlength <- n.season - n.replength 
            
            # array of seasonal matrices
            mpms <- array(dim = c(3,3,2*n.season))
            for(i in 1:(2*n.season)){
              if(i %in% c(1:n.replength, n.season + 1:n.replength)){
                mpms[,,i] <- mpm.r
              }else{
                mpms[,,i] <- mpm.w
              }
            }
            
            # population statistics of the seasonal model
            for(i in 1:n.season){
              mpm.Bi <- mpms[,,i]
              mpm.Di <- mpms[,,i+1]
              for(j in 2:(n.season-1)){
                mpm.Di <- mpms[,,i+j] %*% mpm.Di
              }
              # annual population matrix A
              mpm.Ai <- mpm.Di %*% mpm.Bi
              mpm.A[g,wj,wa,wg1,n,p,,,i] <- mpm.Ai
              # annual population growth rate lambda
              lambda.i <- Re(eigen(mpm.Ai)$values[1])
              lambda[g,wj,wa,wg1,n,p,i] <- lambda.i
              # leading right and left eigenvectors of A 
              w.i <- Re(eigen(mpm.Ai)$vectors[,1])
              v.i <- Re(eigen(t(mpm.Ai))$vectors[,1])
              # sensitivity 
              sensit.i <- t(mpm.Di) %*% matrix(abs(v.i), length(v.i)) %*% matrix(abs(w.i), 1) / 
                sum(abs(w.i)*abs(v.i))
              sensit[g,wj,wa,wg1,n,p,,,i] <- sensit.i
              # elasticity
              elast.i <- t(mpm.Di) %*% matrix(abs(v.i), length(v.i)) %*% matrix(abs(w.i), 1) / 
                sum(abs(w.i)*abs(v.i)) / lambda.i * mpm.Bi
              elast[g,wj,wa,wg1,n,p,,,i] <- elast.i
              elast.4rates[g,wj,wa,wg1,n,p,,i] <- c(sum(elast.i[,1:2]), 
                                                    sum(elast.i * 
                                                          matrix(c(-mpm.Bi[2,1], 0,            0,
                                                                   mpm.Bi[2,1],  -mpm.Bi[3,2], 0,
                                                                   0,            mpm.Bi[3,2],  0), 
                                                                 3, byrow = T) / 
                                                          mpm.Bi, na.rm = T),
                                                    elast.i[1,3],
                                                    elast.i[3,3])
            }
          }
        }
      }
    }
  }
}



# ============================================================================
#  2. Threshold chemical concentration that decreases the annual population 
#       growth rate below 1
# ============================================================================
library(expm)

# ----------------------------------------
#   2.1 Functions used for calculation 
# ----------------------------------------

thres.null.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                      s1*cx*g1,     s2*cx*(1-g2), 0,
                      0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  abs(log(lambda.expy))
}

thres.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                        s1*cx*g1,     s2*cx*(1-g2), 0,
                        0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*cx*(1-g1*wg), 0,                    0,
                        s1*w1*cx*g1*wg,     s2*w2*cx*(1-g2*wg), 0,
                        0,                    s2*w2*cx*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  abs(log(lambda.exp))
}

thres.null.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  abs(log(lambda.expy))
}

thres.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3*cx),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  abs(log(lambda.exp))
}

thres.null.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  abs(log(lambda.expy))
}

thres.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  abs(log(lambda.exp))
}

thres.null.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                      s1*g1*cx,     s2*(1-g2*cx), 0,
                      0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  abs(log(lambda.expy))
}

thres.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                        s1*g1*cx,     s2*(1-g2*cx), 0,
                        0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg*cx), 0,                    0,
                        s1*w1*g1*wg*cx,     s2*w2*(1-g2*wg*cx), 0,
                        0,                    s2*w2*g2*wg*cx,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  abs(log(lambda.exp))
}

# ----------------------------------------
#   2.2 calculation  
# ----------------------------------------

# result storage
lambda.null.threshold.sj <- array(dim=c(2,6))
lambda.threshold.sj <- array(dim = c(2,10,10,10,6,6))
lambda.null.threshold.sa <- array(dim=c(2,6))
lambda.threshold.sa <- array(dim = c(2,10,10,10,6,6))
lambda.null.threshold.f <- array(dim=c(2,6))
lambda.threshold.f <- array(dim = c(2,10,10,10,6,6))
lambda.null.threshold.g <- array(dim=c(2,6))
lambda.threshold.g <- array(dim = c(2,10,10,10,6,6))

# Hill coefficient
h <- 4

pb <- txtProgressBar(min = 1, max = 72000, style = 3)
for(g in 1:2){
  g1 <- c(0.1, 0.9)[g]
  g2 <- c(0.1, 0.9)[g]
  for(n in 1:6){
    n.season <- 10 * n 
    res.null <- optim(0, thres.null.sj, method = "Brent", lower = -5, upper = 5)
    lambda.null.threshold.sj[g,n] <- res.null$par
    res.null <- optim(0, thres.null.sa, method = "Brent", lower = -5, upper = 5)
    lambda.null.threshold.sa[g,n] <- res.null$par
    res.null <- optim(0, thres.null.f, method = "Brent", lower = -5, upper = 5)
    lambda.null.threshold.f[g,n] <- res.null$par
    res.null <- optim(0, thres.null.g, method = "Brent", lower = -5, upper = 5)
    lambda.null.threshold.g[g,n] <- res.null$par
    for(wj in 1:10){
      w1 <- 0.1 * wj 
      w2 <- 0.1 * wj 
      for(wa in 1:10){
        w3 <- 0.1 * wa 
        for(wg1 in 1:10){
          wg <- 0.1 * wg1
          for(p in which(Re(lambda[g,wj,wa,wg1,n,,1])>=1)){
            p.repseason <- 0.1 * p 
            n.replength <- n.season * p.repseason
            res <- optim(0, thres.sj, method = "Brent", lower = -5, upper = 5)
            lambda.threshold.sj[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, thres.sa, method = "Brent", lower = -5, upper = 5)
            lambda.threshold.sa[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, thres.f, method = "Brent", lower = -5, upper = 5)
            lambda.threshold.f[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, thres.g, method = "Brent", lower = -5, upper = 5)
            lambda.threshold.g[g,wj,wa,wg1,n,p] <- res$par
            setTxtProgressBar(pb, (g-1)*36000 + (n-1)*6000 + (wj-1)*600 + 
                                (wa-1)*60 + (wg1-1)*6 + p)
          }
        }
      }
    }
  }
}



# ============================================================================
#  3. Threshold chemical concentration that decreases the annual population 
#       growth rate by five percent
# ============================================================================
library(expm)

# ----------------------------------------
#   3.1 Functions used for calculation 
# ----------------------------------------

ec95.null.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                      s1*cx*g1,     s2*cx*(1-g2), 0,
                      0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  mpm <- matrix(c(s1*(1-g1), 0,         f3,
                  s1*g1,     s2*(1-g2), 0,
                  0,            s2*g2,     s3), 3, byrow = T)
  mpm.y <- mpm %^% n.season
  lambda.y <- Re(eigen(mpm.y)$values[1])
  abs(lambda.expy/lambda.y - 0.95)
}

ec95.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                        s1*cx*g1,     s2*cx*(1-g2), 0,
                        0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*cx*(1-g1*wg), 0,                    0,
                        s1*w1*cx*g1*wg,     s2*w2*cx*(1-g2*wg), 0,
                        0,                    s2*w2*cx*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  mpm.r <- matrix(c(s1*(1-g1), 0,         f3,
                    s1*g1,     s2*(1-g2), 0,
                    0,         s2*g2,     s3), 3, byrow = T)
  mpm.w <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                    s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                    0,                 s2*w2*g2*wg,     s3*w3),
                  3, byrow = T)
  mpm <- (mpm.w %^% (n.season-n.replength)) %*% (mpm.r %^% n.replength)
  lambda <- Re(eigen(mpm)$values[1])
  abs(lambda.exp/lambda - 0.95)
}

ec95.null.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  mpm <- matrix(c(s1*(1-g1), 0,         f3,
                  s1*g1,     s2*(1-g2), 0,
                  0,         s2*g2,     s3), 3, byrow = T)
  mpm.y <- mpm %^% n.season
  lambda.y <- Re(eigen(mpm.y)$values[1])
  abs(lambda.expy/lambda.y - 0.95)
}

ec95.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3*cx),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  mpm.r <- matrix(c(s1*(1-g1), 0,         f3,
                    s1*g1,     s2*(1-g2), 0,
                    0,         s2*g2,     s3), 3, byrow = T)
  mpm.w <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                    s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                    0,                 s2*w2*g2*wg,     s3*w3),
                  3, byrow = T)
  mpm <- (mpm.w %^% (n.season-n.replength)) %*% (mpm.r %^% n.replength)
  lambda <- Re(eigen(mpm)$values[1])
  abs(lambda.exp/lambda - 0.95)
}

ec95.null.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  mpm <- matrix(c(s1*(1-g1), 0,         f3,
                  s1*g1,     s2*(1-g2), 0,
                  0,         s2*g2,     s3), 3, byrow = T)
  mpm.y <- mpm %^% n.season
  lambda.y <- Re(eigen(mpm.y)$values[1])
  abs(lambda.expy/lambda.y - 0.95)
}

ec95.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  mpm.r <- matrix(c(s1*(1-g1), 0,         f3,
                    s1*g1,     s2*(1-g2), 0,
                    0,         s2*g2,     s3), 3, byrow = T)
  mpm.w <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                    s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                    0,                 s2*w2*g2*wg,     s3*w3),
                  3, byrow = T)
  mpm <- (mpm.w %^% (n.season-n.replength)) %*% (mpm.r %^% n.replength)
  lambda <- Re(eigen(mpm)$values[1])
  abs(lambda.exp/lambda - 0.95)
}

ec95.null.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                      s1*g1*cx,     s2*(1-g2*cx), 0,
                      0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  lambda.expy <- Re(eigen(mpm.expy)$values[1])
  mpm <- matrix(c(s1*(1-g1), 0,        f3,
                  s1*g1,     s2*(1-g2), 0,
                  0,         s2*g2,     s3), 3, byrow = T)
  mpm.y <- mpm %^% n.season
  lambda.y <- Re(eigen(mpm.y)$values[1])
  abs(lambda.expy/lambda.y - 0.95)
}

ec95.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                        s1*g1*cx,     s2*(1-g2*cx), 0,
                        0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg*cx), 0,                    0,
                        s1*w1*g1*wg*cx,     s2*w2*(1-g2*wg*cx), 0,
                        0,                    s2*w2*g2*wg*cx,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  lambda.exp <- Re(eigen(mpm.exp)$values[1])
  mpm.r <- matrix(c(s1*(1-g1), 0,         f3,
                    s1*g1,     s2*(1-g2), 0,
                    0,         s2*g2,     s3), 3, byrow = T)
  mpm.w <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                    s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                    0,                 s2*w2*g2*wg,     s3*w3),
                  3, byrow = T)
  mpm <- (mpm.w %^% (n.season-n.replength)) %*% (mpm.r %^% n.replength)
  lambda <- Re(eigen(mpm)$values[1])
  abs(lambda.exp/lambda - 0.95)
}

# ----------------------------------------
#   3.2 Calculation 
# ----------------------------------------

# result storage
lambda.null.ec95.sj <- array(dim=c(2,6))
lambda.ec95.sj <- array(dim = c(2,10,10,10,6,6))
lambda.null.ec95.sa <- array(dim=c(2,6))
lambda.ec95.sa <- array(dim = c(2,10,10,10,6,6))
lambda.null.ec95.f <- array(dim=c(2,6))
lambda.ec95.f <- array(dim = c(2,10,10,10,6,6))
lambda.null.ec95.g <- array(dim=c(2,6))
lambda.ec95.g <- array(dim = c(2,10,10,10,6,6))

# Hill coefficient
h <- 4 

pb <- txtProgressBar(min = 1, max = 72000, style = 3)
for(g in 1:2){
  g1 <- c(0.1, 0.9)[g]
  g2 <- c(0.1, 0.9)[g]
  for(n in 1:6){
    n.season <- 10 * n 
    res.null <- optim(0, ec95.null.sj, method = "Brent", lower = -5, upper = 5)
    lambda.null.ec95.sj[g,n] <- res.null$par
    res.null <- optim(0, ec95.null.sa, method = "Brent", lower = -5, upper = 5)
    lambda.null.ec95.sa[g,n] <- res.null$par
    res.null <- optim(0, ec95.null.f, method = "Brent", lower = -5, upper = 5)
    lambda.null.ec95.f[g,n] <- res.null$par
    res.null <- optim(0, ec95.null.g, method = "Brent", lower = -5, upper = 5)
    lambda.null.ec95.g[g,n] <- res.null$par
    for(wj in 1:10){
      w1 <- 0.1 * wj 
      w2 <- 0.1 * wj 
      for(wa in 1:10){
        w3 <- 0.1 * wa 
        for(wg1 in 1:10){
          wg <- 0.1 * wg1
          for(p in 1:6){
            p.repseason <- 0.1 * p 
            n.replength <- n.season * p.repseason
            res <- optim(0, ec95.sj, method = "Brent", lower = -5, upper = 5)
            lambda.ec95.sj[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, ec95.sa, method = "Brent", lower = -5, upper = 5)
            lambda.ec95.sa[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, ec95.f, method = "Brent", lower = -5, upper = 5)
            lambda.ec95.f[g,wj,wa,wg1,n,p] <- res$par
            res <- optim(0, ec95.g, method = "Brent", lower = -5, upper = 5)
            lambda.ec95.g[g,wj,wa,wg1,n,p] <- res$par
            setTxtProgressBar(pb, (g-1)*36000 + (n-1)*6000 + (wj-1)*600 + 
                                (wa-1)*60 + (wg1-1)*6 + p)
          }
        }
      }
    }
  }
}


# ============================================================================
#  4. Response of the annual population growth rate to chemical exposure
# ============================================================================
library(expm)

# ----------------------------------------
#   4.1 Functions to obtain the annual
#         populatioon growth rate
# ----------------------------------------

prop.dec.null.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                      s1*cx*g1,     s2*cx*(1-g2), 0,
                      0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  Re(eigen(mpm.expy)$values[1])
}

prop.dec.sj <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*cx*(1-g1), 0,            f3,
                        s1*cx*g1,     s2*cx*(1-g2), 0,
                        0,            s2*cx*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*cx*(1-g1*wg), 0,                    0,
                        s1*w1*cx*g1*wg,     s2*w2*cx*(1-g2*wg), 0,
                        0,                    s2*w2*cx*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  Re(eigen(mpm.exp)$values[1])
}

prop.dec.null.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  Re(eigen(mpm.expy)$values[1])
}

prop.dec.sa <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3*cx), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3*cx),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  Re(eigen(mpm.exp)$values[1])
}

prop.dec.null.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                      s1*g1,     s2*(1-g2), 0,
                      0,         s2*g2,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  Re(eigen(mpm.expy)$values[1])
}

prop.dec.f <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1), 0,         f3*cx,
                        s1*g1,     s2*(1-g2), 0,
                        0,         s2*g2,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg), 0,                 0,
                        s1*w1*g1*wg,     s2*w2*(1-g2*wg), 0,
                        0,                 s2*w2*g2*wg,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  Re(eigen(mpm.exp)$values[1])
}

prop.dec.null.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                      s1*g1*cx,     s2*(1-g2*cx), 0,
                      0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.expy <- mpm.exp %^% n.season
  Re(eigen(mpm.expy)$values[1])
}

prop.dec.g <- function(x){
  cx <- 1 / (1 + exp(h*x))
  mpm.r.exp <- matrix(c(s1*(1-g1*cx), 0,            f3,
                        s1*g1*cx,     s2*(1-g2*cx), 0,
                        0,            s2*g2*cx,     s3), 3, byrow = T)
  mpm.w.exp <- matrix(c(s1*w1*(1-g1*wg*cx), 0,                    0,
                        s1*w1*g1*wg*cx,     s2*w2*(1-g2*wg*cx), 0,
                        0,                    s2*w2*g2*wg*cx,     s3*w3),
                      3, byrow = T)
  mpm.exp <- (mpm.w.exp %^% (n.season-n.replength)) %*% (mpm.r.exp %^% n.replength)
  Re(eigen(mpm.exp)$values[1])
}

# ----------------------------------------
#   4.2 Calculation
# ----------------------------------------

# result storage
lambda.null.prop.dec.sj <- array(dim=c(2,6,102))
lambda.prop.dec.sj <- array(dim = c(2,10,10,10,6,6,102))
lambda.null.prop.dec.sa <- array(dim=c(2,6,102))
lambda.prop.dec.sa <- array(dim = c(2,10,10,10,6,6,102))
lambda.null.prop.dec.f <- array(dim=c(2,6,102))
lambda.prop.dec.f <- array(dim = c(2,10,10,10,6,6,102))
lambda.null.prop.dec.g <- array(dim=c(2,6,102))
lambda.prop.dec.g <- array(dim = c(2,10,10,10,6,6,102))

# Hill coefficient
h <- 4 
logconc <- c(-Inf, -50:50/25)
pb <- txtProgressBar(min = 1, max = 72000, style = 3)
for(g in 1:2){
  g1 <- c(0.1, 0.9)[g]
  g2 <- c(0.1, 0.9)[g]
  for(n in 1:6){
    n.season <- 10 * n 
    for(q in 1:102){
      conc <- logconc[q]
      lambda.null.prop.dec.sj[g,n,q] <- prop.dec.null.sj(conc)
      lambda.null.prop.dec.sa[g,n,q] <- prop.dec.null.sa(conc)
      lambda.null.prop.dec.f[g,n,q] <- prop.dec.null.f(conc)
      lambda.null.prop.dec.g[g,n,q] <- prop.dec.null.g(conc)
    }
    for(wj in 1:10){
      w1 <- 0.1 * wj 
      w2 <- 0.1 * wj 
      for(wa in 1:10){
        w3 <- 0.1 * wa 
        for(wg1 in 1:10){
          wg <- 0.1 * wg1
          for(p in 1:6){
            p.repseason <- 0.1 * p 
            n.replength <- n.season * p.repseason
            for(q in 1:102){
              conc <- logconc[q]
              lambda.prop.dec.sj[g,wj,wa,wg1,n,p,q] <- prop.dec.sj(conc)
              lambda.prop.dec.sa[g,wj,wa,wg1,n,p,q] <- prop.dec.sa(conc)
              lambda.prop.dec.f[g,wj,wa,wg1,n,p,q] <- prop.dec.f(conc)
              lambda.prop.dec.g[g,wj,wa,wg1,n,p,q] <- prop.dec.g(conc)
            }
            setTxtProgressBar(pb, (g-1)*36000 + (n-1)*6000 + (wj-1)*600 + 
                                (wa-1)*60 + (wg1-1)*6 + p)
          }
        }
      }
    }
  }
}

save.image("Three-stage_model1.RData")


# ============================================================================
#  5. Figures
# ============================================================================
library(tagcloud)
library(scales)

# ----------------------------------------
#   Fig S8 (a)
# ----------------------------------------
range(log(lambda[,,,,,,1])/log(10))
range(log(lambda.null)/log(10))
freq <- array(0, dim=c(6,6,90))
for(i in 1:2){
  for(j in 1:6){
    for(k in 1:6){
      freq[i,j,] <- freq[i,j,] + hist(log(lambda[i,,,,k,j,1])/log(10) - 
                                        log(lambda.null[i,k])/log(10),
                                      breaks = -80:10*2, plot=F)$counts
    }
  }
}

png("Three-stage1_Fig2.png", width = 3240, height = 3240, res = 300)
par(mfrow=c(1,1), mar=c(5,7,4,2))
plot(NA,NA,bty="L", xlim=c(-160,20), ylim=c(0,max(apply(freq,3,sum))),
     xlab="", ylab="", main="", cex.axis=1.8, las=1)
rect(-80:9*2,0,-79:10*2,apply(freq[,1,],2,sum),
     col = alpha("black", 0))
freq.cumsum <- apply(apply(freq,c(2,3),sum),2,cumsum)
for(k in 2:6){
  rect(-80:9*2, freq.cumsum[k-1,], -79:10*2, freq.cumsum[k,],
       col = alpha("black", 0.2*(k-1)))
}
segments(0,0,0,max(apply(freq,3,sum)), col=1, lwd=5, lty=3)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
dev.off()


# ----------------------------------------
#   Fig S9
# ----------------------------------------
library(scales)

# Preparing subdivided data 
ther.sj <- array(0,dim=c(6,6,53))
for(k in 1:2){
  for(j in 1:6){
    for(i in 1:6){
      exp.cr <- lambda.threshold.sj[k,,,,j,i] 
      cnt.cr <- lambda.null.threshold.sj[k,j]
      exp.cr[exp.cr > 4.99999] <- Inf
      if(cnt.cr > 4.99999) cnt.cr <- Inf
      hist.cr0 <- (exp.cr - cnt.cr)
      hist.cr <- hist.cr0
      hist.cr[is.nan(hist.cr) | hist.cr == Inf | hist.cr == -Inf] <- NA
      ther.sj[j,i,1:50] <- ther.sj[j,i,1:50] + hist((hist.cr)/log(10), breaks = -25:25/5, plot=F)$counts
      ther.sj[j,i,51] <- ther.sj[j,i,51] + length(which(is.nan(hist.cr0)))
      ther.sj[j,i,52] <- ther.sj[j,i,52] + length(which(hist.cr0 == Inf))
      ther.sj[j,i,53] <- ther.sj[j,i,53] + length(which(hist.cr0 == -Inf))
    }
  }
}
ther.sj.cumsum <- apply(apply(ther.sj,c(2,3),sum),2,cumsum)

ther.sa <- array(0,dim=c(6,6,53))
for(k in 1:2){
  for(j in 1:6){
    for(i in 1:6){
      exp.cr <- lambda.threshold.sa[k,,,,j,i] 
      cnt.cr <- lambda.null.threshold.sa[k,j]
      exp.cr[exp.cr>4.99999] <- Inf
      if(cnt.cr>4.99999) cnt.cr <- Inf
      hist.cr0 <- (exp.cr - cnt.cr)
      hist.cr <- hist.cr0
      hist.cr[is.nan(hist.cr) | hist.cr == Inf | hist.cr == -Inf] <- NA
      ther.sa[j,i,1:50] <- ther.sa[j,i,1:50] + hist((hist.cr)/log(10), breaks = -25:25/5, plot=F)$counts
      ther.sa[j,i,51] <- ther.sa[j,i,51] + length(which(is.nan(hist.cr0)))
      ther.sa[j,i,52] <- ther.sa[j,i,52] + length(which(hist.cr0 == Inf))
      ther.sa[j,i,53] <- ther.sa[j,i,53] + length(which(hist.cr0 == -Inf))
    }
  }
}
ther.sa.cumsum <- apply(apply(ther.sa,c(2,3),sum),2,cumsum)

ther.f <- array(0,dim=c(6,6,53))
for(k in 1:2){
  for(j in 1:6){
    for(i in 1:6){
      exp.cr <- lambda.threshold.f[k,,,,j,i] 
      cnt.cr <- lambda.null.threshold.f[k,j]
      exp.cr[exp.cr>4.99999] <- Inf
      if(cnt.cr>4.99999) cnt.cr <- Inf
      hist.cr0 <- (exp.cr - cnt.cr)
      hist.cr <- hist.cr0
      hist.cr[is.nan(hist.cr) | hist.cr == Inf | hist.cr == -Inf] <- NA
      ther.f[j,i,1:50] <- ther.f[j,i,1:50] + hist((hist.cr)/log(10), breaks = -25:25/5, plot=F)$counts
      ther.f[j,i,51] <- ther.f[j,i,51] + length(which(is.nan(hist.cr0)))
      ther.f[j,i,52] <- ther.f[j,i,52] + length(which(hist.cr0 == Inf))
      ther.f[j,i,53] <- ther.f[j,i,53] + length(which(hist.cr0 == -Inf))
    }
  }
}
ther.f.cumsum <- apply(apply(ther.f,c(2,3),sum),2,cumsum)

ther.g <- array(0,dim=c(6,6,53))
for(k in 1:2){
  for(j in 1:6){
    for(i in 1:6){
      exp.cr <- lambda.threshold.g[k,,,,j,i] 
      cnt.cr <- lambda.null.threshold.g[k,j]
      exp.cr[exp.cr>4.99999] <- Inf
      if(cnt.cr>4.99999) cnt.cr <- Inf
      hist.cr0 <- (exp.cr - cnt.cr)
      hist.cr <- hist.cr0
      hist.cr[is.nan(hist.cr) | hist.cr == Inf | hist.cr == -Inf] <- NA
      ther.g[j,i,1:50] <- ther.g[j,i,1:50] + hist((hist.cr)/log(10), breaks = -25:25/5, plot=F)$counts
      ther.g[j,i,51] <- ther.g[j,i,51] + length(which(is.nan(hist.cr0)))
      ther.g[j,i,52] <- ther.g[j,i,52] + length(which(hist.cr0 == Inf))
      ther.g[j,i,53] <- ther.g[j,i,53] + length(which(hist.cr0 == -Inf))
    }
  }
}
ther.g.cumsum <- apply(apply(ther.g,c(2,3),sum),2,cumsum)

# Drawing histograms 

png("Three-stage1_Fig3.png", width = 5040, height = 3240, res = 300)
par(oma=c(0,0,0,0)); layout(matrix(c(rep(1,18),rep(2,3),
                                     rep(3,6),rep(4,5),rep(5,7), rep(6,3)), nrow = 2, ncol = 21, byrow = T))

par(mar=c(9.1,9.1,4.1,2.1))
plot(NA,NA,bty="n", xlim=c(0,sum(ther.sj[,,1:53])), ylim=c(0,4.5),
     xlab="", ylab="", main="", cex.axis=3, las=1, xaxt="n", yaxt="n")
rect(c(0,cumsum(tapply(ther.sj.cumsum[6,-53], c(rep(1,50),2:3), sum))),3,
     cumsum(tapply(ther.sj.cumsum[6,], c(rep(1,50),2:4), sum)),4, density = c(1:4)*5, lwd=1)
rect(c(0,cumsum(tapply(ther.sa.cumsum[6,-53], c(rep(1,50),2:3), sum))),2,
     cumsum(tapply(ther.sa.cumsum[6,], c(rep(1,50),2:4), sum)),3, density = c(1:4)*5, lwd=1)
rect(c(0,cumsum(tapply(ther.sa.cumsum[6,-53], c(rep(1,50),2:3), sum))),2,
     cumsum(tapply(ther.sa.cumsum[6,], c(rep(1,50),2:4), sum)),3, density = c(1:4)*5, angle=135, lwd=1)
rect(c(0,cumsum(tapply(ther.f.cumsum[6,-53], c(rep(1,50),2:3), sum))),1,
     cumsum(tapply(ther.f.cumsum[6,], c(rep(1,50),2:4), sum)),2, density = c(1:4)*5, lwd=1)
rect(c(0,cumsum(tapply(ther.g.cumsum[6,-53], c(rep(1,50),2:3), sum))),0,
     cumsum(tapply(ther.g.cumsum[6,], c(rep(1,50),2:4), sum)),1, density = c(1:4)*5, lwd=1)
text(par()$usr[1], 4.5, "(a)", cex=3.5, adj = 0, xpd = T)
par(mgp=c(3,1.5,0))
axis(1, (seq(0, sum(ther.sj[,,1:53]), length=5)), 
     format(seq(0, sum(ther.sj[,,1:53]), length=5)/sum(ther.sj[,,1:53]),nsmall=2),cex.axis=3)
par(mgp=c(3,1,0))

par(mar=c(9.1,0,4.1,0))
plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")
rect(0.2, 0:3/5,
     0.4, 1:4/5,
     density = 1:4 *5)
rect(0.2, 0:3/5,
     0.4, 1:4/5,
     density = 1:4 *5, angle = c(1,3,1,3)*45)

par(mar=c(5.1,12.1,4.1,2.1),mgp=c(9,1,0))
plot(NA,NA,bty="n", xlim=c(-2,2), ylim=c(0,max(apply(ther.sj[,,1:53],3,sum))),
     xlab="", ylab="Frequency", main="", cex.axis=3, cex.lab=3, las=1, xaxt="n", yaxt="n")
axis(2, 0:10*6000, 0:10*6000, cex.axis=3, las=1)
par(mgp=c(3,1.5,0))
axis(1, -2:2*2, -2:2*2, cex.axis=3)
par(mgp=c(3,1,0))
rect(c(-25:24/5),0,
     c(-24:25/5),apply(ther.sj[,1,1:50],2,sum),
     col = alpha("black", 0))
for(k in 2:6){
  rect(c(-25:24/5), ther.sj.cumsum[k-1,1:50], 
       c(-24:25/5), ther.sj.cumsum[k,1:50],
       col = alpha("black", 0.2*(k-1)))
}
text(-3, max(apply(ther.sj[,,1:53],3,sum))*1.05, "(b)", cex=3.5, adj = 0, xpd=T)
segments(0,0,0,max(apply(ther.sj,3,sum)), col=1, lwd=3, lty=3)

par(mar=c(5.1,5.1,4.1,2.1))
plot(NA,NA,bty="n", xlim=c(-2,2), ylim=c(0,max(apply(ther.f[,,1:53],3,sum))),
     xlab="", ylab="", main="", cex.axis=3, las=1, xaxt="n", yaxt="n")
axis(2, 0:10*6000, 0:10*6000, cex.axis=3, las=1)
par(mgp=c(3,1.5,0))
axis(1, -2:2*2, -2:2*2, cex.axis=3)
par(mgp=c(3,1,0))
rect(c(-25:24/5),0,
     c(-24:25/5),apply(ther.f[,1,1:50],2,sum),
     col = alpha("black", 0))
for(k in 2:6){
  rect(c(-25:24/5), ther.f.cumsum[k-1,1:50], 
       c(-24:25/5), ther.f.cumsum[k,1:50],
       col = alpha("black", 0.2*(k-1)))
}
text(-3, max(apply(ther.f[,,1:53],3,sum))*1.05, "(c)", cex=3.5, adj = 0, xpd=T)
segments(0,0,0,max(apply(ther.f,3,sum)), col=1, lwd=3, lty=3)

plot(NA,NA,bty="n", xlim=c(-4,2), ylim=c(0,max(apply(ther.g[,,1:53],3,sum))),
     xlab="", ylab="", main="", cex.axis=3, las=1, xaxt="n", yaxt="n")
axis(2, 0:10*6000, 0:10*6000, cex.axis=3, las=1)
par(mgp=c(3,1.5,0))
axis(1, -2:2*2, -2:2*2, cex.axis=3)
par(mgp=c(3,1,0))
rect(c(-25:24/5),0,
     c(-24:25/5),apply(ther.g[,1,1:50],2,sum),
     col = alpha("black", 0))
for(k in 2:6){
  rect(c(-25:24/5), ther.g.cumsum[k-1,1:50], 
       c(-24:25/5), ther.g.cumsum[k,1:50],
       col = alpha("black", 0.2*(k-1)))
}
text(-5, max(apply(ther.g[,,1:53],3,sum))*1.05, "(d)", cex=3.5, adj = 0, xpd=T)
segments(0,0,0,max(apply(ther.g,3,sum)), col=1, lwd=3, lty=3)

par(mar=c(5.1,0,4.1,0))
plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")
rect(0.2, 0:5/6,
     0.4, 1:6/6,
     col=alpha("black",0.2*(0:5)))

layout(matrix(1)); par(mar=c(5.1,4.1,4.1,2.1), mfrow=c(1,1), oma=c(0,0,0,0))
dev.off()


# ----------------------------------------
#   Fig S11
# ----------------------------------------
png("Three-stage1_Fig4.png", width = 4000, height = 3200, res = 300)
layout(matrix(c(rep(1,2),rep(2,2),rep(3,1),
                rep(4,2),rep(5,2),rep(6,1)), nrow = 2, ncol = 5, byrow = T))
par(mar=c(3,3,3,3),oma=c(4,4,0,0))

# Elasticity to s1 
x1 <- apply(elast.4rates[1,,,1,6,1,1,1:60],c(1,2),sum) - elast.4rates.null[1,1] * 60
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[x1 < 0] <- cols.b[x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = NA,lwd=3)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

# Elasticity to s2 
x1 <- apply(elast.4rates[1,,,1,6,1,4,1:60],c(1,2),sum) - elast.4rates.null[1,4] * 60
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[x1 < 0] <- cols.b[x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = NA,lwd=3)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

# Legend 
plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")
rect(0.2, 0:17/22.5,
     0.35, 1:18/22.5,
     col=alpha(c(smoothPalette(9:1,pal="Blues"),smoothPalette(1:9,pal="Reds")),0.7),
     border = NA)
segments(0.4,0,0.4,0.8,lwd=2,col="gray25")
segments(0.4,0:4*0.2,0.55,0:4*0.2,lwd=2,col="gray25")

# Elasticity to f 
x1 <- apply(elast.4rates[1,,,1,6,1,3,1:60],c(1,2),sum) - elast.4rates.null[1,3] * 60
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[x1 < 0] <- cols.b[x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = NA,lwd=3)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

# Elasticity to g 
x1 <- apply(elast.4rates[1,,,1,6,1,2,1:60],c(1,2),sum) - elast.4rates.null[1,2] * 60
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),60,0), max = 60, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[x1 < 0] <- cols.b[x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = NA,lwd=3)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

# Legend 
plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")
rect(0.1,6:8/10,0.4,7:9/10,
     col = c("black", NA, NA), 
     density = c(NA,15,0),
     border = "gray25",lwd=2)

par(mar=c(5.1,4.1,4.1,2.1), mfrow=c(1,1), oma=c(0,0,0,0))
layout(matrix(1))
dev.off()


# ----------------------------------------
#   Fig S13
# ----------------------------------------

library(tagcloud)
library(scales)

png("Three-stage1_Fig5.png", width = 4200, height = 6000, res = 300)
layout(matrix(c(rep(1,2),rep(2,1),rep(3,2),
                rep(4,2),rep(5,1),rep(6,2),
                rep(7,2),rep(8,1),rep(9,2),
                rep(10,2),rep(11,1),rep(12,2)), nrow = 4, ncol = 5, byrow = T))
par(mar=c(5,6,4,1),oma=c(5,7,0,5))
par(mgp=c(3,1.5,0))

x1 <- lambda.ec95.sj[1,,,1,6,1] - lambda.null.ec95.sj[1,6]
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[!is.na(x1) & x1 < 0] <- cols.b[!is.na(x1) & x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
par(mgp=c(3,1.5,0))
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

par(mar=c(5,0,3,2))
plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")
rect(0.1, 0:17/36 + 0.4,
     0.25, 1:18/36 + 0.4,
     col=alpha(c(smoothPalette(9:1,pal="Blues"),smoothPalette(1:9,pal="Reds")),0.7),
     border = NA)
segments(0.3,0.4,0.3,0.9,lwd=2,col="gray25")
segments(0.3,seq(0.5,1,length=5)-0.1,0.45,seq(0.5,1,length=5)-0.1,lwd=2,col="gray25")

rect(0.1,1:3/10-0.1,0.25,2:4/10-0.1,
     col = c("black", NA, NA), 
     density = c(NA,15,0),
     border = "gray25",lwd=2)

par(mar=c(5,7,5,1))
plot(logconc, 1-(lambda.null.prop.dec.sj[1,6,]/lambda.null.prop.dec.sj[1,6,1]), type="l", col=1, ylim=c(0,1), lwd=3, 
     ylab="", xlab="", las=1, cex.lab=3, cex.axis=2.5)
points(logconc, 1-(lambda.prop.dec.sj[1,1,10,1,6,1,]/lambda.prop.dec.sj[1,1,10,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.5), lwd=3)
points(logconc, 1-(lambda.prop.dec.sj[1,10,1,1,6,1,]/lambda.prop.dec.sj[1,10,1,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.2), lwd=3)

par(mar=c(5,6,4,1))
x1 <- lambda.ec95.sa[1,,,1,6,1] - lambda.null.ec95.sa[1,6]
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[!is.na(x1) & x1 < 0] <- cols.b[!is.na(x1) & x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
par(mgp=c(3,1.5,0))
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")

par(mar=c(5,7,5,1))
plot(logconc, 1-(lambda.null.prop.dec.sa[1,6,]/lambda.null.prop.dec.sa[1,6,1]), type="l", col=1, ylim=c(0,1), lwd=3, 
     ylab="", xlab="", las=1, cex.lab=3, cex.axis=2.5)
points(logconc, 1-(lambda.prop.dec.sa[1,1,10,1,6,1,]/lambda.prop.dec.sa[1,1,10,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.5), lwd=3)
points(logconc, 1-(lambda.prop.dec.sa[1,10,1,1,6,1,]/lambda.prop.dec.sa[1,10,1,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.2), lwd=3)

par(mar=c(5,6,4,1))
x1 <- lambda.ec95.f[1,,,1,6,1] - lambda.null.ec95.f[1,6]
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[!is.na(x1) & x1 < 0] <- cols.b[!is.na(x1) & x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
par(mgp=c(3,1.5,0))
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")

par(mar=c(5,7,5,1))
plot(logconc, 1-(lambda.null.prop.dec.f[1,6,]/lambda.null.prop.dec.f[1,6,1]), type="l", col=1, ylim=c(0,1), lwd=3, 
     ylab="", xlab="", las=1, cex.lab=3, cex.axis=2.5)
points(logconc, 1-(lambda.prop.dec.f[1,1,10,1,6,1,]/lambda.prop.dec.f[1,1,10,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.5), lwd=3)
points(logconc, 1-(lambda.prop.dec.f[1,10,1,1,6,1,]/lambda.prop.dec.f[1,10,1,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.2), lwd=3)

par(mar=c(5,6,4,1))
x1 <- lambda.ec95.g[1,,,1,6,1] - lambda.null.ec95.g[1,6]
cols.a <- smoothPalette(c(sapply(x1, function(x) max(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Reds")[1:9]))), 
                        na.color = "black")
cols.b <- smoothPalette(c(sapply(x1, function(x) -1*min(x, 0)),0.8,0), max = 0.8, 
                        palfunc = colorRampPalette(smoothPalette(
                          1:9,palfunc = colorRampPalette(smoothPalette(1:9,pal="Blues")[1:9]))), 
                        na.color = "black")
cols.a <- cols.a[-(101:102)]; cols.b <- cols.b[-(101:102)]
cols.a[!is.na(x1) & x1 < 0] <- cols.b[!is.na(x1) & x1 < 0]
plot(1,1,type="n",xlim = c(0.5,10.5),ylim=c(0.5,10.5),las=1,bty="n",xaxt="n",yaxt="n",
     xlab="",
     ylab="")
par(mgp=c(3,1.5,0))
axis(1, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
axis(2, 1+(0:3)*3, format((1+(0:3)*3)/10,nsmall = 1), las=1, cex.axis = 2.5)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = alpha(cols.a[1:100],alpha = 0.7),
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     density = c(0, 15)[as.numeric(lambda[1,,,1,6,1,1]<1)+1],
     border = "gray25",lwd=2)
rect(rep(1:10,10)-0.5, rep(1:10,each=10)-0.5,
     rep(1:10,10)+0.5, rep(1:10,each=10)+0.5,
     col = c(NA, 1)[as.numeric(round(apply(elast[1,,,1,6,1,,,],c(1,2),sum)/60,4)!=1)+1],
     border = "gray25",lwd=2)

plot(NA,NA,bty="n", xlim=c(0,1), ylim=c(0,1),
     xlab="", ylab="", main="", xaxt="n", yaxt="n")

par(mar=c(5,7,5,1))
plot(logconc, 1-(lambda.null.prop.dec.g[1,6,]/lambda.null.prop.dec.g[1,6,1]), type="l", col=1, ylim=c(0,1), lwd=3, 
     ylab="", xlab="", las=1, cex.lab=3, cex.axis=2.5)
points(logconc, 1-(lambda.prop.dec.g[1,1,10,1,6,1,]/lambda.prop.dec.g[1,1,10,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.5), lwd=3)
points(logconc, 1-(lambda.prop.dec.g[1,10,1,1,6,1,]/lambda.prop.dec.g[1,10,1,1,6,1,1])^(1), type="l", 
       col=alpha(1,0.2), lwd=3)

dev.off()

