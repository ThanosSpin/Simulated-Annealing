## Model configurations, energy, magnetization, probability factors and partition function ####

binseq <- function(N) {
  expand.grid(replicate(N, c(-1,1), simplify=FALSE))
}

all.configs <- function(N) {
  X <- binseq(N*N)
  data <- list()
  for (i in 1:(2^{N*N})) {
    data[[i]] <- matrix(X[i,],nrow = N,ncol = N)
  }
  data
}

all.configs.2 <- all.configs(2)
all.configs.3 <- all.configs(3)
all.configs.4 <- all.configs(4)

alldown.config <- function(N) { matrix(-1,N,N) }
allup.config <- function(N) { matrix(1,N,N) }

per.neigh <- function(N,i,j) {
  if (i == 1) im1 <- N else im1 <- i-1
  if (j == 1) jm1 <- N else jm1 <- j-1
  if (i == N) ip1 <- 1 else ip1 <- i+1
  if (j == N) jp1 <- 1 else jp1 <- j+1
  list(c(im1,j),c(i,jm1),c(i,jp1),c(ip1,j))
}

random.config <- function(N) {
  config <- matrix(0,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      config[i,j] <- 2*sample(0:1,1)-1
    }
  }
  config
}

energy <- function(config, J) {
  config <- as.matrix(config)
  N <- nrow(config);
  e <- 0;
  for (i in 1:N) {
    for (j in 1:N) {
      coord <- per.neigh(N,i,j)
      for (b in seq_along(coord)) {
        w <- coord[[b]]; 
        e <- e - (J * as.numeric(config[i,j]) * as.numeric(config[w[1],w[2]]))
      }
    }
  }
  e
}

bolz.fact <- function(energy, temp) { exp(- energy/temp) }

magnetization <- function(config) { sum(as.vector(config)) }

partition <- function(configurations, J, temp) {
  k <- length(configurations);
  probs <- rep(NA,k)
  for (i in 1:k) {
    config <- configurations[[i]];
    probs[i] <- bolz.fact(energy(config,J),temp);
  }
  sum(probs)
}

prob.config <- function(config, part, J, temp) {
  bolz.fact(energy(config, J), temp)/part
}

## Metropolis algorithm ####

flip.config <- function(config) {
  i <- sample(1:N, 1)
  j <- sample(1:N, 1)
  config[i,j] <- -config[i,j]
  config
}

update.config <- function(config, J, temp) {
  new.config <- flip.config(config)
  dE <- energy(new.config, J) - energy(config, J)
  factor <- suppressWarnings(rbinom(1,1,bolz.fact(dE, temp)))
  if (is.na(factor)  & factor == 1) {return(new.config)}
    
    else {
      return(config)}
}

update.config.rep <- function(config, J, temp, L) {
  for (n in 1:L) {
    config <- update.config(config,J,temp)
  }
  config
}
update.config.rep.data <- function(config, J, temp, L) {
  e <- rep(0,L)
  m <- rep(0,L)
  for (n in 1:L) {
    config <- update.config(config,J,temp);
    e[n] <- energy(config,J)
    m[n] <- magnetization(config)
  }
  data <- list(e,m)
}
run.sim <- function(N, J, temp, iterations) {
  config <- random.config(N)
  update.config.rep(config, J, temp, iterations)
}
plot.config <- function(config) {
  image(config, useRaster=TRUE, axes=FALSE, col = grey(seq(0, 1)))
}
  

isinghm=function(niter,n,m=n,beta){
  x=sample(c(0,1),n*m,prob=c(0.5,0.5),replace=TRUE)
  x=matrix(x,n,m)
  for (i in 1:niter){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
      for (l in 1:m){
        n0=xneig4(x,sampl1[k],sampl2[l],x[sampl1[k],sampl2[l]])
        n1=xneig4(x,sampl1[k],sampl2[l],1-x[sampl1[k],sampl2[l]])
        if (runif(1)<exp(beta*(n1-n0)))
          x[sampl1[k],sampl2[l]]=1-x[sampl1[k],sampl2[l]]
      }}}
  x
}

prepa = runif(1, 0, 2)
prop = isinghm(1000, 4, 4, prepa)
image(1:128, 1:128,prop)
