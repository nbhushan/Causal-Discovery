# For Peregrine only:
# module add R/3.4.2-foss-2016a-X11-20160819
# n.bhushan@rug.nl


# TABLE OF CONTENTS
# 0. Prepare environment
# 1. Describe the simulation study
# 2. Set up parameters of the study
#    2.1 Prepare the program for parallel processing
# 3. Run the simulation
# 4. Analyze data after simulation                    
# 

# 0. Prepare environment ----
rm(list=ls())
if (!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])
list.of.packages <- c("simsalapar", "sn", "CompareCausalNetworks", "pcalg", "backShift","SID","doParallel","foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

require(simsalapar)
require(sn)
require(CompareCausalNetworks)
require(pcalg)
require(backShift)
require(SID)
require(doParallel)
require(foreach)


# 1. Describe the simulation study ----
# Theories in psychology are often tentative ideas about causal relationships
# between variables underlying the phenomenon under study and graphical
# causal models are an useful representation of such theories. Historically, the
# gold standard to test the validity of such theories are randomized controlled
# trials (RCTs). But oftentimes, RCTs are not feasible due to various regulatory,
# ethical, or practical constraints. In such cases, causal discovery algorithms help
# discover probabilistic causal relationships between variables of interest from
# observational data. In this talk, we access three such procedures which use
# conditional independence (d-separation) to infer underlying causal structures;
# the PC algorithm (Spirtes et al., 2000), LiNGAM (Shimizu et al., 2006), and
# the FCI algorithm (Zhang, 2008). The PC algorithm assumes a linear model
# with Gaussian errors and no latent variables. The LiNGAM algorithm relaxes
# the Gaussian error assumption and retains assumptions of linearity and absence
# of latent variables. Due to the large number of possible confounders
# in psychology, it is unrealistic to assume the absence of latent variables. In
# such cases, the FCI algorithm allows for latent variables while retaining linear
# Gaussian assumptions. To validate these procedures, we perform a simulation
# study varying the sample size, number of variables, number of latent variables,
# error distribution, and graph density. We score these procedures using three
# metrics; the F score, structural Hamming distance (SHD), the structural intervention
# distance (SID); discuss the results of our study and discuss further
# implications of such procedures for theory development in psychology.
# # END SECTION

# 2. Set up parameters of the study ----
varList <- varlist(
  n.sim = list(type = "N", expr = quote(N[sim]), value = 15),
  n = list(type = "grid", value = c(50,150,300,600)),
  p = list(type = "grid", value = c(4,8,16,32)),
  s = list(type = "grid", value = c(0.1,0.3,0.6, 0.9)),
  alpha = list(type = "grid", value = c(-6,0,6)),
  rho = list(type = "grid", value = c(0,0.4,0.6)),
  algo = list(type = "inner", value = c("pc","LINGAM")))

varList2 <- varlist(
  n.sim = list(type = "N", expr = quote(N[sim]), value = 15),
  n = list(type = "grid", value = c(50,150)),
  p = list(type = "grid", value = c(4,8)),
  s = list(type = "grid", value = c(0.1,0.3)),
  alpha = list(type = "grid", value = c(-6,0,6)),
  rho = list(type = "grid", value = c(0,0.4,0.6)),
  algo = list(type = "inner", value = c("pc","LINGAM")))

# set up the workhorse doOne()
cdSim <- function(n, p, s, alpha, rho, algo){
  # generate A assuming variables are causally ordered
  cont <- TRUE
  while(cont){
    AS <- 0*diag(p)
    for (k in 1:p){
      for (k2 in 1:p){
        if(k2>k) AS[k,k2] <- rbinom(1,1,s)
      }
    }
    if(sum(AS)!=0) cont <- FALSE
  }
  #simulate standardised weights unformly between -1 and 1 and assign to adjacency matrix
  beta <- AS * matrix( runif(p^2,-1,1),nrow=p)
  ###### simulate noise
  Sigma <- matrix(rho,nrow=p,ncol=p)
  diag(Sigma) <- 1
  Alpha <- rep(alpha, p)
  #simualte from multivariate skewed normal distribution
  noise <- rmsn(n=n, xi=rep(0,p), Sigma, Alpha,  tau=0, dp=NULL)
  #generate data
  inv <- solve(diag(p) - beta)
  X <- noise%*%inv
  ## loop over all methods and compute and print/plot estimate
  for (method in algo){
    Ahat <- getParents(X, environment = NULL, method=method, alpha=0.1, 
                       onlyObservationalData=TRUE, directed = TRUE)
    assign(paste("wmat", method, sep = "_"), Ahat)
  }
  getAdj <- function(G){
    return (1L * (G!=0))
  }
  normalizedHammingdist <- function(distance){
    return (distance / ((p * (p-1))/2))
  }  
  normalizedSiddist <- function(distance){
    return (distance / (p * (p-1)))
  }  
  #compute metrics
  h1 <- hammingDist(getAdj(AS), wmat_pc)
  h2 <- hammingDist(getAdj(AS), wmat_LINGAM)
  #h3 <- hammingDist(getAdj(AS), getAdj(wmat_fci))
  s1 <- structIntervDist(getAdj(AS), wmat_pc)$sid
  s2 <- structIntervDist(getAdj(AS), wmat_LINGAM)$sid
 #s3 <- structIntervDist(getAdj(AS), getAdj(wmat_fci))$sid)
  #c(h1,h3,si1,si3)
  hamming <- c(normalizedHammingdist(h1),normalizedHammingdist(h2))
  sid <- c(normalizedSiddist(s1),normalizedSiddist(s2))
  temp <- rbind(hamming,sid)
  dimnames(temp)[[2]] <- c("pc","LINGAM")
  temp
}

###------- test code -------------------------------------------
#automatic test
#test <- doCheck(cdSim, varList, nChks = 5)

#manual test
testM2 <- cdSim(n = 100, p = 4, s = 0.4, alpha = 0, rho = 0, algo = c("pc","LINGAM"))

###------- The workhorse -------------------------------------------

#set seed
seedy <- c("seq")

#run smulation
res <- doForeach(varList, sfile = "cdSimNoFCI.rds", doOne = cdSim, keepSeed=TRUE,
                 seed=seedy, extraPkgs = c("sn", "CompareCausalNetworks", "SID"))

###------- Save results -------------------------------------------

val  <- getArray(res, "value")
err  <- getArray(res, "error")
warn <- getArray(res, "warning")
time <- getArray(res, "time")

###################################################
### results-as-df
###################################################
dval <- array2df(val)
derror <- array2df(err)
dwarn <- array2df(warn)
dtime <- array2df(time)

###################################################
### write tables to csv
###################################################
write.csv(dval, "valuesNew.csv")# % errors or warnings
write.csv(derror, "errors.csv")# % errors or warnings
write.csv(dwarn, "warnings.csv")# % errors or warnings
write.csv(dtime, "times.csv")# % errors or warnings