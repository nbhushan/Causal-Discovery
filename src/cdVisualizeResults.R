rm(list=ls())

require(simsalapar)
require(abind)
require(reshape2)
require(tidyr)
require(wesanderson)
require(xtable)
require(sjstats)
require(car)

## Make function nestedloop available to generate data for
## nested loop plot
source("scripts/nestedloop.r")
source("scripts/cdNestedViz.R")

#BEST DAG(sid lower bound)
resBest <- array(numeric(), c(2, 2, 4, 4, 4, 3, 3, 0))
resBest <- getArray(readRDS("Results/cdSim00.rds"), "value")
names(dimnames(resBest))[[1]]<-"metric"


#SID
dfBestSID <- getDf(resBest, metric = "sid")
dfBestseSID <- getDf_se(resBest, metric = "sid")
dfBestSID_N <- getDf_N(resBest, metric = "sid")

#SHD
dfBestSHD <- getDf(resBest, metric = "hamming")
dfBestseSHD <- getDf_se(resBest, metric = "hamming")



#Worst DAG(sid upper bound)
resWorst <- array(numeric(), c(2, 2, 4, 4, 4, 3, 3, 0))
resWorst <- getArray(readRDS("Results/cdSim01.rds"), "value")
names(dimnames(resWorst))[[1]]<-"metric"

#SID
dfWorstSID <- getDf(resWorst, metric = "sid")
dfWorstSIDse <- getDf_se(resWorst, metric = "sid")
dfWorstSID_n <- getDf_N(resWorst, metric = "sid")

#Merge SID
SID_n <- merge(dfBestSID_N,dfWorstSID_n, by=c("p","s", "alpha","rho","LINGAM"),suffixes = c(".lb",".ub"))
SID_n <- SID_n[with(SID_n, order(p, s, alpha, rho)), ]


SID <- merge(dfBestSID,dfWorstSID,by=c("n","p","s", "alpha","rho","LINGAM"),suffixes = c(".lb",".ub"))
SID <- SID[with(SID, order(n, p, s, alpha, rho)), ]

SID.se <- merge(dfBestseSID,dfWorstSIDse,by=c("n","p","s", "alpha","rho","LINGAM"),suffixes = c(".lb",".ub"))
SID.se <- SID.se[with(SID.se, order(n, p, s, alpha, rho)), ]

#SHD
SHD <- getDf(resWorst, metric = "hamming")
SHD.se <- getDf_se(resWorst, metric = "hamming")
SHD_N <- getDf_N(resWorst, metric = "hamming")



#Save graphs
nldata <- nestedloop(SHD.se,
                     varnames=c("n", "p", "s", "alpha", "rho"),
                     varlabels=
                       c("Sample size",
                         "No. variables",
                         "sparsity",
                         "non-normality", 
                         "confounding"))

nldata2 <- nestedloop(SHD_N,
                     varnames=c("p", "s", "alpha", "rho"),
                     varlabels=
                       c("No. variables",
                         "sparsity",
                         "non-normality", 
                         "confounding"))


createGraph(pd=nldata2, metric = "hamming", xlabel = "4 X 4 X 4 X 3 X 3 conditions", figName = "images/SHD_n.pdf")

anova.SID.LINGAM <- anovaSID(resBest)
anova.SID.PC <- anovaSID(resBest, algo = "pc")

anova.SHD.LINGAM <- anovaSID(resBest, metric = "hamming")
anova.SHD.PC <- anovaSID(resBest, metric="hamming", algo = "pc")

#save tables
#Obtain tables

#Run ANOVA
anovaSID <- function(x, metric = "sid", algo="LINGAM"){
  SIDRes <- x[metric = metric,algo=algo , , , , , ,]
  dat <- array2df(SIDRes)
  mod <- aov(value ~ n * p * s * alpha * rho, data = dat)
  carMod <- anova_stats(car::Anova(mod, type = 3, p.adjust.method=TRUE))
  carMod[,'etasq']<-NULL
  carMod[,'omegasq']<-NULL
  return(xtable(carMod))
}


