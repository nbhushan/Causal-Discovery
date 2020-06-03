if (!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])
list.of.packages <- c("simsalapar", "abind", "limmma", "pcalg", "xtable","SID","sjstats","car")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")



require(simsalapar)
require(abind)
require(limmma)
require(xtable)
require(sjstats)
require(car)
require(papeR)
rm(list=ls())

#Run ANOVA
anva <- function(x, metric="hamming"){
    dat <- x[metric=metric,,,,,,,]
    dat <- array2df(dat)
    if (metric=="hamming"){
      p <- as.numeric(levels(dat[,'p']))[dat[,'p']]
      dist <- dat[, 'value'] * ((p * (p-1))/2)
      dat[,'value']<-dist
    }
    dat[,'n.sim']<-NULL
    mod <- aov(value ~ algo * n * p * s * alpha * rho, data = dat)
    carMod <- anova_stats(car::Anova(mod, type = 3, p.adjust.method=TRUE))
    carMod[,'etasq']<-NULL
    carMod[,'omegasq']<-NULL
    return(xtable(carMod))
}

#Obtain tables
textab <- function(x, v, metric="hamming", a="0", r="0.0"){
  SIDRes <- x[metric=metric,,,,,alpha=a,rho=r,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  Var.mad <- apply(SIDRes, non.sim.margins, sd)
  fval <- formatC(Var, digits=1, format="f")
  fmad <- paste0("(", format(round(Var.mad, 1), scientific=FALSE, trim=TRUE), ")")
  nc <- nchar(fmad)
  sm <- nc == min(nc)
  fmad[sm] <- paste0("\\ \\,", fmad[sm])
  fres <- array(paste(fval, fmad),
                dim=dim(fval), dimnames=dimnames(fval))
  ft <- ftable(fres, row.vars=c("n","p"), col.vars=c("s", "algo"))
  return(toLatex(ft, vList = v))
}

#Obtain tables
graphs <- function(x, v, metric="hamming", a="0", r="0.0"){
  #plot SID
  if (metric=="hamming"){
    lab <- "SHD"
  }else (lab<-"SID")
  SIDRes <- res[metric=metric,,,,,alpha=a,rho=r,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  dimnames(Var)[["n"]] <- paste0("n==", dimnames(Var)[["n"]])
  dimnames(Var)[["p"]] <- paste0("p==", dimnames(Var)[["p"]])
  #dimnames(Var)[["alpha"]] <- paste0("alpha==", dimnames(Var)[["alpha"]])
  dimnames(Var)[["s"]] <- paste0("s==", dimnames(Var)[["s"]])
  return(mayplot(Var, v, row.vars = "p", col.vars = "s", xvar = "n",
          type="b", ylab = lab, ylim="local", do.n.sim=TRUE,
          main="test"))
}

#p on x axis
graphs2 <- function(x, v, metric="hamming", a="0", r="0.0"){
  #plot SID
  if (metric=="hamming"){
    lab <- "SHD"
  }else (lab<-"SID")
  SIDRes <- res[metric=metric,,,,,alpha=a,rho=r,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  dimnames(Var)[["n"]] <- paste0("n==", dimnames(Var)[["n"]])
  dimnames(Var)[["p"]] <- paste0("p==", dimnames(Var)[["p"]])
  #dimnames(Var)[["alpha"]] <- paste0("alpha==", dimnames(Var)[["alpha"]])
  dimnames(Var)[["s"]] <- paste0("s==", dimnames(Var)[["s"]])
  return(mayplot(Var, v, row.vars = "s", col.vars = "n", xvar = "p",
                 type="b", ylab = lab, ylim="local", do.n.sim=TRUE,
                 main="test"))
}

#keep p and alpha fixed, vary n,rho,s
graphs3 <- function(x, v, metric="hamming", pp="4", a="0"){
  #plot SID
  if (metric=="hamming"){
    lab <- "SHD"
  }else (lab<-"SID")
  SIDRes <- res[metric=metric,,,p=pp,,alpha=a,,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  dimnames(Var)[["n"]] <- paste0("n==", dimnames(Var)[["n"]])
  dimnames(Var)[["rho"]] <- paste0("rho==", dimnames(Var)[["rho"]])
  #dimnames(Var)[["alpha"]] <- paste0("alpha==", dimnames(Var)[["alpha"]])
  dimnames(Var)[["s"]] <- paste0("s==", dimnames(Var)[["s"]])
  return(mayplot(Var, v, row.vars = "s", col.vars = "n", xvar = "rho",
                 type="b", ylab = lab, ylim="local", do.n.sim=TRUE,
                 main="test"))
}


#keep p and rho fixed, vary n,alpha,s
graphs4 <- function(x, v, metric="hamming", pp="4", r="0.0"){
  #plot SID
  if (metric=="hamming"){
    lab <- "SHD"
  }else (lab<-"SID")
  SIDRes <- res[metric=metric,,,p=pp,,,rho=r,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  dimnames(Var)[["n"]] <- paste0("n==", dimnames(Var)[["n"]])
  #dimnames(Var)[["rho"]] <- paste0("rho==", dimnames(Var)[["rho"]])
  dimnames(Var)[["alpha"]] <- paste0("alpha==", dimnames(Var)[["alpha"]])
  dimnames(Var)[["s"]] <- paste0("s==", dimnames(Var)[["s"]])
  return(mayplot(Var, v, row.vars = "s", col.vars = "n", xvar = "alpha",
                 type="b", ylab = lab, ylim="local", do.n.sim=TRUE,
                 main="test"))
}


#Obtain boxplots
boxgraphs <- function(x, v, metric="hamming", a="0", r="0.0"){
  #plot SID
  if (metric=="hamming"){
    lab <- "SHD"
  }else (lab<-"SID")
  Var <- res[metric=metric,,,,,alpha=a,rho=r,]
  dimnames(Var)[["n"]] <- paste0("n==", dimnames(Var)[["n"]])
  dimnames(Var)[["p"]] <- paste0("p==", dimnames(Var)[["p"]])
  #dimnames(Var)[["alpha"]] <- paste0("alpha==", dimnames(Var)[["alpha"]])
  dimnames(Var)[["s"]] <- paste0("s==", dimnames(Var)[["s"]])
  return(mayplot(Var, v, row.vars = "p", col.vars = "n", xvar = "s",
                ylab = lab))
}


###------- Load Results -------------------------------------------
#first set source file
files <- list.files(path = '.', pattern = '^cdSim[0-9]\\.rds$')
FUNC <- function (x) getArray(readRDS(x), "value")
res <- array(numeric(),c(2,3,4,4,4,3,3,0)) 
for (f in files) {
  temp <- FUNC(f)
  res <- abind(res, temp, use.dnns = TRUE)
  rm(temp)
}
rm(f, files)

names(dimnames(res))[[1]]<-"metric"

varList <- varlist(
  n.sim = list(type = "N", expr = quote(N[sim]), value = 105),
  n = list(type = "grid", value = c(50,150,300,600)),
  p = list(type = "grid", value = c(4,8,16,32)),
  s = list(type = "grid", value = c(0.1,0.3,0.6, 0.9)),
  alpha = list(type = "grid", value = c(-6,0,6)),
  rho = list(type = "grid", value = c(0,0.4,0.6)),
  algo = list(type = "inner", value = c("pc","LINGAM", "fci")))

toLatex(varList)

###------- Save graphs -------------------------------------------

### impact on n


pdf(file <- "Images//SHDn.pdf", width=15, paper = "a4r") 
graphs(res, varList)
dev.off.pdf(file)

pdf(file <- "Images//SIDp.pdf", width=15, paper = "a4r") 
graphs2(res, varList, metric="sidd")
dev.off.pdf(file)



pdf(file <- "Images//HammingAlphaZeroRhoZero.pdf", width=15, paper = "a4r") 
graphs3(res, varList)
dev.off.pdf(file)

pdf(file <- "Images//SIDAlphaZeropFour.pdf", width=15, paper = "a4r") 
graphs3(res, varList, "sidd")
dev.off.pdf(file)

pdf(file <- "Images//HammingAlphaNegSixRhoZero.pdf", width=15, paper = "a4r") 
graphs(res, varList, a = "-6")
dev.off.pdf(file)

pdf(file <- "Images//SIDAlphaNegSixZeroRhoZero.pdf", width=15, paper = "a4r") 
graphs(res, varList, "sidd", a = "-6")
dev.off.pdf(file)

pdf(file <- "Images//HammingAlphaNegSixZeroRhoFour.pdf", width=15, paper = "a4r") 
graphs(res, varList, a="-6", r="0.4")
dev.off.pdf(file)

pdf(file <- "Images//SIDAlphaNegSixZeroRhoFour.pdf", width=15, paper = "a4r") 
graphs(res, varList, "sidd", a="-6", r="0.4")
dev.off.pdf(file)



###------- Output tables -------------------------------------------
shdTab <- textab(res, varList)
sidTab <- textab(res, varList, "sidd")

###Box-plots
pdf(file <- "Images//boxplot_SHDAlphaNegSixZeroRhoFour.pdf", width=15, paper = "a4r") 
boxgraphs(res, metric = "hamming", varList, a="0")
dev.off.pdf(file)

pdf(file <- "Images//boxplot_SIDAlphaNegSixZeroRhoFour.pdf", width=15, paper = "a4r") 
boxgraphs(res, metric = "sidd", varList, a="0")
dev.off.pdf(file)


#impact of latent confounders
pdf(file <- "Images//SHDAlphaZeroPFour.pdf", width=15, paper = "a4r") 
graphs3(res, varList, pp = "32")
dev.off.pdf(file)

#impact of non-normality
pdf(file <- "Images//SHDrhoZeroPFour.pdf", width=15, paper = "a4r") 
graphs4(res, varList)
dev.off.pdf(file)

pdf(file <- "Images//SIDRhoZeroPFour.pdf", width=15, paper = "a4r") 
graphs4(res, varList,metric = "sidd")
dev.off.pdf(file)

#Run ANOVASID
sidAnova <- anva(res, "sidd")
shdAnova <-  anva(res, "hamming")



