### Linear SEM
require(pcalg)
require(lavaan)
require(qgraph)
require(RColorBrewer)
require(CompareCausalNetworks)
require(eeptools)
require(MVN)
require(naniar)
require(mice)
#for limgam exploration
require(doSNOW)
require(rlecuyer)



rm(list = ls())


#Buurkracht
load("BuurkrachtScales.Rdata")

# some variables are reported as factors, but cor() only works for numericals
for (i in 1:ncol(Scales)) {
  if (class(Scales[[i]]) == "character") {
    Scales[[i]] <- as.numeric(Scales[[i]])
  }
}


#convert birthyear to age
t <- strptime(Scales$Birthyear, format = "%Y")
for (i in 1:length(t)) {
  if (is.na(t[i])) {
    t[i] <- as.POSIXct.Date(Sys.Date())
  }
}

age <- floor(age_calc(as.Date(t), units = "years"))
Scales$Birthyear <- age
#Scales$Birthyear<-NULL
rm(age, t, Membership, i)

Scales2 <- na.omit(Scales)

#temp
tScales <- mice(Scales)
completeScales <- complete(tScales,1)

#test for multivariate normality
mvn.result.mardia <- mvn(data = Scales2, mvnTest = "mardia")
mvn.result.energy <- mvn(data = Scales2, mvnTest = "energy")
mvn.result.hz <- mvn(data = Scales2, mvnTest = "hz")
mvn.result.royston  <- mvn(data = Scales2, mvnTest = "royston")



#groups scales by category
buur_groups <- list(
  `Personal factors` = c(1:8),
  `Factors related to the social context` = c(9:15),
  `Evaluations of energy companies and the government` = c(16:17),
  `Sustainable energy intentions and behaviours` = c(18:27),
  `Socio-demographics` = c(28:31),
  `Membership` <- c(32)
)

scales.personal <- Scales2[,c(1:8,23)]


# apply PC algorithm
buur.pc.fit <- getParents(
    Scales2,
    environment=NULL,
    method = c("pc"),
    onlyObservationalData=TRUE,
    pointConf = TRUE
  )
saveRDS(buur.pc.fit,"buur.pc.fit.RDS")

# apply LiNGAM
buur.lingam.fit <- getParents(
  Scales2,
  environment=NULL,
  method = c("LINGAM"),
  onlyObservationalData=TRUE,
  pointConf = TRUE
)
qgraph(
  buur.lingam.fit,
  layout = "spring",
  labels = TRUE,
  nodeNames = gsub("_", " ", colnames(Scales2)),
  legend.cex = 0.3,
  groups = buur_groups,
  palette = "ggplot2",
  legend = TRUE,
  legend.mode = "style1",
  legend.cex = 0.32,
  vsize = 3.8,
  layoutOffset = c(-.2, 0),
  #GLratio=1.0,
  minimum = 0.01,
  weighted = FALSE
)
saveRDS(buur.lingam.fit, "buur.lingam.fit.RDS")
#plot Lingam
buur.lingam.qgraph <-
  qgraph(
    buur.lingam.fit,
    layout = "spring",
    labels = TRUE,
    nodeNames = gsub("_", " ", colnames(Scales2)),
    legend.cex = 0.3,
    groups = buur_groups,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    layoutOffset = c(-.2, 0),
    #GLratio=1.0,
    filename = "buurLINGAM",
    filetype = "pdf",
    minimum = 0.01,
    weighted = FALSE
  )
saveRDS(buur.lingam.qgraph, "buur.lingam.qgraph.RDS")
#plot PC

buur.pc.qgraph <-
  qgraph(
    buur.pc.fit,
    layout = "spring",
    labels = TRUE,
    nodeNames = gsub("_", " ", colnames(Scales2)),
    legend.cex = 0.3,
    groups = buur_groups,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    layoutOffset = c(-.2, 0),
    #GLratio=1.0,
    filename = "buurPC",
    filetype = "pdf",
    minimum = 0.01
  )

saveRDS(buur.pc.qgraph,"buur.pc.qgraph.RDS")

#compare causal networks stability section
#LINGAM

#plot Lingam
# lingam test
thresh <- seq(0.6, 0.9, 0.1)
ev <- seq(2, 30, 2)
sim <- expand.grid(thresh, ev)
sim$edges <- NA
names(sim) <- c("Threshhold", "EV", "no. edges")
sim <- sim[,c(2,1,3)]


#sim function
buur.lingam.stable.sim <- function(ev, thresh, nsim){
  res <- getParentsStable(
    Scales2,
    threshold = thresh,
    environment = NULL,
    nodewise = T,
    sampleObservations = 0.7,
    EV = ev,
    nsim = nsim,
    onlyObservationalData = TRUE,
    method = c("pc"),
    verbose = F
  )
  return(sum(as.logical(res)))
}

# 4.b Run the analysis on all availible Cores
cluster<-makeCluster(4, type = "SOCK")
clusterSetupRNG(cluster, seed = 29012001) 
registerDoSNOW(cluster)

foreach(i= 1:nrow(sim),
        .packages ='CompareCausalNetworks') %do% {sim[i,3] <- buur.lingam.stable.sim(ev = sim[i,1], thresh = sim[i,2], nsim = 10)}

simPC <- sim

####
simLingam$Threshhold <- as.factor(simLingam$Threshhold)
simPC$Threshhold <- as.factor(simPC$Threshhold)

lin <- subset(simLingam, simLingam$Threshhold =="0.7")
pc <- subset(simPC, simPC$Threshhold =="0.7")

ggplot(simLingam, aes(x=EV, y=`no. edges`, main="LiNGAM")) + 
                                geom_line()  +
                                facet_grid(cols = vars(Threshhold))

ggplot(simPC, aes(x=EV, y=`no. edges`, main="PC")) + 
  geom_line()  +
  facet_grid(cols = vars(Threshhold))

###

buur.lingam.stable <-
  getParentsStable(
    Scales2,
    threshold = 0.75,
    environment = NULL,
    nodewise = T,
    sampleObservations = 0.7,
    EV = 2,
    nsim = 100,
    onlyObservationalData = TRUE,
    method = c("LINGAM"),
    verbose = F,
    setOptions=list(pointConf = TRUE)
  )
qgraph(buur.lingam.stable)
range(buur.lingam.stable)
sum(as.logical(buur.lingam.stable))
saveRDS(buur.lingam.stable,"buur.lingam.stable.RDS")


#PC
buur.pc.stable <-
  getParentsStable(
    Scales2,
    threshold = 0.75,
    environment = NULL,
    sampleObservations = 0.7,
    nodewise = T,
    EV = 2,
    nsim = 100,
    onlyObservationalData = TRUE,
    method = c("pc"),
    verbose = FALSE
  )
range(buur.pc.stable)
sum(as.logical(buur.pc.stable))
saveRDS(buur.pc.stable,"buur.pc.stable.RDS")

#plot stability Lingam
buur.CCN.lingam.stable.qgraph <-
  qgraph(
    buur.lingam.stable,
    layout = "spring",
    labels = TRUE,
    nodeNames = gsub("_", " ", colnames(Scales)),
    legend.cex = 0.3,
    groups = buur_groups,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    filename = "buurLINGAM_stabilitySelection",
    filetype = "pdf"
    #minimum = 0.9,
  )

buur.pc.stable <- readRDS("buur.pc.stable.RDS")
buur.CCN.PC.stable.qgraph <-
  qgraph(
    buur.pc.stable,
    layout = "spring",
    labels = TRUE,
    nodeNames = gsub("_", " ", colnames(Scales)),
    legend.cex = 0.3,
    groups = buur_groups,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    #GLratio=1.0,
    filename = "buurPC_stabilitySelection",
    filetype = "pdf"
    #minimum = 0.9
  )

##################################################################################
#ESS Data
#### Load standardized data ####
MyData         <- readRDS("HandledData.rds")
longnames      <- readRDS("longnames.rds")
gr             <- readRDS("gr.rds")
namesRelabeled <- readRDS("namesRelabeled.rds")
MyData$countryLong <- NULL
MyData$country <- NULL

MyData2 <- MyData[!is.na(MyData)]
MyData2 <- na.omit(MyData)


# apply LiNGAM
ESS.lingam.fit <- getParents(
  X = MyData2,
  method = c("LINGAM"),
  mode = c("raw"),
  directed = TRUE,
  pointConf = TRUE
)

saveRDS(ESS.lingam.fit,"ESS.lingam.fit.RDS")

#Apply PC
# apply PC algorithm
ESS.pc.fit <- pc(
  suffStat = list(
    C = cor_auto(MyData2,
                 ordinalLevelMax = 5),
    n = nrow(MyData2)
  ),
  indepTest = gaussCItest,
  labels = longnames,
  alpha = 0.0000001
)

saveRDS(ESS.pc.fit,"ESS.pc.fit.RDS")

ESS.lingam.stable <-
  getParentsStable(
    MyData2,
    threshold = 0.75,
    environment = NULL,
    sampleObservations = 0.7,
    nodewise = T,
    EV = 2,
    nsim = 100,
    onlyObservationalData = TRUE,
    method = c("LINGAM"),
    verbose = FALSE
  )
range(ESS.lingam.stable)
sum(as.logical(ESS.lingam.stable))
saveRDS(ESS.lingam.stable,"ESS.lingam.stable.RDS")

ESS.PC.stable <-
  getParentsStable(
    MyData2,
    threshold = 0.75,
    environment = NULL,
    sampleObservations = 0.7,
    nodewise = T,
    EV = 2,
    nsim = 100,
    onlyObservationalData = TRUE,
    method = c("pc"),
    verbose = FALSE
  )
range(ESS.PC.stable)
sum(as.logical(ESS.PC.stable))
saveRDS(ESS.PC.stable,"ESS.PC.stable.RDS")

#plot ESS Lingam
ESS.lingam.qgraph <-
  qgraph(
    ESS.lingam.fit,
    layout = "spring",
    labels = namesRelabeled,
    nodeNames = longnames,
    legend.cex = 0.3,
    groups = gr,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    filename = "ESS_LINGAM",
    filetype = "pdf",
    minimum=0.1
  )

#plot ESS PC
ESS.PC.qgraph <-
  qgraph(
    ESS.pc.fit,
    layout = "spring",
    labels = namesRelabeled,
    nodeNames = longnames,
    legend.cex = 0.3,
    groups = gr,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "style1",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    filename = "ESS_PC",
    filetype = "pdf",
    minimum = 0.1
  )

#plot stability Lingam
ESS.lingam.stable.qgraph <-
  qgraph(
    ESS.lingam.stable,
    layout = "groups",
    labels = namesRelabeled,
    nodeNames = longnames,
    legend.cex = 0.3,
    groups = gr,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "groups",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    filename = "ESSLINGAM_stabilitySelectionCircle",
    filetype = "pdf",
    minimum = 75
  )

ESS.PC.stable.qgraph <-
  qgraph(
    ESS.PC.stable,
    layout = "groups",
    labels = namesRelabeled,
    nodeNames = longnames,
    legend.cex = 0.3,
    groups = gr,
    palette = "ggplot2",
    legend = TRUE,
    legend.mode = "groups",
    legend.cex = 0.32,
    vsize = 3.8,
    edge.labels = FALSE,
    edge.label.cex = 0.4,
    edge.label.bg = TRUE,
    edge.label.position = 0.5,
    weighted = FALSE,
    layoutOffset = c(-.2, 0),
    #GLratio=1.0,
    filename = "ESSPC_stabilitySelectionCircle",
    filetype = "pdf",
    minimum=75
  )
