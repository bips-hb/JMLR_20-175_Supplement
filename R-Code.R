#THIS SOURCE CODE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR AND THE JOURNAL
#OF MACHINE LEARNING RESEARCH (JMLR) AND JMLR'S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL 
#WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
#A PARTICULAR PURPOSE, AND ANYWARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND 
#RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR JMLR, NOR JMLR'S PUBLISHERS 
#AND DISTRIBUTORS, WILL BE LIABLE FOR DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting 
#the generality of the foregoing, neither the author, nor JMLR, nor JMLR's publishers and 
#distributors, warrant that the Source Code will be error-free, will operate without interruption, 
#or will meet the needs of the user.



####################################################################################################
### Code for Figure 5 in:                                                                        ###
### Witte J, Henckel L, Maathuis MH, Didelez V (2020): On efficient adjustment in causal graphs  ###
####################################################################################################

### R version 3.6.2
library(graph) # version 1.62.0
library(pcalg) # version 2.6-3

### create the DAG from Figure 4(b) with specified edge weights
trueDAG <- graphNEL(
   nodes=c("V_1", "X", "V_4", "V_3", "V_2", "V_5", "Y"),
   edgeL=list("V_1" = list(edges=c("X", "V_2", "V_3"), weights=c(1, 1, 1)),
              "X"   = list(edges=c("V_3", "V_4"),     weights=c(1, 1)),
              "V_4" = list(edges=NULL,               weights=NULL),
              "V_3" = list(edges=c("V_2", "V_5", "Y"), weights=c(0.5, 0.5, 1)),
              "V_2" = list(edges="Y",                weights=0.5),
              "V_5" = list(edges="Y",                weights=0.5),
              "Y"   = list(edges=NULL,               weights=NULL)),
   edgemode="directed")

### the true causal effect of X on Y is 1.5
ida(2, 7, trueCov(trueDAG), trueDAG, method="local", type="pdag")

### corresponding CPDAG (bi-directed edges in the plot are interpreted as undirected edges)
trueCPDAG <- dag2cpdag(trueDAG)
trueCPDAGamat <- as(trueCPDAG, "matrix")
dimnames(trueCPDAGamat) <- list(1:7, 1:7)

### the multiset of possible causal effects is (0, 0, 1.5, 2.5, 2.5)
round(ida(2, 7, trueCov(trueDAG), trueCPDAG, method="local", type="pdag"), 10)

### simulate data from trueDAG and estimate possible causal effects using local vs. optimal IDA
sim <- function(DAG, CPDAGamat) {
   simdata <- rmvDAG(40, DAG, "normal")
   local_result <- ida(2, 7, cov(simdata), CPDAGamat, method="local")
   O_result <- c(ida(2, 7, cov(simdata), CPDAGamat, method="optimal"))
   return(cbind(local_result, O_result))
}
set.seed(28359)
result <- replicate(1000, sim(trueDAG, trueCPDAGamat))

### create density plot in the style of Maathuis et al. (2009)
plot(density(result[ , 2, ]), xlab="Estimated causal effect", lty="dotted", main="")
lines(density(result[ , 1, ]))
lines(x=c(0, 0), y=c(0, 0.6), lwd=2)
lines(x=rep(2.5, 2), y=c(0, 0.6), lwd=2)
lines(x=rep(1.5, 2), y=c(0, 0.3), lwd=2)




####################################################################################################
### Code for the simulation study in:                                                            ###
### Witte J, Henckel L, Maathuis MH, Didelez V (2020): On efficient adjustment in causal graphs  ###
####################################################################################################

### R version 3.6.0 / 3.6.1
library(pcalg) # version 2.6-3
library(MASS) # version 1.4-5
library(parallel) # version 3.6.2
library(abind) # version 1.4-5
library(ggplot2) # version 3.1.1

### We manually deleted the isValidGraph check in ida() to speed up computation

####################################################################################################
### define functions                                                                             ###
####################################################################################################

sim2 <- function(n, p, truecovm, trueCPDAGamat, x, y) {
   # generate data according to trueDAG
   dat <- mvrnorm(n=n, mu=rep(0, p), Sigma=truecovm)
   colnames(dat) <- colnames(trueCPDAGamat)
   
   # estimate CPDAG
   score <- new("GaussL0penObsScore", dat)
   estCPDAG <- ges(score)
   estCPDAGamat <- 1*wgtMatrix(estCPDAG$essgraph)
   dimnames(estCPDAGamat) <- list(1:p, 1:p)
   
   # estimate possible causal effects using different variants of IDA
   # true CPDAG
   trueLocal <- abs(ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(trueCPDAGamat)))
   trueOptimal <- abs(ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(trueCPDAGamat), 
                          method="optimal"))
   
   # estimated CPDAG
   estLocal <- abs(ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(estCPDAGamat)))
   estOptimal <- abs(ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(estCPDAGamat), 
                         method="optimal"))
   
   # minimum absolute values
   tLocal <- min(trueLocal)
   tOptimal <- min(trueOptimal)
   eLocal <- min(estLocal)
   eOptimal <- min(estOptimal)

   return(c(tLocal, tOptimal, eLocal, eOptimal))
}

sim <- function(p, d, n){

   nonamenable <- FALSE
   nonnull <- FALSE
   while (!nonamenable | !nonnull) {
      
      wfun <- function(m) {-1^rbinom(m, 1, 0.5) * runif(m, 0.1, 1)}
      trueDAG <- randDAG(p, d, method="er", wFUN=wfun)
      truecovm <- trueCov(trueDAG, back.compatible=TRUE)
      DAGamat <- t(as(trueDAG, "matrix"))
      
      # draw x and y
      ok <- FALSE
      while (!ok) {
         x <- sample(1:p, 1)
         desX <- setdiff(possDe(DAGamat, x, possible=FALSE, type="dag"), x)
         if (length(desX) < 1) {next}
         if (length(desX)==1) {y <- desX} else {y <- sample(desX, 1)}
         ok <- TRUE
      }
      
      # true CPDAG
      trueCPDAG <- dag2cpdag(trueDAG)
      trueCPDAGamat <- t(as(trueCPDAG, "matrix"))
      dimnames(trueCPDAGamat) <- list(1:p, 1:p)
      
      nonamenable <- !pcalg:::isAmenable(trueCPDAGamat, x, y, type="cpdag")
      
      if (nonamenable) {
      # determine true smallest possible effect in CPDAG
      alltrue <- ida(x, y, trueCov(trueDAG, back.compatible=TRUE), t(trueCPDAGamat), 
                     type="pdag")
      true <- min(abs(alltrue))
      if (true>(10^(-7))) {nonnull <- TRUE}
      }
   }
   
   res2 <- replicate(100, sim2(n, p, truecovm, trueCPDAGamat, x, y))
   diffres2 <- ( res2 - true )^2
   MSE <- apply(diffres2, 1, mean)
   
   true_o_l <- MSE[2] / MSE[1]
   est_o_l <- MSE[4] / MSE[3]
   
   ratios <- c(true_o_l, est_o_l)
   names(ratios) <- c("o:l true", "o:l est")
   
   return(ratios)
}


####################################################################################################
### run simulation                                                                               ###
####################################################################################################

nrep <- 1000

set.seed(90641)
res_10_2_100 <- replicate(nrep, sim(p=10, d=2, n=100))
res_20_2_100 <- replicate(nrep, sim(p=20, d=2, n=100))
res_50_2_100 <- replicate(nrep, sim(p=50, d=2, n=100))
set.seed(33516)
res_10_2_1000 <- replicate(nrep, sim(p=10, d=2, n=1000))
res_20_2_1000 <- replicate(nrep, sim(p=20, d=2, n=1000))
set.seed(51901)
res_50_2_1000 <- replicate(nrep, sim(p=50, d=2, n=1000))
res_100_2_1000 <- replicate(nrep, sim(p=100, d=2, n=1000))
set.seed(00217)
res_10_3_100 <- replicate(nrep, sim(p=10, d=3, n=100))
res_20_3_100 <- replicate(nrep, sim(p=20, d=3, n=100))
res_50_3_100 <- replicate(nrep, sim(p=50, d=3, n=100))
res_10_3_1000 <- replicate(nrep, sim(p=10, d=3, n=1000))
set.seed(08745)
res_20_3_1000 <- replicate(nrep, sim(p=20, d=3, n=1000))
set.seed(24570)
res_50_3_1000 <- replicate(nrep, sim(p=50, d=3, n=1000))
set.seed(91555)
res_100_3_1000 <- replicate(nrep, sim(p=100, d=3, n=1000))
set.seed(23779)
res_10_4_100 <- replicate(nrep, sim(p=10, d=4, n=100))
res_20_4_100 <- replicate(nrep, sim(p=20, d=4, n=100))
set.seed(98426)
res_50_4_100 <- replicate(nrep, sim(p=50, d=4, n=100))
set.seed(46902)
res_10_4_1000 <- replicate(nrep, sim(p=10, d=4, n=1000))
res_20_4_1000 <- replicate(nrep, sim(p=20, d=4, n=1000))
set.seed(71798)
res_50_4_1000 <- replicate(nrep, sim(p=50, d=4, n=1000))
set.seed(32899)
res_100_4_1000 <- replicate(nrep, sim(p=100, d=4, n=1000))
set.seed(68470)
res_100_2_100 <- replicate(nrep, sim(p=100, d=2, n=100))

cl<-makeCluster(12)
clusterCall(cl, function() library(pcalg))
clusterCall(cl, function() library(MASS))
clusterExport(cl, c("sim", "sim2", "p", "n", "d"))
clusterSetRNGStream(cl, 98759)
pres_100_3_100 <- parLapply(cl, 1:nrep, function(j){
   sim(p=100, d=3, n=100)
})
res_100_3_100 <- do.call(cbind, pres_100_3_100)
clusterSetRNGStream(cl, 07963)
pres_100_4_100 <- parLapply(cl, 1:nrep, function(j){
   sim(p=100, d=4, n=100)
})
res_100_4_100 <- do.call(cbind, pres_100_4_100)
stopCluster(cl)


####################################################################################################
### plots                                                                                        ###
####################################################################################################

### Figure 6

# geometric mean and median of MSE ratios (scenario A)
geomean <- apply(res_100_4_1000, 1, function(x) {exp(mean(log(x)))} )
media <- apply(res_100_4_1000, 1, median, na.rm=TRUE)

# data for plotting (scenario A)
ggdata <- data.frame( x=rep(c("true CPDAG","estimated CPDAG"), each=1000),
                      y=c(res_100_4_1000[1, ], res_100_4_1000[2, ]) )
ggdata$x <- factor(ggdata$x, levels=c("true CPDAG","estimated CPDAG"))

# plot Figure 6 A
ggplot(data=ggdata, aes(x=x, y=y)) +
   geom_violin() + 
   geom_point(aes(x=1, y=geomean[1]), shape=20) + 
   geom_point(aes(x=2, y=geomean[2]), shape=20) + 
   geom_point(aes(x=1, y=media[1]), shape=3) + 
   geom_point(aes(x=2, y=media[2]), shape=3) +
   geom_hline(yintercept=1, lty=2) +
   ylim(0, 2.75) +
   annotate("text", label="A", x=0.6, y=2.6, size=8) +
   labs(x="", y="RMSE optimal IDA / local IDA") +
   theme(text=element_text(size=13, face="plain", color="black"),
         axis.text.x=element_text(c("true CPDAG","estimated CPDAG", color="black")),
         axis.title.y=element_text(size=12, face="plain"), axis.ticks.x=element_blank()) 

# geometric mean and median of MSE ratios (scenario B)
geomean <- apply(res_10_4_100, 1, function(x) {exp(mean(log(x)))} )
media <- apply(res_10_4_100, 1, median, na.rm=TRUE)

# data for plotting (scenario B)
ggdata <- data.frame( x=rep(c("true CPDAG","estimated CPDAG"), each=1000),
                      y=c(res_10_4_100[1, ], res_10_4_100[2, ]) )
ggdata$x <- factor(ggdata$x, levels=c("true CPDAG","estimated CPDAG"))

# plot Figure 6 B
ggplot(data=ggdata, aes(x=x, y=y)) +
   geom_violin() + 
   geom_point(aes(x=1, y=geomean[1]), shape=20) + 
   geom_point(aes(x=2, y=geomean[2]), shape=20) + 
   geom_point(aes(x=1, y=media[1]), shape=3) + 
   geom_point(aes(x=2, y=media[2]), shape=3) +
   geom_hline(yintercept=1, lty=2) +
   ylim(0, 2.75) +
   annotate("text", label="B", x=0.6, y=2.6, size=8) +
   labs(x="", y="RMSE optimal IDA / local IDA") +
   theme(text=element_text(size=13, face="plain", color="black"),
         axis.text.x=element_text(c("true CPDAG","estimated CPDAG", color="black")),
         axis.title.y=element_text(size=12, face="plain"), axis.ticks.x=element_blank())

### Figures 9 and 10

res <- abind(res_10_2_100, res_10_3_100, res_10_4_100, res_20_2_100, res_20_3_100, res_20_4_100,
             res_50_2_100, res_50_3_100, res_50_4_100, res_100_2_100, res_100_3_100, res_100_4_100,
             res_10_2_1000, res_10_3_1000, res_10_4_1000, res_20_2_1000, res_20_3_1000,
             res_20_4_1000, res_50_2_1000, res_50_3_1000, res_50_4_1000, res_100_2_1000,
             res_100_3_1000, res_100_4_1000, along=3)

# geometric mean and median of MSE ratios
geomean <- apply(res[ , , ], c(1,3), function(x) {exp(mean(log(x)))} )
media <- apply(res[ , , ], c(1,3), median, na.rm=TRUE)

# data for plotting - true CPDAG
p <- factor(rep(c("p=10","p=20","p=50","p=100"), times=2, each=3000),
            levels=c("p=10","p=20","p=50","p=100"))
d <- factor(rep(c("d=2","d=3","d=4"), times=8, each=1000))
n <- factor(rep(c("n=100","n=1000"), each=12000))
y <- c(res[1, , ])

ggdata <- data.frame(p=p, d=d, n=n, y=y)

p2 <- factor(rep(c("p=10","p=20","p=50","p=100"), times=2, each=3),
             levels=c("p=10","p=20","p=50","p=100"))
d2 <- factor(rep(c("d=2","d=3","d=4"), times=8))
n2 <- factor(rep(c("n=100","n=1000"), each=12))
yg <- geomean[1, ]
ym <- media[1, ]

gdata <- data.frame(p=p2, d=d2, n=n2, y=yg)
mdata <- data.frame(p=p2, d=d2, n=n2, y=ym)

# plot Figure 9
ggplot(data=ggdata, aes(x="", y=y)) + geom_violin() +
   geom_hline(yintercept=1, lty=2)  +
   ylim(0, 3.8) +
   labs(x="", y="RMSE optimal IDA / local IDA") +
   theme(text=element_text(size=13, face="plain", color="black"),
         axis.title.y=element_text(size=12, face="plain"), axis.ticks.x=element_blank(),
         axis.text.x=element_blank(), panel.grid.major.x=element_blank()) +
   geom_point(data=gdata, aes(x="", y=yg), shape=20) +
   geom_point(data=mdata, aes(x="", y=ym), shape=3) +
   facet_grid(p ~ n + d)

# data for plotting - estimated CPDAG
p <- factor(rep(c("p=10","p=20","p=50","p=100"), times=2, each=3000),
            levels=c("p=10","p=20","p=50","p=100"))
d <- factor(rep(c("d=2","d=3","d=4"), times=8, each=1000))
n <- factor(rep(c("n=100","n=1000"), each=12000))
y <- c(res[3, , ])

ggdata <- data.frame(p=p, d=d, n=n, y=y)

p2 <- factor(rep(c("p=10","p=20","p=50","p=100"), times=2, each=3),
             levels=c("p=10","p=20","p=50","p=100"))
d2 <- factor(rep(c("d=2","d=3","d=4"), times=8))
n2 <- factor(rep(c("n=100","n=1000"), each=12))
yg <- geomean[3, ]
ym <- media[3, ]

gdata <- data.frame(p=p2, d=d2, n=n2, y=yg)
mdata <- data.frame(p=p2, d=d2, n=n2, y=ym)

# plot Figure 10
ggplot(data=ggdata, aes(x="", y=y)) + geom_violin() +
   geom_hline(yintercept=1, lty=2)  +
   ylim(0, 3.8) +
   labs(x="", y="RMSE optimal IDA / local IDA") +
   theme(text=element_text(size=13, face="plain", color="black"),
         axis.title.y=element_text(size=12, face="plain"), axis.ticks.x=element_blank(),
         axis.text.x=element_blank(), panel.grid.major.x=element_blank()) +
   geom_point(data=gdata, aes(x="", y=yg), shape=20) +
   geom_point(data=mdata, aes(x="", y=ym), shape=3) +
   facet_grid(p ~ n + d)



####################################################################################################
### Additional code for the high-dimensional setting in:                                         ###
### Witte J, Henckel L, Maathuis MH, Didelez V (2020): On efficient adjustment in causal graphs  ###
####################################################################################################

### PART I #########################################################################################
### searching for graphs that are non-amenable and have a non-null minimum absolute causal effect###
####################################################################################################

library(pcalg) # version 2.6-3
library(parallel) # version 3.6.2

### We manually deleted the isValidGraph check in ida() and optAdjSet() and set checkInput=FALSE
### in hasExtension() to speed up computation

####################################################################################################
### PART I define functions                                                                      ###
####################################################################################################

sim <- function(p, d){
  counter <- 0
  
  nonamenable <- FALSE
  nonnull <- FALSE
  while (!nonamenable | !nonnull) {
    counter <- counter + 1
    cat(counter)
    
    wfun <- function(m) {-1^rbinom(m, 1, 0.5) * runif(m, 0.1, 1)}
    trueDAG <- randDAG(p, d, method="er", wFUN=wfun)
    truecovm <- trueCov(trueDAG, back.compatible=TRUE)
    DAGamat <- t(as(trueDAG, "matrix"))
    
    # draw x and y
    ok <- FALSE
    ok_counter <- 1
    while (!ok) {
      ok_counter <- ok_counter + 1
      if (ok_counter > 1500) {break}
      x <- sample(1:p, 1)
      desX <- setdiff(possDe(DAGamat, x, possible=FALSE, type="dag"), x)
      if (length(desX) < 1) {next}
      if (length(desX)==1) {y <- desX} else {y <- sample(desX, 1)}
      ok <- TRUE
    }
    
    if (!ok) {next}
    
    # true CPDAG
    trueCPDAG <- dag2cpdag(trueDAG)
    trueCPDAGamat <- t(as(trueCPDAG, "matrix"))
    dimnames(trueCPDAGamat) <- list(1:p, 1:p)
    
    nonamenable <- !pcalg:::isAmenable(trueCPDAGamat, x, y, type="cpdag")
    
    cat("h")
    
    if (nonamenable) {
      # reduce CPDAG to connected component of x and y
      ma <- ( ( trueCPDAGamat + t(trueCPDAGamat) )>0 )*1
      cc <- NULL
      cc_new <- c(x,y)
      while(!identical(cc,cc_new)) {
        cc <- cc_new
        cc_add <- unlist(sapply(cc, function(z){which(ma[z, ]==1)}))
        cc_new <- sort(unique(union(cc,cc_add)))
      }
      newCPDAGamat <- trueCPDAGamat[cc_new, cc_new]
      newcovm <- truecovm[cc_new, cc_new]
      newx <- which(colnames(newCPDAGamat)==x)
      newy <- which(colnames(newCPDAGamat)==y)
      newp <- length(cc_new)
      
      # determine true smallest possible effect in CPDAG
      alltrue <- ida(newx, newy, newcovm, t(newCPDAGamat), 
                      type="cpdag")
      true <- min(abs(alltrue))
      if (true>(10^(-7))) {nonnull <- TRUE}
    }
  }
  
  dimnames(newCPDAGamat) <- list(1:newp, 1:newp)
  dimnames(newcovm) <- list(1:newp, 1:newp)
  
  return( list(trueCPDAGamat, truecovm, x, y, newCPDAGamat, newcovm, newx, newy) )
}


####################################################################################################
### PART I run simulation                                                                        ###
####################################################################################################

nrep <- 100

cl<-makeCluster(100)
clusterCall(cl, function() library(pcalg, lib="~/Rpackages"))
clusterExport(cl, c("sim", "ida", "hasExtension", "optAdjSet"))

clusterSetRNGStream(cl, 47863) #d=1, p=1000
res_1000_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=1000, d=1)
})

clusterSetRNGStream(cl, 11435) #d=1, p=250, 750
res_250_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=250, d=1)
})
res_750_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=750, d=1)
})

clusterSetRNGStream(cl, 93196) #d=1, p=100, 50, 20, 10
res_100_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=100, d=1)
})
res_50_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=50, d=1)
})
res_20_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=20, d=1)
})
res_10_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=10, d=1)
})

clusterSetRNGStream(cl, 19478) #d=1, p=500
res_500_1_100 <- parLapply(cl, 1:nrep, function(j){
  sim(p=500, d=1)
})

stopCluster(cl)


### PART II ########################################################################################
### Data generation, graph search and causal effect estimation                                   ###
####################################################################################################

### R version 3.6.0 / 3.6.1
library(pcalg) # version 2.6-3
library(MASS) # version 1.4-5
library(abind) # version 1.4-5
library(ggplot2) # version 3.1.1
library(rcausal)

### We manually deleted the isValidGraph check in ida() to speed up computation

####################################################################################################
### PART II define functions                                                                     ###
####################################################################################################

sim2 <- function(re, n, p, trueCPDAGamat, truecovm, x, y, newCPDAGamat, newcovm, newx, newy) {
  
  ma <- ( ( trueCPDAGamat + t(trueCPDAGamat) )>0 )*1
  cc <- NULL
  cc_new <- c(x,y)
  while(!identical(cc,cc_new)) {
    cc <- cc_new
    cc_add <- unlist(sapply(cc, function(z){which(ma[z, ]==1)}))
    cc_new <- sort(unique(union(cc,cc_add)))
  }
  
  # generate data according to trueDAG
  dat <- mvrnorm(n=n, mu=rep(0, p), Sigma=truecovm)
  
  # estimate possible causal effects using different variants of IDA
  # true CPDAG
  trueLocal <- abs(ida(x.pos=newx, y.pos=newy, mcov=cov(dat)[cc_new, cc_new], graphEst=t(newCPDAGamat)))
  trueOptimal <- abs(ida(x.pos=newx, y.pos=newy, mcov=cov(dat)[cc_new, cc_new], graphEst=t(newCPDAGamat), 
                          method="optimal"))
  
  # estimate CPDAG
  tetradest <- tetradrunner(df=dat, algoId="fges", scoreId="sem-bic", dataType="continuous",
                            numCategoriesToDiscretize=NULL, verbose=FALSE, penaltyDiscount=0.5)
  
  estCPDAGamat <- matrix(0, nrow=p, ncol=p)
  for (i in tetradest$edges) {
    A <- as.numeric(unlist(regmatches(i,gregexpr("[[:digit:]]+\\.*[[:digit:]]*", i)))[1])
    B <- as.numeric(unlist(regmatches(i,gregexpr("[[:digit:]]+\\.*[[:digit:]]*", i)))[2])
    arrow <- gsub('[0-9]+', '', i)
    estCPDAGamat[B,A] <- 1
    if (arrow==" --- ") {
      estCPDAGamat[A,B] <- 1 
    }
  }
  
  dimnames(estCPDAGamat) <- list(1:p, 1:p)
  
  # estimated CPDAG
  ma <- ( ( estCPDAGamat + t(estCPDAGamat) )>0 )*1
  cc <- NULL
  cc_nw <- c(x,y)
  while(!identical(cc,cc_nw)) {
    cc <- cc_nw
    cc_add <- unlist(sapply(cc, function(z){which(ma[z, ]==1)}))
    cc_nw <- sort(unique(union(cc,cc_add)))
  }
  nwCPDAGamat <- estCPDAGamat[cc_nw, cc_nw]
  nwx <- which(colnames(nwCPDAGamat)==x)
  nwy <- which(colnames(nwCPDAGamat)==y)
  nwp <- length(cc_nw)
  dimnames(nwCPDAGamat) <- list(1:nwp, 1:nwp)
  
  
  # estimate possible causal effects using different variants of IDA
  # estimated CPDAG
  estLocal <- abs(ida(x.pos=nwx, y.pos=nwy, mcov=cov(dat)[cc_nw,cc_nw], graphEst=t(nwCPDAGamat)))
  estOptimal <- abs(ida(x.pos=nwx, y.pos=nwy, mcov=cov(dat)[cc_nw,cc_nw], graphEst=t(nwCPDAGamat), 
                         method="optimal"))
  
  # minimum absolute values
  tLocal <- min(trueLocal)
  tOptimal <- min(trueOptimal)
  eLocal <- min(estLocal)
  eOptimal <- min(estOptimal)
  
  return(c(tLocal, tOptimal, eLocal, eOptimal))
}


sim <- function(re, n, p){
  
  trueCPDAGamat <- re[[1]]
  truecovm <- re[[2]]
  x <- re[[3]]
  y <- re[[4]]
  newCPDAGamat <- re[[5]]
  newcovm <- re[[6]]
  newx <- re[[7]]
  newy <- re[[8]]
  
  # true causal effect
  
  alltrue <- ida(newx, newy, newcovm, t(newCPDAGamat), type="cpdag")
  true <- min(abs(alltrue))
  
  res2 <- tryCatch( replicate(10, sim2(re, n, p, trueCPDAGamat, truecovm, x, y, newCPDAGamat, newcovm, newx, newy)), error=function(e){paste("Error")} )
  if (res2 %in% "Error") {
    ratios <- c(NA, NA)
    names(ratios) <- c("o:l true", "o:l est")
    return(ratios)
  }
  
  diffres2 <- ( res2 - true )^2
  MSE <- apply(diffres2, 1, mean)
  
  true_o_l <- MSE[2] / MSE[1]
  est_o_l <- MSE[4] / MSE[3]
  
  ratios <- c(true_o_l, est_o_l)
  names(ratios) <- c("o:l true", "o:l est")
  
  return(c(ratios, MSE))
}


getEndResult <- function(finBIC) {
  fin2 <- matrix(unlist(finBIC), nrow=6)
  NAs <- apply(fin2, 2, anyNA)
  fin3 <- fin2[ ,!NAs]
  
  geomean <- apply(fin3, 1, function(x) {exp(mean(log(x)))} )
  media <- apply(fin3, 1, median)
  
  return(rbind(geomean, media))
}


####################################################################################################
### PART II run simulation                                                                       ###
####################################################################################################

set.seed(76770)
finBIC <- lapply(res_1000_1_100, sim, n=100, p=1000)
# 4 errors
getEndResult(finBIC)
# TRUE: geo 0.6876592, med 0.7456218
# EST: geo 0.7942288, med 1.0004899

set.seed(05297)
finBIC <- lapply(res_750_1_100, sim, n=100, p=750)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.673056816, med 0.767127376
# EST: geo 0.936857903, med 1.036208148

set.seed(94024)
finBIC <- lapply(res_500_1_100, sim, n=100, p=500)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.6654397, med 0.6909629
# EST: geo 1.1851532, med 1.0937967

set.seed(03300)
finBIC <- lapply(res_250_1_100, sim, n=100, p=250)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.599493051, med 0.68380938
# EST: geo 0.734648597, med 1.02317966

set.seed(97255)
finBIC <- lapply(res_100_1_100, sim, n=100, p=100)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.66125061, med 0.711398606
# EST: geo 0.92739776, med 1.000153340

set.seed(90399)
finBIC <- lapply(res_50_1_100, sim, n=100, p=50)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.658235141, med 0.691471457
# EST: geo 0.986352904, med 1.015066597

set.seed(24755)
finBIC <- lapply(res_20_1_100, sim, n=100, p=20)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.691222492, med 0.759835200 
# EST: geo 1.043744036, med 1.006750876

set.seed(49635)
finBIC <- lapply(res_10_1_100, sim, n=100, p=10)
# no errors
getEndResult(finBIC)
# TRUE: geo 0.714272213, med 0.766230052 
# EST: geo 1.035017691, med 1.058864116

