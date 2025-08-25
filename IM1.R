rm(list=ls())

library(splines)
library(brinla)
library(INLA)
library(Metrics)
library(verification)
INLA:::inla.dynload.workaround() 

# January data 
load("C:/Users/mapso/Box/2. Project/171019 - NSF Transmission planning/West Texas Mesonet/Final_code_share/Data_load.RData")

#########################################################
# make a ws.mat data file for INLA
#########################################################
data = as.matrix( dat[,c(1:no.station)] ) 
ws.mat = matrix(c( data ), length(c(data)),1)
ws.mat1 = matrix(rep(1:dim(data)[2], each=dim(data)[1]), length(c(data)),1) #location number
ws.mat = cbind(ws.mat, ws.mat1)
ws.mat = as.data.frame(ws.mat)
ws.mat$time <- rep( seq(1,24), length.out=(dim(data)[1]*dim(data)[2]) )
colnames(ws.mat) = c("ws","location", "time")

t =1:nrow(dat)
freq = 24

## Diurnal pattern: Use trigonometric funnction for diurnal effects
d_s1=as.matrix(sin(2*pi*t/freq))
d_c1=as.matrix(cos(2*pi*t/freq))
d_s2=as.matrix(sin(4*pi*t/freq))
d_c2=as.matrix(cos(4*pi*t/freq))
d_s3=as.matrix(sin(6*pi*t/freq))
d_c3=as.matrix(cos(6*pi*t/freq))
d_s4=as.matrix(sin(8*pi*t/freq))
d_c4=as.matrix(cos(8*pi*t/freq))
d_s5=as.matrix(sin(10*pi*t/freq))
d_c5=as.matrix(cos(10*pi*t/freq))

ws.mat$b_0  = rep(rep(1,max(t)), no.station)
ws.mat$b_1  = rep(d_s1, no.station)
ws.mat$b_2  = rep(d_c1, no.station)
ws.mat$b_3  = rep(d_s2, no.station)
ws.mat$b_4  = rep(d_c2, no.station)
ws.mat$b_5  = rep(d_s3, no.station)
ws.mat$b_6  = rep(d_c3, no.station)
ws.mat$b_7  = rep(d_s4, no.station)
ws.mat$b_8  = rep(d_c4, no.station)
ws.mat$b_9  = rep(d_s5, no.station)
ws.mat$b_10 = rep(d_c5, no.station)

n.hour = length(unique(ws.mat$time)) 
idx.hour = unique(ws.mat$time)
n.time = nrow(data)

# coordinate for INLA dataset
for (i in 1:no.station){
  if (i==1){
    mat.coord = matrix( rep(coord[i,c(3:4)], each=dim(data)[1]), ncol=2)
  }
  if (i>1){
    mat.coord = rbind(mat.coord, matrix( rep(coord[i,c(3:4)], each=dim(data)[1]), ncol=2))
  }
}
mat.coord[,1] = as.numeric(as.character(mat.coord[,1]))
mat.coord[,2] = as.numeric(as.character(mat.coord[,2]))

ws.mat$lat = mat.coord[,1]
ws.mat$lon = mat.coord[,2]
location = cbind(ws.mat$lat, ws.mat$lon)

ws.mat$day = rep(rep(1:31, each = 24), no.station)


####################################################################################################
pred.p1 = matrix(NA, n.time, no.station)
pred.p1.sd = matrix(NA, n.time, no.station)
pred.p1.95ub = matrix(NA, n.time, no.station)
pred.p1.95lb = matrix(NA, n.time, no.station)
pred.p1.90ub = matrix(NA, n.time, no.station)
pred.p1.90lb = matrix(NA, n.time, no.station)
pred.p1.mod  = matrix(NA, n.time, no.station)
pred.p1.fixed= matrix(NA, n.time, no.station)
pred.p1.spatial.random= matrix(NA, n.time, no.station)
pred.p1.daily.random= matrix(NA, n.time, no.station)

for(ii in 1:no.station){
  
  station.te = ii 
  cat("-------Testing station-------\n",c(station.te),"\n------------------------------------------------------------\n")
  
  library(INLA)
  idx.te.loc = which(ws.mat$location == station.te)
  idx.tr.loc = setdiff(c(1:dim(ws.mat)[1]),idx.te.loc)
  
  ws.mat.tr = ws.mat[idx.tr.loc, ]
  ws.mat.te = ws.mat[idx.te.loc, ]
  ws.mat.tr$location.1 = ws.mat.tr$location
  ws.mat.tr$location.2 = ws.mat.tr$location
  ws.mat.tr$location.3 = ws.mat.tr$location
  ws.mat.tr$location.4 = ws.mat.tr$location
  ws.mat.tr$location.5 = ws.mat.tr$location
  ws.mat.tr$location.6 = ws.mat.tr$location
  ws.mat.tr$location.7 = ws.mat.tr$location
  ws.mat.tr$location.8 = ws.mat.tr$location
  ws.mat.tr$location.9 = ws.mat.tr$location
  ws.mat.tr$location.10 = ws.mat.tr$location
  
  ws.mat.tr$day.1 = ws.mat.tr$day
  ws.mat.tr$day.2 = ws.mat.tr$day
  ws.mat.tr$day.3 = ws.mat.tr$day
  ws.mat.tr$day.4 = ws.mat.tr$day
  ws.mat.tr$day.5 = ws.mat.tr$day
  ws.mat.tr$day.6 = ws.mat.tr$day
  ws.mat.tr$day.7 = ws.mat.tr$day
  ws.mat.tr$day.8 = ws.mat.tr$day
  ws.mat.tr$day.9 = ws.mat.tr$day
  ws.mat.tr$day.10 = ws.mat.tr$day
  
  location.tr = cbind( as.numeric(ws.mat[idx.tr.loc,]$lon), as.numeric(ws.mat[idx.tr.loc,]$lat) )
  location.te = cbind( as.numeric(ws.mat[idx.te.loc,]$lon), as.numeric(ws.mat[idx.te.loc,]$lat) )
  
  # m1 <- inla.mesh.2d(loc=location.tr, max.edge=c(0.002, 0.003), offset = c(0.001,0.003))
  m1 <- inla.mesh.2d(loc=location.tr, max.edge=c(1, 3), offset = c(1,3))
  par(mfrow=c(1,1))
  m1 <- inla.mesh.2d(loc=location.tr, max.edge=c(0.1, 0.3), offset = c(0.1,0.3))
  m1$n
  plot(m1, main="")
  points(location.tr[,1], location.tr[,2], pch=19, col=2)
  
  ## create SPDE objects: This produces a stationary field model 
  spde.0  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.1  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.2  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.3  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.4  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.5  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.6  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.7  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.8  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.9  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  spde.10  <- inla.spde2.matern(mesh=m1, alpha=2, constr=TRUE)
  
  A0  <- inla.spde.make.A(m1, loc=location.tr)
  A1  <- inla.spde.make.A(m1, loc=location.tr)
  A2  <- inla.spde.make.A(m1, loc=location.tr)
  A3  <- inla.spde.make.A(m1, loc=location.tr)
  A4  <- inla.spde.make.A(m1, loc=location.tr)
  A5  <- inla.spde.make.A(m1, loc=location.tr)
  A6  <- inla.spde.make.A(m1, loc=location.tr)
  A7  <- inla.spde.make.A(m1, loc=location.tr)
  A8  <- inla.spde.make.A(m1, loc=location.tr)
  A9  <- inla.spde.make.A(m1, loc=location.tr)
  A10  <- inla.spde.make.A(m1, loc=location.tr)
  
  s.index.0  <- inla.spde.make.index(name="spatial.field.0",  n.spde=spde.0$n.spde)
  s.index.1  <- inla.spde.make.index(name="spatial.field.1",  n.spde=spde.1$n.spde)
  s.index.2  <- inla.spde.make.index(name="spatial.field.2",  n.spde=spde.2$n.spde)
  s.index.3  <- inla.spde.make.index(name="spatial.field.3",  n.spde=spde.3$n.spde)
  s.index.4  <- inla.spde.make.index(name="spatial.field.4",  n.spde=spde.4$n.spde)
  s.index.5  <- inla.spde.make.index(name="spatial.field.5",  n.spde=spde.5$n.spde)
  s.index.6  <- inla.spde.make.index(name="spatial.field.6",  n.spde=spde.6$n.spde)
  s.index.7  <- inla.spde.make.index(name="spatial.field.7",  n.spde=spde.7$n.spde)
  s.index.8  <- inla.spde.make.index(name="spatial.field.8",  n.spde=spde.8$n.spde)
  s.index.9  <- inla.spde.make.index(name="spatial.field.9",  n.spde=spde.9$n.spde)
  s.index.10  <- inla.spde.make.index(name="spatial.field.10",  n.spde=spde.10$n.spde)
  
  stack_est <- inla.stack(data=list(y=ws.mat.tr$ws), 
                          A=list(Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_0) %*% A0, 
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_1) %*% A1,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_2) %*% A2,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_3) %*% A3,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_4) %*% A4,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_5) %*% A5,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_6) %*% A6,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_7) %*% A7,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_8) %*% A8,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_9) %*% A9,
                                 Diagonal(length(ws.mat.tr$b_0), ws.mat.tr$b_10) %*% A10, 1,1,1,1,1,1,1,1,1,1,1,1),
                          tag='est',
                          effects=list(s.index.0,
                                       s.index.1,
                                       s.index.2,
                                       s.index.3,
                                       s.index.4,
                                       s.index.5,
                                       s.index.6,
                                       s.index.7,
                                       s.index.8,
                                       s.index.9,
                                       s.index.10,
                                       list(j.0 = ws.mat.tr$day,   weight.0 = ws.mat.tr$b_0), 
                                       list(j.1 = ws.mat.tr$day.1, weight.1 = ws.mat.tr$b_1), 
                                       list(j.2 = ws.mat.tr$day.2, weight.2 = ws.mat.tr$b_2), 
                                       list(j.3 = ws.mat.tr$day.3, weight.3 = ws.mat.tr$b_3), 
                                       list(j.4 = ws.mat.tr$day.4, weight.4 = ws.mat.tr$b_4), 
                                       list(j.5 = ws.mat.tr$day.5, weight.5 = ws.mat.tr$b_5), 
                                       list(j.6 = ws.mat.tr$day.6, weight.6 = ws.mat.tr$b_6), 
                                       list(j.7 = ws.mat.tr$day.7, weight.7 = ws.mat.tr$b_7), 
                                       list(j.8 = ws.mat.tr$day.8, weight.8 = ws.mat.tr$b_8), 
                                       list(j.9 = ws.mat.tr$day.9, weight.9 = ws.mat.tr$b_9), 
                                       list(j.10 = ws.mat.tr$day.10, weight.10 = ws.mat.tr$b_10), 
                                       data.frame(b.0 = ws.mat.tr$b_0,
                                                  b.1 = ws.mat.tr$b_1,
                                                  b.2 = ws.mat.tr$b_2,
                                                  b.3 = ws.mat.tr$b_3,
                                                  b.4 = ws.mat.tr$b_4,
                                                  b.5 = ws.mat.tr$b_5,
                                                  b.6 = ws.mat.tr$b_6,
                                                  b.7 = ws.mat.tr$b_7,
                                                  b.8 = ws.mat.tr$b_8,
                                                  b.9 = ws.mat.tr$b_9,
                                                  b.10 = ws.mat.tr$b_10)))
  
  formula = y ~ 0 + b.0 + b.1 + b.2 + b.3 + b.4 + b.5 + b.6 + b.7 + b.8 + b.9 + b.10 + 
    f(spatial.field.0,  model=spde.0) + f(j.0,  model = "iid") + 
    f(spatial.field.1,  model=spde.1) + f(j.1,  weight.1,  model = "iid") + 
    f(spatial.field.2,  model=spde.2) + f(j.2,  weight.2,  model = "iid") + 
    f(spatial.field.3,  model=spde.3) + f(j.3,  weight.3,  model = "iid") + 
    f(spatial.field.4,  model=spde.4) + f(j.4,  weight.4,  model = "iid") + 
    f(spatial.field.5,  model=spde.5) + f(j.5,  weight.5,  model = "iid") + 
    f(spatial.field.6,  model=spde.6) + f(j.6,  weight.6,  model = "iid") + 
    f(spatial.field.7,  model=spde.7) + f(j.7,  weight.7,  model = "iid") + 
    f(spatial.field.8,  model=spde.8) + f(j.8,  weight.8,  model = "iid") + 
    f(spatial.field.9,  model=spde.9) + f(j.9,  weight.9,  model = "iid") + 
    f(spatial.field.10, model=spde.10)+ f(j.10, weight.10, model = "iid")  
  
  result.diurnal <- inla(formula,
                         data=inla.stack.data(stack_est, 
                                              spde0 = spde.0, 
                                              spde1 = spde.1, 
                                              spde2 = spde.2, 
                                              spde3 = spde.3, 
                                              spde4 = spde.4, 
                                              spde5 = spde.5, 
                                              spde6 = spde.6, 
                                              spde7 = spde.7, 
                                              spde8 = spde.8, 
                                              spde9 = spde.9, 
                                              spde10 = spde.10), 
                         family = "gaussian",
                         control.predictor=list(A=inla.stack.A(stack_est), compute=TRUE),
                         control.compute=list(config = TRUE),
                         control.inla= list(strategy = "gaussian", int.strategy = "eb"),
                         #control.mode=list(theta=theta.ini, restart = TRUE), 
                         num.threads = 1,
                         verbose= TRUE)
  
  theta.ini = result.diurnal$mode$theta
  names(theta.ini)=NULL
  theta.ini
  
  est.idx = inla.stack.index(stack_est, tag = "est")$data
  diurnal.est = result.diurnal$summary.fitted.values$mean[est.idx]
  diurnal.est.mat = matrix(diurnal.est, ncol=(no.station-1))
  
  result.diurnal$summary.fixed
  result.diurnal$summary.fitted.values$mean
  
  library(brinla)
  bri.hyperpar.summary(result.diurnal)
  # for hyperparameters
  bri.hyperpar.plot(result.diurnal)
  # for fixed parameters
  bri.fixed.plot(result.diurnal)
  
  ###########################################################################
  #EXTRACT POSTERIOR SUMMARY STATISTICS
  ###########################################################################
  #-- Extract results for fixed effects - covariate coeffs --#
  beta = result.diurnal$summary.fixed[,"mean"]
  beta_sd = result.diurnal$summary.fixed[,"sd"]
  
  #-- Extract results for sigma2eps (1/precision)
  sigma2e_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for the Gaussian observations")
  sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
  sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
  sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
  sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
  cat("-----Results for sigma2--------------------------------\n",c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),"\n----------------------------------------------------------\n")
  # sigma2
  var.error = sigma2e_m1
  
  ####################################################################################################
  ####################################################################################################
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.0", spde.0)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_0------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = var.nom.m1
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.1", spde.1)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_1------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.2", spde.2)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_2------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.3", spde.3)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_3------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.4", spde.4)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_4------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.5", spde.5)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_5------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.6", spde.6)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_6------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.7", spde.7)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_7------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.8", spde.8)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_8------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.9", spde.9)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_9------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  
  #-- Extract results for the spatial covariance function parameters --#
  result.diurnal.field = inla.spde2.result(result.diurnal, name="spatial.field.10", spde.10)
  var.nom.marg = result.diurnal.field$marginals.variance.nominal[[1]]
  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
  cat("-----Results for tau^2_10------------------------------\n",c(var.nom.m1, var.nom.stdev, var.nom.quantiles),"\n----------------------------------------------------------\n")
  var.space = c(var.space, var.nom.m1)
  ###############################################################################################################
  
  #-- Extract results for delta2 (1/precision): day-to-day variability
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.0")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = sigma2d_m1
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.1")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.2")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.3")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.4")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.5")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.6")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.7")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.8")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x, result.diurnal$marginals.hyperpar$"Precision for j.9")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  sigma2d_marg = inla.tmarginal(function(x) 1/x,  result.diurnal$marginals.hyperpar$"Precision for j.10")
  sigma2d_m1 = inla.emarginal(function(x) x, sigma2d_marg)
  sigma2d_m2 = inla.emarginal(function(x) x^2, sigma2d_marg)
  sigma2d_stdev = sqrt(sigma2d_m2 - sigma2d_m1^2)
  sigma2d_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2d_marg)
  cat("-----Results for delta2day--------------------------------\n",c(sigma2d_m1, sigma2d_stdev, sigma2d_quantiles),"\n----------------------------------------------------------\n")
  var.day = c(var.day, sigma2d_m1)
  
  var.error 
  sum(var.space)
  sum(var.day)
  
  var.error / (var.error +sum(var.space)+sum(var.day))*100
  sum(var.space) / (var.error +sum(var.space)+sum(var.day))*100
  sum(var.day) / (var.error +sum(var.space)+sum(var.day))*100
  
  ####################################################################################################################################################
  ####################################################################################################################################################
  # range parameter: unit - km
  range.nom.marg = result.diurnal.field$marginals.range.nominal[[1]]
  range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
  range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
  range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
  range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
  cat("-----Results for the range (km) ------------------------------\n", c(range.nom.m1, range.nom.stdev, range.nom.quantiles), "\n----------------------------------------------------------\n")
  
  ####################################################################################################################################################
  A0p  <- inla.spde.make.A(m1, loc=location.te)
  A1p  <- inla.spde.make.A(m1, loc=location.te)
  A2p  <- inla.spde.make.A(m1, loc=location.te)
  A3p  <- inla.spde.make.A(m1, loc=location.te)
  A4p  <- inla.spde.make.A(m1, loc=location.te)
  A5p  <- inla.spde.make.A(m1, loc=location.te)
  A6p  <- inla.spde.make.A(m1, loc=location.te)
  A7p  <- inla.spde.make.A(m1, loc=location.te)
  A8p  <- inla.spde.make.A(m1, loc=location.te)
  A9p  <- inla.spde.make.A(m1, loc=location.te)
  A10p  <- inla.spde.make.A(m1, loc=location.te)
  
  ####################################################################################################################################################
  # sampling (MCMC) starts for posterior predictive distribution
  N = 10^4 /2
  samples = inla.posterior.sample(N, result.diurnal) 
  contents = result.diurnal$misc$configs$contents
  
  effect = "spatial.field.0"
  id.effect = which(contents$tag==effect)
  ind.effect.s0 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.1"
  id.effect = which(contents$tag==effect)
  ind.effect.s1 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.2"
  id.effect = which(contents$tag==effect)
  ind.effect.s2 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.3"
  id.effect = which(contents$tag==effect)
  ind.effect.s3 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.4"
  id.effect = which(contents$tag==effect)
  ind.effect.s4 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.5"
  id.effect = which(contents$tag==effect)
  ind.effect.s5 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.6"
  id.effect = which(contents$tag==effect)
  ind.effect.s6 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.7"
  id.effect = which(contents$tag==effect)
  ind.effect.s7 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.8"
  id.effect = which(contents$tag==effect)
  ind.effect.s8 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.9"
  id.effect = which(contents$tag==effect)
  ind.effect.s9 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "spatial.field.10"
  id.effect = which(contents$tag==effect)
  ind.effect.s10 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.0"
  id.effect = which(contents$tag==effect)
  ind.effect.j0 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.1"
  id.effect = which(contents$tag==effect)
  ind.effect.j1 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.2"
  id.effect = which(contents$tag==effect)
  ind.effect.j2 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.3"
  id.effect = which(contents$tag==effect)
  ind.effect.j3 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.4"
  id.effect = which(contents$tag==effect)
  ind.effect.j4 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.5"
  id.effect = which(contents$tag==effect)
  ind.effect.j5 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.6"
  id.effect = which(contents$tag==effect)
  ind.effect.j6 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.7"
  id.effect = which(contents$tag==effect)
  ind.effect.j7 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.8"
  id.effect = which(contents$tag==effect)
  ind.effect.j8 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.9"
  id.effect = which(contents$tag==effect)
  ind.effect.j9 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  effect = "j.10"
  id.effect = which(contents$tag==effect)
  ind.effect.j10 = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  # look at one sample
  samples[[1]]$latent[ind.effect.s0, , drop=F]
  samples[[1]]$hyperpar
  tail(samples[[1]]$latent, n=50)
  
  # fixed effect #########
  b.0 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+1]))
  b.1 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+2]))
  b.2 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+3]))
  b.3 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+4]))
  b.4 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+5]))
  b.5 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+6]))
  b.6 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+7]))
  b.7 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+8]))
  b.8 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+9]))
  b.9 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+10]))
  b.10 <- unlist(lapply(samples, function(x) x$latent[max(ind.effect.j10)+11]))
  
  # spatial random effect #########
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s0])
  s.eff.0 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.0) = rownames(samples[[1]]$latent)[ind.effect.s0]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s1])
  s.eff.1 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.1) = rownames(samples[[1]]$latent)[ind.effect.s1]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s2])
  s.eff.2 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.2) = rownames(samples[[1]]$latent)[ind.effect.s2]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s3])
  s.eff.3 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.3) = rownames(samples[[1]]$latent)[ind.effect.s3]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s4])
  s.eff.4 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.4) = rownames(samples[[1]]$latent)[ind.effect.s4]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s5])
  s.eff.5 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.5) = rownames(samples[[1]]$latent)[ind.effect.s5]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s6])
  s.eff.6 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.6) = rownames(samples[[1]]$latent)[ind.effect.s6]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s7])
  s.eff.7 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.7) = rownames(samples[[1]]$latent)[ind.effect.s7]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s8])
  s.eff.8 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.8) = rownames(samples[[1]]$latent)[ind.effect.s8]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s9])
  s.eff.9 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.9) = rownames(samples[[1]]$latent)[ind.effect.s9]
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect.s10])
  s.eff.10 = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff.10) = rownames(samples[[1]]$latent)[ind.effect.s10]
  
  # spatial random effect for unobserved location #########
  library(Matrix)
  s.eff.0 = Matrix(s.eff.0, sparse=TRUE)
  b.s0.hat = Matrix(crossprod(t(A0p), t(s.eff.0)), sparse = TRUE)
  
  s.eff.1 = Matrix(s.eff.1, sparse=TRUE)
  b.s1.hat = Matrix(crossprod(t(A1p), t(s.eff.1)), sparse = TRUE)
  
  s.eff.2 = Matrix(s.eff.2, sparse=TRUE)
  b.s2.hat = Matrix(crossprod(t(A2p), t(s.eff.2)), sparse = TRUE)
  
  s.eff.3 = Matrix(s.eff.3, sparse=TRUE)
  b.s3.hat = Matrix(crossprod(t(A3p), t(s.eff.3)), sparse = TRUE)
  
  s.eff.4 = Matrix(s.eff.4, sparse=TRUE)
  b.s4.hat = Matrix(crossprod(t(A4p), t(s.eff.4)), sparse = TRUE)
  
  s.eff.5 = Matrix(s.eff.5, sparse=TRUE)
  b.s5.hat = Matrix(crossprod(t(A5p), t(s.eff.5)), sparse = TRUE)
  
  s.eff.6 = Matrix(s.eff.6, sparse=TRUE)
  b.s6.hat = Matrix(crossprod(t(A6p), t(s.eff.6)), sparse = TRUE)
  
  s.eff.7 = Matrix(s.eff.7, sparse=TRUE)
  b.s7.hat = Matrix(crossprod(t(A7p), t(s.eff.7)), sparse = TRUE)
  
  s.eff.8 = Matrix(s.eff.8, sparse=TRUE)
  b.s8.hat = Matrix(crossprod(t(A8p), t(s.eff.8)), sparse = TRUE)
  
  s.eff.9 = Matrix(s.eff.9, sparse=TRUE)
  b.s9.hat = Matrix(crossprod(t(A9p), t(s.eff.9)), sparse = TRUE)
  
  s.eff.10 = Matrix(s.eff.10, sparse=TRUE)
  b.s10.hat = Matrix(crossprod(t(A10p), t(s.eff.10)), sparse = TRUE)
  
  # day-to-day random effect #########
  samples.effect.j0 = lapply(samples, function(x) x$latent[ind.effect.j0])
  samples.effect.j1 = lapply(samples, function(x) x$latent[ind.effect.j1])
  samples.effect.j2 = lapply(samples, function(x) x$latent[ind.effect.j2])
  samples.effect.j3 = lapply(samples, function(x) x$latent[ind.effect.j3])
  samples.effect.j4 = lapply(samples, function(x) x$latent[ind.effect.j4])
  samples.effect.j5 = lapply(samples, function(x) x$latent[ind.effect.j5])
  samples.effect.j6 = lapply(samples, function(x) x$latent[ind.effect.j6])
  samples.effect.j7 = lapply(samples, function(x) x$latent[ind.effect.j7])
  samples.effect.j8 = lapply(samples, function(x) x$latent[ind.effect.j8])
  samples.effect.j9 = lapply(samples, function(x) x$latent[ind.effect.j9])
  samples.effect.j10 = lapply(samples, function(x) x$latent[ind.effect.j10])
  
  d.s0.hat = matrix(NA, n.time, N)
  d.s1.hat = matrix(NA, n.time, N)
  d.s2.hat = matrix(NA, n.time, N)
  d.s3.hat = matrix(NA, n.time, N)
  d.s4.hat = matrix(NA, n.time, N)
  d.s5.hat = matrix(NA, n.time, N)
  d.s6.hat = matrix(NA, n.time, N)
  d.s7.hat = matrix(NA, n.time, N)
  d.s8.hat = matrix(NA, n.time, N)
  d.s9.hat = matrix(NA, n.time, N)
  d.s10.hat = matrix(NA, n.time, N)
  
  for(n in 1:N){
    d.s0.hat[,n] = rep(samples.effect.j0[[n]], each=24)
    d.s1.hat[,n] = rep(samples.effect.j1[[n]], each=24)
    d.s2.hat[,n] = rep(samples.effect.j2[[n]], each=24)
    d.s3.hat[,n] = rep(samples.effect.j3[[n]], each=24)
    d.s4.hat[,n] = rep(samples.effect.j4[[n]], each=24)
    d.s5.hat[,n] = rep(samples.effect.j5[[n]], each=24)
    d.s6.hat[,n] = rep(samples.effect.j6[[n]], each=24)
    d.s7.hat[,n] = rep(samples.effect.j7[[n]], each=24)
    d.s8.hat[,n] = rep(samples.effect.j8[[n]], each=24)
    d.s9.hat[,n] = rep(samples.effect.j9[[n]], each=24)
    d.s10.hat[,n] = rep(samples.effect.j10[[n]], each=24)
  }
  
  # standard deviation for the Gaussian observations
  sigma <- unlist(lapply(samples, function(x) 1/sqrt(x$hyperpar["Precision for the Gaussian observations"])))
  eps = rnorm(N, sd=sigma)
  
  yhat.sim = matrix(NA, length(ws.mat.te$ws), N)
  yhat.sim.fixed = matrix(NA, length(ws.mat.te$ws), N)
  yhat.sim.spatial.random = matrix(NA, length(ws.mat.te$ws), N)
  yhat.sim.daily.random = matrix(NA, length(ws.mat.te$ws), N)
  
  for(k in 1:N){
    if ( k == 1000 || k == 2000 || k == 3000 || k == 4000 || k == 5000 ){ print(k) }
    
    # fixed effect
    yhat.sim.fixed[,k] =  b.0[k]*ws.mat.te$b_0 + 
      b.1[k]*ws.mat.te$b_1 + 
      b.2[k]*ws.mat.te$b_2 + 
      b.3[k]*ws.mat.te$b_3 + 
      b.4[k]*ws.mat.te$b_4 + 
      b.5[k]*ws.mat.te$b_5 +
      b.6[k]*ws.mat.te$b_6 + 
      b.7[k]*ws.mat.te$b_7 + 
      b.8[k]*ws.mat.te$b_8 + 
      b.9[k]*ws.mat.te$b_9 + 
      b.10[k]*ws.mat.te$b_10
    
    # spatial random effect
    yhat.sim.spatial.random[,k] = b.s0.hat[,k]*ws.mat.te$b_0 + 
      b.s1.hat[,k]*ws.mat.te$b_1 + 
      b.s2.hat[,k]*ws.mat.te$b_2 + 
      b.s3.hat[,k]*ws.mat.te$b_3 + 
      b.s4.hat[,k]*ws.mat.te$b_4 + 
      b.s5.hat[,k]*ws.mat.te$b_5 +
      b.s6.hat[,k]*ws.mat.te$b_6 + 
      b.s7.hat[,k]*ws.mat.te$b_7 + 
      b.s8.hat[,k]*ws.mat.te$b_8 + 
      b.s9.hat[,k]*ws.mat.te$b_9 + 
      b.s10.hat[,k]*ws.mat.te$b_10
    
    # daily random effect
    yhat.sim.daily.random[,k] = d.s0.hat[,k]*ws.mat.te$b_0 + 
      d.s1.hat[,k]*ws.mat.te$b_1 + 
      d.s2.hat[,k]*ws.mat.te$b_2 + 
      d.s3.hat[,k]*ws.mat.te$b_3 + 
      d.s4.hat[,k]*ws.mat.te$b_4 + 
      d.s5.hat[,k]*ws.mat.te$b_5 +
      d.s6.hat[,k]*ws.mat.te$b_6 + 
      d.s7.hat[,k]*ws.mat.te$b_7 + 
      d.s8.hat[,k]*ws.mat.te$b_8 + 
      d.s9.hat[,k]*ws.mat.te$b_9 + 
      d.s10.hat[,k]*ws.mat.te$b_10
    
    yhat.sim[,k] = yhat.sim.fixed[,k] + yhat.sim.spatial.random[,k] + yhat.sim.daily.random[,k] + eps[k] 
  }
  
  post.pred.yhat = matrix(NA, nrow(yhat.sim), 7)
  post.pred.yhat[, 1] = rowMeans(yhat.sim, na.rm = TRUE)
  
  for(l in 1:nrow(yhat.sim)){
    post.pred.yhat[l, 2] = sd(yhat.sim[l,], na.rm=TRUE)
    post.pred.yhat[l, 3:7] = quantile(yhat.sim[l,], probs = c(0.025, 0.05, 0.5, 0.95, 0.975), na.rm = TRUE)
  }
  
  colnames(post.pred.yhat) = c("mean", "sd", "2.5%","0.5%", "50%", "95%", "97.5%")
  
  post.pred.fixed = rowMeans(yhat.sim.fixed, na.rm = TRUE)
  post.pred.spatial.random = rowMeans(yhat.sim.spatial.random, na.rm = TRUE)
  post.pred.daily.random = rowMeans(yhat.sim.daily.random, na.rm = TRUE)
  
  ## 
  pred.p1[,ii] = post.pred.yhat[,1]
  pred.p1.sd[,ii] = post.pred.yhat[,2]
  pred.p1.95lb[,ii] = post.pred.yhat[,3]
  pred.p1.90lb[,ii] = post.pred.yhat[,4]
  pred.p1.mod[,ii]  = post.pred.yhat[,5]
  pred.p1.90ub[,ii] = post.pred.yhat[,6]
  pred.p1.95ub[,ii] = post.pred.yhat[,7]
  pred.p1.fixed[,ii] = post.pred.fixed
  pred.p1.spatial.random[,ii] = post.pred.spatial.random
  pred.p1.daily.random[,ii] = post.pred.daily.random
  
  plot(dat[,1])
  lines(pred.p1[,1], col=2)
  
  rm(samples)
}

plot(dat[c(1:168),1]*10, ylim=c(0,20), xlab="Time index (1:168)", ylab="Wind speed (m/s)", main="Prediction: The first week of April")
lines(pred.p1[c(1:168),1]*10, col=2, lwd=2)
lines(pred.p1.95ub[c(1:168),1]*10, col=2, lty=2, lwd=1)
lines(pred.p1.95lb[c(1:168),1]*10, col=2, lty=2, lwd=1)
