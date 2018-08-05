#---- Fitting a univariate dual change score model in dynr---- 

options(scipen=999) #Suppresses scientific notation
library(dynr)

#Read in data and set up dynr data structure
thedat = read.table('./Data/LDS.dat') 
np=500
T = 10

thedat$ID = rep(1:np,each=T)
thedat$Time = rep(1:T,np)
View(thedat)
data <- dynr.data(thedat, id="ID", time="Time", 
                  observed=c("V1"))

#Measurement model
meas <- prep.measurement(
    values.load=matrix(c(1,0),ncol=2,byrow=TRUE),
    params.load=matrix(rep("fixed",2),ncol=2),
    state.names=c("readLevel","readSlope"),
    obs.names=c("V1")
  )

#Initial mean and covariance structure
initial <- prep.initial(
    values.inistate=c(-1,.5),
    params.inistate=c('mu_readLevel', 
                      'mu_readSlope'), 
    values.inicov=matrix(c(.2,.02,
                           .02,.1),byrow=TRUE,ncol=2),
    params.inicov=matrix(c("v_11","c_12",
                           "c_12","v_22"),byrow=TRUE,ncol=2)) #initial covariance matrix is freely estimated.

  #Process (latent) and measurement error covariance matrices  
  mdcov <- prep.noise(
    values.latent= diag(rep(0,2)), 
    params.latent=diag(rep("fixed",2)), 
    values.observed=.5, 
    params.observed='readErrorV')
  
  #Formulas for dynamic model
  formula =
    list(readLevel~ (1+beta.read)*readLevel + readSlope,
         readSlope~ readSlope
    )
  
  #Dynamic model
  dynm  <- prep.formulaDynamics(formula=formula,
                                startval=c(beta.read = 1.62
                                ), isContinuousTime=FALSE) 
  
  
  #Put all recipes together into a dynr model
  model <- dynr.model(dynamics=dynm, measurement=meas,
                      noise=mdcov, initial=initial, data=data,
                      outfile="uDCS.c")
  model$'param.names'
  
  uDCS <- dynr.cook(model)
  summary(uDCS)
  

#Plot the key modeling equations for model inspection
p3 <-plotFormula(model,ParameterAs=model$param.names)
print(p3)


