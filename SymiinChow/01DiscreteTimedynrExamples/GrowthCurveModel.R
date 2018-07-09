# Fitting a growth curve model in dynr
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This demo illustrates how to fit a linear
# growth curve model with random intercept and
# slope using dynr and extract intercept and 
# slope estimates.
#
# Author: Sy-Miin Chow
# Last modified: 7/8/18
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ---- dynr code starts here ---- 
library(dynr)
d = read.table("../Data/GrowthCurveExample.csv",
               sep=",", header = TRUE)
View(d)

# Set up data. Specify names of variables to be used for modeling
data <- dynr.data(d, id="subject", time="time", 
                  observed=c("y1"), covariates="Deltat")

# Measurement model
# The values provided in values.XXX are starting values for the corresponding
# parameters. If an entry is specified as "fixed" under params.XXX, then it is
# fixed at the value specified.
meas <- prep.measurement(
  values.load=matrix(c(1,0),ncol=2,byrow=TRUE), 
  params.load=matrix(rep("fixed",2),ncol=2),
  state.names=c("Level","Slope"),
  obs.names=c("y1")
)

# Formulas for dynamic model
formula2 =list(
  Level~ Level + Deltat*Slope,
  Slope~ Slope
)

# Dynamic model
dynm  <- prep.formulaDynamics(formula=formula2,
                              isContinuousTime=FALSE)

#Here I am giving an example of how to specify the dynamic model using prep.matrixDynamics.
#Note however that in the current version of dynr, you cannot incorporate covariate Deltat
#into the transition matrix (values.dyn and params.dyn) so the specification below
#assumes that you successive measurement occasions are equally spaced.The alternative
#prep.formulaDynamics is more flexible and allows
#for the inclusion of covariates in the dynamic formulas.
dynm2 <- prep.matrixDynamics(
  values.dyn=matrix(c(1,1,
                      0,1),byrow=TRUE,ncol=2),
  params.dyn=matrix(c('fixed','fixed',
                      'fixed','fixed'),byrow=TRUE,ncol=2),
  isContinuousTime=FALSE)

initial <- prep.initial(
  values.inistate=matrix(c(-1,
                           .5),ncol=1,byrow=TRUE),
  params.inistate=c('mu_Level0',
                    'mu_Slope0'), 
  values.inicov=matrix(c(5,1,
                         1,5),byrow=TRUE,ncol=2),
  params.inicov=matrix(c("v_11","c_12",
                         "c_12","v_22"),
                       byrow=TRUE,ncol=2)) 

#A hypothetical example of how to incorporate a predictor, u1, of
#the initial condition of the latent variables.
#initial2 <- prep.initial(
#  values.inistate=matrix(
#    c(-1, .1,
#      .5,  .1), byrow=TRUE,ncol=2), #nrow = numLatentState, ncol = numCovariates + 1
#  params.inistate=matrix(
#    c('mu_Level0', 'mu_Level1',
#      'mu_Slope0', 'mu_Slope1'), byrow=TRUE,
#    ncol=2), #Level has free intercept, mu_Level0, and free u1 effect, mu_Level1, 
#             #Slope has free intercept, mu_Slope0, and free u1 effect mu_Slope1.
#  values.inicov=matrix(c(5,1,
#                         1,5),byrow=TRUE,ncol=2),
#  params.inicov=matrix(c("v_11","c_12",
#                         "c_12","v_22"),
#                       byrow=TRUE,ncol=2),  
#  covariates='u1')


# Process noise and measurement error cov matrices
mdcov <- prep.noise(
  values.latent=diag(rep(0,2)), 
  params.latent=diag(rep("fixed",2)), 
  values.observed=.5, 
  params.observed='ErrorV')

# Put all the recipes together into a dynr model
GCmodel <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data,#transform=trans,
                    outfile="GrowthCurve.c")

# Create a LaTeX file showing all the equations. Don't run if you don't
# already use LaTeX on your machine and has all the dependencies set up.
# Use plotFormula instead (see below).
printex(GCmodel, ParameterAs=GCmodel$'param.names', show=FALSE, printInit=TRUE, 
        printRS=FALSE,
        outFile="GrowthCurve.tex")
tools::texi2pdf("GrowthCurve.tex")
system(paste(getOption("pdfviewer"), "GrowthCurve.pdf"))

plotFormula(GCmodel, ParameterAs=GCmodel$'param.names')

GrowthCurveResults <- dynr.cook(GCmodel,debug_flag=TRUE)
#debug_flag=FALSE gives you limited by-products
#Tip: if you want to perform Kalman filtering and smoothing at the 
#parameter values provided without estimating the parameters
#then set optimization_flag = FALSE, hessian_flag = FALSE.

# true values: readErrorV=0.25; mu_level=0.5; mu_slope=0.25; v_11=10, c_12=1, v_22=2
summary(GrowthCurveResults)
plotFormula(GCmodel, ParameterAs=coef(GrowthCurveResults))

dynr.ggplot(GrowthCurveResults,GCmodel,style=2)
plot(GrowthCurveResults,GCmodel)


#You can also extract other by-products from the cooked object
#by using the ``@'' symbol. For instance, to get smoothed latent variables
#estimates of the latent level and slope (factor score estimates), we can 
#extract eta_smooth_final as:
GrowthCurveResults@eta_smooth_final 
#These smoothed latent variable estimates are
#equivalent to factor scores from the regression approach

#Filtered latent variable estimates
#Under some conditions, equivalent to factor scores from the Barlett approach
GrowthCurveResults@eta_filtered

LV = GrowthCurveResults@eta_smooth_final
dim(LV) #This is of dimension number of latent variables x (total # of subjects and time points)

LV = t(LV) #Transposing the LV data set and put things into a data frame.
LV = data.frame(ID = d$subject, time = d$time, levelEst = LV[,1], slopeEst = LV[,2])
boxplot(LV$slopeEst) #Create boxplot of the person-specific linear slope estimates
