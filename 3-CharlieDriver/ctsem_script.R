####to update ctsem via github, uncomment the following lines (some small upgrades since CRAN release)
# install.packages("devtools")
# library(devtools)
# install_github("cdriveraus/ctsem")
# install.packages('bridgesampling')

library(ctsem)
set.seed(3)

#recommended to adjust these if laptop does not have 4 or more cores (min 3 chains recommended, but not necessary here)
setcores <- 7
setchains <- 7

#### Data generation (this needs to be run, but not necessary to understand!)
Tpoints <- 20
nmanifest <- 4
nlatent <- 2
nsubjects <- 20

#random effects
age <- rnorm(nsubjects) #standardised
cint1 <- rnorm(nsubjects, 2, .3) + age * .5
cint2 <- cint1 * .5 + rnorm(nsubjects, 1, .2) + age * .5
tdpredeffect <- rnorm(nsubjects, 5, .3) + age * .5

for (i in 1:nsubjects) {
  #generating model
  gm <- ctModel(Tpoints = Tpoints, n.manifest = nmanifest, n.latent = nlatent,
                n.TDpred = 1,
    LAMBDA = matrix(c(1, 0, 0, 0, 0, 1, .8, 1.3),
                    nrow = nmanifest, ncol = nlatent),
    DRIFT = matrix(c(-.3, .1, 0, -.5), nlatent, nlatent),
    TDPREDEFFECT = matrix(c(tdpredeffect[i], 0), nrow = nlatent),
    TDPREDMEANS = matrix(c(rep(0, Tpoints - 10), 1, rep(0, 9)), ncol = 1),
    DIFFUSION = matrix(c(.5, 0, 0, .5), 2, 2),
    CINT = matrix(c(cint1[i], cint2[i]), ncol = 1),
    T0VAR = diag(2, nlatent, nlatent),
    MANIFESTVAR = diag(.5, nmanifest))

  #generate data
  newdat <- ctGenerate(ctmodelobj = gm, n.subjects = 1, burnin = 2,
                       dtmat <- rbind(c(rep(.5, 8), 3, rep(.5, Tpoints - 9))),
                       wide = FALSE)
  newdat[,'id'] <- i #set id for each subject
  newdat <- cbind(newdat, age[i]) #include time independent predictor
  if(i == 1) {
    dat <- newdat[1:(Tpoints - 10), ] #pre intervention data
    dat2 <- newdat #including post intervention data
  }
  if(i > 1) {
    dat <- rbind(dat, newdat[1:(Tpoints - 10), ])
    dat2 <- rbind(dat2, newdat)
  }
}
colnames(dat)[ncol(dat)] <- 'age'
colnames(dat2)[ncol(dat)] <- 'age'

#plot generated data for sanity
plot(age)
matplot(dat[, gm$manifestNames], type = 'l', pch = 1)
plotvar <- 'Y1'
plot(dat[dat[, 'id'] == 1, 'time'], dat[dat[, 'id'] == 1, plotvar], type = 'l',
     ylim = range(dat[, plotvar], na.rm = TRUE))
for(i in 2:nsubjects){
  points(dat[dat[, 'id'] == i, 'time'], dat[dat[, 'id'] == i, plotvar],
         type = 'l', col = i)
}

#### Model fitting (from here it is good to understand!)

#Specify univariate linear growth curve
#page 5 of https://cran.r-project.org/web/packages/ctsem/vignettes/hierarchical.pdf documents these arguments (or use ?ctModel )
m1 <- ctModel(
  n.manifest = 1, n.latent = 1, n.TIpred = 1, type = 'stanct',
  manifestNames = c('Y1'), latentNames = c('L1'), TIpredNames = 'age',
  DRIFT = matrix(-1e-5, nrow = 1, ncol = 1),
  DIFFUSION = matrix(1e-5, nrow = 1, ncol = 1),
  CINT = matrix(c('cint1'), ncol = 1),
  T0MEANS = matrix(c('t0m1'), ncol = 1),
  T0VAR = matrix(0, nrow = 1, ncol = 1),
  LAMBDA = diag(1),
  MANIFESTMEANS = matrix(0, ncol = 1),
  MANIFESTVAR = matrix(c('merror1'), nrow = 1, ncol = 1))

#modify between subject aspects -- alternatively, run: edit(m1$pars)
m1$pars$indvarying[-which(m1$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
m1$pars$age_effect[-which(m1$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE

plot(m1) #plot prior distributions

#fit
f1 <- ctStanFit(datalong = dat, ctstanmodel = m1, cores = setcores,
                chains = setchains, plot = TRUE,
                control = list(max_treedepth = 10), iter = 150)
summary(f1)
plot(f1)

#plots of individual subject models v data
ctKalman(f1, timestep = .01, plot = TRUE, subjects = 1:4,
         kalmanvec = c('y', 'etasmooth'))
ctKalman(f1, timestep = .01, plot = TRUE, subjects = 1,
         kalmanvec = c('y', 'ysmooth'))

ctStanPlotPost(f1) #compare prior to posterior distributions
ctStanPlotPost(f1, priorwidth = FALSE) #rescale to width of posterior

ctStanPostPredict(f1) #compare randomly generated data from posterior to observed data

cf <- ctCheckFit(f1) #compare covariance of randomly generated data to observed cov
plot(cf)

#accessing the stan object directly
library(rstan)
postsamples <- extract(f1$stanfit, pars = 'Ygen') #extract data generated from posterior
plot(f1$data$time, postsamples$Ygen[1, 1, , 1]) #1st iteration (though randomly shuffled by extract), 1st chain (merged by default), all measurement occasions, 1st manifest variable
points(f1$data$time, f1$data$Y[, 1], col = 'red') #1st manifest variable



#Specify model including dynamics
m2 <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
  manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
  DRIFT=matrix('drift11',nrow=1,ncol=1),
  DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
  CINT=matrix(c('cint1'),ncol=1),
  T0MEANS=matrix(c('t0m1'),ncol=1),
  T0VAR=matrix('t0var11',nrow=1,ncol=1),
  LAMBDA = diag(1),
  MANIFESTMEANS=matrix(0,ncol=1),
  MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))

m2$pars$indvarying[-which(m2$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
m2$pars$age_effect[-which(m2$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE


f2 <- ctStanFit(datalong = dat, ctstanmodel = m2, cores = setcores,chains = setchains,plot=TRUE,
  control=list(max_treedepth=7),iter=150)

summary(f2,parmatrices=TRUE,timeinterval=1)

ctKalman(f2,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
ctKalman(f2,timestep=.01,plot=TRUE,subjects=1:4,kalmanvec=c('y','etasmooth'))
ctKalman(f2,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))

ctStanPlotPost(f2)

ctStanPostPredict(f2)

#model comparison
library(bridgesampling)
ml1<-bridge_sampler(samples = f1$stanfit)
ml2<-bridge_sampler(samples = f2$stanfit)
bf(ml1, ml2)







#Include intervention
m3 <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
  manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
  n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
  TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #and sets intervention effect matrix
  DRIFT=matrix('drift11',nrow=1,ncol=1),
  DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
  CINT=matrix(c('cint1'),ncol=1),
  T0MEANS=matrix(c('t0m1'),ncol=1),
  T0VAR=matrix('t0var11',nrow=1,ncol=1),
  LAMBDA = diag(1),
  MANIFESTMEANS=matrix(0,ncol=1),
  MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))

m3$pars$indvarying[-which(m3$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
m3$pars$age_effect[-which(m3$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE

f3 <- ctStanFit(datalong = dat2, ctstanmodel = m3, cores = setcores,chains = setchains,plot=TRUE,
  control=list(max_treedepth=7),iter=150)

summary(f3,parmatrices=TRUE)

ctKalman(f3,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
ctKalman(f3,timestep=.01,plot=TRUE,subjects=1:4,kalmanvec=c('y','etasmooth'))
ctKalman(f3,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))

ctStanPlotPost(f3)

ctStanPostPredict(f3, datarows=0:100)






#include 2nd latent process

#use either full explicit specification
m4 <- ctModel(n.manifest = 2,n.latent = 2,n.TIpred = 1, type = 'stanct', #no of vars updated
  manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),TIpredNames = 'age', #new variable names added on this line
  n.TDpred=1,TDpredNames = 'TD1',
  TDPREDEFFECT=matrix(c('tdpredeffect1','tdpredeffect2'),nrow=2,ncol=1),
  DRIFT=matrix(c('drift11','drift21','drift12','drift22'),nrow=2,ncol=2),
  DIFFUSION=matrix(c('diffusion11','diffusion21',0,'diffusion22'),nrow=2,ncol=2),
  CINT=matrix(c('cint1','cint2'),nrow=2,ncol=1),
  T0MEANS=matrix(c('t0m1','t0m2'),nrow=2,ncol=1),
  T0VAR=matrix(c('t0var11','t0var21',0,'t0var22'),nrow=2,ncol=2),
  LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2),
  MANIFESTMEANS=matrix(c(0,0),nrow=2,ncol=1),
  MANIFESTVAR=matrix(c('merror1',0,0,'merror2'),nrow=2,ncol=2))

#restrict between subjects variation / covariate effects
m4$pars$indvarying[-which(m4$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
m4$pars$age_effect[-which(m4$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE

#or rely on defaults (MANIFESTMEANS now free instead of CINT -- no substantive difference for one indicator factors)
m4 <- ctModel(n.manifest = 2,n.latent = 2,n.TIpred = 1, type = 'stanct', #no of vars updated
  manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),TIpredNames = 'age', #new variable names added on this line
  n.TDpred=1,TDpredNames = 'TD1',
  LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2))

#restrict between subjects variation / covariate effects
m4$pars$indvarying[-which(m4$pars$matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT'))] <- FALSE
m4$pars$age_effect[-which(m4$pars$matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT'))] <- FALSE

f4 <- ctStanFit(datalong = dat2, ctstanmodel = m4, cores = setcores,chains = setchains,plot=TRUE,
  control=list(max_treedepth=7),iter=150)

summary(f4,parmatrices=TRUE)

ctStanDiscretePars(f4,plot=TRUE) #auto and cross regressive plots over time

ctKalman(f4,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
ctKalman(f4,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','etasmooth'))
ctKalman(f4,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))

ctStanPlotPost(f4)

ctStanPostPredict(f4, datarows=0:100)



#non-linear dedpendencies - based on m3 model (including intervention)
m3nl <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
  manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
  n.TDpred=1,TDpredNames = 'TD1',
  TDPREDEFFECT=matrix(c('PARS[1,1] * state[1]'),nrow=1,ncol=1), #specify intervention as dependent on extra parameter in PARS matrix, and latent process 1
  PARS=matrix(c('tdpredeffect'),1,1),
  DRIFT=matrix('drift11',nrow=1,ncol=1),
  DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
  CINT=matrix(c('cint1'),ncol=1),
  T0MEANS=matrix(c('t0m1'),ncol=1),
  T0VAR=matrix('t0var11',nrow=1,ncol=1),
  LAMBDA = diag(1),
  MANIFESTMEANS=matrix(0,ncol=1),
  MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))

m3nl$pars$indvarying[-which(m3nl$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
m3nl$pars$age_effect[-which(m3nl$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE

f3nl <- ctStanFit(datalong = dat2, ctstanmodel = m3nl, cores = setcores, chains = setchains,plot=TRUE,
  control=list(max_treedepth=7),iter=150)

summary(f3nl)

#ctKalman and ctStanDiscr plotting functions still need updating for non-linearities! (as of ctsem v 2.6.5)

ctStanPostPredict(f3nl, datarows=1:100)
