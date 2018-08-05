#To simulate data for dynamic factor analysis model with auto- and cross-regressions for 2 factors, 2lags (0 and 1 lag) each
rm(list=ls(all=TRUE))
nt=1 	# number of time points
ne=2   	# number of states
ny=6	# number of observed 
nx=0    # number of fixed regressors
np=500	# number of subjects
filey=paste0('CFA.dat')  # output file for obs y
npad=1 # start up
ist=npad+1   
ntt=nt+npad
# S
S=matrix(c(
  1,0,
  1.2,0,
  .8,0,
  0,1,
  0,.9,
  0,1.1
),ny,ne,byrow=TRUE)
# Q
Q=matrix(c(
  2.5,.6,
  .6,2.5
),ne,ne,byrow=TRUE)
# H
#H=matrix(c(
#  .5,-.3,
#  -.2,.4),ne,ne,byrow=TRUE)
H = matrix(rep(0,4),ncol=ne)

# R
R=diag(c(.8,.6,2,1,1.5,2))
# c
c=matrix(c(0,0),ne,1,byrow=TRUE)
# d
d=matrix(c(3,2,4,5,3,4),ny,1,byrow=TRUE)
## Z
# states a t=0
a0=matrix(c(0,0),ne,1,byrow=TRUE)
# cholesky of Q & R (assumes these are positive definite)
Qs = Q
Rs = R
if (sum(diag(Q))> 0) Qs = chol(Q)
if (sum(diag(R))> 0) Rs = chol(R)
# innov z residuals e
a=matrix(0,ntt,ne)
y=matrix(0,ntt,ny)
x=matrix(0,ntt,nx)
yall=matrix(0,nt*np,ny)
all = matrix(0,nt*np,ne)


for (j in 1:np){
  a[1,1:ne] = a0
  for (i in 2:ntt) 
  {
    ztmp=t(rnorm(ne)%*%Qs)
    etmp=t(rnorm(ny)%*%Rs)
    
    atmp=as.matrix(a[i-1,1:ne])
    atmp=H%*%atmp+ztmp+c
    a[i,1:ne]=t(atmp)
    ytmp=S%*%atmp+etmp+d
    y[i,1:ny]=ytmp
  }
  yall[ (1+(j-1)*nt):(j*nt),1:ny] = y[(ist:ntt),1:ny]
  all[ (1+(j-1)*nt):(j*nt),1:ne] = a[(ist:ntt),1:ne]
  
}
#
yx=yall
nyx=nx+ny
if (nx>0) { yx=cbind(y,x) }
write.table(yx,file=filey,append=FALSE,col.names = FALSE,row.names = FALSE)
