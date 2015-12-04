# MSQC -  Multivariate Statistical Quality Control
# Primary Source - Multivariate Statistical Quality Control Using R - By Edgar Santos-Fernández
# https://books.google.co.in/books?id=SZfRrDr7VHYC&printsec=frontcover#v=onepage&q&f=false


library(MSQC)
library("MSQC", lib.loc="~/R/win-library/3.1")
data(package="MSQC")
library(qcc)
citation("qcc")
#Scrucca, L. (2004). qcc: an R package for quality control
#charting and statistical process control. R News 4/1, 11-17.
set.seed(20)

x <-round(rnorm(120,20,2),2)
length <-matrix(x, ncol=4, byrow=TRUE)
par(mfrow=c(1,2))
qcc(length, type="xbar", std.dev="RMSDF"); qcc(length, type="R")
#
cap<-qcc(length, type="xbar", nsigmas=3, plot=T)
process.capability(cap, spec.limits=c(14,26))
#

mu <-c(0,0) # The Mean vector - "mu" of a Bivariate Normal Distribution 
sigma <- matrix(c(10,3,3,6),2,2) # Sigma the Covariance Matrix 
rho <- sigma[1,2] / (sqrt(sigma[1,1] * sigma[2,2])) # https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html
var1<-seq(12,12,.7)
var2<-var1
f<-matrix(0, length(var1), length(var1))
for(i in 1:length(var1))
  {
for(j in 1:length(var1))
  {
  f[i,j]<-1/(2*pi*sqrt(sigma[1,1]*sigma[2,2]*(1-rho^2)))*
    exp(-1/(2*(1-rho^2)) * ((var1[i] - mu[1]) ^ 2 / sigma[1,1] + 
    (var2[j] - mu[2]) ^ 2 /sigma[2,2]-2 * rho * ((var1[i] - mu[1]) * 
  (var2[j] - mu[2])) / (sqrt(sigma[1,1]) * sqrt(sigma[2,2]))))}}
# 
str(f)
#
persp(var1,var2,f,xlab="Var-1",ylab="Var-2",zlab="f(var1,var2)", 
      theta=30,phi=30,r=500,d=0.02,expand=0.6,ltheta=90,lphi=180,
      nticks=4)
# In the Above seen Perspective Plot - persp {graphics}
# f == the matrix containing values to be plotted.

