rm(list=ls())
df = read.csv('~/Documents/wolpertlab/CurlFieldFamilies/dat.csv')
#df$Subject = as.factor(df$Subject)
#df$Session = as.factor(df$Session)
#df$Target = as.factor(df$Target)
#df$Force = scale(df$Force)
#df$Object = df$Object-3

df2 = subset(df,Target==2 & Session==3 & !IsOutlier)[c('Object','Force','Subject')]
xName = "Object" ; yName = "Force" ; sName="Subject"
fileNameRoot = "HierLinRegressData-Jags-"
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("~/Documents/R/DBDA2Eprograms/Jags-Ymet-XmetSsubj-MrobustHier.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=df2 , xName=xName , yName=yName , sName=sName ,
                    numSavedSteps=50000 , thinSteps=7 , saveName=fileNameRoot, nChains = 8 )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
# #------------------------------------------------------------------------------- 
# # Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in c("beta0mu","beta1mu")){#},"nu","sigma","beta0[1]","beta1[1]") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=df2 , xName=xName , yName=yName , sName=sName ,
          compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

library(rmcorr)
for (ti in 1:3){
  df2 = subset(df,Session==3 & Target==ti & !IsOutlier)
  rmc = rmcorr(Subject,Object,Force,df2)
  print(rmc)
  plot(rmc)
}

library(lmerTest)
library(sjPlot)
df2 = subset(df,Session==3 & !IsOutlier)
df2$Object = df2$Object-1
mdl = lmerTest::lmer(Force ~ Object + Target + (1 + Object | Subject), df2)
summary(mdl)
coef(mdl)
plot_model(mdl,type='est')
plot_model(mdl,type='re')

library(rjags)
library(runjags)
source("~/Documents/R/DBDA2Eprograms/DBDA2E-utilities.R")
