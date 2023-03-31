# CultusBassASMR.R
# Monthly time-step age-structured mark-recapture model for Cultus Lake bass fit to 
#   data collected by Wendy Margetts (TRU)
# Estimates will be very uncertain because only two years are availble (difficult to separate
#   recruitment from selectivity, etc.)
# January 15, 2023

library(R2jags)
#library(jagstools)
library(coda)
library(rjags)
library(runjags)
library(here)

# flag to start estimation (=1) or not (=0)
EST <- TRUE

# set path
dummy <- function(){}
path <- getSrcDirectory(dummy)
setwd(path)

# load data (data file already created)
load(file=paste0(path,"/CultusBassASMR.Rdata"))

# break the list apart and create each object
nam <- names(data)
for(i in 1:length(nam)){
  assign( nam[i], data[[i]] )
}

AR     <- 1                  # age at recruitment
A      <- 5                  # maximum age
nyr    <- 2                  # number of years
nT     <- nrow( Marked )     # number of time-steps
nny    <- 9                  # time-step when new year happens
nF     <- 6                  # number of estimated fishing mortality rates (May-Sept + one for other months)

Linf   <- 462
vbK    <- 0.3676
iF     <- c(1:(nF-1),rep(nF,12-(nF-1)),1:(nF-1))   # indexes of F's to use

# Specify model in JAGS
sink(paste0(path,"/ASMR.bug"))
cat("
    model{
    # Priors
    M ~ dlnorm(-1.6,pow(0.8,-2))        # natural mortality rate
    for( a in AR:(A-1) ){
      lNT[a] ~ dunif(0,10)              # terminal abundance by age
      NT[a] <- exp( lNT[a] )
    } #t
    for( a in AR:A ){
      sel[a] ~ dbeta(1,1)               # age-specific selectivity
    } #a
    for( i in 1:nF ){
      lFunobs[i] ~ dunif(-10,1)         # unobserved fishing mortality in some months
    } #t
    for( t in 1:nT ){
      Funobs[t] <- exp(lFunobs[iF[t]])  # fishing mortality put in appropriate months
    } #t

    # Model
    for( a in AR:A ){
      Lt[a] <- Linf*(1-exp(-K*a))       # length at age
    } #a

    # Forecast marked abundance
    for( a in AR:A ){
      Mhat[1,a] <- 0                    # initial abundance of marked fish (=0)
    } #a
    Mhat[nny,1] <- 0                    # initial marked abundance in January, year-2

    # estimated marked abundance in first calendar year
    for( t in 2:(nny-1) ){              # marked abundance progressing within an age across months 
                                        #   in first year
      for( a in AR:A ){
        Mhat[t,a] <- ( Mhat[t-1,a] + MarkedMort[t-1,a] - Marked[t-1,a] ) * 
                        exp( -M/12 - Funobs[t-1]*sel[a] )
      } #a
    } #t
    # estimated marked abundance in first month of new year
    for( a in (AR+1):A ){
      Mhat[nny,a] <- ( Mhat[nny-1,a-1] + MarkedMort[nny-1,a-1] - Marked[nny-1,a-1] ) *
                        exp( -M/12 - Funobs[nny-1]*sel[a-1] )
    } #a
    # estimated marked abundance in second calendar year
    for( t in (nny+1):nT ){
      for( a in AR:A ){
        Mhat[t,a] <- ( Mhat[t-1,a] + MarkedMort[t-1,a] - Marked[t-1,a] ) *
                        exp( -M/12 - Funobs[t-1]*sel[a] )
      } #a
    } #t
    
    # Hindcast unmarked abundance
    for(a in AR:(A-1)){
      Uhat[nT,a] <- NT[a]               # terminal year unmarked fish
    } #a
    for(t in 1:nT){
      Uhat[t,A] <- 0                    # terminal age unmarked fish
    } #t
    # estimated marked abundance in first calendar year
    for(t in 1:Tback[nny]){             # loops only increment up, so use vectors Aback and Tback
      for(a in AR:(A-1)){               # to result in a backwards progression
        Uhat[Tback[t],Aback[a]] <- Uhat[Tback[t]+1,Aback[a]]/
          exp(-M/12-Funobs[Tback[t]]*sel[Aback[a]]) + 
          Marked[Tback[t],Aback[a]] + UnmarkedMort[Tback[t],Aback[a]]   # unmarked abundance
      } #a
    } #t
    # estimated marked abundance in first month of second calendar year
    for(a in AR:(A-1)){
      Uhat[Tback[Tback[nny-1]],Aback[a]] <- Uhat[Tback[Tback[nny-1]]+1,Aback[a]+1]/
        exp(-M/12-Funobs[Tback[Tback[nny-1]]]*sel[Aback[a]]) + 
        Marked[Tback[Tback[nny-1]],Aback[a]] + UnmarkedMort[Tback[Tback[nny-1]],Aback[a]]   # unmarked abundance
    } #a
    # estimated marked abundance in second calendar year
    for(t in Tback[nny-2]:(nT-1)){
      for(a in AR:(A-1)){
        Uhat[Tback[t],Aback[a]] <- Uhat[Tback[t]+1,Aback[a]]/
          exp(-M/12-Funobs[Tback[t]]*sel[Aback[a]]) + 
          Marked[Tback[t],Aback[a]] + UnmarkedMort[Tback[t],Aback[a]]   # unmarked abundance
      } #a
    } #t
  
    # MLE of exploitation rate
    for(t in 1:nT){
      V[t] <- sum( sel * (Uhat[t,]+Mhat[t,]))
      Fobs[t] <- -log( 1 - csum[t]/V[t])  # observed exploitation rate 
    } #t

    # Predict catches by gear
    for(t in 1:nT){
      for(a in AR:A){
        Z[t,a] <- M + (Fobs[t] + Funobs[t])*sel[a]        # total mortality
        # predicted unmarked catch
        uchat[t,a] <- Uhat[t,a] * Fobs[t]*sel[a]/Z[t,a] * (1-exp(-Z[t,a]))
        # predicted marked catch
        mchat[t,a] <- Uhat[t,a] * Fobs[t]*sel[a]/Z[t,a] * (1-exp(-Z[t,a]))   
      } #a
    } #t

    # Likelihood
    Nhat <- Uhat + Mhat
    for(t in 1:nT){
      for(a in AR:A){
        mobs[t,a] ~ dbin( Mhat[t,a]/(Uhat[t,a]+Mhat[t,a]),cobs[t,a])
        #uobs[s,a] ~ dpois(uchat[g,s,a]+1e-10)
        #mobs[s,a] ~ dpois(mchat[g,s,a]+1E-10)
      } #a
    } #t

  } #Model
    ",fill = TRUE)
sink()

# backwards index over ages
Aback <- (A-1):AR
# backwards index over time
Tback <- (nT-1):1
# sum of catch each gear, year
usum <- rowSums(UnmarkedCatch)
msum <- rowSums(MarkedCatch)
jags.data <- list( AR=AR, A=A, nF=nF, nT=nT, iF=iF, Aback=Aback, Tback=Tback,
                  nny=nny, Linf=Linf, K=vbK, usum=usum, msum=msum, csum=usum+msum, Marked=Marked,
                  UnmarkedMort=UnmarkedMort, MarkedMort=MarkedMort,
                  uobs=UnmarkedCatch, mobs=MarkedCatch, cobs=UnmarkedCatch+MarkedCatch )

jags.inits <- function(){
  list( M=rnorm(1,0.4,0.004), lNT=rnorm(A-AR,4,0.1), 
        lsel=rbeta(A-AR+1,1,1), lFunobs=rnorm(nF,-3,0.03) )
}

parameters <- c("M","sel","NT","Funobs","V")

# MCMC settings
ni <- 20000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R
if(EST){
  #  asmr <- run.jags(model=paste0(mod.dir,"/ASMR.bug"), 
  #                   data=jags.data,
  #                   monitor=parameters,inits=jags.inits,
  #                   method="parallel",modules = c("glm","bugs"),
  #                   n.chains=3,sample=1000,
  #                   adapt=1000,burnin=10000,thin=10)
  asmr <- jags.parallel(data=jags.data, parameters.to.save=parameters,
                        inits=jags.inits, n.thin=1, n.cluster = 3,
                        n.iter=200, n.chains=3,
                        model.file=paste0(path,"/ASMR.bug"), 
                        n.burnin=100, DIC=FALSE)
  save(asmr,file="SMB posterior.RData")
}

