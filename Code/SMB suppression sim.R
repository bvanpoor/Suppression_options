# SMB suppression sim.R
# Brett van Poorten
# January 19, 2023
# This individual-based simulation model will simulate different suppression plans for Cultus SMB, 
#   for use with an estimation model to evaluate bias and time to reasonable estimates

# Notes:
#  In future, it would be good to fit effort each year to that estimated from the stated choice 
#    survey. That way, if creel surveys were shitty or missing in some years, effort could still be 
#    imputed from sequential improvements in priors provided by the stated choice survey

# Libraries
library(here)
library(coda)
library(TMB)
library(reshape2)
library(ggplot2)

# file paths
path    <- here()
setwd( paste0( path,"/Code"))

SIM     <- TRUE                    # flag to determine if data should be simulated
EST     <- TRUE                    # flag to determine if parameters should be estimated
TMB     <- FALSE                   # flag to determine if using fast (frequentist) estimation

# Indices
AR      <- 1                       # age at recruitment to gear
A       <- 6                       # maximum age (plus group)
init_yr <- 7                       # number of years between establishment and now
sim_yr  <- 10                      # number of years to project forward
nT      <- init_yr + sim_yr        # simulation years (2017 to now; 20 years forward)
nA      <- A-AR+1                  # number of age-classes to simulate
treg    <- nT

# Control parameters
# exploitation rate of marking method
UMark   <- c( rep( 0, init_yr ), rep(0.2, sim_yr ))#rlnorm( sim_yr, log(0.3), 0.1 )) #0.1          
# effort (days) devoted to nest destruction (helps account for nest re-creation) - in integers
ENest   <- c( rep( 0, init_yr ),rep(c(2,0,2,0),sim_yr))#rep( c(0,0,0,0), length = sim_yr )) #rep(2,nT)          
# exploitation rate of male removal
UMale   <- c( rep( 0, init_yr ), rep (c(0.3,0.3,0.0,0.0), length.out=sim_yr )) #0.3          
# years when regulations are tightened to mandate harvest; pHarv becomes pHarv2
ASyr    <- rep( FALSE, nT )
if ( treg <= nT )
  ASyr[treg:nT]  <- TRUE
# CV of acoustic estimate of fishing mortality when acoustics is used
AcousCV <- c( rep( 0,init_yr ), rep( 0.1, sim_yr ))

# Monitoring parameters
# monitoring fishing effort
FMon    <- c(rep(FALSE,init_yr),rep(TRUE,sim_yr))
F_Cover <- 1/(1+exp(-rnorm(nT,-2, .25)))   # fishery coverage as a proportion of catch enumerated
# CV in estimated fishing effort (can vary based on number of surveys and stratification, etc.)
F_CV    <- rep( 0.2, nT )

# Parameters
R0      <- 2500                    # unfished equilibrium recruitment
Reck    <- 9.2                     # recruitment compensation ratio
M       <- 0.25                    # natural mortality
Linf    <- 50                      # asymptotic length
K       <- 0.35                    # von Bertalanfy growth parameter
CV      <- 0.1                     # coefficient of variation in length at age
Sel     <- plogis((1:nA-2)/.35)      # age-specific selectivity
EFishSel<- c(.4,plogis((1:(nA-1)-2)/-1.5))  # age-specific selectivity to electrofishing (marking)
Mat50   <- 2.5                     # age-at-50% maturity
MatSig  <- 0.5                     # standard deviation in logistic maturity
EggMass <- 10000                   # eggs per kg of female
pFem    <- 0.5                     # sex ratio (50%)
UNest   <- 0.4                     # instantaneous proportion of nests destroyed per unit effort
CPUE0   <- 5                       # mean CPUE at equilibrium maximum density
E_CV    <- 0.2                     # CV in annual effort
N2017   <- c(rep(0,nA-1),20)       # numbers at age in 2017
pHarv   <- 0.3                     # proportion of legal fish that are voluntarily harvested
pHarv2  <- 0.8                     # proportion of legal harvest if regulations strengthened
DisMort <- 0.1                     # release mortality
procSig <- 0.5                     # process error in recruitment
Omega   <- rep(0,nT)               # process error on initial 
Nang    <- 500                     # maximum number of angler days
pTotE   <- 0.9                     # prop max angler days if fishing happened on unfished pop
UtBassCa<- 0.5                     # maximum utility for bass catch rate (fish/day)
UtBassCb<- -0.5                    # loss in utility at 1 CPUE
UtBassCc<- -1                      # exponent of bass catch rate
UtBassWa<- 0                       # median utility for bass catch size at 1 kg
UtBassWb<- 0.05                    # increase in utility with each 100g
UtNatC  <- 0.01                    # utility for catch rates of native fish (fish/day)
UtNatN  <- 0.01                    # utility for proportional change in native fish abundance
UtSupp  <- c( 0, -0.1, -0.2, -0.4 ) # utility for status quo, fishing with reg changes, nest 
                                   # destruction and nest destruction with male removal
pRebuild<- 0.3                     # proportion of males that rebuild nests after they are destroyed
Egrow   <- 1                       # rate at which effort grows per year
Sel     <- Sel/max(Sel)*.95
EFishSel<- EFishSel/max(EFishSel)*.95

pars    <- list( R0=R0, 
                 Reck=Reck, 
                 M=M, 
                 Linf=Linf, 
                 K=K, 
                 CV=CV, 
                 Sel=Sel, 
                 EFishSel=EFishSel,
                 Mat50=Mat50, 
                 MatSig=MatSig, 
                 EggMass=EggMass,
                 pFem=pFem, 
                 UNest=UNest,
                 N2017=N2017,
                 pHarv=pHarv, 
                 pHarv2=pHarv2,
                 DisMort=DisMort, 
                 procSig=procSig,
                 Omega=Omega, 
                 Nang=Nang, 
                 pTotE=pTotE,
                 E_CV=E_CV,
                 CPUE0=CPUE0,
                 UtBassCa=UtBassCa, 
                 UtBassCb=UtBassCb, 
                 UtBassCc=UtBassCc, 
                 UtBassWa=UtBassWa, 
                 UtBassWb=UtBassWb, 
                 UtNatC=UtNatC, 
                 UtNatN=UtNatN, 
                 UtSupp=UtSupp,
                 UMark=UMark, 
                 ENest=ENest, 
                 UMale=UMale, 
                 pRebuild=pRebuild,
                 Egrow=Egrow )

##################################################################################################
#########################################   PROCEDURES   #########################################
##################################################################################################

"RandomParameters" <- function( ParsIn ){
  # randomize parameters and processes
  R0         <- rlnorm( 1, log(ParsIn$R0), 0.2 )
  Reck       <- rlnorm( 1, log(ParsIn$Reck), 0.2 )
  M          <- rlnorm( 1, log(ParsIn$M), 0.2 )
  Linf       <- rlnorm( 1, log(ParsIn$Linf), 0.2 )
  K          <- rlnorm( 1, log(ParsIn$K), 0.1 )
  CV         <- rlnorm( 1, log(ParsIn$CV), 0.1 )
  Sel        <- 1 / ( 1 + exp(-rnorm( nA, log(ParsIn$Sel/(1-ParsIn$Sel)), 0.2 ) ))
  EFishSel   <- 1 / ( 1 + exp(-rnorm( nA, log(ParsIn$EFishSel/(1-ParsIn$EFishSel)), 0.2 ) ))
  Mat50      <- rlnorm( 1, log(ParsIn$Mat50), 0.2 )
  MatSig     <- rlnorm( 1, log(ParsIn$MatSig), 0.2 )
  EggMass    <- rlnorm( 1, log(ParsIn$EggMass), 0.2 )
  pFem       <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$pFem/(1-ParsIn$pFem)), 0.2 ) ))
  UNest      <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$UNest/(1-ParsIn$UNest)), 0.2 ) ))
  CPUE0      <- rlnorm( 1, log(ParsIn$CPUE0), 0.2 )
  E_CV       <- rlnorm( 1, log(ParsIn$E_CV), 0.2 )
  N2017      <- rpois( nA, N2017 )
  pHarv      <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$pHarv/(1-ParsIn$pHarv)), 0.2 ) ))
  pHarv2      <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$pHarv2/(1-ParsIn$pHarv2)), 0.2 ) ))
  DisMort    <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$DisMort/(1-ParsIn$DisMort)), 0.2 ) ))
  Omega      <- rnorm( nT, ParsIn$Omega, ParsIn$procSig )
  Nang       <- rlnorm( 1, log(ParsIn$Nang), 0.2 )
  pTotE      <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$pTotE/(1-ParsIn$pTotE)), 0.2 ) ))
  UtBassCa   <- rnorm( 1, ParsIn$UtBassCa, abs( ParsIn$UtBassCa*0.1 ))
  UtBassCb   <- rnorm( 1, ParsIn$UtBassCb, abs( ParsIn$UtBassCb*0.1 ))
  UtBassCc   <- rnorm( 1, ParsIn$UtBassCc, abs( ParsIn$UtBassCc*0.1 ))
  UtBassWa   <- rnorm( 1, ParsIn$UtBassWa, abs( ParsIn$UtBassWa*0.1 ))
  UtBassWb   <- rnorm( 1, ParsIn$UtBassWb, abs( ParsIn$UtBassWb*0.1 ))
  UtNatC     <- rnorm( 1, ParsIn$UtNatC, abs( ParsIn$UtNatC*0.2 ))
  UtNatN     <- rnorm( 1, ParsIn$UtNatN, abs( ParsIn$UtNatN*0.2 ))
  UtSupp     <- rnorm( 4, ParsIn$UtSupp, abs( ParsIn$UtSupp*0.2 ))
  UMark      <- rlnorm( nT, log(ParsIn$UMark), 0.2 )
  ENest      <- rpois( nT, ParsIn$ENest )
  UMale      <- rlnorm( nT, log(ParsIn$UMale), 0.2 )
  pRebuild   <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$pRebuild/(1-ParsIn$pRebuild)), 0.2 ) ))
  Egrow      <- 1 / ( 1 + exp(-rnorm( 1, log(ParsIn$Egrow/(1-ParsIn$Egrow)), 0.2 ) ))
   
  ParsOut   <- list( R0=R0, 
                     Reck=Reck, 
                     M=M, 
                     Linf=Linf, 
                     K=K, 
                     CV=CV, 
                     Sel=Sel, 
                     EFishSel=EFishSel,
                     Mat50=Mat50, 
                     MatSig=MatSig, 
                     EggMass=EggMass,
                     pFem=pFem, 
                     UNest=UNest,
                     N2017=N2017,
                     CPUE0=CPUE0,
                     E_CV=E_CV,
                     pHarv=pHarv, 
                     pHarv2=pHarv2,
                     DisMort=DisMort, 
                     procSig=procSig,
                     Omega=Omega,
                     Nang=Nang, 
                     pTotE=pTotE,
                     UtBassCa=UtBassCa, 
                     UtBassCb=UtBassCb, 
                     UtBassCc=UtBassCc, 
                     UtBassWa=UtBassWa, 
                     UtBassWb=UtBassWb, 
                     UtNatC=UtNatC, 
                     UtNatN=UtNatN, 
                     UtSupp=UtSupp,
                     UMark=UMark, 
                     ENest=ENest, 
                     UMale=UMale, 
                     pRebuild=pRebuild,
                     Egrow=Egrow )
  return( ParsOut )
} # RandomParameters()

"Initialize" <- function( pars, RAND=TRUE ){
  # read in parameters
  if( RAND )
    randpars   <- RandomParameters( ParsIn=pars )
  else
    randpars   <- pars
  R0         <- randpars$R0
  Reck       <- randpars$Reck
  M          <- randpars$M
  Linf       <- randpars$Linf
  K          <- randpars$K
  Sel        <- randpars$Sel
  EFishSel   <- randpars$EFishSel
  Mat50      <- randpars$Mat50
  MatSig     <- randpars$MatSig
  EggMass    <- randpars$EggMass
  pFem       <- randpars$pFem
  UNest      <- randpars$UNest
  N2017      <- randpars$N2017
  pTotE      <- randpars$pTotE
  CPUE0      <- randpars$CPUE0
  E_CV       <- randpars$E_CV
  pHarv      <- randpars$pHarv
  pHarv2     <- randpars$pHarv2
  DisMort    <- randpars$DisMort
  Omega      <- randpars$Omega
  Nang       <- randpars$Nang
  UtBassCa   <- randpars$UtBassCa
  UtBassCb   <- randpars$UtBassCb
  UtBassCc   <- randpars$UtBassCc
  UtBassWa   <- randpars$UtBassWa
  UtBassWb   <- randpars$UtBassWb
  
  # set equilibrium conditions
  Age        <- AR:A
  # length-at-age
  Lt         <- Linf*( 1 - exp( -K*Age ))
  # weight-at-age
  Wt         <- 1E-5 * Lt^3
  # probability of maturity
  Fec        <- Wt * EggMass
  # eggs per female
  Mat        <- 1 / ( 1 + exp( -( Age - Mat50 )/MatSig ))
  # survivorship
  lx         <- exp( -M*(Age-1) )
  lx[A]      <- lx[A] / ( 1 - exp( -M ))
  lxF        <- rep( 1, nA )
  Sel        <- Sel
  EFishSel   <- EFishSel
  
  # unfished abundance
  N0         <- R0 * lx
  # unfished vulnerable abundance
  V0         <- sum( N0 * Sel )
  # expected weight at unexploited
  W0         <- sum( N0 * Sel * Wt )/V0
  # catchability of bass ( calculated: CPUE0/sum(R0*lx*Sel) where CPUE0 = 10 fish/day)
  q          <- CPUE0 / V0
  # exploitation rate on if fishing on unfished population
  F0         <- q * Nang * pTotE
  
  # total utility that drives effort (no suppression is done in the optimal 
  #   situation, so are avoided here)
  UtotHat    <- UtBassCa + UtBassCb*CPUE0^UtBassCc + UtBassWa + (W0-1)*10 * UtBassWb
  #pfish      <- exp( UtotHat ) / ( exp( UtOth ) + exp( UtotHat ))
  # rearrange to solve for UtOth
  # utility of not fishing for bass when bass abundance highest (assuming 95% max effort)
  UtOth      <- log( exp( UtotHat )/pTotE - exp( UtotHat ) )
  
  
  for( a in (AR+1):A )
    lxF[a]   <- lxF[a-1] * exp( -M - F0*Sel[a]*( pHarv + (1-pHarv)*DisMort ))
  lxF[A]     <- lxF[A] / ( 1 - exp( -M -F0*Sel[A]*( pHarv + (1-pHarv)*DisMort )))
  
  # calculate Beverton-Holt recruitment parameters
  # eggs per recruit
  phi0       <- sum( Fec * Mat * lx ) * pFem            # phi0 adjusted by the sex ratio 
                                                        # so any unpaired females don't spawn
  # recruitment function
  RecAlpha     <- Reck / phi0
  RecBeta      <- log( Reck )/( R0 * phi0 )
  
  out          <- list()
  out$R0       <- R0
  out$Reck     <- Reck
  out$M        <- M
  out$Linf     <- Linf
  out$K        <- K
  out$CV       <- CV
  out$Sel      <- Sel
  out$EFishSel <- EFishSel
  out$Mat50    <- Mat50
  out$MatSig   <- MatSig
  out$EggMass  <- EggMass
  out$pFem     <- pFem
  out$UNest    <- UNest
  out$F0       <- F0
  out$CPUE0    <- CPUE0
  out$E_CV     <- E_CV
  out$pHarv    <- pHarv
  out$pHarv2   <- pHarv2
  out$DisMort  <- DisMort
  out$procSig  <- procSig
  out$Omega    <- Omega
  out$Nang     <- randpars$Nang
  out$pTotE    <- pTotE
  out$UtBassCa <- UtBassCa
  out$UtBassCb <- UtBassCb
  out$UtBassCc <- UtBassCc
  out$UtBassWa <- UtBassWa
  out$UtBassWb <- UtBassWb
  out$UtNatC   <- randpars$UtNatC
  out$UtNatN   <- randpars$UtNatN
  out$UtSupp   <- randpars$UtSupp
  out$UMark    <- randpars$UMark
  out$ENest    <- randpars$ENest
  out$UMale    <- randpars$UMale
  out$pRebuild <- randpars$pRebuild
  out$Egrow    <- randpars$Egrow

  out$Age      <- Age
  out$Lt       <- Lt
  out$Wt       <- Wt
  out$Mat      <- Mat
  out$Fec      <- Fec
  out$lx       <- lx
  out$lxF      <- lxF
  out$phi0     <- phi0
  out$RecA     <- RecAlpha
  out$RecB     <- RecBeta
  out$N0       <- N0
  out$N2017    <- N2017
  out$q        <- q
  out$UtOth    <- UtOth
  return( out )
} # Initialize()

"Dynamics" <- function( init, itreg=treg, FCover=F_Cover ){
  iASyr       <- rep( FALSE, nT )
  if ( itreg <= nT )
    iASyr[itreg:nT]  <- TRUE
#  UMark      <- iUMark
#  ENest      <- iENest
#  UMale      <- iUMale
  
  # Myr is the years when fish are marked
  # NSyr is the years when nest suppression happens
  # MSyr is the years when male removal happens
  # ASyr is the years when regulations are tightened to mandate harvest pHarv becomes pHarv2
  
  # parameters and initial rates
  M          <- init$M
  CV         <- init$CV
  Sel        <- init$Sel
  EFishSel   <- init$EFishSel
  pFem       <- init$pFem
  UNest      <- init$UNest
  E_CV       <- init$E_CV
  pHarv      <- init$pHarv
  pHarv2     <- init$pHarv2
  DisMort    <- init$DisMort
  Omega      <- init$Omega
  Nang       <- init$Nang
  UtBassCa   <- init$UtBassCa
  UtBassCb   <- init$UtBassCb
  UtBassCc   <- init$UtBassCc
  UtBassWa   <- init$UtBassWa
  UtBassWb   <- init$UtBassWb
  UtSupp     <- init$UtSupp
  UMark      <- init$UMark
  ENest      <- init$ENest
  UMale      <- init$UMale
  pRebuild   <- init$pRebuild
  Egrow      <- init$Egrow
  Age        <- init$Age
  Lt         <- init$Lt
  Wt         <- init$Wt
  Mat        <- init$Mat
  Fec        <- init$Fec
  Alpha      <- init$RecA
  Beta       <- init$RecB
  N2017      <- init$N2017
  q          <- init$q
  UtOth      <- init$UtOth

  iMyr       <- UMark > 0
  iNSyr      <- ENest > 0
  iMSyr      <- UMale > 0
  # establish objects to be filled
  UtF <- MtF <- UtM <- MtM <- matrix( nrow=nT+1, ncol=nA )      # numbers at large over time
  FMarks <- MMarks <- FRecs <- MRecs <- matrix( nrow=nT, ncol=nA ) # numbers marked each year
  FCap <- MCap <- FRecap <- MRecap <- rep( 0,nT )               # sum of numbers capture and
                                                                # recaptured by efishing
  CUMal <- CMMal <- matrix( nrow=nT, ncol=nA )                  # catch of unmarked and marked males off nests
  Cr_CU <- Cr_CM <- vector( )                                   # catch rates observed in creel

  Fem0       <- rbinom( nA, N2017, pFem )
  FMarks[1,] <- rbinom( nA, Fem0, EFishSel*UMark[1] )                # number of marked females
  MMarks[1,] <- rbinom( nA, N2017 - Fem0, EFishSel*UMark[1] )        # number of marked males
  FRecs[1,]  <- rep( 0, nA )
  MRecs[1,]  <- rep( 0, nA )
  FCap[1]    <- sum( FMarks[1,] )
  MCap[1]    <- sum( MMarks[1,] )
  MtF[1, ]   <- FMarks[1,]
  MtM[1, ]   <- MMarks[1,]
  UtF[1, ]   <- Fem0 - MtF[1,]                              # unmarked females in first year
  UtM[1, ]   <- N2017 - Fem0 - MtM[1,]                      # unmarked males in first year
  
  # Mature unmarked and marked fish
  UMatMal    <- rbinom( nA, UtM[1,], Mat )                  # unmarked mature males 
  UMatFem    <- rbinom( nA, UtF[1,], Mat )                  # unmarked mature females 
  MMatMal    <- rbinom( nA, MtM[1,], Mat )                  # marked mature males 
  MMatFem    <- rbinom( nA, MtF[1,], Mat )                  # marked mature females 
  Eggs       <- vector()                                    # actual eggs laid
  predEggs   <- vector()                                    # predicted eggs laid ignoring male removal and nest destruction
  ipHarv     <- vector()
  CUtF <- CUtM <- CMtF <- CMtM <- matrix( nrow=nT, ncol=nA )     # catch of unmarked and marked fish in fishery
  Et <- Ft <- Zt <- ZEst <- Catch <- vector()                      # fishing effort, fishing and total mortality, estimated total Z and catch each year
  CNest      <- matrix( data=0, nrow=nT, ncol=max(ENest,1) ) # number of nests destroyed each year
  CMale      <- matrix( nrow=nT, ncol=nA )        # number of males removed each year
  # eggs calculated as the number of eggs produced by females; assumes the number of spawning males 
  #    and females is proportional (1:1 pairing; no more)
  MalSpawn   <- UMatMal + MMatMal
  FemSpawn   <- UMatFem + MMatFem
  predEggs[1]<- sum( FemSpawn*Fec )
  Nests      <- MalSpawn
  # nests destroyed
  tmp_CU <- tmp_CM <- matrix(0,nrow=max(ENest,1),ncol=nA)
  if(ENest[1]>0){
    for( i in 1:ENest[1] ){
      NestsDest  <- rbinom( nA, Nests, UNest )
      tmp_CU[i,] <- rbinom( nA, UMatMal, UMale[1] )             # unmarked male removals off nests
      tmp_CM[i,] <- rbinom( nA, MMatMal, UMale[1])              # marked male removals off nests
      UMatMal    <- UMatMal - tmp_CU[i,]
      MMatMal    <- MMatMal - tmp_CM[i,]
      rebuilds   <- rbinom( nA, NestsDest, pRebuild )
      CNest[1,i] <- sum( NestsDest )
      Nests      <- Nests - NestsDest + rebuilds
    }
  }
  CUMal[1,]  <- colSums( tmp_CU )
  CMMal[1,]  <- colSums( tmp_CM )
  # proportional reduction in spawning - assume 1:1 pairing
  pReduc     <- min( 1, sum( Nests ) / sum( FemSpawn ) )  
  Eggs[1]    <- as.integer( sum( FemSpawn*Fec )*pReduc )
  
  # fishery effort and catch
  ECPUE      <- sum( N2017*Sel ) * q                      # expected CPUE
  EWt        <- sum( N2017*Sel*Wt ) / sum( N2017*Sel )    # expected mean weight
  # total utility that drives effort
  UtotHat    <- UtBassCa + UtBassCb*ECPUE^UtBassCc + UtBassWa + (EWt-1)*10 * UtBassWb +
                  UtSupp[2]*iASyr[1] + UtSupp[3]*iNSyr[1] + UtSupp[4]*iMSyr[1]
  ipHarv[1]  <- pHarv * (!iASyr[1]) + pHarv2 * iASyr[1]
  pfish      <- exp( UtotHat ) / ( exp( UtOth ) + exp( UtotHat ))   # proportion of total angler days 
  Et[1]      <- pfish * Nang * ( 1 + rnorm(1,0,1)*E_CV )            # Effort in year-1
  Ft[1]      <- q*Et[1]
  Zt         <- M + Ft[1] * Sel * ( ipHarv[1] + ( 1-ipHarv[1])*DisMort )
  
  # total catch
  NTemp      <- UtF[1,] + UtM[1,] - CUMal[1,] + MtF[1,] + MtM[1,] - CMMal[1,]
  Catch[1]   <- sum( NTemp * q * Sel / Zt * ( 1-exp(-Zt)) )
  
  # predicted F from acoustics
  ZEst[1]    <- M + Ft[1] * ( ipHarv[1] + ( 1-ipHarv[1] ) * DisMort )
  ZEst[1]    <- ZEst[1] * ( 1 + rnorm( 1, 0, AcousCV[1] ))
  
  # harvest
  pharvest   <- ( Ft[1]*Sel*ipHarv[1] ) / Zt * ( 1-exp(-Zt) )
  CUtF[1,]   <- as.integer( UtF[1,] * pharvest )               # fishery catch of unmarked females
  CUtM[1,]   <- as.integer( (UtM[1,]-CUMal[1,]) * pharvest )   # fishery catch of unmarked males
  CMtF[1,]   <- as.integer( MtF[1,] * pharvest )               # fishery catch of marked females
  CMtM[1,]   <- as.integer( (MtM[1,]-CMMal[1,]) * pharvest )   # fishery catch of marked males

  # CPUE estimates in creel
  Cr_CU[1]  <- sum( ( UtF[1,] + UtM[1,] - CUMal[1,] ) * q * ipHarv[1] * Sel / Zt * ( 1-exp(-Zt)) )
#    rbinom( nA, CUtF[1,] + CUtM[1,], FCover[1] ) / ( Et[1] * FCover[1]) * FMon[1]
  Cr_CM[1]  <- sum( ( MtF[1,] + MtM[1,] - CMMal[1,] ) * q * ipHarv[1] * Sel / Zt * ( 1-exp(-Zt)) )
#    rbinom( nA, CMtF[1,] + CMtM[1,], FCover[1] ) / ( Et[1] * FCover[1]) * FMon[1]

  # project forward through time
  for( t in 2:nT ){
    # calculate recruits and add them to the unmarked population
    Rt         <- rbinom( 1, Eggs[t-1], Alpha * exp( -Beta * Eggs[t-1] + Omega[t] ))  # recruits
    UtF[t,1]   <- rbinom( 1, Rt, pFem )                            # females among recruits
    UtM[t,1]   <- Rt - UtF[t,1]                                    # males among recruits
    # unmarked females
    UtF[t,2:nA]   <- rbinom( nA-1, UtF[t-1,1:(nA-1)], exp(-Zt[1:(nA-1)]) ) 
    UtF[t,nA]     <- UtF[t,nA] + rbinom( 1, UtF[t-1,nA], exp(-Zt[nA]) )
    # unmarked males
    UtM[t,2:nA]   <- rbinom( nA-1, UtM[t-1,1:(nA-1)], exp(-Zt[1:(nA-1)]) ) 
    UtM[t,nA]     <- UtM[t,nA] + rbinom( 1, UtM[t-1,nA], exp(-Zt[nA]) )
    # marked females
    MtF[t,]       <- c( 0, rbinom( nA-1, MtF[t-1,1:(nA-1)], exp(-Zt[1:(nA-1)]) ))
    MtF[t,nA]     <- MtF[t,nA] + rbinom( 1, MtF[t-1,nA], exp(-Zt[nA]) )
    # marked males
    MtM[t,]       <- c( 0, rbinom( nA-1, MtM[t-1,1:(nA-1)], exp(-Zt[1:(nA-1)]) ))
    MtM[t,nA]     <- MtM[t,nA] + rbinom( 1, MtM[t-1,nA], exp(-Zt[nA]) )
    
    # calculate number of fish marked and move them from unmarked to marked population
    FMarks[t,] <- rbinom( nA, UtF[t,], EFishSel*UMark[t] )         # number of marked females
    MMarks[t,] <- rbinom( nA, UtM[t,], EFishSel*UMark[t] )         # number of marked males
    FRecs[t,]  <- rbinom( nA, MtF[t,], EFishSel*UMark[t] )         # number of recaped females
    MRecs[t,]  <- rbinom( nA, MtM[t,], EFishSel*UMark[t] )         # number of recaped males
    FCap[t]    <- sum( FMarks[t,] )                                # sum of marked females
    MCap[t]    <- sum( MMarks[t,] )                                # sum of marked males
    FRecap[t]  <- sum( FRecs[t,] )                                 # sum of recaptured females
    MRecap[t]  <- sum( MRecs[t,] )                                 # sum of recaptured males
    UtF[t, ]   <- UtF[t,] - FMarks[t,]                             # unmarked females
    UtM[t, ]   <- UtM[t,] - MMarks[t,]                             # unmarked males
    MtF[t, ]   <- MtF[t,] + FMarks[t,]                             # total number of marked females
    MtM[t, ]   <- MtM[t,] + MMarks[t,]                             # total number of marked males
    
    # number of unmarked and marked males caught off nests
    UMatMal    <- rbinom( nA, UtM[t,], Mat )                      # unmarked mature males 
    UMatFem    <- rbinom( nA, UtF[t,], Mat )                      # unmarked mature females 
    MMatMal    <- rbinom( nA, MtM[t,], Mat )                      # marked mature males 
    MMatFem    <- rbinom( nA, MtF[t,], Mat )                      # marked mature females 
    CUMal[t, ] <- rbinom( nA, UtM[t,], Mat*UMale[t] )            # unmarked male removals off nests
    CMMal[t, ] <- rbinom( nA, MtM[t,], Mat*UMale[t])             # marked male removals off nests
    
    # number of remaining spawners
    # assume number of nests is equal to male spawners before removals
    MalSpawn   <- UMatMal + MMatMal
    FemSpawn   <- UMatFem + MMatFem
    predEggs[t]<- sum( FemSpawn*Fec )
    Nests      <- as.integer( ( UtM[t,] + MtM[t,] ) * Mat )# MalSpawn
    # nests destroyed in an initial run
    tmp_CU <- tmp_CM <- matrix(0,nrow=max(ENest,1),ncol=nA)
    if(ENest[t]>0){
      for( i in 1:ENest[t] ){
        NestsDest  <- rbinom( nA, Nests, UNest )
        tmp_CU[i,] <- rbinom( nA, UtM[t,], Mat*UMale[t] )             # unmarked male removals off nests
        tmp_CM[i,] <- rbinom( nA, MtM[t,], Mat*UMale[t])              # marked male removals off nests
        UtM[t,]    <- UtM[t,] - tmp_CU[i,]
        MtM[t,]    <- MtM[t,] - tmp_CM[i,]
        rebuilds   <- rbinom( nA, NestsDest, pRebuild )
        CNest[t,i] <- sum( NestsDest )
        Nests      <- Nests - NestsDest + rebuilds
      }
      # this overwrites  previous calculation if nest destruction happens more than once per year
      CUMal[t,]  <- colSums( tmp_CU )
      CMMal[t,]  <- colSums( tmp_CM )
    }

    # proportional reduction in spawning - assume 1:1 pairing
    pReduc     <- min( 1, sum( Nests ) / sum( FemSpawn ) )  
    Eggs[t]    <- as.integer( sum( FemSpawn*Fec )*pReduc )
    
    # fishery effort and catch
    # proportion of legal fish harvested (changes if regulations become stricter) 
    ipHarv[t]     <- pHarv * (!iASyr[t]) + pHarv2 * iASyr[t]            
    Nt         <- UtM[t,] + UtF[t,] + MtM[t,] + MtF[t,] - CUMal[t,] - CMMal[t,]
    ECPUE      <- sum( Nt*Sel ) * q                      # expected CPUE
    EWt        <- sum( Nt*Sel*Wt ) / sum( Nt*Sel )    # expected mean weight
    UtotHat    <- UtBassCa + UtBassCb*ECPUE^UtBassCc + UtBassWa + (EWt-1)*10 * UtBassWb +
      UtSupp[2]*iASyr[t] + UtSupp[3]*iNSyr[t] + UtSupp[4]*iMSyr[t]# total utility that drives effort
    pfish      <- exp( UtotHat ) / ( exp( UtOth ) + exp( UtotHat )) # proportion of total angler days 
    # Effort in year-1
    Et[t]      <- pfish * Nang * ( 1 + rnorm(1,0,1)*E_CV ) * Egrow + (1-Egrow) * Et[t-1]        
    Ft[t]      <- q*Et[t]
    Zt         <- M + Ft[t] * Sel * ( ipHarv[t] + (1-ipHarv[t])*DisMort )
    
    # total catch per unit effort
    NTemp      <- UtF[t,] + UtM[t,] + MtF[t,] + MtM[t,] - CUMal[t,] - CMMal[t,]
    Catch[t]   <- sum( NTemp * q * Sel / Zt * ( 1-exp(-Zt)) )
    
    # predicted F from acoustics
    ZEst[t]    <- M + Ft[t] * ( ipHarv[t] + ( 1-ipHarv[t] ) * DisMort )
    ZEst[t]    <- ZEst[t] * ( 1 + rnorm( 1, 0, AcousCV[t] ))
    
    # harvest in the fishery
    pharvest   <- ( Ft[t]*Sel*ipHarv[t] ) / Zt * ( 1-exp(-Zt) )
    CUtF[t,]   <- rbinom( nA, UtF[t,], pharvest )                     # fishery catch of unmarked females
    CUtM[t,]   <- rbinom( nA, UtM[t,], pharvest )         # fishery catch of unmarked males
    CMtF[t,]   <- rbinom( nA, MtF[t,], pharvest )                     # fishery catch of marked females
    CMtM[t,]   <- rbinom( nA, MtM[t,], pharvest )         # fishery catch of marked males
 
    # HPUE estimates in creel
    Cr_CU[t]  <- sum( ( UtF[t,] + UtM[t,] - CUMal[t,] ) * q * ipHarv[t] * Sel / Zt * ( 1-exp(-Zt)) )
    #    rbinom( nA, CUtF[t,] + CUtM[t,], FCover[t] ) / ( Et[t] * FCover[t]) * FMon[t]
    Cr_CM[t]  <- sum( ( MtF[t,] + MtM[t,] - CMMal[t,] ) * q * ipHarv[t] * Sel / Zt * ( 1-exp(-Zt)) )
    #    rbinom( nA, CMtF[t,] + CMtM[t,], FCover[t] ) / ( Et[t] * FCover[t]) * FMon[t]
  } # t
  
  Nt         <- UtF + UtM + MtF + MtM
  Vt         <- t( apply( Nt, 1, '*', Sel ) )
  Ct         <- CUtF + CUtM + CMtF + CMtM
  Marks      <- (FMarks+MMarks)
  Recaps     <- (FRecs+MRecs)
  sim_yrs    <- (init_yr+1):nT
  tmp.yrs    <- sim_yrs[-sim_yr]
  Nest.out  <- glm( (Eggs/predEggs)[tmp.yrs] ~ 0 + ENest[tmp.yrs] + UMale[tmp.yrs], 
                     family=gaussian(link='log'))$coefficients
  
  out        <- list()   # elements to return from function
  out$FMarks <- FMarks
  out$MMarks <- MMarks
  out$FRecs  <- FRecs
  out$MRecs  <- MRecs
  out$UtM    <- UtM
  out$UtF    <- UtF
  out$MtM    <- MtM
  out$MtF    <- MtF
  out$FCap   <- FCap
  out$MCap   <- MCap
  out$FRecap <- FRecap
  out$MRecap <- MRecap
  out$Marks  <- Marks
  out$Recaps <- Recaps
  out$CUMal  <- CUMal
  out$CMMal  <- CMMal
  out$CNest  <- CNest
  out$Eggs   <- Eggs
  out$Et     <- Et
  out$Ft     <- Ft
  out$Nt     <- Nt
  out$Vt     <- Vt
  out$Ct     <- Ct
  out$pHarv  <- ipHarv
  out$CUtF   <- CUtF
  out$CUtM   <- CUtM
  out$CMtF   <- CMtF
  out$CMtM   <- CMtM
  out$pHarv  <- ipHarv
  out$Catch  <- Catch
  out$ZEst   <- ZEst
  out$aNest  <- Nest.out[1]
  out$bNest  <- Nest.out[2]
  out$Cr_CU  <- Cr_CU[FMon]
  out$Cr_CM  <- Cr_CM[FMon]

  return( out )
  
} # Dynamics()

"ErrorEval" <- function(n=1, sim=SIM, est=EST, iUMark=UMark,iENest=ENest,iUMale=UMale,itreg=treg){
#  n=1
#  sim=SIM
#  est=EST
#  iUMark=UMark
#  iENest=ENest
#  iUMale=UMale
#  itreg=treg
  sim_pars <- est_pars <- matrix(nrow=n, ncol=6+sim_yr+nA*3)
  sim_info <- est_info <- matrix(nrow=n, ncol=nA+sim_yr*2+sim_yr-1+4)
  sim_N    <- est_N    <- matrix(nrow=n, ncol=nA*sim_yr)
  
  sim_yrs          <- (init_yr+1):nT
  EMon_yrs         <- which( FMon )
  UMark            <- iUMark
  ENest            <- iENest
  UMale            <- iUMale
  regyr            <- itreg
  if( treg == nT )
    regyr <- 4
  ipars            <- pars
  ipars$UMark      <- UMark
  ipars$ENest      <- ENest
  ipars$UMale      <- UMale

  for( i in 1:n ){
    ii               <<- i
    init             <- Initialize( pars=ipars )
    dyn              <- Dynamics( init=init, itreg=treg )
    Mon_E            <- dyn$Et * ( 1 + rnorm( n=nT, 0, 1) * F_CV ) * FMon
    init_N           <- (dyn$UtF+dyn$UtM)[init_yr+1,]
    aNest <- bNest   <- 0
    if( sum(ENest) > 0 ){
      aNest          <- dyn$aNest
    }
    if( sum(UMale) > 0 ){
      bNest          <- dyn$bNest
    }
    
    DATA <- list( 
      AR=AR, 
      A=A, 
      nT=sim_yr, 
      regyr=regyr-1,                       # year when regs are tightened (minus 1 to account for
                                             # zero-based counting in cpp)
      Linf=init$Linf,                      # von B asymptotic size 
      K=init$K,                            # von B metabolic growth parameter
      LWa=1e-5,                            # scalar of length-weight relationship
      LWb=3,                               # exponent of length-weight relationship
      EggMass=init$EggMass,                # eggs per gram of female
      pFem=init$pFem,                      # equilibrium and initial proportion of females
      Mat50=init$Mat50,                    # age at 50% maturity
      MatSig=init$MatSig,                  # variation in logistic maturity-at-age
      DisMort=init$DisMort,                # discard mortality
      eps=0.001,                           # lower limit at which penalties start applying 
      proc_sig=init$procSig,               # process error around recruitment
      proc_sig2=max(F_CV),                 # process error around effort
      pHarv=dyn$pHarv[sim_yrs],                     # proportion harvested
      Age=init$Age,                        # Ages included (from AR:A)
      ENest=ENest[sim_yrs],                # number of days devoted to destroying nests
      # estimated fishing effort
      E_Hat=Mon_E[EMon_yrs],               # estimated fishing effort
      # marks applied each year
      Cap=(dyn$FCap+dyn$MCap)[sim_yrs],    # total annual capture of unmarked fish
      Recap=(dyn$FRecap+dyn$MRecap)[sim_yrs], # total annual recapture of marked fish
      Catch=dyn$Catch[sim_yrs],            # total catch of fish (whether harvested or not)
      ZEst=dyn$ZEst[sim_yrs],              # total mortality of vulnerable fish (M+F+discard mort)
      Cr_CU=dyn$Cr_CU,                     # CPUE of unmarked fish in creel
      Cr_CM=dyn$Cr_CM,                     # CPUE of marked fish in creel
      Marks=t( dyn$Marks[sim_yrs,] ),
      Recs=t( dyn$Recaps[sim_yrs,] ),
      CUMal=t(dyn$CUMal[sim_yrs,]),        # removal of unmarked males
      CMMal=t(dyn$CMMal[sim_yrs,])         # removal of marked males
    )
    
    PARAMETERS <- list(
      logR0=log(init$R0),
      logRecK=log(init$Reck-1),
      logM=log(init$M),
      logq=log(init$q),
      aNest=0,
      bNest=0,
      logSigma=2,
      logSigma_CPUE=0,
      logSigma_Z=0,
      logitSel=qlogis(init$Sel),
      logitEFishSel=qlogis(init$EFishSel),
      logInitN=log( pmax( init_N, 1)),
      Eps=rnorm(sim_yr-1,0,0.001),
      Zeta=rnorm(sim_yr,0,0.001)
    )
    
    if(est){
      compile( "SMB_suppression.cpp")
      dyn.load( dynlib( "SMB_suppression" ))
      # PHASE 1
      map <- list(
        aNest=factor(NA),
        bNest=factor(NA),
        Eps=factor(rep(NA,sim_yr-1)),
        Zeta=factor(rep(NA,sim_yr))
      )
      if( sum(ENest) > 0 ){
        map$aNest <- NULL
      } 
      if( sum(UMale) > 0 ){
        map$bNest <- NULL
      }
      obj <- MakeADFun( data=DATA, parameters=PARAMETERS, map=map, 
                        DLL="SMB_suppression", silent = TRUE )
      
      opt <- nlm( f=obj$fn, p=obj$par, gradient=obj$gr, hessian=TRUE, print.level = 1,
                  iterlim=500 )#, steptol=1e-7)
      # PHASE 2
      p1  <- obj$env$parList(opt$estimate)          # Parameter estimate after phase 1
      map$Zeta <- NULL                              # remove Zeta from map, so it is now estimated
      obj <- MakeADFun( data=DATA, parameters=p1, random=c("Eps","Zeta"), map=map, 
                        DLL="SMB_suppression", silent = TRUE )
      
      opt <- nlm( f=obj$fn, p=obj$par, gradient=obj$gr, hessian=TRUE, print.level = 1,
                  iterlim=500 )#, steptol=1e-7)
      reprt <- obj$report
      sdreprt <- sdreport(obj)
      if( is.finite( sum(sdreprt$sd ) )){
        sim_pars[i,] <- c(
          init$R0,
          init$Reck,
          init$M,
          init$q,
          dyn$aNest,
          dyn$bNest,
          dyn$pHarv[sim_yrs],
          init$Sel,
          init$EFishSel,
          init_N
        )
        sim_info[i,] <- c(
          dyn$Nt[nT,],                   # final abundance
          dyn$Nt[sim_yrs,1],             # recruitment = used to evaluate effectiveness of nest destruction
          rowMeans(outer(init$UMark[sim_yrs[-1]],init$EFishSel)),       # electrofishing sampling rate
          init$UMale[sim_yrs],           # sampling rate of males off nests
          dyn$pHarv[c(init_yr+1,nT)],    # proportional harvest before and after reg changes
          dyn$aNest,                     # effect of nests destruction on egg deposition
          dyn$bNest                      # effect of male removal on egg deposition
        )
        sim_N[i,]    <- as.vector( t( dyn$Nt[sim_yrs,] ))
          
        est_pars[i,] <- sdreprt$value
        est_info[i,] <- c(
          reprt()$Nt[,sim_yr],           # length: nA
          reprt()$Nt[1,],                # length: sim_yr
          rowMeans(outer(reprt()$EfishU[2:sim_yr],
                         sdreprt$value[6+sim_yr+nA+1:nA])),# length: sim_yr-1
          reprt()$MaleU,                 # length: sim_yr
          sdreprt$value[c(7,6+sim_yr)],  # length: 2
          sdreprt$value[5],              # length: 1
          sdreprt$value[6]               # length: 1
        )
        est_N[i,] <- as.vector( reprt()$Nt )
      }
    }
  }
  out <- list()
  out$init      <- init
  out$dyn       <- dyn
  out$obj       <- obj
  out$opt       <- opt
  out$reprt     <- reprt()
  out$sdreprt   <- sdreprt
  out$sim       <- sim_pars
  out$est       <- est_pars
  out$diff      <- (est_pars-sim_pars)/sim_pars
  out$sim_info  <- sim_info
  out$est_info  <- est_info
  out$diff_info <- (est_info-sim_info)/sim_info
  out$sim_N     <- sim_N
  out$est_N     <- est_N
  out$diff_N    <- (est_N-sim_N)/sim_N
  return(out)
} # ErrorEval()

"decision.decisions.decisions" <- function( load.prev=TRUE ){
  n.opt    <- 10
  # effort time series options for nest destruction
  Nest.opt <- list(
    scen1 <- rep( 0, nT ),
    scen2 <- c( rep( 0, init_yr ), rep( 2, sim_yr ) ),
    scen3 <- c( rep( 0, init_yr ), rep( 2, sim_yr ) ),
    scen4 <- c( rep( 0, init_yr ), rep( 2, sim_yr ) ),
    scen5 <- c( rep( 0, init_yr ), rep( 2, sim_yr ) ),
    scen6 <- c( rep( 0, init_yr ), rep( 0, 2 ), rep( 2, sim_yr-2 ) ),
    scen7 <- c( rep( 0, init_yr ), rep( 2, 6 ), rep( 0, 2 ), rep( 2, sim_yr-8 ) ),
    scen8 <- c( rep( 0, init_yr ), rep( 2, 2 ), rep( 0, 2 ), rep( 2, 2 ), rep( 0, 2 ), rep( 2, sim_yr-8 ) ),
    scen9 <- c( rep( 0, init_yr ), rep( 2, 3 ), rep( 0, 3 ), rep( 2, sim_yr-6 ) ),
    scen10 <- c( rep( 0, init_yr ), rep( 2, 3 ), 0, rep( 2, 2 ), 0, rep( 2, sim_yr-7 ) )
  )
  # exploitation rate time series options for male removals off nests
  Male.opt <- list(
    scen1 <- rep( 0, nT ),
    scen2 <- c( rep( 0, init_yr ), rep( 0, 3 ), rep( 0.3, 2 ), rep( 0, 2 ), rep( 0.3, sim_yr-7 ) ),
    scen3 <- c( rep( 0, init_yr ), rep( 0.3, 3 ), rep( 0, 6 ), rep( 0.3, sim_yr-9 ) ),
    scen4 <- c( rep( 0, init_yr ), rep( 0, 2 ), rep( 0.3, sim_yr-2 ) ),
    scen5 <- c( rep( 0, init_yr ), rep( 0, 2 ), rep( 0.3, 1 ), rep( 0, 2 ), rep( 0.3, sim_yr-5 ) ),
    scen6 <- c( rep( 0, init_yr ), rep( 0, 2 ), rep( 0.3, sim_yr-2 ) ),
    scen7 <- c( rep( 0, init_yr ), rep( 0, 2 ), rep( 0.3, 2 ), rep( 0, 4 ), rep( 0.3, sim_yr-8 ) ),
    scen8 <- c( rep( 0, init_yr ), rep( 0.3, 2 ), rep( 0, 6 ), rep( 0.3, sim_yr-8 ) ),
    scen9 <- c( rep( 0, init_yr ), rep( 0.3, 3 ), rep( 0, 3 ), rep( 0.3, sim_yr-6 ) ),
    scen10 <- c( rep( 0, init_yr ), rep( 0.3, 3 ), 0, rep( 0.3, 2 ), rep( 0, sim_yr-6 ) )
  )
  # years when regs are tightened
  reg.opt <- list(
    scen1 <- nT+1,
    scen2 <- init_yr+6,
    scen3 <- init_yr+4,
    scen4 <- init_yr+5,
    scen5 <- init_yr+4,
    scen6 <- init_yr+1,
    scen7 <- init_yr+5,
    scen8 <- init_yr+3,
    scen9 <- init_yr+7,
    scen10 <- init_yr+1
  )
  
  # create a decision table to be filled
  dec.table <- data.frame( "Scenario"=1:n.opt, 
                           "RMSE" = NA,
                           "RMSE.Male"=NA,
                           "RMSE.Nest"=NA,
                           "RMSE.NestImpact"=NA,
                           "RMSE.MaleImpact"=NA)
  for( i in 1:n.opt ){
    if( load.prev ){
      load( paste0("Scenario",i,".RData"))
    } else {
      eval <- ErrorEval( n=100, sim=SIM, est=EST, iUMark=UMark,iENest=Nest.opt[[i]],
                         iUMale=Male.opt[[i]],itreg=reg.opt[[i]])
      save( eval, file=paste0( "Scenario",i,".RData"))
    }
    
    # effectiveness of different treatments based on RMSE
    dec.table$RMSE[i]          <- sqrt( mean( (log(eval$est_N+1e-5)-log(eval$sim_N+1e-5))^2, na.rm=TRUE))
      # Male effectiveness estimated from UMale
    iMale <- nA+2*sim_yr-1+1:sim_yr
    dec.table$RMSE.Male[i]        <- sqrt(mean((eval$est_info[,iMale]-
                                          eval$sim_info[,iMale])^2, na.rm=TRUE))
      # Nest effectiveness estimated from recruitment
    iNest <- nA+(1:sim_yr)
    dec.table$RMSE.Nest[i]        <- sqrt(mean((log(eval$est_info[,iNest]+1e-5)-
                                           log(eval$sim_info[,iNest]+1e-5))^2, na.rm=TRUE))
    
    iImpact <- nA+(sim_yr*3-1)+2+1:2
      # Nest impact on recruitment
    dec.table$RMSE.NestImpact[i]  <- sqrt( mean( (eval$est_info[,iImpact[1]]-
                                            eval$sim_info[,iImpact[1]])^2, na.rm=TRUE))
      # Male removal impact on recruitment
    dec.table$RMSE.MaleImpact[i]  <- sqrt( mean( (eval$est_info[,iImpact[2]]-
                                                 eval$sim_info[,iImpact[2]])^2, na.rm=TRUE))
    
  }
  save( dec.table, file="Decision table.RData")
  dec.df <- as.data.frame( dec.table )
}

##################################################################################################
############################################   PLOT   ############################################
##################################################################################################

"Proportional.Error" <- function( scenario, Info=FALSE ){
  par.nam <- c(
    expression(paste("R"[0])),
    expression(kappa),
    "M","q",
    expression(paste(alpha)["Nest"]),
    expression(paste(beta)["Nest"]),
    expression(paste("s"["a,1"])),
    expression(paste("s"["a,2"])),
    expression(paste("s"["a,3"])),
    expression(paste("s"["a,4"])),
    expression(paste("s"["a,5"])),
    expression(paste("s"["a,6+"])),
    expression(paste("s"["EF,1"])),
    expression(paste("s"["EF,2"])),
    expression(paste("s"["EF,3"])),
    expression(paste("s"["EF,4"])),
    expression(paste("s"["EF,5"])),
    expression(paste("s"["EF,6+"])),
    expression(paste("N0"["1"])),
    expression(paste("N0"["2"])),
    expression(paste("N0"["3"])),
    expression(paste("N0"["4"])),
    expression(paste("N0"["5"])),
    expression(paste("N0"["6+"]))
  )
  
  par2.nam <- c(
    expression(paste("NT"["1"])),
    expression(paste("NT"["2"])),
    expression(paste("NT"["3"])),
    expression(paste("NT"["4"])),
    expression(paste("NT"["5"])),
    expression(paste("NT"["6+"])),
    expression(paste("R"["1"])),
    expression(paste("R"["2"])),
    expression(paste("R"["3"])),
    expression(paste("R"["4"])),
    expression(paste("R"["5"])),
    expression(paste("R"["6"])),
    expression(paste("R"["7"])),
    expression(paste("R"["8"])),
    expression(paste("R"["9"])),
    expression(paste("R"["10"])),
    expression(paste("U"["2"]^"EF")),
    expression(paste("U"["3"]^"EF")),
    expression(paste("U"["4"]^"EF")),
    expression(paste("U"["5"]^"EF")),
    expression(paste("U"["6"]^"EF")),
    expression(paste("U"["7"]^"EF")),
    expression(paste("U"["8"]^"EF")),
    expression(paste("U"["9"]^"EF")),
    expression(paste("U"["10"]^"M")),
    expression(paste("U"["1"]^"M")),
    expression(paste("U"["2"]^"M")),
    expression(paste("U"["3"]^"M")),
    expression(paste("U"["4"]^"M")),
    expression(paste("U"["5"]^"M")),
    expression(paste("U"["6"]^"M")),
    expression(paste("U"["7"]^"M")),
    expression(paste("U"["8"]^"M")),
    expression(paste("U"["9"]^"M")),
    expression(paste("U"["10"]^"M")),
    expression(paste(alpha["Nest"])),
    expression(paste(beta["Nest"]))
  )
  
  load( file=paste0("Scenario",scenario,".RData" ))
  diff <- as.data.frame(eval$diff)
  diff <- diff[,-(6+(1:sim_yr))]
  names(diff) <- par.nam
  melt.diff <- melt( diff )
  
  if( !Info ){
    plot.out <- ggplot( data=melt.diff, aes( x=variable, 
                                             y=value) ) +
      geom_boxplot() + 
      theme_classic( base_size=14) +
      xlab("Parameter") + ylab("Proportional error") +
      scale_x_discrete( labels=par.nam) +
      ylim( -1,10 ) +
      geom_hline(yintercept=0, linetype="dashed", 
                 colour="red",size=2)
    pdf( file=paste0(path,"/Figures/Parameter error.pdf"))
      plot.out
    dev.off()
    jpeg( filename = paste0(path,"/Figures/Parameter error.jpg"),
          width = 960, height = 960 )
      plot.out
    dev.off()
  } else {  
    diff2 <- eval$diff_info
    diff2 <- as.data.frame(diff2[,-(nA+sim_yr*3-1+1:2)])
    melt.diff2 <- melt( diff2 )
    plot.out <- ggplot( data=melt.diff2, aes( x=variable, y=value) ) +
      geom_boxplot() +
      theme_classic( base_size=14) +
      xlab("Parameter") + ylab("Proportional error") +
      scale_x_discrete( labels=par2.nam) +
      ylim(-1,5) +
      geom_hline(yintercept=0, linetype="dashed", colour="red",size=2)
    pdf( file=paste0(path,"/Figures/Predictive error.pdf"), width=9)
      plot.out
    dev.off()
    jpeg( filename = paste0(path,"/Figures/Predictive error.jpg"),
          width = 960, height = 960 )
      plot.out
    dev.off()
  }  
}

#eval <- ErrorEval(n=1)
#eval <- ErrorEval(n=100)
#decision.decisions.decisions()