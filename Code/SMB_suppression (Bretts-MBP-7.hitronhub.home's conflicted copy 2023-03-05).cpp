#include <TMB.hpp>
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // READ IN THE DATA
  DATA_INTEGER(AR);                    // age at recruitment
  DATA_INTEGER(A);                     // maximum age (a plus-group)
  DATA_INTEGER(nT);                    // number of years
  DATA_SCALAR(Linf);                   // von Bertalanffy asymptotic length
  DATA_SCALAR(K);                      // von Bertalanffy metabolic growth parameter
//  DATA_SCALAR(CV);                     // cv in length-at-age
  DATA_SCALAR(LWa);                    // scalar of length-weight relationship
  DATA_SCALAR(LWb);                    // exponent of length-weight relationship
  DATA_SCALAR(EggMass);                // eggs per gram of fish
  DATA_SCALAR(pFem);                   // equilibrium and initial proportion females
  DATA_SCALAR(Mat50);                  // age at 50% maturity
  DATA_SCALAR(MatSig);                 // logistic slope in age-at-maturity
  DATA_SCALAR(DisMort);                // discard mortality (CONSIDER ESTIMATING LATER)
  DATA_SCALAR(eps);                    // lower limit at which penalties start applying 
  DATA_SCALAR(proc_sig);               // standard deviation in recruitment process error
  DATA_SCALAR(proc_sig2);              // standard deviation in effort process error
  DATA_VECTOR(pHarv);                  // proportion harvested in status-quo management
  DATA_VECTOR(Age);                    // ages included in the model
  DATA_VECTOR(ENest);                  // number of days devoted to destroying nests each year
  DATA_VECTOR(E_Hat);                  // effort estimated in creel surveys
  DATA_VECTOR(FCap);                   // sum of unmarked females captured in e-fishing
  DATA_VECTOR(MCap);                   // sum of unmarked males captured in e-fishing
  DATA_VECTOR(FRecap);                 // sum of marked females recaptured in e-fishing
  DATA_VECTOR(MRecap);                 // sum of marked males recaptured in e-fishing
  DATA_VECTOR(Catch);                  // total catch of fish (whether harvested or not)
  DATA_VECTOR(ZEst);                   // total mortality of vulnerable fish (M+F+discard mort)
  DATA_ARRAY(FMarks);                  // marked females at start of each year
  DATA_ARRAY(MMarks);                  // marked males at start of each year
  DATA_ARRAY(FRecs);                   // recaped females at start of each year
  DATA_ARRAY(MRecs);                   // recaped males at start of each year
  DATA_ARRAY(CUMal);                   // unmarked males removed off nests
  DATA_ARRAY(CMMal);                   // marked males removed off nests
  DATA_ARRAY(Cr_CU);                   // CPUE of unmarked fish in creel
  DATA_ARRAY(Cr_CM);                   // CPUE of marked fish in creel

  // DECLARE PARAMETERS TO BE ESTIMATED
  PARAMETER(logR0);                    // log of unfished equilibrium recruits
  PARAMETER(logRecK);                  // log of recruitment compensation ratio
  PARAMETER(logM);                     // log of natural mortality rate
  PARAMETER(logq);                     // log of fishery catchability
  PARAMETER(aNest);                    // coefficient of impact of nest destruction
  PARAMETER(bNest);                    // coefficient of impact of male removal
  PARAMETER(logSigma);                 // log of standard deviation around creel catch rates
  PARAMETER_VECTOR(logSel);            // selectivity-at-age on logit-scale
  PARAMETER_VECTOR(logEFishSel);       // elecrofishing selectivity-at-age on logit-scale
  PARAMETER_VECTOR(logInitN);          // log of initial abundance of unmarked fish
  PARAMETER_VECTOR(Eps);               // process error in annual recruitment (random variable)
  PARAMETER_VECTOR(Zeta);              // process error in annual fishing effort (random variable)

  // DECLARE VARIABLES TO BE CALCULATED
  int nA=A-AR+1;
  Type R0=exp(logR0);
  Type RecK=exp(logRecK);
  Type M=exp(logM);
  Type q=exp(logq);

  Type Tiny = 0.00001;
  Type pen=0.0;
  Type phi0;                    // spawners (eggs) per recruit
  Type RecA;                    // Beverton-Holt A parameter
  Type RecB;                    // Beverton-Holt B parameter
  Type FEfishU;                 // exploitation rate of females to electrofishing
  Type MEfishU;                 // exploitation rate of males to electrofishing
  Type Sigma = exp(logSigma);   // standard deviation around creel catch rates
  vector<Type> Lt(nA);          // length-at-age
  vector<Type> Wt(nA);          // weight-at-age
  vector<Type> Mat(nA);         // maturity-at-age
  vector<Type> Fec(nA);         // fecundity-at-age
  vector<Type> lx(nA);          // survivorship-to-age
  vector<Type> Sel(nA);         // selectivity-at-age
  vector<Type> EFishSel(nA);    // electrofishing selectivity-at-age
  vector<Type> InitN(nA);       // abundance in first year
  vector<Type> Eggs(nT);        // eggs created each year
  vector<Type> Nests(nT);       // nests each year 
  vector<Type> ExpMal(nA);      // exploitation rate of males on nests
  vector<Type> St(nT);          // spawner index - predicted eggs with adjustment for nest destruction
  vector<Type> Et_pred(nT);     // predicted annual effort
  vector<Type> Zt(nA);          // total mortality each year
  vector<Type> FVt(nA);         // vulnerable female population to electrofishing
  vector<Type> MVt(nA);         // vulnerable male population to electrofishing
  vector<Type> pcatch(nA);      // proportion of fish caught and harvested in fishery
  vector<Type> tempNt(nA);      // placeholder for total abundance
  vector<Type> Catch_hat(nT);   // total catch (release + harvest ) per year
  vector<Type> Acoust_hat(nT);  // total mortality estimated from acoustic fish
  array<Type> UtF(nA,nT);       // unmarked females in population
  array<Type> UtM(nA,nT);       // unmarked males in population
  array<Type> MtF(nA,nT);       // marked females in population
  array<Type> MtM(nA,nT);       // marked males in population
  array<Type> FMark_hat(nA,nT); // predicted females marked at the start of each year
  array<Type> MMark_hat(nA,nT); // predicted males marked at the start of each year
  array<Type> FRec_hat(nA,nT);  // predicted females recaped at the start of each year
  array<Type> MRec_hat(nA,nT);  // predicted males recaped at the start of each year
  array<Type> CPE_U(nA,nT);     // creel CPUE of unmarked fish
  array<Type> CPE_M(nA,nT);     // creel CPUE of marked fish

  // DECLARE OBJECTIVE FUNCTION
  Type nll = 0.0;
  // INITIALIZE OJIVES AND CALCULATE RECRUITMENT PARAMETERS
  Lt = Linf * ( Type(1.) - exp( -K*Age ));
  Wt = LWa * pow( Lt, LWb );
  Mat = Type(1.) / ( Type(1.) + exp( -( Age - Mat50 )/MatSig ));
  Fec = Wt * EggMass;
  lx = exp( -M*(Age-AR) );
  lx(nA-AR) /= Type(1.)-exp(-M);
//  std::cout <<"Length\n"<<Lt<<"\n";
//  std::cout <<"Weight\n"<<Wt<<"\n";
//  std::cout <<"Mat\n"<<Mat<<"\n";
//  std::cout <<"Fec\n"<<Fec<<"\n";
//  std::cout <<"lx\n"<<lx<<"\n";
  for( int a=0; a<(nA-1); a++){
    Sel(a)=exp(logSel(a));
    EFishSel(a)=exp(logEFishSel(a));
  }
  Sel(nA-1)=Type(1.);
  EFishSel(nA-1)=Type(1.);
//  std::cout <<"Sel: "<<Sel<<"\n"<<"EFishSel\t"<<EFishSel<<"\n";
  
  phi0 = ( Mat * Fec * lx ).sum() * pFem;
  RecA = RecK / phi0;
  RecB = log( RecK )/( R0 * phi0 );
//  std::cout <<"phi0: "<<phi0<<"\tR0: "<<R0<<"\tRecK: "<<RecK<<"\n";
//  std::cout <<"RecA: "<<RecA<<"\tRecB: "<<RecB<<"\n";
  
  // FIRST-YEAR CONDITIONS
  // initialize abundance
  InitN = exp( logInitN );
  MtF.col(0) = FMarks.col(0);
  MtM.col(0) = MMarks.col(0);
  // InitN is estimated as # + 1, which should keep the estimated # positive
  for( int a=0; a<nA; a++ ){
    UtF(a,0) = posfun( InitN(a) - Type(1.), eps, pen ) * pFem;             
    UtM(a,0) = posfun( InitN(a) - Type(1.), eps, pen ) * (Type(1.)-pFem);
  } // a
//  std::cout<<"Marks\n"<<FMarks.col(0)<<"\n"<<MMarks.col(0)<<"\n"<<"abund\n"<<UtF.col(0)<<"\n"<<UtM.col(0)<<"\n"<<MtF.col(0)<<"\n"<<MtM.col(0)<<"\n";
  // note: don't predict catches in first year, because you aren't estimating initial pop abund.

  // index of nest removals is calculated as a measure of total potential impact on eggs and recruitment
  Nests(0) = ( (CUMal.col(0) + CMMal.col(0)) / (UtM.col(0)+MtM.col(0)) * Mat ).sum();  // assume nests are proportional to the number of males

  // calculate reproductive output
  Eggs(0) = ( ( UtF.col(0)+MtF.col(0) ) * Fec * Mat ).sum();
  St.col(0) = Eggs(0) * exp( aNest*ENest(0) + bNest*Nests(0) );
//  ExpMal = ( CUMal.col(0) + CMMal.col(0) ) / ( UtM.col(0) + MtM.col(0) + Tiny );
  //  std::cout<<"init\t"<<Eggs(0)<<"\n";

  // calculate fishing impact (predicted HPUE)
//  std::cout<<"E_Hat\n"<<E_Hat<<"\n"<<"FCap\n"<<FCap<<"\n";
  Et_pred(0) = E_Hat(0) * (Type(1.) + Zeta(0));
  Zt = M + q*Et_pred(0) * Sel * ( pHarv(0) + ( Type(1.) - pHarv(0) )*DisMort ); // total mortality
  pcatch = q * Sel * pHarv(0) / Zt * ( Type(1.) - exp(-Zt) );
//  std::cout<<"Sel\n"<<Sel<<"\nZ\n"<<Zt<<"\nE_Hat\n"<<E_Hat(0)<<"\npcatch\n"<<pcatch<<"\nEFishSel\n"<<EFishSel<<"\n";

  // total catch
  tempNt = ( UtF.col(0) + UtM.col(0) + MtF.col(0) + MtM.col(0) - CUMal.col(0) - CMMal.col(0) );
  Catch_hat(0) = ( tempNt * q * Sel / Zt * ( Type(1.) - exp(-Zt) )).sum();
//  std::cout<<"time 1: Nt="<<tempNt<<"\tCatch="<<Catch_hat(0)<<"\n";
  
  // total mortality of vulnerable fish as estimated from acoustically tagged fish
  Acoust_hat(0) = M + q*Et_pred(0) * ( pHarv(0) + ( Type(1.) - pHarv(0) )*DisMort );
//  std::cout<<"t: "<<0<<" M: "<<M<<" q: "<<q<<" Et: "<<E_Hat(0)<<" Zeta: "<<Zeta(0)<<" Epred: "<<Et_pred(0)<<" pHarv: "<<pHarv(0)<<"\n";
  
  // harvest per unit effort
  CPE_U.col(0) = ( UtF.col(0) + UtM.col(0) - CUMal.col(0) ) * pcatch;  // Baranov eq'n for HPUE
  CPE_M.col(0) = ( MtF.col(0) + MtM.col(0) - CMMal.col(0) ) * pcatch;  // Baranov eq'n for HPUE
//  std::cout<<"CPUE\n"<<CPE_UF.col(0)<<"\n"<<CPE_UM.col(0)<<"\n"<<CPE_MF.col(0)<<"\n"<<CPE_MM.col(0)<<"\n";

  // PROJECT FOR FUTURE YEARS
  for( int t=1; t<nT; t++ ){
    // establish recruitment into the different sex/marking groups
    UtF(AR-1,t) = posfun( St(t-1)*RecA * exp(-St(t-1)*RecB + Eps(t-1) ) * pFem - FMarks(AR-1,t), eps, pen );
    UtM(AR-1,t) = posfun( St(t-1)*RecA * exp(-St(t-1)*RecB + Eps(t-1) ) * ( Type(1.)-pFem ) - MMarks(AR-1,t), eps, pen);
    MtF(AR-1,t) = FMarks(AR-1,t);
    MtM(AR-1,t) = MMarks(AR-1,t);
    for( int a=1; a<nA; a++ ){  // loop from age-2 to oldest age
      UtF(a,t) = UtF(a-1,t-1) * exp( -Zt(a-1) );
      UtM(a,t) = posfun( UtM(a-1,t-1) - CUMal(a-1,t-1), eps, pen ) * exp( -Zt(a-1) );
      MtF(a,t) = MtF(a-1,t-1) * exp( -Zt(a-1) );
      MtM(a,t) = posfun( MtM(a-1,t-1) - CMMal(a-1,t-1), eps, pen ) * exp( -Zt(a-1) );
    } // a
     // adjust the plus-group
    UtF(nA-1,t) += UtF(nA-1,t-1)*exp( -Zt[nA-1] );
    UtM(nA-1,t) += posfun( UtM(nA-1,t-1) - CUMal(nA-1,t-1), eps, pen )*exp(-Zt[nA-1]);
    MtF(nA-1,t) += MtF(nA-1,t-1)*exp( -Zt[nA-1] );
    MtM(nA-1,t) += posfun( MtM(nA-1,t-1) - CMMal(nA-1,t-1), eps, pen )*exp(-Zt[nA-1]);

    for( int a=1; a<nA; a++ ){
      // redistribute marks at start of year
      UtF(a,t) = posfun( UtF(a,t) - FMarks(a,t), eps, pen );
      UtM(a,t) = posfun( UtM(a,t) - MMarks(a,t), eps, pen );
      MtF(a,t) += FMarks(a,t);
      MtM(a,t) += MMarks(a,t);
    }
//    std::cout<<"UtF\n"<<UtF.col(t)<<"\nUtM\n"<<UtM.col(t)<<"\nMtF\n"<<MtF.col(t)<<"\nMtM\n"<<MtM.col(t)<<"\n";

    // predicted e-fishing catch
    FVt = ( UtF.col(t) + MtF.col(t) ) * EFishSel;
    MVt = ( UtM.col(t) + MtM.col(t) ) * EFishSel;
    FEfishU = ( FCap(t) + FRecap(t) ) / ( FVt ).sum();
    MEfishU = ( MCap(t) + MRecap(t) ) / ( MVt ).sum();
//    std::cout<<"FEfishU: "<<FEfishU<<"\nMEfishU: "<<MEfishU<<"\n";
     
    FMark_hat.col(t) = UtF.col(t) * EFishSel * FEfishU;
    MMark_hat.col(t) = UtM.col(t) * EFishSel * MEfishU;
    FRec_hat.col(t) = MtF.col(t) * EFishSel * FEfishU;
    MRec_hat.col(t) = MtM.col(t) * EFishSel * MEfishU;
//    std::cout<<"UtF\n"<<UtF.col(t)<<"\nEFishSel\n"<<EFishSel<<"\nFEfishU\n"<<FEfishU<<"\n";
//    std::cout<<"FMark\n"<<FMark_hat.col(t)<<"\nMMark\n"<<MMark_hat.col(t)<<"\nFRec\n"<<FRec_hat.col(t)<<"\nMRec\n"<<MRec_hat.col(t)<<"\n";
//    std::cout<<"t\t"<<t<<"\nUtF\n"<<UtF.col(t)<<"\nUtM\n"<<UtM.col(t)<<"\nMtF\n"<<MtF.col(t)<<"\nMtM\n"<<MtM.col(t)<<"\n";
    
    // index of nests is calculated prior to destruction as a measure of total potential to be destroyed
    Nests(t) = ( (CUMal.col(t) + CMMal.col(t)) / (UtM.col(t)+MtM.col(t)) * Mat ).sum();  // assume nests are proportional to the number of males
    
    // calculate reproductive output
    Eggs(t) = ( ( UtF.col(t)+MtF.col(t) ) * Fec * Mat ).sum();
    St(t) = Eggs(t) * exp( aNest*ENest(t) + bNest*Nests(t) );
//    ExpMal = ( CUMal.col(t) + CMMal.col(t) ) / ( UtM.col(t) + MtM.col(t) + Tiny );
    //    std::cout<<"Nests: "<<Nests(t)<<"\nEggs: "<<Eggs(t)<<"\nSt: "<<St(t)<<"\n";
    
    // calculate fishing impact (predicted HPUE)
    Et_pred(t) = E_Hat(t) * (Type(1.) + Zeta(t));
    Zt = M + q*Et_pred(t) * ( pHarv(t) + ( Type(1.) - pHarv(t) )*DisMort ); // total mortality
    pcatch = q*Sel*pHarv(t) / Zt * ( Type(1.)-exp(-Zt) );  
    // pcatch omits effort, so when multiplied by abundance, gives HPUE

    // total catch
    tempNt = ( UtF.col(t) + UtM.col(t) + MtF.col(t) + MtM.col(t) - CUMal.col(t) - CMMal.col(t) );
    Catch_hat(t) = ( tempNt * q * Sel / Zt * ( Type(1.) - exp(-Zt) )).sum();
    
    // total mortality of vulnerable fish as estimated from acoustically tagged fish
    Acoust_hat(t) = M + q*Et_pred(t) * ( pHarv(t) + ( Type(1.) - pHarv(t) )*DisMort );
//    std::cout<<"t: "<<t<<" M: "<<M<<" q: "<<q<<" Et: "<<E_Hat(t)<<" Zeta: "<<Zeta(t)<<" Epred: "<<Et_pred(t)<<" pHarv: "<<pHarv(t)<<"\n";
    
    CPE_U.col(t) = ( UtF.col(t) + UtM.col(t) - CUMal.col(t) ) * pcatch;  // Baranov eq'n for HPUE
    CPE_M.col(t) = ( MtF.col(t) + MtM.col(t) - CMMal.col(t) ) * pcatch;
//    std::cout<<"Zt\n"<<Zt<<"\npcatch\n"<<pcatch<<"\nUF\n"<<CPE_UF.col(t)<<"\nUM\n"<<CPE_UM.col(t)<<"\nMF\n"<<CPE_MF.col(t)<<"\nMM\n"<<CPE_MM.col(t)<<"\n";
  } // t

  // CALCULATE LIKELIHOODS
  Type EFnll = 0.0;      // neg log likelihood for electrofishing data
  Type CRnll = 0.0;      // neg log likelihood for creel estimates of HPUE
  Type Cnll  = 0.0;      // neg log likelihood for creel estimates of CPUE
  Type Znll  = 0.0;      // neg log likelihood for acoustic estimates of Z
  Type REnll = 0.0;      // neg log likelihood for random effect

  // fit to fishery CPUE
  Cnll -= dnorm( log( Catch(0)+Tiny ), log( Catch_hat(0)+Tiny ), Sigma_CPUE, true );
  
  // fit to acoustic estimates of ZEst
  Znll -= dnorm( log( ZEst(0) ), log( Acoust_hat(0) ), Sigma_Z, true );
  
  for( int t=1; t<nT; t++ ){  // start at second time-step since no FMark_hat calc'd in 1st
    for( int a=0; a<nA; a++ ){
      // fit to electrofishing catch
      EFnll -= dpois( FMarks(a,t), FMark_hat(a,t)+Tiny, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tFMarks\t"<<FMarks(a,t)<<"\t"<<FMark_hat(a,t)<<"\t"<<EFnll<<"\n";
      EFnll -= dpois( MMarks(a,t), MMark_hat(a,t)+Tiny, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tMMarks\t"<<MMarks(a,t)<<"\t"<<MMark_hat(a,t)<<"\t"<<EFnll<<"\n";
      EFnll -= dpois( FRecs(a,t), FRec_hat(a,t)+Tiny, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tFRecs2\t"<<FRecs(a,t)<<"\t"<<FRec_hat(a,t)<<"\t"<<EFnll<<"\n";
      EFnll -= dpois( MRecs(a,t), MRec_hat(a,t)+Tiny, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tMRecs\t"<<MRecs(a,t)<<"\t"<<MRec_hat(a,t)<<"\t"<<EFnll<<"\n";
      
      // fit to fishery HPUE
      CRnll -= dnorm( log( Cr_CU(a,t)+Tiny ), log( CPE_U(a,t)+Tiny ), Sigma, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tCRU HPUE obs: "<<Cr_CU(a,t)<<"\tCRU HPUE pred: "<<CPE_U(a,t)<<"\n";
      CRnll -= dnorm( log( Cr_CM(a,t)+Tiny ), log( CPE_M(a,t)+Tiny ), Sigma, true );
//      std::cout<<"a "<<a<<" t "<<t<<"\tCRM HPUE obs: "<<Cr_CM(a,t)<<"\tCRU HPUE pred: "<<CPE_M(a,t)<<"\n";
    }
    // fit to fishery CPUE
    Cnll -= dnorm( log( Catch(t)+Tiny ), log( Catch_hat(t)+Tiny ), Sigma_CPUE, true );
//    std::cout<<"CRU CPUE obs: "<<Catch(t)<<"\tCRM CPUE pred: "<<Catch_hat(t)<<"\n";
    
    // fit to acoustic estimates of ZEst
    Znll -= dnorm( log( ZEst(t) ), log( Acoust_hat(t) ), Sigma_Z, true );
//    std::cout<<"obs Z: "<<ZEst(t)<<"\test Z: "<<Acoust_hat(t)<<"\n";
  }

  for( int t=0; t<(nT-1); t++ ){
    REnll -= ( dnorm( Eps(t), Type(0.0), proc_sig+Tiny, true ) );
//    std::cout<<"t "<<t<<"\t"<<Eps(t)<<"\t"<<REnll<<"\n";
  }
  
  nll = EFnll + REnll + CRnll + Cnll + Znll;
  nll += pen;
//  std::cout<<"nll\t"<<EFnll<<"\t"<<REnll<<"\t"<<pen<<"\t"<<nll<<"\n";

  Type meanM=log(0.25);
  Type sdM=0.25;
//  nll -= dnorm( meanM, logM, sdM, true );

  // SPIT OUT VARIABLES
  ADREPORT(R0);                    // unfished equilibrium recruits
  ADREPORT(RecK);                  // recruitment compensation ratio
  ADREPORT(M);                     // natural mortality rate
  ADREPORT(q);                     // fishery catchability
  ADREPORT(Sigma);                 // standard deviation around creel catch rates
  ADREPORT(Sel);                   // selectivity-at-age
  ADREPORT(EFishSel);              // elecrofishing selectivity-at-age
  ADREPORT(InitN);                 // initial abundance of unmarked fish

  REPORT(FMarks);
  REPORT(FMark_hat);
  REPORT(MMarks);
  REPORT(MMark_hat);
  REPORT(FRecs);
  REPORT(FRec_hat);
  REPORT(MRecs);
  REPORT(MRec_hat);

  return(nll);
}
