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
  DATA_INTEGER(regyr);                 // year when regs are tightened
  DATA_SCALAR(Linf);                   // von Bertalanffy asymptotic length
  DATA_SCALAR(K);                      // von Bertalanffy metabolic growth parameter
  DATA_SCALAR(LWa);                    // scalar of length-weight relationship
  DATA_SCALAR(LWb);                    // exponent of length-weight relationship
  DATA_SCALAR(EggMass);                // eggs per gram of fish
  DATA_SCALAR(pFem);                   // equilibrium and initial proportion females
  DATA_SCALAR(Mat50);                  // age at 50% maturity
  DATA_SCALAR(MatSig);                 // logistic slope in age-at-maturity
  DATA_SCALAR(eps);                    // lower limit at which penalties start applying 
  DATA_SCALAR(DisMort);                // discard mortality
  DATA_SCALAR(proc_sig);               // standard deviation in recruitment process error
  DATA_SCALAR(proc_sig2);              // standard deviation in effort process error
  DATA_VECTOR(pHarv);                  // proportion harvested in status-quo management
  DATA_VECTOR(Age);                    // ages included in the model
  DATA_VECTOR(ENest);                  // number of days devoted to destroying nests each year
  DATA_VECTOR(E_Hat);                  // effort estimated in creel surveys
  DATA_VECTOR(Cap);                    // sum of unmarked fish captured in e-fishing
  DATA_VECTOR(Recap);                  // sum of unmarked fish captured in e-fishing
  DATA_VECTOR(Catch);                  // total catch of fish (whether harvested or not)
  DATA_VECTOR(ZEst);                   // total mortality of vulnerable fish (M+F+discard mort)
  DATA_VECTOR(Cr_CU);                  // CPUE of unmarked fish in creel
  DATA_VECTOR(Cr_CM);                  // CPUE of marked fish in creel
  DATA_ARRAY(Marks);                   // marked fish at start of each year
  DATA_ARRAY(Recs);                    // recaped fish at start of each year
  DATA_ARRAY(CUMal);                   // unmarked males removed off nests
  DATA_ARRAY(CMMal);                   // marked males removed off nests

  // DECLARE PARAMETERS TO BE ESTIMATED
  PARAMETER(logR0);                    // log of unfished equilibrium recruits
  PARAMETER(logRecK);                  // log of recruitment compensation ratio
  PARAMETER(logM);                     // log of natural mortality rate
  PARAMETER(logq);                     // log of fishery catchability
  PARAMETER(aNest);                    // coefficient of impact of nest destruction
  PARAMETER(bNest);                    // coefficient of impact of male removal
  PARAMETER(logSigma);                 // log of standard deviation around creel harvest rates
  PARAMETER(logSigma_CPUE);            // log of standard deviation around creel catch rates
  PARAMETER(logSigma_Z);               // log of standard deviation around acoustic tag total mortality
  PARAMETER_VECTOR(logitSel);          // selectivity-at-age on logit-scale
  PARAMETER_VECTOR(logitEFishSel);     // elecrofishing selectivity-at-age on logit-scale
  PARAMETER_VECTOR(logInitN);          // log of initial abundance of unmarked fish
  PARAMETER_VECTOR(Eps);               // process error in annual recruitment (random variable)
  PARAMETER_VECTOR(Zeta);              // process error in annual fishing effort (random variable)

  // DECLARE VARIABLES TO BE CALCULATED
  int nA=A-AR+1;
  Type R0=exp(logR0);
  Type RecK=exp(logRecK)+Type(1.);
  Type M=exp(logM);
  Type q=exp(logq);

  Type Tiny = 0.00001;
  Type pen=0.0;
  Type phi0;                    // spawners (eggs) per recruit
  Type RecA;                    // Beverton-Holt A parameter
  Type RecB;                    // Beverton-Holt B parameter
  Type Sigma = exp(logSigma);   // standard deviation around creel catch rates
  Type Sigma_CPUE = exp(logSigma_CPUE);
  Type Sigma_Z = exp(logSigma_Z);
  vector<Type> Lt(nA);          // length-at-age
  vector<Type> Wt(nA);          // weight-at-age
  vector<Type> Mat(nA);         // maturity-at-age
  vector<Type> Fec(nA);         // fecundity-at-age
  vector<Type> lx(nA);          // survivorship-to-age
  vector<Type> Sel(nA);         // selectivity-at-age
  vector<Type> EfishU(nT);      // exploitation rate of fish to electrofishing
  vector<Type> MaleU(nT);       // exploitation rate of males off nests
  vector<Type> EFishSel(nA);    // electrofishing selectivity-at-age
  vector<Type> InitN(nA);       // abundance in first year
  vector<Type> Eggs(nT);        // eggs created each year
  vector<Type> St(nT);          // spawner index - predicted eggs with adjustment for nest destruction
  vector<Type> Et_pred(nT);     // predicted annual effort
  vector<Type> Zt(nA);          // total mortality each year
  vector<Type> Vt(nA);          // vulnerable population to electrofishing
  vector<Type> pcatch(nA);      // proportion of fish caught and harvested in fishery
  vector<Type> tempNt(nA);      // placeholder for total abundance
  vector<Type> Catch_hat(nT);   // total catch (release + harvest ) per year
  vector<Type> Acoust_hat(nT);  // total mortality estimated from acoustic fish
  vector<Type> FMark_hat(nA);   // predicted females marked at the start of each year
  vector<Type> MMark_hat(nA);   // predicted males marked at the start of each year
  array<Type> UtF(nA,nT);       // unmarked females in population
  array<Type> UtM(nA,nT);       // unmarked males in population
  array<Type> MtF(nA,nT);       // marked females in population
  array<Type> MtM(nA,nT);       // marked males in population
  array<Type> Nt(nA,nT);        // total population for reporting
  array<Type> Mark_hat(nA,nT); // predicted fish marked at the start of each year
  array<Type> Rec_hat(nA,nT);  // predicted fish recaped at the start of each year
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
    //std::cout <<"Length\n"<<Lt<<"\n";
    //std::cout <<"Weight\n"<<Wt<<"\n";
    //std::cout <<"Mat\n"<<Mat<<"\n";
    //std::cout <<"Fec\n"<<Fec<<"\n";
    //std::cout <<"lx\n"<<lx<<"\n";

  Sel = Type(1.)/(Type(1.)+exp( -logitSel ));
//  Sel /= max( Sel );
  EFishSel = Type(1.)/(Type(1.)+exp( -logitEFishSel ));
//  EFishSel /= max( EFishSel );
    //std::cout<<"Sel\n"<<Sel<<"\nEFishSel\n"<<EFishSel<<"\n";

  phi0 = ( Mat * Fec * lx ).sum() * pFem;
  RecA = RecK / phi0;
  RecB = log( RecK )/( R0 * phi0 );
    //std::cout <<"phi0: "<<phi0<<"\tR0: "<<R0<<"\tRecK: "<<RecK<<"\n";
    //std::cout <<"RecA: "<<RecA<<"\tRecB: "<<RecB<<"\n";
  
  // FIRST-YEAR CONDITIONS
  // initialize abundance
  InitN = exp( logInitN );
  MtF.col(0) = Marks.col(0) * pFem;
  MtM.col(0) = Marks.col(0) * (Type(1.)-pFem);
  // InitN is estimated as # + 1, which should keep the estimated # positive
  for( int a=0; a<nA; a++ ){
    UtF(a,0) = InitN(a) * pFem;             
    UtM(a,0) = InitN(a) * (Type(1.)-pFem);
  } // a
    //std::cout<<"t\t"<<Type(0)<<"\nUtF\n"<<UtF.col(0)<<"\nUtM\n"<<UtM.col(0)<<"\nMtF\n"<<MtF.col(0)<<"\nMtM\n"<<MtM.col(0)<<"\n";
    //std::cout<<"Marks\n"<<FMarks.col(0)<<"\n"<<MMarks.col(0)<<"\n"<<"abund\n"<<UtF.col(0)<<"\n"<<UtM.col(0)<<"\n"<<MtF.col(0)<<"\n"<<MtM.col(0)<<"\n";
  // note: don't predict catches in first year, because you aren't estimating initial pop abund.

  // index of nest removals is calculated as a measure of total potential impact on eggs and recruitment
  MaleU(0) = ( CUMal.col(0) + CMMal.col(0) ).sum();
  MaleU(0) /= (( UtM.col(0) + MtM.col(0) ) * Mat ).sum();

  // calculate reproductive output
  Eggs(0) = ( ( UtF.col(0)+MtF.col(0) ) * Fec * Mat ).sum();
  St.col(0) = Eggs(0) * exp( -(aNest*ENest(0) + bNest*MaleU(0)) );
    //std::cout<<"init\t"<<Eggs(0)<<"\n";

  // calculate fishing impact (predicted HPUE)
    //std::cout<<"E_Hat\n"<<E_Hat<<"\n"<<"FCap\n"<<FCap<<"\n";
  Et_pred(0) = E_Hat(0) * (Type(1.) + Zeta(0));
  Zt = M + q*Et_pred(0) * Sel * ( pHarv(0) + ( Type(1.) - pHarv(0) )*DisMort ); // total mortality
  pcatch = q * Sel * pHarv(0) / Zt * ( Type(1.) - exp(-Zt) );
    //std::cout<<"Sel\n"<<Sel<<"\nZ\n"<<Zt<<"\nE_Hat\n"<<E_Hat(0)<<"\npcatch\n"<<pcatch<<"\nEFishSel\n"<<EFishSel<<"\n";

  // total catch
  for( int a=0; a<nA; a++ ){
    tempNt(a) = posfun( UtF(a,0) + UtM(a,0) + MtF(a,0) + MtM(a,0) - CUMal(a,0) - CMMal(a,0), eps, pen );
  }
  Catch_hat(0) = ( tempNt * q * Sel / Zt * ( Type(1.) - exp(-Zt) )).sum();
    //std::cout<<"t: "<<0<<"\t Catch_hat\n"<<Catch_hat(0)<<"\nCatch\n"<<Catch(0)<<"\ntempNt: "<<tempNt<<"\nq: "<<q<<"\nZt: "<<Zt<<"\n";
  
  // total mortality of vulnerable fish as estimated from acoustically tagged fish
  Acoust_hat(0) = M + q*Et_pred(0) * ( pHarv(0) + ( Type(1.) - pHarv(0) )*DisMort );
    //std::cout<<"t: "<<0<<" M: "<<M<<" q: "<<q<<" Et: "<<E_Hat(0)<<" Zeta: "<<Zeta(0)<<" Epred: "<<Et_pred(0)<<" pHarv: "<<pHarv(0)<<"\n";
  
  // harvest per unit effort
  for( int a=0; a<nA; a++ ){
    CPE_U(a,0) = posfun( UtF(a,0) + UtM(a,0) - CUMal(a,0), eps, pen ) * pcatch(a);  // Baranov eq'n for HPUE
    CPE_M(a,0) = posfun( MtF(a,0) + MtM(a,0) - CMMal(a,0), eps, pen ) * pcatch(a);  // Baranov eq'n for HPUE
  }
    //std::cout<<"CPUE\n"<<CPE_UF.col(0)<<"\n"<<CPE_UM.col(0)<<"\n"<<CPE_MF.col(0)<<"\n"<<CPE_MM.col(0)<<"\n";

  // PROJECT FOR FUTURE YEARS
  for( int t=1; t<nT; t++ ){
    // establish recruitment into the different sex/marking groups
      //std::cout<<"t: "<<t<<" Eggs: "<<St(t-1)<<"\n";
    UtF(0,t) = (St(t-1)*RecA * exp(-St(t-1)*RecB + Eps(t-1))) * pFem;
    UtM(0,t) = (St(t-1)*RecA * exp(-St(t-1)*RecB + Eps(t-1))) * ( Type(1.)-pFem );
    MtF(0,t) = Type(0.);
    MtM(0,t) = Type(0.);
    for( int a=1; a<nA; a++ ){  // loop from age-2 to oldest age
      UtF(a,t) = UtF(a-1,t-1) * exp( -Zt(a-1) );
      UtM(a,t) = posfun( UtM(a-1,t-1) - CUMal(a-1,t-1), eps, pen ) * exp( -Zt(a-1) );
      MtF(a,t) = MtF(a-1,t-1) * exp( -Zt(a-1) );
      MtM(a,t) = posfun( MtM(a-1,t-1) - CMMal(a-1,t-1), eps, pen ) * exp( -Zt(a-1) );
    } // a
    // adjust the plus-group
    UtF(nA-1,t) += UtF(nA-1,t-1) * exp( -Zt[nA-1] );
    UtM(nA-1,t) += posfun( UtM(nA-1,t-1) - CUMal(nA-1,t-1), eps, pen ) * exp(-Zt[nA-1]);
    MtF(nA-1,t) += MtF(nA-1,t-1) * exp( -Zt[nA-1] );
    MtM(nA-1,t) += posfun( MtM(nA-1,t-1) - CMMal(nA-1,t-1), eps, pen ) * exp(-Zt[nA-1]);

    // predicted e-fishing catch
    Vt = ( UtF.col(t) + UtM.col(t) + MtF.col(t) + MtM.col(t) ) * EFishSel;
    EfishU(t) = ( Cap(t) + Recap(t) ) / ( Vt ).sum();
      //std::cout<<"t: "<<t<<" Vt: "<<Vt<<"\n"<<UtF.col(t)<<"\n"<<UtM.col(t)<<"\n";
    
    FMark_hat = UtF.col(t) * EFishSel * EfishU(t);
    MMark_hat = UtM.col(t) * EFishSel * EfishU(t);
      //std::cout<<"FMark\n"<<FMark_hat<<"\nMMark\n"<<MMark_hat<<"\n";

    Mark_hat.col(t) = FMark_hat + MMark_hat;
    Rec_hat.col(t) = ( MtF.col(t) + MtM.col(t) ) * EFishSel * EfishU(t);
      //std::cout<<"Mark\n"<<Mark_hat.col(t)<<"\nRec\n"<<Rec_hat.col(t)<<"\nEfishU\t"<<EfishU(t)<<"\n";
      //std::cout<<"t\t"<<t<<"\nUtF\n"<<UtF.col(t)<<"\nUtM\n"<<UtM.col(t)<<"\nMtF\n"<<MtF.col(t)<<"\nMtM\n"<<MtM.col(t)<<"\n";
    
    for( int a=0; a<nA; a++ ){
      UtF(a,t) = posfun( UtF(a,t)-FMark_hat(a), eps, pen );
      UtM(a,t) = posfun( UtM(a,t)-MMark_hat(a), eps, pen );
      MtF(a,t) = posfun( MtF(a,t)+FMark_hat(a), eps, pen );
      MtM(a,t) = posfun( MtM(a,t)+MMark_hat(a), eps, pen );
    }
      //std::cout<<"try 2\n"<<Type(1.) - EFishSel * EfishU(t)<<"\n"<<UtM.col(t)<<"\nZ\n"<<Zt<<"\n";
    
    // index of nests is calculated prior to destruction as a measure of total potential to be destroyed
    MaleU(t) = ( CUMal.col(t) + CMMal.col(t) ).sum();
    MaleU(t) /= (( UtM.col(t) + MtM.col(t) ) * Mat ).sum();

    // calculate reproductive output
    Eggs(t) = ( ( UtF.col(t)+MtF.col(t) ) * Fec * Mat ).sum();
    St(t) = Eggs(t) * exp( -(aNest*ENest(t) + bNest*MaleU(t)) );
      //std::cout<<"Nests: "<<Nests(t)<<"\nEggs: "<<Eggs(t)<<"\nSt: "<<St(t)<<"\n";
    
    // calculate fishing impact (predicted HPUE)
    Et_pred(t) = E_Hat(t) * (Type(1.) + Zeta(t));
    Zt = M + q*Et_pred(t) * Sel * ( pHarv(t) + ( Type(1.) - pHarv(t) )*DisMort ); // total mortality
    pcatch = q*Sel*pHarv(t) / Zt * ( Type(1.)-exp(-Zt) );  
    // pcatch omits effort, so when multiplied by abundance, gives HPUE

    // total catch rate (harvest and release)
    for( int a=0; a<nA; a++ ){
      tempNt(a) = posfun( UtF(a,t) + UtM(a,t) + MtF(a,t) + MtM(a,t) - CUMal(a,t) - CMMal(a,t) , eps, pen );
    }
    Catch_hat(t) = ( tempNt * q * Sel / Zt * ( Type(1.) - exp(-Zt) )).sum();
      //std::cout<<"t: "<<t<<"\t Catch_hat\n"<<Catch_hat(t)<<"\nCatch\n"<<Catch(t)<<"\ntempNt: "<<tempNt<<"\nq: "<<q<<"\nSel\n"<<Sel<<"\nZt: "<<Zt<<"\n";
    
    // total mortality of vulnerable fish as estimated from acoustically tagged fish
    Acoust_hat(t) = M + q*Et_pred(t) * ( pHarv(t) + ( Type(1.) - pHarv(t) )*DisMort );
      //std::cout<<"t: "<<t<<" M: "<<M<<" q: "<<q<<" Et: "<<E_Hat(t)<<" Zeta: "<<Zeta(t)<<" Epred: "<<Et_pred(t)<<" pHarv: "<<pHarv(t)<<"\n";
    
    for( int a=0; a<nA; a++ ){
      CPE_U(a,t) = posfun( UtF(a,t) + UtM(a,t) - CUMal(a,t), eps, pen ) * pcatch(a);  // Baranov eq'n for HPUE
      CPE_M(a,t) = posfun( MtF(a,t) + MtM(a,t) - CMMal(a,t), eps, pen ) * pcatch(a);
    }
      //std::cout<<"t\t"<<t<<"\tpcatch\n"<<pcatch<<"\nCPUEU\n"<<CPE_U.col(t)<<"\nUtF\n"<<UtF.col(t)<<"\nUtM\n"<<UtM.col(t)<<"\n";
      //std::cout<<"CPUEU\n"<<CPE_M.col(t)<<"\nMtF\n"<<MtF.col(t)<<"\nMtM\n"<<MtM.col(t)<<"\nEggs\n"<<Eggs(t)<<"\n";
      //std::cout<<"M: "<<M<<"\tq: "<<q<<"Et_pred: "<<Et_pred(t)<<"\nSel\n"<<Sel<<"\n";
  } // t
  Nt = UtM + UtF + MtM + MtF;
    //std::cout<<"Nt\n"<<Nt<<"\n";

  // CALCULATE LIKELIHOODS
  Type EFnll = 0.0;      // neg log likelihood for electrofishing data
  Type CRnll = 0.0;      // neg log likelihood for creel estimates of HPUE
  Type Cnll  = 0.0;      // neg log likelihood for creel estimates of CPUE
  Type Znll  = 0.0;      // neg log likelihood for acoustic estimates of Z
  Type REnll = 0.0;      // neg log likelihood for random effect

  for( int t=1; t<nT; t++ ){  // start at second time-step since no FMark_hat calc'd in 1st
    for( int a=0; a<nA; a++ ){
      // fit to electrofishing catch
      if( !(Marks(a,t) == 0 & Mark_hat(a,t) == 0 )){
        EFnll -= dpois( Marks(a,t), Mark_hat(a,t), true );
      }
        //std::cout<<"a "<<a<<" t "<<t<<"\tMarks\t"<<Marks(a,t)<<"\t"<<Mark_hat(a,t)<<"\t"<<EFnll<<"\n";
      if( !(Recs(a,t) == 0 & Rec_hat(a,t) == 0 )){
        EFnll -= dpois( Recs(a,t), Rec_hat(a,t), true );
      }
        //std::cout<<"a "<<a<<" t "<<t<<"\tRecs2\t"<<Recs(a,t)<<"\t"<<Rec_hat(a,t)<<"\t"<<EFnll<<"\n";
    }
  }

  for( int t=0; t<nT; t++ ){
    // fit to fishery HPUE
    CRnll -= dnorm( log( Cr_CU(t)+Tiny ), log( CPE_U(t)+Tiny ), Sigma, true );
      //std::cout<<" t "<<t<<"\tCRU HPUE obs: "<<Cr_CU(t)<<"\tCRU HPUE pred: "<<CPE_U(t)<<"\n";
    CRnll -= dnorm( log( Cr_CM(t)+Tiny ), log( CPE_M(t)+Tiny ), Sigma, true );
      //std::cout<<" t "<<t<<"\tCRM HPUE obs: "<<Cr_CM(t)<<"\tCRU HPUE pred: "<<CPE_M(t)<<"\n";

    // fit to fishery CPUE
    Cnll -= dnorm( log( Catch(t)+Tiny ), log( Catch_hat(t)+Tiny ), Sigma_CPUE, true );
      //std::cout<<"CRU CPUE obs: "<<Catch(t)<<"\tCRM CPUE pred: "<<Catch_hat(t)<<"\n";
    
    // fit to acoustic estimates of Z
    Znll -= dnorm( log( ZEst(t) ), log( Acoust_hat(t) ), Sigma_Z, true );
      //std::cout<<"obs Z: "<<ZEst(t)<<"\test Z: "<<Acoust_hat(t)<<"\n";

    // prior on random effect around fishing mortality
    REnll -= dnorm( Zeta(t), Type(0.), proc_sig2, true);
  } // t
  
  for( int t=0; t<(nT-1); t++ ){
      // prior on random effect around recruitment
    REnll -= dnorm( Eps(t), Type(0.), proc_sig, true );
      //std::cout<<"t "<<t<<"\t"<<Eps(t)<<"\t"<<Zeta(t)<<"\n"<<REnll<<"\n";
  } // t
  
  nll = EFnll + CRnll + Cnll + Znll + REnll;
  nll += pen;

  // PRIORS
  Type meanM=log(0.25);
  Type sdM=0.25;
  Type meanR0=log(8.);
  Type sdR0=2.;
  Type meanRecK=2.;
  Type sdRecK=0.5;
  // prior on M
  nll -= dnorm( logM, meanM, sdM, true );
  // vague prior on R0
  nll -= dnorm( logR0, meanR0, sdR0, true );
  // informative prior on reck
  nll -= dnorm( logRecK, meanRecK, sdRecK, true );
  // vague priors on selectivity - uninformative, but keeps them away from large logistic tails
  nll -= dnorm( logitSel, Type(0.), Type(1.7), true).sum();
  nll -= dnorm( logitEFishSel, Type(0.), Type(1.7), true).sum();
    //std::cout<<"nll\t"<<EFnll<<"\t"<<CRnll<<"\t"<<Znll<<"\t"<<REnll<<"\t"<<pen<<"\t"<<nll<<"\n";
  
  // SPIT OUT VARIABLES
  ADREPORT(R0);                    // unfished equilibrium recruits
  ADREPORT(RecK);                  // recruitment compensation ratio
  ADREPORT(M);                     // natural mortality rate
  ADREPORT(q);                     // fishery catchability
  ADREPORT(aNest);                 // coefficient of impact of nest destruction
  ADREPORT(bNest);                 // coefficient of impact of male removal
  ADREPORT(pHarv);                // proportion harvested
  ADREPORT(Sel);                   // selectivity-at-age
  ADREPORT(EFishSel);              // elecrofishing selectivity-at-age
  ADREPORT(InitN);                 // initial abundance of unmarked fish

  REPORT(Marks);
  REPORT(Mark_hat);
  REPORT(Recs);
  REPORT(Rec_hat);
  REPORT(Cr_CU);
  REPORT(Cr_CM);
  REPORT(Catch);
  REPORT(Catch_hat);
  REPORT(ZEst);
  REPORT(Acoust_hat);
  REPORT(UtM);
  REPORT(UtF);
  REPORT(MtM);
  REPORT(MtF);
  REPORT(Nt);
  REPORT(EfishU);
  REPORT(MaleU);
  
  return(nll);
}
