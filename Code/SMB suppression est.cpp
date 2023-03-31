#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in the data
  DATA_INTEGER(AR);                    // age at recruitment
  DATA_INTEGER(A);                     // maximum age (a plus-group)
  DATA_INTEGER(nT);                    // number of years
  DATA_SCALAR(Linf);                   // von Bertalanffy asymptotic length
  DATA_SCALAR(K);                      // von Bertalanffy metabolic growth parameter
  DATA_SCALAR(CV);                     // cv in length-at-age
  DATA_SCALAR(LWa);                    // scalar of length-weight relationship
  DATA_SCALAR(LWb);                    // exponent of length-weight relationship
  DATA_SCALAR(EggMass);                // eggs per gram of fish
  DATA_SCALAR(Mat50);                  // age at 50% maturity
  DATA_SCALAR(MatSig);                 // logistic slope in age-at-maturity
  DATA_SCALAR(DisMort);                // discard mortality (CONSIDER ESTIMATING LATER)
  DATA_SCALAR(pHarv);                  // proportion harvested in status-quo management
  DATA_SCALAR(pHarv2);                 // proportion harvested if regulations clarified
  DATA_VECTOR(Age);                    // ages included in the model
  DATA_ARRAY(FMarks);                  // marked females at start of each year
  DATA_ARRAY(MMarks);                  // marked males at start of each year
  DATA_ARRAY(CUtF);                    // actual fishery harvest of unmarked females (not needed?)
  DATA_ARRAY(CUtM);                    // actual fishery harvest of unmarked males (not needed?)
  DATA_ARRAY(CMtF);                    // actual fishery harvest of marked females (not needed?)
  DATA_ARRAY(CMtM);                    // actual fishery harvest of marked males (not needed?)
  DATA_ARRAY(CUMal);                   // unmarked males removed off nests
  DATA_ARRAY(CMMal);                   // marked males removed off nests
  DATA_ARRAY(Cr_CUF);                  // CPUE of unmarked females in creel
  DATA_ARRAY(Cr_CUM);                  // CPUE of unmarked males in creel
  DATA_ARRAY(Cr_CMF);                  // CPUE of marked females in creel
  DATA_ARRAY(Cr_CMM);                  // CPUE of marked males in creel
  
  // declare parameters to be estimated 
  PARAMETER(logR0);                    // log of unfished equilibrium recruits
  PARAMETER(logRecK);                  // log of recruitment compensation ratio
  PARAMETER(logq);                     // log of fishery catchability
  PARAMETER_VECTOR(logitSel);          // selectivity-at-age on logit-scale
  PARAMETER_VECTOR(logInitN);          // log of initial abundance
  PARAMETER_VECTOR(Eps);               // process error in annual recruitment (random variable)
  
  // declare variables to be calculated
  Type R0=exp(logR0);
  Type RecK=exp(logRecK);
  Type q=exp(logq);
  vector<Type> Sel(nA);
  Type InitN
  Type phi0;                // spawners (eggs) per recruit
  Type RecA;                // Beverton-Holt A parameter
  Type RecB;                // Beverton-Holt B parameter
  vector<Type> Lt(nA);      // length-at-age
  vector<Type> Wt(nA);      // weight-at-age
  vector<Type> Mat(nA);     // maturity-at-age
  vector<Type> Fec(nA);     // fecundity-at-age
  vector<Type> lx(nA);      // survivorship-to-age
  
  // objective function
  Type nll = 0.0;
  
  // INITIALIZE OJIVES AND CALCULATE RECRUITMENT PARAMETERS
  Lt = Linf * ( Type(1.0) - exp( -K*Age ));
  Wt = LWa * Lt ^ LWb;
  Mat = Type(1.0)/(Type(1.0) - exp( -( Lt - Mat50 )/MatSig ));
  Fec = Wt * EggMass;
  lx = exp( -M*(Age-AR) );
  lx(nA-AR) /= Type(1.0)-exp(-M);
  Sel=Type(1.0)/(Type(1.0)-exp(-logitSel));
  
  phi0 = sum( Mat * Fec * lx ) * pFem;
  RecA = RecK / phi0;
  RecB = ( RecK-Type(1.0) )/( R0 * phi0 );
}