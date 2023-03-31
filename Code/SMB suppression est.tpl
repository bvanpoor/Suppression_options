<<<<<<< Local Changes
GLOBALS_SECTION
  #include "stats.cxx"
  #include <fstream>

DATA_SECTION
  int seed; 														// random number seed
  !! ifstream ifs("seed.txt");
  !! ifs>>seed;
  init_int debug;                 									// flag indicating whether model should be debugged
  init_int nT;														// number of years of data
  init_int AR;														// age of recruitment
  init_int A;														// plus age group
  init_int nA;														// number of age classes
  init_number Linf;    												// von Bertalanffy asymptotic length
  init_number K;													// von Bertalanffy metabolic growth parameter
  init_number CV;													// coefficient of variation in length at age
  init_number Mat50;												// length at 50% maturity
  init_number MatSig;												// logistic slope of length at maturity
  init_number EggMass;												// eggs per gram of fish
  init_number pFem;													// initial sex ratio
  init_number pHarv;												// initial proportion voluntarily harvested
  init_number pHarv2;												// proportion voluntarily harvested if regs tightened
  init_matrix FMarks(1,nT,AR,A);									// marked females each year
  init_matrix MMarks(1,nT,AR,A);									// marked males each year
  init_matrix CUtF(1,nT,AR,A);										// harvest of unmarked females in the fishery
  init_matrix CUtM(1,nT,AR,A);										// harvest of unmarked males in the fishery
  init_matrix CMtF(1,nT,AR,A);										// harvest of marked females in the fishery
  init_matrix CMtM(1,nT,AR,A);										// harvest of marked males in the fishery
  init_matrix CUMal(1,nT,AR,A);										// removal of unmarked males from nests
  init_matrix CMMal(1,nT,AR,A);										// removal of marked males from nests
  init_matrix Nest(1,nT,AR,A);										// nests destroyed each year
  init_matrix Cr_CUF(1,nT,AR,A);									// CPUE of unmarked females in creel
  init_matrix Cr_CUM(1,nT,AR,A);									// CPUE of unmarked males in creel
  init_matrix Cr_CMF(1,nT,AR,A);									// CPUE of marked females in creel
  init_matrix Cr_CMM(1,nT,AR,A);									// CPUE of marked males in creel
  init_number datacheck;
  !! if(datacheck!=999)	{
	  cout<<"data input error! datacheck = "<<datacheck<<endl;	
	  exit(1);
  }
     
PARAMETER_SECTION
  init_number logR0;
  init_number logReck;
  init_number logM;
  init_number logq;
  init_number logDisMort;
  init_vector logitSel(AR,A);
  init_bounded_dev_vector proc_error(1,nT);
  init_vector logNinit(AR,A);
  init_number dend2;

  //VARIABLES
  number R0;
  number Reck;
  number M;
  number q;
  number Sela;
  number Selb;
  number DisMort;
  number RecA;
  number RecB;
  vector Lt(AR,A);
  vector Wt(AR,A);
  vector Fec(AR,A);
  vector Mat(AR,A);
  vector Sel(AR,A);
  vector lx(AR,A);
  
//  4darray Pcap(1,na,1,ny,1,ny,1,na);								// capture probability array
//  matrix Pnocap(1,ny,1,na);											// capture probability matrix for tags never recaptured
//  vector temp(1,na);
//  number N0;														// initial abundance
//  sdreport_number M;
  objective_function_value ff;										// objective function
  number tiny;														// wee little number (prevents numbers from going to zero)
 
  LOCAL_CALCS
//    if(debug==0){
//    }

PROCEDURE_SECTION

  ff.initialize();
//  Pnocap.initialize();
//  nodes.initialize();
//  MR.initialize();
//  Rels.initialize();
//  catpred.initialize();
//  double pi=3.141593;

  if(debug==1) cout<<"Start init\n"<<endl;
//  init();
//  if(debug==1) cout<<"Start Dyn_Model\n"<<endl;
//  Dyn_Model();
//  if(debug==1) cout<<"Start Likelihoods\n"<<endl;
//  Likelihoods();
//  if(debug==1) cout<<"Model proceeding\n"<<endl;
  
//_________________________________________
FUNCTION init
  dvar_variable PhiE;
  R0 = exp( logR0 );
  RecK = exp( logRecK );
  
  for( int i=AR; i<=A; i++){
    Lt(i) = Linf * ( 1 - exp( -K * i ));
    Wt(i) = 1e-5 * Lt(i) ^ 3;
    Fec(i) = Wt(i) * EggMass;
    Mat(i) = 1./( 1. + exp( -( Lt(i)-Mat50)/MatSig));
    Sel(i) = 1./(1. + exp( logitSel(i) ));
	lx(i) = exp( -M * (i-AR));
  }
  lx(A) = lx(A)/(1-exp(-M));
  
  PhiE = sum( Fec * Mat * lx );
  RecA = log(RecK)/PhiE;
  RecB = (RecK-1)/(R0*PhiE);
   
TOP_OF_MAIN_SECTION

  arrmblsize =25000000;            
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(22000);  
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(32000000);  
  gradient_structure::set_CMPDIF_BUFFER_SIZE(560000000); 
  gradient_structure::set_MAX_NVAR_OFFSET(1800);


REPORT_SECTION
  //report<<munofish<<endl;
  //report<<munofish<<endl;

//  ofstream ofsm2("est pars1.txt",ios::app);
//  ofsm2<<objective_function_value::pobjfun->gmax<<"\t"<<seed<<"\t"<<ff<<"\t"<<exp(lnN0)<<"\t"<<exp(lnroe)<<"\t"<<exp(lnM)<<"\t"<<exp(lnF)<<"\t"<<movest<<endl;
=======
# simulated data for smallmouth bass suppression
written on2023-02-19 11:11:49
# debug
0
# nT
27
# AR
1
# A
5
# nA
5
# Linf
50
# K
0.35
# CV
0.1
# Mat50
2.5
# MatSig
0.5
# EggMass
10000
# pFem
0.5>>>>>>> External Changes
