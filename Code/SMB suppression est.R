GLOBALS_SECTION
  #include "stats.cxx"
  #include <fstream>

DATA_SECTION
  int seed; 														// random number seed
  !! ifstream ifs("seed.txt");
  !! ifs>>seed;
  init_int debug;                 									// flag indicating whether model should be debugged
  init_int sim;														// flag indicating whether data should be simulated (=1) or not (=0)
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
  init_matrix Nests(1,nT,AR,A);										// nests destroyed each year
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
  init_number logSela;
  init_number logSelb;
  init_number logDisMort;
  init_bounded_dev_vector proc_error(1,nT);
  init_vector logNinit(AR,A);
  init_number dend2;

  //VARIABLES
  4darray Surv(1,na,1,ny,1,ny,1,na);								// survival probability array
  4darray Pcap(1,na,1,ny,1,ny,1,na);								// capture probability array
  matrix Pnocap(1,ny,1,na);											// capture probability matrix for tags never recaptured
  vector temp(1,na);
  vector idist(1,na);												// asymptotic initial distribution
  number N0;														// initial abundance
  number res;														// residency - NOT USED
  number pp;														// capture probability
  number roe;														// overdispersion parameter
  number tempLHF;
  matrix Rels(1,na,1,ny);											// number of releases
  vector MR(1,na);													// mark rate by area
  sdreport_number M;
  matrix F(1,ny,1,na);												// sampling rate
  matrix N(1,ny,1,na);												// number of fish by area and year
  matrix NT(1,ny,1,na);												// number of tagged fish by area and year (after sampling)
  matrix NM(1,ny,1,na);												// number of tagged fish by area and year (before sampling)
  matrix catpred(1,ny,1,na);										// predicted recaptures
  matrix nocappred(1,ny,1,na);										// predicted new captures
  4darray cappred(1,na,1,ny,1,ny,1,na);								// predicted total captures
  matrix mov(1,na,1,na);											// gravity matrix
  //vector nodes(1,na*ny*2+na*na+1);
  vector nodes(1,na*ny+nF+nMP+1+na);								// something to do with reporting
  vector sim_pars(1,3+nF+nMP);										// simulated parameters; used to compare with estimated parameters
  vector est_pars(1,3+nF+nMP);										// estimated parameters
  vector par_err(1,3+nF+nMP);										// proportional error in parameters (est-act)/act
  objective_function_value ff;										// objective function
  number tiny;														// wee little number (prevents numbers from going to zero)
 
  LOCAL_CALCS
    if(debug==0){
       for(int i =1;i<=nF;i++){
         lnF(i)=log(0.01);
       }
       for(int i =1;i<=nMP;i++){
         movest(i)=0.;
       } 
       lnM=Mmu;
       lnN0=log(N0ini);
    }

PROCEDURE_SECTION

  ff.initialize();
  Pnocap.initialize();
  nodes.initialize();
  MR.initialize();
  Rels.initialize();
  catpred.initialize();
  double pi=3.141593;

//  if( sim==0 ){
//	for( int i=1; i<=ntag; i++ ){
//		sim_tagdat(i)=tagdat(i);
//	}
//  }
  if(debug==1) cout<<"Start init\n"<<endl;
  init();
  if(debug==1) cout<<"Start Dyn_Model\n"<<endl;
  Dyn_Model();
  if(debug==1) cout<<"Start Likelihoods\n"<<endl;
  Likelihoods();
  if(debug==1) cout<<"Model proceeding\n"<<endl;
  
//_________________________________________
FUNCTION init
  // initialize variables in the model
  N0=mfexp(lnN0);
  F=0.;
  N=0.;
  roe=mfexp(lnroe);
  //cout<<F<<endl;
  tiny=1E-20;
  M=mfexp(lnM);
  
  for(int i =1;i<=nF;i++){
    int y=Find(i,1);
    int a=Find(i,2);
    //if(debug==1)  cout<<"y "<<y<<"\ta "<<a<<endl;
    F(y,a)=mfexp(lnF(i)); 
  }
  
  if(fullmov==0){ 
	// Gravity model calc
    if(debug==1) cout<<"gravity\n"<<endl;
    //cout<<nMP<<endl;
    //cout<<movest<<endl;
    for(int a =1;a<=na;a++){
      mov(a,1)=0.;													// first row is set to zero
    }  

    for(int a=1;a<=na;a++){
      for(int a2=2;a2<=na;a2++){
        mov(a,a2)=movest(a2-1);     								// gravity weight for that receiving area
      }
      mov(a,a)+=mfexp(movest(nMP));  									 	// add on the viscosity parameter to staying in the same area
    }
    
  }else{ 
	// fully prescribed Markov matrix
    if(debug==1) cout<<"markov\n"<<endl;
    //cout<<nMP<<endl;
    //cout<<movest<<endl;
    for(int a =1;a<=na;a++){
      mov(a,1)=0.;													// first row is set to zero
    } 
    int i = 0;
    
    for(int a2=2;a2<=na;a2++){
      for(int a =1;a<=na;a++){
        i+=1;
        mov(a,a2)=movest(i);										// probability of moving to any cell is proportional to its estimated G
      }
    }
  }
    
  if(debug==1) cout<<"normalize\n"<<endl;
  mov=mfexp(mov);													// all values must be exponentiated
  for(int a = 1; a<=na;a++){
    mov(a)=mov(a)/sum(mov(a));										// make movement probabilities sum to 1
  }
  if(debug==1)	cout<<"mov\n"<<mov<<endl;
    
  // set the asymptotic population distribution by numerical approximation
  for(int a = 1; a<=na;a++){										// initialize initial distribution to uniform across areas
   idist(a)=1./na;
  }
  for(int i=1; i<=100; i++){										// numerically solve for initial distribution (it should take much fewer than 100 iterations)
    idist=idist*mov;
  }
  if(debug==1) cout<<"parms\nN0 "<<N0<<"\nF "<<F<<"\nM "<<M<<"\nmov "<<mov<<endl;
  if(debug==1) cout<<"idist\t"<<idist<<endl;
  

//_________________________________________
FUNCTION Dyn_Model
  // run the actual dynamics of the model

  // initialize numbers of the different combinations of fish
  N(1)=N0*idist;													// total number of fish
  NM(1)=N0*idist;													// number of unmarked fish availble to be tagged
  NT(1)=N0*idist;													// number of unmarked fish immediately after each years' tagging (i.e. includes newly marked fish)
  
  // redistribute fish through years and areas and remove fish due to mortality
  for(int y = 2; y<=ny; y++){
    N(y)=N(y-1)*mov*exp(-M); // 
    NM(y)=NT(y-1)*mov*exp(-M); 										// # before F for catch pred calculation
    NT(y)=elem_prod(NT(y-1),exp(-M-F(y-1)))*mov; 					// now accounting for F (N is a neg binomial model)
  }
  
  for(int a = 1; a<=na; a++){										// initialize probability of capture and survival
    for(int y = 1; y<=ny; y++){
      //Pnocap(y,a)=1.;
      for(int y2 = 1; y2<=ny; y2++){
        for(int a2 = 1; a2<=na; a2++){
          Pcap(a,y,y2,a2)=0.;
          Surv(a,y,y2,a2)=0.;
        }
      }
    }
  }
 
  for(int a = 1; a<=na;a++){										// probability of recapture and survival immediately after capture
    for(int y =1; y<=ny;y++){
      Surv(a,y,y,a)=1.; //mfexp(-M/2);
      Pcap(a,y,y,a)=0.; //mfexp(-(M+F(y,a))/2.)*(1.-mfexp(-F(y,a)/2.));		// Brett added the half fishing mortality option to be consistent with the paper - could be changed, but inconsistent with simulation code
    } 
  }  
  
  for(int a = 1; a<=na; a++){										// probability of surviving and being in state a
    for(int y = 1; y<=(ny-1); y++){
      for(int y2 = (y+1); y2<=ny; y2++){
        //cout<<y2<<endl;
        temp=elem_prod(Surv(a)(y)(y2-1),exp(-F(y2-1)-M));
        Surv(a)(y)(y2)=temp*mov;
        Pcap(a)(y)(y2)=elem_prod(Surv(a)(y)(y2),(1-exp(-F(y2))));
        //cout<<Surv(a)(y)(y2)<<endl;
        //cout<<Pcap(a)(y)(y2)<<endl;
        //cout<<"-----"<<endl;
      }
    }
  }
  //cout<<"F\n"<<F<<endl;
   
//_________________________________________
FUNCTION Likelihoods
  // Calculate likelihoods of data given the model and prior probabilies
  //cout<<Pnocap<<endl;
  //cout<<"------"<<endl;
  
  if(sim==1){
  	ntag = sim_ntag;
  }
  
  if(debug==1) cout<<"theta[U] "<<ff<<endl;
  for(int a = 1; a<=na; a++){										// predicted catch of fish never recaptured (theta[U}])
    for(int y = 1; y<=(ny-1); y++){
      Pnocap(y,a)=1-sum(Pcap(a)(y));
   
      //if(cattot(y,a)>0){
      nocappred(y,a)=cattot(y,a)*Pnocap(y,a); 
      if(debug==1)	cout<<y<<" "<<a<<" nocap: "<<nocappred(y,a)<<" "<<nocap(y,a)<<" "<<cattot(y,a)<<" "<<Pnocap(y,a)<<endl;      
      switch(LHF_U){
          	  		    
        case 1:  // SSQ
          ff+=pow(nocap(y,a)-nocappred(y,a),2);    
          break;
          	  		   
        case 2:  // Log-normal
          ff+=0.5*log(2.*pi)+log(nocap(y,a)*CV+tiny)+(0.5/(CV*CV))*log((nocap(y,a)/nocappred(y,a))+tiny)*log((nocap(y,a)/nocappred(y,a))+tiny);
          break;
          	     
        case 3:  // Multinomial / binomial
          ff-=nocap(y,a)*log(Pnocap(y,a));
		  //cout<<nocap(y,a)<<" "<<log(Pnocap(y,a))<<endl;
          break;
            
        case 4:  // Poisson
          //ff-=nocap(y,a)*log(nocappred(y,a))-nocappred(y,a);
          tempLHF=(pow(nocappred(y,a),nocap(y,a))*mfexp(-nocappred(y,a)))/exp(factln(nocap(y,a)));
          ff-=log(tempLHF);
          break;
              
        case 5:  // Negative binomial w O.D.
      	  tempLHF=(exp(gammln(roe+nocap(y,a)))/(exp(gammln(roe))*exp(factln(nocap(y,a))))) * pow(roe/(roe+nocappred(y,a)),roe) * pow(nocappred(y,a)/(nocappred(y,a)+roe),nocap(y,a)); // negbin
          ff-=log(tempLHF); 
      	  break;
          
      }
      //}
      //cout<<ff<<" nocappred"<<endl;    
      Rels(a,y)+=nocap(y,a);
      
    }
  }

  if(debug==1) cout<<"theta[C] "<<ff<<endl;
  for(int i = 1; i<=nF;i++){
    
    int y = Find(i,1);
    int a = Find(i,2);
        
    pp=(1-exp(-F(y,a)));
    catpred(y,a)=NM(y,a)*pp;										// predicted capture of previously untagged fish (theta[C])
    
    if(debug==1)	cout<<y<<" "<<a<<" Ct:"<<NM(y,a)<<" "<<catpred(y,a)<<" "<<cat(y,a)<<endl;
    switch(LHF_C){
    	  		    
      case 1:  // SSQ
    	ff+=pow(cat(y,a)-catpred(y,a),2);    
     	break;
    	  		   
      case 2:  // Log-normal
    	ff+=0.5*log(2.*pi)+log(cat(y,a)*CV+tiny)+(0.5/(CV*CV))*log((cat(y,a)/catpred(y,a))+tiny)*log((cat(y,a)/catpred(y,a))+tiny);
		//cout<<cat(y,a)<<" "<<catpred(y,a)<<endl;
    	break;
    	     
      case 3:  // Binomial
    	ff-=cat(y,a)*log(pp)+(NM(y,a)-cat(y,a))*log(1-pp);  
        break;
      
      case 4:  // Poisson
       // ff-=cat(y,a)*log(catpred(y,a)+tiny)-catpred(y,a);  
        tempLHF=(pow(catpred(y,a),cat(y,a))*mfexp(-catpred(y,a)))/exp(factln(cat(y,a)));
        ff-=log(tempLHF);
        break;
        
      case 5:  // Negative binomial w O.D.
		tempLHF=(exp(gammln(roe+cat(y,a)))/(exp(gammln(roe))*exp(factln(cat(y,a))))) * pow(roe/(roe+catpred(y,a)),roe) * pow(catpred(y,a)/(catpred(y,a)+roe),cat(y,a)); // negbin
        ff-=log(tempLHF+tiny); 
	break;
    
    }
    
    //cout<<ff<<" "<<pp<<" "<<catpred(y,a)<<" "<<NM(y,a)<<" catpred"<<endl;
    
  }
  
  if(debug==1) cout<<"theta[R] "<<ff<<endl;
  for(int i = 1; i<=ntag; i++){										// predicted recapture of marked fish (theta[R])
    
    int a=tagdat(i,1);
    int y=tagdat(i,2);
    int y2=tagdat(i,3);
    int a2=tagdat(i,4);
	//cout<<i<<" "<<tagdat(i)<<" "<<Pcap(a,y,y2,a2)<<Surv(a,y,y2,a2)<<endl;//cappred<<endl;
    
    cappred(a,y,y2,a2)=Pcap(a,y,y2,a2)*cattot(y,a);
    if(debug==1)	cout<<a<<" "<<y<<" "<<y2<<" "<<a2<<" tagdat:"<<cappred(a,y,y2,a2)<<" "<<tagdat(i,5)<<endl;
    switch(LHF_R){
              	  		    
      case 1:  // SSQ
        ff+=pow(tagdat(i,5)- cappred(a,y,y2,a2),2);    
        break;
              	  		   
      case 2:  // Log-normal
        ff+=0.5*log(2.*pi)+log(tagdat(i,5)*CV+tiny)+(0.5/(CV*CV))*log((tagdat(i,5)/cappred(a,y,y2,a2))+tiny)*log((tagdat(i,5)/cappred(a,y,y2,a2))+tiny);
        break;
              	     
      case 3:  // Multinomial / binomial
        ff-=tagdat(i,5)*log(Pcap(a,y,y2,a2));
		//cout<<a<<" "<<y<<" "<<y2<<" "<<a2<<" "<<tagdat(i,5)<<" "<<log(Pcap(a,y,y2,a2))<<endl;
        break;
                
      case 4:  // Poisson
        //ff-=tagdat2(i,5)*log(cappred(a,y,y2,a2))-cappred(a,y,y2,a2);  
        tempLHF=(pow(cappred(a,y,y2,a2),tagdat(i,5))*mfexp(-cappred(a,y,y2,a2)))/exp(factln(tagdat(i,5)));
		ff-=log(tempLHF);
      
        break;
                  
      case 5:  // Negative binomial w O.D.
        tempLHF=(exp(gammln(roe+tagdat(i,5)))/(exp(gammln(roe))*exp(factln(tagdat(i,5))))) * pow(roe/(roe+cappred(a,y,y2,a2)),roe) * pow(cappred(a,y,y2,a2)/(cappred(a,y,y2,a2)+roe),tagdat(i,5)); // negbin
        ff-=log(tempLHF); 
        break;
              
    }
    
   // cout<<ff<<" cappred"<<endl;  
    
    Rels(a,y)+=tagdat(i,5);
    
  }
  //cout<<Rels<<endl;
  
  if(debug==1) cout<<"theta[Y] "<<ff<<endl;
  for(int a=1;a<=na;a++){											// telemetry data (mmovement matrix)
    for(int a2=1;a2<=na;a2++){
      ff-=tel(a,a2)*log(mov(a,a2));
      //cout<<tel(a,a2)<<" "<<log(mov(a,a2)+tiny)<<endl;
    }
  }
  //exit(1);
  
  ff+=dlnorm(M,Mmu,Msd);											// prior on natural mortality
  ff+=dlnorm(roe,roemu,roesd);										// prior on overdispersion parameter (if used)
  
  if(debug==1){	
	cout<<"ff "<<ff<<endl; 
	exit(1);
  }

  MR=0;
  for(int a = 1; a<=na; a++){
    for(int y = 1; y<=(ny-1); y++){
      for(int a2 = 1; a2<=na; a2++){
        MR(a2)+=Rels(a,y)*Surv(a)(y)(ny-1)(a2)/N(ny-1)(a);
      }
    }
  }
   
  //cout<<N<<endl;
  //cout<<cat<<endl;
  for(int a=1;a<=na;a++){
    for(int y=1;y<=ny;y++){
      nodes((a-1)*ny+y)=N(y,a);
    }
  }
  for(int i=1;i<=nF;i++){
    nodes(na*ny+i)=mfexp(lnF(i));
  }    
  for(int i=1;i<=nMP;i++){
    nodes(na*ny+nF+i)=mfexp(movest(i));
  }    
  nodes(na*ny+nF+nMP+1)=M;
  for(int a=1;a<=na;a++){
    nodes(na*ny+nF+nMP+1+a)=MR(a);
  }  
  
  if(mceval_phase()){
    ofstream nodesout("nodes.cha", ios::app);
    nodesout<<nodes<<"\n";
  }
 
   
TOP_OF_MAIN_SECTION

  arrmblsize =25000000;            
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(22000);  
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(32000000);  
  gradient_structure::set_CMPDIF_BUFFER_SIZE(560000000); 
  gradient_structure::set_MAX_NVAR_OFFSET(1800);


REPORT_SECTION
  //report<<munofish<<endl;
  //report<<munofish<<endl;
  report <<"ny, number of years"<<endl;
  report <<ny<<endl;
  report <<"na, number of areas"<<endl;
  report <<na<<endl;
  report <<"nF, number of estimated Fs"<<endl;
  report <<nF<<endl;
  report <<"nMP, number of estimated movement params"<<endl;
  report <<nMP<<endl;
  report <<"M, natural mortality rate"<<endl;
  report <<M<<endl;
  report <<"roe, overdispersion parameter"<<endl;
  report <<roe<<endl;
  report <<"N, number of individuals"<<endl;
  report <<N<<endl;
  report <<"F, fishing mortality rate"<<endl;
  report <<F<<endl;
  report <<"mov, movement probability matrix"<<endl;
  report <<mov<<endl;
  report <<"idist, equilibrium spatial distribution"<<endl;
  report <<idist<<endl;
  report <<"Find, estimated fishing mortality rates"<<endl;
  report <<Find<<endl;
  report <<"Current mark rate by area"<<endl;
  report <<MR<<endl;
  report <<"observed catches"<<endl;
  report<<cat<<endl;
  report<<"predicted catches"<<endl;
  report<<catpred<<endl;
  report<<"OBJ"<<endl;
//  if(LHF==5){
//    report<<ff<<endl;
//  }else{
    report<<ff-dlnorm(roe,roemu,roesd)<<endl; // without roe overdispersion parameter
//  } 
  report<<"npar"<<endl;
//  if(LHF==5){
//    report<<nF+nMP+3<<endl; //with roe Overdispersion parameter
//  }else{
    report<<nF+nMP+2<<endl;
//  } 
  report <<"data check"<<endl;
  report <<datacheck<<endl;    

  ofstream ofsm2("est pars1.txt",ios::app);
  ofsm2<<objective_function_value::pobjfun->gmax<<"\t"<<seed<<"\t"<<ff<<"\t"<<exp(lnN0)<<"\t"<<exp(lnroe)<<"\t"<<exp(lnM)<<"\t"<<exp(lnF)<<"\t"<<movest<<endl;
