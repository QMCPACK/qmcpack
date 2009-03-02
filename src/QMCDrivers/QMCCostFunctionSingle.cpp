//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus {

  QMCCostFunctionSingle::QMCCostFunctionSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCCostFunctionBase(w,psi,h)
  {
  	CSWeight=(1.0); 
  	cout<<" Using QMCCostFunctionSingle::QMCCostFunctionSingle"<<endl;
  }

  /** Clean up the vector */
  QMCCostFunctionSingle::~QMCCostFunctionSingle() {
  }

  /**  Perform the correlated sampling algorthim.
   */
  QMCCostFunctionSingle::Return_t QMCCostFunctionSingle::correlatedSampling() {

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    //Return_t eloc_new;
    //int totalElements=W.getTotalNum()*OHMMS_DIM;

    Return_t wgt_tot=0.0;
    //bool usingWeight=UseWeight;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    int iw=0;
    for(; it!=it_end;++it,++iw)
    {
      Walker_t& thisWalker(**it);
      Return_t* restrict saved = Records[iw];
      W.R=(*it)->R;
      W.update();
      Return_t logpsi=Psi.evaluateDeltaLog(W);
      W.G += *dLogPsi[iw];
      W.L += *d2LogPsi[iw];
      
			RealType KEtemp = H_KE.evaluate(W); 
			Return_t eloc_new = KEtemp + saved[ENERGY_FIXED];
          Return_t weight;
          if (samplePsi2) weight = std::min( std::exp(2.0*(logpsi-saved[LOGPSI_FREE])),MaxWeight) ;
          else weight = std::min( std::exp( logpsi-saved[LOGPSI_FREE] ),MaxWeight) ;
          if (KEtemp<MinKE) weight=0.0;

      //Return_t eloc_new=+saved[ENERGY_FIXED];
      //Return_t weight = UseWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
      
			vector<Return_t>* Dsaved= &(TempDerivRecords[iw]) ;
			vector<Return_t>* HDsaved= &(TempHDerivRecords[iw]) ;
			Psi.evaluateDerivatives(W,KEtemp,OptVariables,*Dsaved,*HDsaved);
// 			for(int l=0; l<NumOptimizables; l++) (*HDsaved)[l] += saved[ENERGY_FIXED] * (*Dsaved)[l];


      saved[ENERGY_NEW]=eloc_new;
      saved[REWEIGHT]=weight;
      wgt_tot+=weight;
    }

    //collect the total weight for normalization and apply maximum weight
    myComm->allreduce(wgt_tot);

    for(int i=0; i<SumValue.size(); i++) SumValue[i]=0.0;
    CSWeight=wgt_tot=1.0/wgt_tot;

    int nw=W.getActiveWalkers();
    for(iw=0; iw<nw;iw++) {
      Return_t* restrict saved = Records[iw];
				Return_t weight=saved[REWEIGHT]*wgt_tot;
				Return_t eloc_new=saved[ENERGY_NEW];
				Return_t delE=std::pow(abs(eloc_new-EtargetEff),PowerE);
				SumValue[SUM_E_BARE] += eloc_new;
				SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
				SumValue[SUM_ABSE_BARE] += delE;
				SumValue[SUM_E_WGT] += eloc_new*saved[REWEIGHT];
				SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*saved[REWEIGHT];
				SumValue[SUM_ABSE_WGT] += delE*saved[REWEIGHT];
				SumValue[SUM_WGT] += saved[REWEIGHT];
				SumValue[SUM_WGTSQ] += saved[REWEIGHT]*saved[REWEIGHT];
        SumValue[SUM_WGTSQ] += weight*weight;
    }

    //collect everything
    myComm->allreduce(SumValue);

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

  void 
  QMCCostFunctionSingle::getConfigurations(const string& aroot) {
    app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << endl;
    app_log() << "    number of walkers before load " << W.getActiveWalkers() << endl;
    Timer t1;
    W.loadEnsemble();
    app_log() << "    Loading time = " << t1.elapsed() << endl;
    app_log() << "    number of walkers after load " << W.getActiveWalkers() << endl;
    //if(aroot.size() && aroot != "invalid") {
    //  app_log() << "  Reading configurations from the previous qmc block" << endl;
    //  HDFWalkerInputCollect wReader(aroot);
    //  wReader.putSingle(W);
    //}
    //
    if(dLogPsi.size() != W.getActiveWalkers())
    {
      delete_iter(dLogPsi.begin(),dLogPsi.end());
      delete_iter(d2LogPsi.begin(),d2LogPsi.end());
      int nptcl=W.getTotalNum();
      int nwtot=W.getActiveWalkers();
      dLogPsi.resize(nwtot,0);
      d2LogPsi.resize(nwtot,0);
      for(int i=0; i<nwtot; ++i) dLogPsi[i]=new ParticleGradient_t(nptcl);
      for(int i=0; i<nwtot; ++i) d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
    }

#if defined(QMCCOSTFUNCTION_DEBUG)
    if(debug_stream) delete debug_stream;
    char fname[16];
    sprintf(fname,"optdebug.p%d",OHMMS::Controller->mycontext());
    debug_stream = new ofstream(fname);
    debug_stream->setf(ios::scientific, ios::floatfield);
    debug_stream->precision(8);

    *debug_stream << "Initial : " << endl;
    for(int i=0; i<OptParams.size(); i++) 
      *debug_stream << " " << IDtag[i] << "=" << OptParams[i];
    *debug_stream << endl;
#endif
  }

 // /** evaluate everything before optimization */
 // void 
 // QMCCostFunctionSingle::getConfigurations(vector<string>& ConfigFile, 
 //     int partid, int nparts) {
 //   if(ConfigFile.size()) {

 //     app_log() << "  Reading configurations from mcwalkerset " << endl;

 //     W.destroyWalkers(W.begin(),W.end());
 //     for(int i=0; i<ConfigFile.size(); i++) {
 //       //JNKIM: needs to change to HDFWalkerInputCollect
 //       //HDFWalkerInput0 wReader(ConfigFile[i],partid,nparts);
 //       HDFWalkerInputCollect wReader(ConfigFile[i]);
 //       wReader.putSingle(W);
 //       //wReader.put(W,-1);
 //     }

 //     //remove input files
 //     ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
 //   }
 // }

  /** evaluate everything before optimization */
  void 
  QMCCostFunctionSingle::checkConfigurations() {

    //dG.resize(W.getTotalNum());
    //dL.resize(W.getTotalNum());
    int numLocWalkers=W.getActiveWalkers();
    Records.resize(numLocWalkers,6);
    
		TempHDerivRecords.resize(numLocWalkers,vector<Return_t>(NumOptimizables,0));
		TempDerivRecords.resize(numLocWalkers,vector<Return_t>(NumOptimizables,0));

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    //int totalElements=W.getTotalNum()*OHMMS_DIM;
    Etarget=0.0;
    Return_t e2sum=0.0;
    for(int iw=0; it!=it_end; ++it,++iw)
    {
      Walker_t& thisWalker(**it);
      W.R=thisWalker.R;
      W.update();
      Return_t*  saved=Records[iw];
      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iw], *d2LogPsi[iw]);
      Return_t e=H.evaluate(W);
      e2sum += e*e;
      Etarget += saved[ENERGY_TOT] = e;
      saved[ENERGY_FIXED] = H.getLocalPotential();
      
      vector<Return_t>* Dsaved=  &(TempDerivRecords[iw]) ;
      vector<Return_t>* HDsaved=  &(TempHDerivRecords[iw]) ;
      Psi.evaluateDerivatives(W,e-saved[ENERGY_FIXED],OptVariables,*Dsaved,*HDsaved);
    }

    //Need to sum over the processors
    vector<Return_t> etemp(3);
    etemp[0]=Etarget;
    etemp[1]=static_cast<Return_t>(numLocWalkers);
    etemp[2]=e2sum;

    myComm->allreduce(etemp);
    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
    NumSamples = static_cast<int>(etemp[1]);

    app_log() << "  VMC Eavg = " << Etarget << endl;
    app_log() << "  VMC Evar = " << etemp[2]/etemp[1]-Etarget*Etarget << endl;
    app_log() << "  Total weights = " << etemp[1] << endl;

    setTargetEnergy(Etarget);

    ReportCounter=0;
  }

  void QMCCostFunctionSingle::resetPsi()
  {
    if(OptVariables.size() < OptVariablesForPsi.size())
    {
      for(int i=0; i<equalVarMap.size(); ++i)
        OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
    }
    else
      for(int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];
      
      Psi.resetParameters(OptVariablesForPsi);
  }
  
  
    void QMCCostFunctionSingle::GradCost(vector<Return_t>& PGradient, vector<Return_t> PM, Return_t FiniteDiff)
  {
  	if (FiniteDiff>0)
  	{
  		
			Return_t dh=1.0/(2.0*FiniteDiff);
			for(int i=0; i<NumOptimizables ; i++) {
				for(int j=0; j<NumOptimizables; j++) OptVariables[j]=PM[j];
				OptVariables[i] = PM[i]+ FiniteDiff;
				Return_t CostPlus = this->Cost(); 
				OptVariables[i] = PM[i]- FiniteDiff;
				Return_t CostMinus = this->Cost(); 
				PGradient[i]= (CostPlus-CostMinus)*dh;
			}
    }
    else
    {
    	for(int j=0; j<NumOptimizables; j++) OptVariables[j]=PM[j];
    	 resetPsi();
 
			 //evaluate new local energies and derivatives
			 NumWalkersEff=correlatedSampling( );
			  //Estimators::accumulate has been called by correlatedSampling


			 curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
			 Return_t curAvg2_w = curAvg_w*curAvg_w;
				
        vector<Return_t> EDtotals(NumOptimizables,0.0);
        vector<Return_t> E2Dtotals(NumOptimizables,0.0);
        vector<Return_t> EDtotals_w(NumOptimizables,0.0);
        vector<Return_t> E2Dtotals_w(NumOptimizables,0.0);
        vector<Return_t> URV(NumOptimizables,0.0);
        vector<Return_t> HD_avg(NumOptimizables,0.0);
        Return_t e2(0.0);


 
				int nw=W.getActiveWalkers();
				for(int iw=0; iw<nw;iw++) 
				{
					Return_t*  saved=Records[iw];
					Return_t weight=saved[REWEIGHT]/SumValue[SUM_WGT];
					Return_t eloc_new=saved[ENERGY_NEW];
     e2 += eloc_new*eloc_new*weight;
					Return_t delta_eloc_new=eloc_new-curAvg_w;
					vector<Return_t>  Dsaved= (TempDerivRecords[iw]);
					vector<Return_t>  HDsaved= (TempHDerivRecords[iw]);
					for(int pm=0; pm<NumOptimizables;pm++) 
					{
                    HD_avg[pm]+= HDsaved[pm];
                    Return_t val;
                    if (samplePsi2) val = (HDsaved[pm] + 2.0*(eloc_new - curAvg_w)*Dsaved[pm]);
                    else val = (HDsaved[pm] + (eloc_new - curAvg_w)*Dsaved[pm]);
                    EDtotals[pm] +=  val;
                    EDtotals_w[pm] += weight*val;
                    
					}
				} 
        myComm->allreduce(EDtotals);
        myComm->allreduce(EDtotals_w);
        myComm->allreduce(HD_avg);
 
        Return_t wgtinv=1.0/static_cast<Return_t>(NumSamples);
				for(int iw=0; iw<nw;iw++) 
				{
					Return_t*  saved=Records[iw];
					Return_t weight=saved[REWEIGHT]*CSWeight;
					Return_t eloc_new=saved[ENERGY_NEW];
					Return_t delta_eloc_new=eloc_new-curAvg_w;
					vector<Return_t> Dsaved= TempDerivRecords[iw];
					vector<Return_t> HDsaved= TempHDerivRecords[iw];
					for(int pm=0; pm<NumOptimizables;pm++) 
					{
                    URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
                    Return_t val;
                    if (samplePsi2) val=  2.0*(eloc_new*eloc_new - e2)*Dsaved[pm]
                                    + 2.0*eloc_new*(HDsaved[pm])
                                    - 2.0*curAvg_w*EDtotals_w[pm];
                    else val = (eloc_new*eloc_new - e2)*Dsaved[pm]
                                    + 2.0*eloc_new*(HDsaved[pm])
                                    - 2.0*curAvg_w*EDtotals_w[pm];
                    E2Dtotals[pm] += val ;
                    E2Dtotals_w[pm] += val * weight;
					}
				}

        myComm->allreduce(E2Dtotals);
        myComm->allreduce(E2Dtotals_w);
        myComm->allreduce(URV);
			
        for (int pm=0; pm<NumOptimizables;pm++)  URV[pm] *=wgtinv;
        for (int j=0; j<NumOptimizables; j++) PGradient[j] = (w_var*E2Dtotals_w[j] + w_en*EDtotals_w[j] + w_w*URV[j]) ;


			 IsValid=true;

			 if (NumWalkersEff < MinNumWalkers*NumSamples)
    {
                  ERRORMSG("CostFunction-> Number of Effective Walkers is too small " << NumWalkersEff<< "Minimum required"<<MinNumWalkers*NumSamples)
// 				 ERRORMSG("Going to stop now.")
				 IsValid=false;
			 }
    }
  }

  ///** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
  // *@return true always
  // *
  // * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
  // * set \f$ {{\bf \alpha_{i+1}}}\f$  
  // */
  //bool
  //QMCCostFunctionSingle::resetWaveFunctions() {

  //  resetPsi();
  //  return true;
  //}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
