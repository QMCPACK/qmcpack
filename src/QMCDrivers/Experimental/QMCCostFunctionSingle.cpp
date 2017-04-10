//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
//#define QMCCOSTFUNCTION_DEBUG

#include <Platforms/sysutil.h>
#include "Utilities/PooledData.h"


namespace qmcplusplus
{

QMCCostFunctionSingle::QMCCostFunctionSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCCostFunctionBase(w,psi,h)
{
  CSWeight=(1.0);
  APP_ABORT("QMCCostFunctionSingle is deprecated");
  app_log()<<" Using QMCCostFunctionSingle::QMCCostFunctionSingle"<< std::endl;
}

/** Clean up the vector */
QMCCostFunctionSingle::~QMCCostFunctionSingle()
{
}

/**  Perform the correlated sampling algorthim.
 */
QMCCostFunctionSingle::Return_t QMCCostFunctionSingle::correlatedSampling(bool needGrad)
{
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  //Return_t eloc_new;
  //int totalElements=W.getTotalNum()*OHMMS_DIM;
  Return_t wgt_tot=0.0;
  //bool usingWeight=UseWeight;
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int iw=0;
  for (; it!=it_end; ++it,++iw)
  {
    Walker_t& thisWalker(**it);
    Return_t* restrict saved = Records[iw];
    W.R=(*it)->R;
    W.update();
// buffer for MultiSlaterDet data
    Return_t logpsi;
    if(usebuffer=="yes")
    {
      Walker_t::Buffer_t& tbuffer=thisWalker.DataSetForDerivatives;
      logpsi=Psi.evaluateDeltaLog(W,tbuffer);
    }
    else
    {
      logpsi=Psi.evaluateDeltaLog(W);
    }
    //Return_t logpsi=Psi.evaluateDeltaLog(W);
    W.G += *dLogPsi[iw];
    W.L += *d2LogPsi[iw];
    Return_t eloc_new = H_KE.evaluate(W) + saved[ENERGY_FIXED];
    Return_t weight = std::exp(2.0*(logpsi-saved[LOGPSI_FREE])) ;
    //Return_t eloc_new=+saved[ENERGY_FIXED];
    //Return_t weight = UseWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
    std::vector<Return_t>* Dsaved= &(TempDerivRecords[iw]) ;
    std::vector<Return_t>* HDsaved= &(TempHDerivRecords[iw]) ;
    if (needGrad)
      Psi.evaluateDerivatives(W, OptVariables,*Dsaved,*HDsaved);
//    for(int l=0; l<NumOptimizables; l++) (*HDsaved)[l] += saved[ENERGY_FIXED] * (*Dsaved)[l];
    saved[ENERGY_NEW]=eloc_new;
    saved[REWEIGHT]=weight;
    wgt_tot+=weight;
  }
  //collect the total weight for normalization and apply maximum weight
  myComm->allreduce(wgt_tot);
  Return_t wgtnorm = (1.0*NumSamples)/wgt_tot;
  wgt_tot=0.0;
  int nw=W.getActiveWalkers();
  for (int iw=0; iw<nw; iw++)
  {
    Return_t* restrict saved = Records[iw];
    saved[REWEIGHT] = std::min(saved[REWEIGHT]*wgtnorm,MaxWeight) ;
    wgt_tot+= saved[REWEIGHT];
  }
  myComm->allreduce(wgt_tot);
  for (int i=0; i<SumValue.size(); i++)
    SumValue[i]=0.0;
  CSWeight=wgt_tot=1.0/wgt_tot;
  for (iw=0; iw<nw; iw++)
  {
    Return_t* restrict saved = Records[iw];
    Return_t weight=saved[REWEIGHT]*wgt_tot;
    Return_t eloc_new=saved[ENERGY_NEW];
    Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
    SumValue[SUM_E_BARE] += eloc_new;
    SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
    SumValue[SUM_ABSE_BARE] += delE;
    SumValue[SUM_E_WGT] += eloc_new*saved[REWEIGHT];
    SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*saved[REWEIGHT];
    SumValue[SUM_ABSE_WGT] += delE*saved[REWEIGHT];
    SumValue[SUM_WGT] += saved[REWEIGHT];
    SumValue[SUM_WGTSQ] += saved[REWEIGHT]*saved[REWEIGHT];
  }
  //collect everything
  myComm->allreduce(SumValue);
  return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
}

void
QMCCostFunctionSingle::getConfigurations(const std::string& aroot)
{
  app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << std::endl;
  app_log() << "    number of walkers before load " << W.getActiveWalkers() << std::endl;
  Timer t1;
  W.loadEnsemble();
  app_log() << "    Loading time = " << t1.elapsed() << std::endl;
  app_log() << "    number of walkers after load " << W.getActiveWalkers() << std::endl;
  //if(aroot.size() && aroot != "invalid") {
  //  app_log() << "  Reading configurations from the previous qmc block" << std::endl;
  //  HDFWalkerInputCollect wReader(aroot);
  //  wReader.putSingle(W);
  //}
  //
  if (dLogPsi.size() != W.getActiveWalkers())
  {
    delete_iter(dLogPsi.begin(),dLogPsi.end());
    delete_iter(d2LogPsi.begin(),d2LogPsi.end());
    int nptcl=W.getTotalNum();
    int nwtot=W.getActiveWalkers();
    dLogPsi.resize(nwtot,0);
    d2LogPsi.resize(nwtot,0);
    for (int i=0; i<nwtot; ++i)
      dLogPsi[i]=new ParticleGradient_t(nptcl);
    for (int i=0; i<nwtot; ++i)
      d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
  }
#if defined(QMCCOSTFUNCTION_DEBUG)
  if (debug_stream)
    delete debug_stream;
  char fname[16];
  sprintf(fname,"optdebug.p%d",OHMMS::Controller->mycontext());
  debug_stream = new std::ofstream(fname);
  debug_stream->setf(std::ios::scientific, std::ios::floatfield);
  debug_stream->precision(8);
  *debug_stream << "Initial : " << std::endl;
  for (int i=0; i<OptParams.size(); i++)
    *debug_stream << " " << IDtag[i] << "=" << OptParams[i];
  *debug_stream << std::endl;
#endif
}

// /** evaluate everything before optimization */
// void
// QMCCostFunctionSingle::getConfigurations(std::vector<std::string>& ConfigFile,
//     int partid, int nparts) {
//   if(ConfigFile.size()) {

//     app_log() << "  Reading configurations from mcwalkerset " << std::endl;

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
QMCCostFunctionSingle::checkConfigurations()
{
  //dG.resize(W.getTotalNum());
  //dL.resize(W.getTotalNum());
  int numLocWalkers=W.getActiveWalkers();
  Records.resize(numLocWalkers,6);
  if(usebuffer == "yes" )
    app_log() <<"Using buffers for temporary storage in QMCCostFunction.\n" << std::endl;
  TempHDerivRecords.resize(numLocWalkers,std::vector<Return_t>(NumOptimizables,0));
  TempDerivRecords.resize(numLocWalkers,std::vector<Return_t>(NumOptimizables,0));
  if ( !H_KE.getHamiltonian("Kinetic") )
  {
    H_KE.addOperator(H.getHamiltonian("Kinetic"),"Kinetic");
//       if ( includeNonlocalH=="yes" )
//       {
//         H_KE.addOperator(H.getHamiltonian("NonLocalECP"),"NonLocalECP");
//         app_log()<<" Adding Non local PP component for optimization"<< std::endl;
//       }
    H_KE.addObservables(W);
  }
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  //int totalElements=W.getTotalNum()*OHMMS_DIM;
  Etarget=0.0;
  Return_t e2sum=0.0;
  for (int iw=0; it!=it_end; ++it,++iw)
  {
    Walker_t& thisWalker(**it);
    W.R=thisWalker.R;
    W.update();
// buffer for MultiSlaterDet data
    Return_t*  saved=Records[iw];
    if(usebuffer=="yes")
    {
      Walker_t::Buffer_t& tbuffer=thisWalker.DataSetForDerivatives;
      Psi.registerDataForDerivatives(W,tbuffer);
      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iw], *d2LogPsi[iw],tbuffer);
    }
    else
    {
      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iw], *d2LogPsi[iw]);
    }
    Return_t e=H.evaluate(W);
    e2sum += e*e;
    Etarget += saved[ENERGY_TOT] = e;
//         if (includeNonlocalH=="yes")
//           saved[ENERGY_FIXED] = 0;
//         else
    saved[ENERGY_FIXED] = H.getLocalPotential();
    saved[REWEIGHT]=thisWalker.Weight=1.0;
    std::vector<Return_t>* Dsaved=  &(TempDerivRecords[iw]) ;
    std::vector<Return_t>* HDsaved=  &(TempHDerivRecords[iw]) ;
    Psi.evaluateDerivatives(W, OptVariables,*Dsaved,*HDsaved);
  }
  //Need to sum over the processors
  std::vector<Return_t> etemp(3);
  etemp[0]=Etarget;
  etemp[1]=static_cast<Return_t>(numLocWalkers);
  etemp[2]=e2sum;
  myComm->allreduce(etemp);
  Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << std::endl;
  app_log() << "  VMC Evar = " << etemp[2]/etemp[1]-Etarget*Etarget << std::endl;
  app_log() << "  Total weights = " << etemp[1] << std::endl;
  setTargetEnergy(Etarget);
  ReportCounter=0;
}

void QMCCostFunctionSingle::resetPsi(bool final_reset)
{
  if (OptVariables.size() < OptVariablesForPsi.size())
  {
    for (int i=0; i<equalVarMap.size(); ++i)
      OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
  }
  else
    for (int i=0; i<OptVariables.size(); ++i)
      OptVariablesForPsi[i]=OptVariables[i];
  if (final_reset)
  {
    Psi.stopOptimization();
    MCWalkerConfiguration::iterator it(W.begin());
    MCWalkerConfiguration::iterator it_end(W.end());
    for (; it!=it_end; ++it)
      (**it).DataSetForDerivatives.clear();
  }
  Psi.resetParameters(OptVariablesForPsi);
}


void QMCCostFunctionSingle::GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff)
{
  if (FiniteDiff>0)
  {
    Return_t dh=1.0/(2.0*FiniteDiff);
    for (int i=0; i<NumOptimizables ; i++)
    {
      for (int j=0; j<NumOptimizables; j++)
        OptVariables[j]=PM[j];
      OptVariables[i] = PM[i]+ FiniteDiff;
      Return_t CostPlus = this->Cost();
      OptVariables[i] = PM[i]- FiniteDiff;
      Return_t CostMinus = this->Cost();
      PGradient[i]= (CostPlus-CostMinus)*dh;
    }
  }
  else
  {
    for (int j=0; j<NumOptimizables; j++)
      OptVariables[j]=PM[j];
    resetPsi();
    //evaluate new local energies and derivatives
    NumWalkersEff=correlatedSampling();
    //Estimators::accumulate has been called by correlatedSampling
    curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
    Return_t curAvg2_w = curAvg_w*curAvg_w;
    curVar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-curAvg_w*curAvg_w;
    std::vector<Return_t> EDtotals(NumOptimizables,0.0);
    std::vector<Return_t> EDtotals_w(NumOptimizables,0.0);
    std::vector<Return_t> E2Dtotals_w(NumOptimizables,0.0);
    std::vector<Return_t> URV(NumOptimizables,0.0);
    std::vector<Return_t> HD_avg(NumOptimizables,0.0);
    Return_t wgtinv = 1.0/SumValue[SUM_WGT];
    Return_t delE_bar;
    int nw=W.getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      Return_t*  saved=Records[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      delE_bar += weight*std::pow(std::abs(eloc_new-EtargetEff),PowerE);
      std::vector<Return_t>  Dsaved= (TempDerivRecords[iw]);
      std::vector<Return_t>  HDsaved= (TempHDerivRecords[iw]);
      for (int pm=0; pm<NumOptimizables; pm++)
      {
        HD_avg[pm]+= HDsaved[pm];
      }
    }
    myComm->allreduce(HD_avg);
    myComm->allreduce(delE_bar);
    for (int pm=0; pm<NumOptimizables; pm++)
      HD_avg[pm] *= 1.0/static_cast<Return_t>(NumSamples);
    for (int iw=0; iw<nw; iw++)
    {
      Return_t*  saved=Records[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delta_l = (eloc_new-curAvg_w);
      bool ltz(true);
      if (eloc_new-EtargetEff<0)
        ltz=false;
      Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
      Return_t ddelE = PowerE*std::pow(std::abs(eloc_new-EtargetEff),PowerE-1);
      std::vector<Return_t> Dsaved= TempDerivRecords[iw];
      std::vector<Return_t> HDsaved= TempHDerivRecords[iw];
      for (int pm=0; pm<NumOptimizables; pm++)
      {
        EDtotals_w[pm] += weight*(HDsaved[pm] + 2.0*Dsaved[pm]*delta_l);
        URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
        if (ltz)
          EDtotals[pm]+= weight*(2.0*Dsaved[pm]*(delE-delE_bar) + ddelE*HDsaved[pm]);
        else
          EDtotals[pm] += weight*(2.0*Dsaved[pm]*(delE-delE_bar) - ddelE*HDsaved[pm]);
      }
    }
    myComm->allreduce(EDtotals_w);
    myComm->allreduce(EDtotals);
    myComm->allreduce(URV);
    Return_t smpinv=1.0/static_cast<Return_t>(NumSamples);
    for (int iw=0; iw<nw; iw++)
    {
      Return_t*  saved=Records[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delta_l = (eloc_new-curAvg_w);
      Return_t sigma_l = delta_l*delta_l;
      std::vector<Return_t> Dsaved= TempDerivRecords[iw];
      std::vector<Return_t> HDsaved= TempHDerivRecords[iw];
      for (int pm=0; pm<NumOptimizables; pm++)
      {
        E2Dtotals_w[pm] += weight*2.0*(Dsaved[pm]*(sigma_l-curVar_w) + delta_l*(HDsaved[pm]-EDtotals_w[pm]));
      }
    }
    myComm->allreduce(E2Dtotals_w);
    for (int pm=0; pm<NumOptimizables; pm++)
      URV[pm] *=smpinv;
    // app_log() << "Gradient:\n";
    for (int j=0; j<NumOptimizables; j++)
    {
      PGradient[j] = 0.0;
      if (std::abs(w_var) > 1.0e-10)
        PGradient[j] += w_var*E2Dtotals_w[j];
      if (std::abs(w_en)  > 1.0e-10)
        PGradient[j] += w_en*EDtotals_w[j];
      if (std::abs(w_w)   > 1.0e-10)
        PGradient[j] += w_w*URV[j];
      if (std::abs(w_abs) > 1.0e-10)
        PGradient[j] += w_abs*EDtotals[j];
    }
    IsValid=true;
    if (NumWalkersEff < MinNumWalkers*NumSamples)
    {
      ERRORMSG("CostFunction-> Number of Effective Walkers is too small " << NumWalkersEff<< "Minimum required"<<MinNumWalkers*NumSamples)
//      ERRORMSG("Going to stop now.")
      IsValid=false;
    }
  }
}

QMCCostFunctionSingle::Return_t QMCCostFunctionSingle::fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap)
{
//     Return_t NWE = NumWalkersEff=correlatedSampling();
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  Return_t curAvg2_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT];
  Return_t wgtinv = 1.0/SumValue[SUM_WGT];
  std::vector<Return_t> D_avg(NumParams(),0);
  int nw=W.getActiveWalkers();
  for (int iw=0; iw<nw; iw++)
  {
    Return_t*  saved=Records[iw];
    Return_t weight=saved[REWEIGHT]*wgtinv;
    std::vector<Return_t> Dsaved= TempDerivRecords[iw];
    for (int pm=0; pm<NumParams(); pm++)
    {
      D_avg[pm]+= Dsaved[pm]*weight;
    }
  }
  myComm->allreduce(D_avg);
  ///zero out matrices before we start
  for (int pm=0; pm<NumParams()+1; pm++)
  {
    for (int pm2=0; pm2<NumParams()+1; pm2++)
    {
      Overlap(pm,pm2)=0;
      Hamiltonian(pm,pm2)=0;
      Variance(pm,pm2)=0;
      H2(pm,pm2) = 0;
    }
  }
  for (int iw=0; iw<nw; iw++)
  {
    Return_t*  saved=Records[iw];
    Return_t weight=saved[REWEIGHT]*wgtinv;
    Return_t eloc_new=saved[ENERGY_NEW];
    std::vector<Return_t> Dsaved= TempDerivRecords[iw];
    std::vector<Return_t> HDsaved= TempHDerivRecords[iw];
    for (int pm=0; pm<NumParams(); pm++)
    {
      Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
      Return_t wfm = (HDsaved[pm] - 2.0*Dsaved[pm]*(eloc_new - curAvg_w) )*weight;
      Return_t wfe = (HDsaved[pm] + Dsaved[pm]*(eloc_new-curAvg_w) )*weight;
      H2(0,pm+1) += wfe*(eloc_new );
      H2(pm+1,0) += wfe*(eloc_new );
      Return_t vterm = HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
      Variance(0,pm+1) += vterm*weight;
      Variance(pm+1,0) += vterm*weight;
      Hamiltonian(0,pm+1) += wfe;
      Hamiltonian(pm+1,0) += wfd*(eloc_new-curAvg_w);
      for (int pm2=0; pm2<NumParams(); pm2++)
      {
        H2(pm+1,pm2+1) += wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
        Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
        Variance(pm+1,pm2+1) += wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
        Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);
      }
    }
  }
  myComm->allreduce(Hamiltonian);
  myComm->allreduce(Overlap);
  myComm->allreduce(Variance);
  myComm->allreduce(H2);
  Overlap(0,0) = 1;
  Hamiltonian(0,0) = curAvg_w ;
  H2(0,0) = curAvg2_w;
  Variance(0,0) = curAvg2_w - curAvg_w*curAvg_w;
  for (int pm=1; pm<NumParams()+1; pm++)
    for (int pm2=1; pm2<NumParams()+1; pm2++)
      Variance(pm,pm2) += Variance(0,0)*Overlap(pm,pm2);
  return 1.0;
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
