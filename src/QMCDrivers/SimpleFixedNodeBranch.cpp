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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include <numeric>
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "QMCDrivers/WalkerControlBase.h"
#include "Estimators/ScalarEstimatorManager.h"
#include "Estimators/DMCEnergyEstimator.h"

namespace qmcplusplus {

SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
FixedNumWalkers(false), QMCCounter(-1), SwapMode(0), Counter(0), 
  Nideal(nideal), NumGeneration(50), 
  MaxCopy(-1), Tau(tau), E_T(0.0), EavgSum(0.0), WgtSum(0.0), PopControl(0.1), 
  WalkerController(0), MyEstimator(0), SwapWalkers("yes")
{
  registerParameters();
  reset();
}

/** copy constructor
 *
 * Copy only selected data members and WalkerController is never copied.
 */
SimpleFixedNodeBranch::SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch):
  FixedNumWalkers(false), SwapMode(0), Counter(0), 
  Nideal(abranch.Nideal), NumGeneration(abranch.NumGeneration), 
  MaxCopy(-1), Tau(abranch.Tau), E_T(abranch.E_T), 
  EavgSum(0.0), WgtSum(0.0), PopControl(abranch.PopControl), 
  WalkerController(0)
{
  registerParameters();
  reset();
}

void SimpleFixedNodeBranch::registerParameters() {
  m_param.add(E_T,"refEnergy","AU");
  m_param.add(NumGeneration,"popControl","int"); 
  m_param.add(MaxCopy,"maxCopy","int"); 
  m_param.add(Nideal,"targetWalkers","int"); 
  m_param.add(EavgSum,"energySum","AU"); 
  m_param.add(WgtSum,"weightSum","none");
  m_param.add(PopControl,"swapTrigger","none");
  m_param.add(SwapWalkers,"swapWalkers","string");
  m_param.add(SwapWalkers,"collect","string");

  //backward compatability
  m_param.add(E_T,"ref_energy","AU"); m_param.add(E_T,"en_ref","AU");
  m_param.add(NumGeneration,"pop_control","int"); 
  m_param.add(NumGeneration,"num_gen","int"); 
  m_param.add(Nideal,"target_walkers","int");
  m_param.add(PopControl,"swap_trigger","none");
  m_param.add(SwapWalkers,"swap_walkers","string");
}

void SimpleFixedNodeBranch::setEstimatorManager(ScalarEstimatorManager* est) {
  MyEstimator=est;
}

void SimpleFixedNodeBranch::initWalkerController(RealType tau, bool fixW) {
  Tau=tau;
  reset();
  if(WalkerController == 0) {
    FixedNumWalkers=fixW;
    if(fixW) {Feed=0.0;logN=0.0;}
    WalkerController = CreateWalkerController(FixedNumWalkers, SwapMode, Nideal, Nmax, Nmin, WalkerController);
    Nmax=WalkerController->Nmax;
    Nmin=WalkerController->Nmin;
  }

  if(fixW) {
    ETrialIndex=-1;
    E_T=0.0;
  } else {
    ETrialIndex = MyEstimator->addColumn("Etrial");
  }

  DMCEnergyEstimator* dmcE= dynamic_cast<DMCEnergyEstimator*>(MyEstimator->getMainEstimator());
  if(dmcE) {
    dmcE->setWalkerControl(WalkerController);
  } else {
    OHMMS::Controller->abort("SimpleFixedNodeBranch::initWalkerController failed."); 
  }

  app_log() << "  QMC counter      = " << QMCCounter << endl;
  app_log() << "  reference energy = " << E_T << endl;
  if(!fixW) {
    app_log() << "  target_walkers = " << Nideal << endl;
    app_log() << "  Max and mimum walkers per node= " << Nmax << " " << Nmin << endl;
    app_log() << "  number of generations (feedback) = " << NumGeneration << " ("<< Feed << ")"<< endl;
  }

}

void SimpleFixedNodeBranch::flush(int counter) 
{
  if(counter==0) 
    if(WalkerController) WalkerController->reset();
  Counter=counter;
}

void
SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& w) {
  int pop_now = WalkerController->branch(iter,w,PopControl);
  EavgSum+=WalkerController->getCurrentValue(WalkerControlBase::EREF_INDEX);
  WgtSum+=WalkerController->getCurrentValue(WalkerControlBase::WALKERSIZE_INDEX);
  if(ETrialIndex>0) {
    E_T = (EavgSum/WgtSum+E_T)*0.5-Feed*log(static_cast<RealType>(pop_now))+logN;
    MyEstimator->setColumn(ETrialIndex,E_T);
  }
}

int SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones) {
  return WalkerController->branch(iter,w,PopControl);
}

void SimpleFixedNodeBranch::reset() {
  Feed = 1.0/(static_cast<RealType>(NumGeneration)*Tau);
  logN = Feed*log(static_cast<RealType>(Nideal));
  app_log() << "  Current Counter = " << Counter << "\n  Trial Energy = " << E_T << endl;
  app_log() << "  Feedback parameter = " << Feed <<endl;
}

void SimpleFixedNodeBranch::finalize() {
  MyEstimator->finalize();
  if(!WalkerController || ETrialIndex<0) {
    MyEstimator->getEnergyAndWeight(EavgSum,WgtSum);
    E_T=EavgSum/WgtSum;
  }
}

/**  Parse the xml file for parameters
 *@param cur current xmlNode 
 *@param LogOut ostream to which the run-time report is sent
 *
 * Few important parameters are:
 * <ul>
 * <li> en_ref: a reference energy
 * <li> num_gen: number of generations \f$N_G\f$ to reach  equilibrium, used in the feedback parameter
 * \f$ feed = \frac{1}{N_G \tau} \f$ 
 * </ul>
 */
bool SimpleFixedNodeBranch::put(xmlNodePtr cur){

  //check dmc/vmc and decide to create WalkerControllerBase
  m_param.put(cur);

  //set the SwapMode
  SwapMode = (SwapWalkers=="yes");

  /* \warning{backward compatability qmc/@collect="yes|no"}
   */
  const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"collect");
  if(t != NULL) {
    SwapMode = xmlStrEqual(t,(const xmlChar*)"yes");
  } 

  //overwrite the SwapMode with the number of contexts
  SwapMode = (SwapMode && (OHMMS::Controller->ncontexts()>1));

  if(SwapMode) {
    app_log() << "  Collective handling of Average/walkers" << endl;
  } else {
    app_log() << "  Node-parallel handling of Average/walkers" << endl;
  }

  reset();

  return true;
}

#if defined(HAVE_LIBHDF5)
void SimpleFixedNodeBranch::write(hid_t grp, bool append) {

  ////This is the case with the VMC.
  //if(!WalkerController) {
  //  MyEstimator->getEnergyAndWeight(EavgSum,WgtSum);
  //  E_T=EavgSum/WgtSum;
  //}
  hsize_t dim=3;
  vector<RealType> esave(dim);
  /** stupid gnu compiler bug */
  //esave[0] = (fabs(E_T) < numeric_limits<RealType>::epsilon())?  EavgSum/WgtSum:E_T;
  esave[0] = E_T;
  esave[1] = EavgSum; 
  esave[2] = WgtSum;
  
  if(append) {
    hid_t dataset = H5Dopen(grp,"Summary");
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&esave[0]);
    H5Dclose(dataset);

    if(LogNorm.size()){
      //lognorm
      dataset = H5Dopen(grp,"LogNorm");
      ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&LogNorm[0]);
      H5Dclose(dataset);
    }
  } else {
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  H5Dcreate(grp, "Summary", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&esave[0]);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    //create lognrom
    if(LogNorm.size()){
      dim=LogNorm.size();
      dataspace  = H5Screate_simple(1, &dim, NULL);
      dataset =  H5Dcreate(grp, "LogNorm", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&LogNorm[0]);
      H5Dclose(dataset);
      H5Sclose(dataspace);
    }
  }
}

void SimpleFixedNodeBranch::read(hid_t grp) {
  herr_t status = H5Eset_auto(NULL, NULL);
  status = H5Gget_objinfo (grp, "Summary", 0, NULL);
  if(status == 0) {
    hid_t dataset = H5Dopen(grp, "Summary");
    vector<RealType> esave(3);
    hsize_t ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &esave[0]);
    H5Dclose(dataset);

    E_T=esave[0]; EavgSum=esave[1]; WgtSum=esave[2];

    //Will add reading LogNorm
    bool normFound(0);
    if(LogNorm.size()) {
      dataset = H5Dopen(grp, "LogNorm");
      ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &LogNorm[0]);
      H5Dclose(dataset);
      normFound=1;
    }
    app_log() << "  Summary is found. BranchEngine is initialized"
      << "\n    E_T     = " << E_T 
      << "\n    EavgSum = " << EavgSum 
      << "\n    WgtSum  = " << WgtSum << endl;
    if(normFound)app_log() << " Normalization factor defined from previous calculation " << endl;
  } else {
    app_log() << "  Summary is not found. Starting from scratch" << endl;
  }
  
}

void SimpleFixedNodeBranch::read(const string& fname) {
  string h5file = fname;
  string ext=getExtension(h5file);
  if(ext != "h5") { //if the filename does not h5 extension, add the extension
    h5file.append(".config.h5");
  }
  hid_t h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  hid_t h_config = H5Gopen(h_file,"config_collection");
  read(h_config);
  H5Gclose(h_config);
  H5Fclose(h_file);
}
#else
void SimpleFixedNodeBranch::write(hid_t grp, bool append) { }
void SimpleFixedNodeBranch::read(hid_t grp) {}
void SimpleFixedNodeBranch::read(const string& fname) {}
#endif
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

