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
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
FixedNumWalkers(false), SwapMode(0), Counter(0), 
  Nideal(nideal), NumGeneration(50), 
  MaxCopy(-1), Tau(tau), E_T(0.0), EavgSum(0.0), WgtSum(0.0), PopControl(0.1), 
  WalkerController(0)
{
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
  reset();
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

  app_log() << "  reference energy = " << E_T << endl;
  if(!fixW) {
    app_log() << "  target_walkers = " << Nideal << endl;
    app_log() << "  Max and mimum walkers per node= " << Nmax << " " << Nmin << endl;
    app_log() << "  number of generations (feedback) = " << NumGeneration << " ("<< Feed << ")"<< endl;
  }
}

/** update the weights and multiplicity of all the walkers
 * @param it starting walker
 * @param it_end ending walker
 * @return new trial energy
void
SimpleFixedNodeBranch::setWeights(MCWalkerConfiguration::iterator it, 
    MCWalkerConfiguration::iterator it_end) {
  while(it != it_end) {
    MCWalkerConfiguration::Walker_t& thisWalker(**it);
    RealType M=thisWalker.Weight;
    if(thisWalker.Age > MaxAge) M = std::min(0.5,M);
    else if(thisWalker.Age > 0) M = std::min(1.0,M);
    thisWalker.Multiplicity = M + Random();
    accumulate(thisWalker.Properties(LOCALENERGY),M);
    ++it; 
  }
  //return E_T = WalkerController->average(EavgSum,WgtSum)-Feed*log(static_cast<RealType>(nw))+logN;
}
*/

int SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones) {

  for(int i=0; i<clones.size(); i++) {
    Counter += clones[i]->Counter;
    EavgSum += clones[i]->EavgSum;
    WgtSum += clones[i]->WgtSum;
  }
  return WalkerController->branch(iter,w,PopControl);
}

void SimpleFixedNodeBranch::reset() {
  //WalkerController = CreateWalkerController(FixedNumWalkers, SwapMode, Nideal, Nmax, Nmin, WalkerController);
  //Nmax=WalkerController->Nmax;
  //Nmin=WalkerController->Nmin;
  //Feed=WalkerController->getFeedBackParameter(NumGeneration,Tau);
  Feed = 1.0/(static_cast<RealType>(NumGeneration)*Tau);
  logN = Feed*log(static_cast<RealType>(Nideal));
  app_log() << "  Current Counter = " << Counter << "\n  Trial Energy = " << E_T << endl;
  app_log() << "  Feedback parameter = " << Feed <<endl;
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
  string toswap("no");
  ParameterSet m_param;
  m_param.add(E_T,"ref_energy","AU"); m_param.add(E_T,"en_ref","AU");
  m_param.add(NumGeneration,"pop_control","int"); m_param.add(NumGeneration,"num_gen","int");
  m_param.add(MaxCopy,"max_branch","int"); m_param.add(MaxCopy,"max_copy","int");
  m_param.add(Nideal,"target_walkers","int");
  m_param.add(EavgSum,"energy_sum","AU");
  m_param.add(WgtSum,"weight_sum","none");
  m_param.add(toswap,"swap_walkers","string");
  m_param.add(PopControl,"swap_trigger","none");
  m_param.put(cur);

  if(toswap == "yes") {
    SwapMode = 1;
  } else {//check qmc/@collect
    const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"collect");
    if(t != NULL) {
      toswap=(const char*)t;
      SwapMode = (toswap == "yes");
    }
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
  hsize_t dim=3;
  vector<RealType> esave(3);
  /** stupid gnu compiler bug */
  esave[0] = (fabs(E_T) < numeric_limits<RealType>::epsilon())? EavgSum/WgtSum:E_T;
  esave[1]=EavgSum; esave[2]=WgtSum;
  
  if(append) {
    hid_t dataset = H5Dopen(grp,"Summary");
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&esave[0]);
    H5Dclose(dataset);
  } else {
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  H5Dcreate(grp, "Summary", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&esave[0]);
    H5Dclose(dataset);
    H5Sclose(dataspace);
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
    app_log() << "  Summary is found. BranchEngine is initialized"
      << "\n    E_T     = " << E_T 
      << "\n    EavgSum = " << EavgSum 
      << "\n    WgtSum  = " << WgtSum << endl;
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

  /*
int SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& W) {

  return WalkerController->branch(iter,W);
  int maxcopy=10;

  typedef MCWalkerConfiguration::Walker_t Walker_t;

  MCWalkerConfiguration::iterator it(W.begin());
  int iw=0, nw = W.getActiveWalkers();

  vector<Walker_t*> good, bad;
  vector<int> ncopy;
  ncopy.reserve(nw);

  int num_walkers=0;
  MCWalkerConfiguration::iterator it_end(W.end());
  while(it != it_end) {
    int nc = std::min(static_cast<int>((*it)->Multiplicity),maxcopy);
    if(nc) {
      num_walkers += nc;
      good.push_back(*it);
      ncopy.push_back(nc-1);
    } else {
      bad.push_back(*it);
    }
    iw++;it++;
  }

  //remove bad walkers
  for(int i=0; i<bad.size(); i++) delete bad[i];

  if(good.empty()) {
    ERRORMSG("All the walkers have died. Abort. ")
    OHMMS::Controller->abort();
  }

  //check if the projected number of walkers is too small or too large
  if(num_walkers>Nmax) {
    int nsub=0;
    int nsub_target=num_walkers-static_cast<int>(0.9*Nmax);
    int i=0;
    while(i<ncopy.size() && nsub<nsub_target) {
      if(ncopy[i]) {ncopy[i]--; nsub++;}
      i++;
    }
    num_walkers -= nsub;
  } else  if(num_walkers < Nmin) {
    int nadd=0;
    int nadd_target = static_cast<int>(Nmin*1.1)-num_walkers;
    if(nadd_target>good.size()) {
      WARNMSG("The number of walkers is running low. Requested walkers " << nadd_target << " good walkers = " << good.size())
    }
    int i=0;
    while(i<ncopy.size() && nadd<nadd_target) {
      ncopy[i]++; nadd++;i++;
    }
    num_walkers +=  nadd;
  }

  //WalkerControl
  //MPI Send to the master, MPI Irecv by the master
  //send the total number of walkers to the master
 
  //clear the WalkerList to populate them with the good walkers
  W.clear();
  W.insert(W.begin(), good.begin(), good.end());

  int cur_walker = good.size();
  for(int i=0; i<good.size(); i++) { //,ie+=ncols) {
    for(int j=0; j<ncopy[i]; j++, cur_walker++) {
      W.push_back(new Walker_t(*(good[i])));
    }
  }

  int nw_tot = W.getActiveWalkers();

  //WalkerControl
  //Master check if criteria is met and send back 0/1, total walkers, max, min

  //if(swap) nw_tot= swapWalkers();
  if(SwapWalkers) gsum(nw_tot,0);

  //set Weight and Multiplicity to default values
  iw=0;
  it=W.begin();
  it_end=W.end();
  while(it != it_end) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }

  return nw_tot;
   //return w.branch(10,Nmax,Nmin,SwapWalkers);
}
  */
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

