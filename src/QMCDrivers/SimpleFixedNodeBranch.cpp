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
#include <numeric>
#include "OhmmsData/ParameterSet.h"
using namespace qmcplusplus;

SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
SwapMode(0), Counter(0), Nideal(nideal), NumGeneration(50), 
  MaxCopy(-1), Tau(tau), E_T(0.0), EavgSum(0.0), WgtSum(0.0), PopControl(0.1),
WalkerController(0)
{
  Feed = 1.0/(static_cast<RealType>(NumGeneration)*Tau);
  logN = Feed*log(static_cast<RealType>(Nideal));
}

void SimpleFixedNodeBranch::reset() {

  WalkerController = CreateWalkerController(SwapMode, Nideal, Nmax, Nmin, WalkerController);
  Nmax=WalkerController->Nmax;
  Nmin=WalkerController->Nmin;
  Feed = 1.0/(static_cast<RealType>(NumGeneration)*Tau);
  logN = Feed*log(static_cast<RealType>(Nideal));
  LOGMSG("Current Counter = " << Counter << " Trial Energy = " << E_T)
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
bool SimpleFixedNodeBranch::put(xmlNodePtr cur, OhmmsInform *LogOut){
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
    LogOut->getStream() << "# Walkers are swapped to balance the load " << endl;
    SwapMode = 1;
  }

  reset();
  LogOut->getStream() << "#target_walkers = " << Nideal << endl;
  LogOut->getStream() << "#Max and mimum walkers per node= " << Nmax << " " << Nmin << endl;
  LogOut->getStream() << "#reference energy = " << E_T << endl;
  LogOut->getStream() << "#number of generations (feedbac) = " << NumGeneration 
    << " ("<< Feed << ")"<< endl;
  return true;
}


void SimpleFixedNodeBranch::write(hid_t grp) {
  hsize_t dim=3;
  vector<RealType> esave(3);
  /** stupid gnu compiler bug */
  esave[0] = (fabs(E_T) < numeric_limits<RealType>::epsilon())? EavgSum/WgtSum:E_T;
  esave[1]=EavgSum; esave[2]=WgtSum;
  hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
  hid_t dataset =  
    H5Dcreate(grp, "Summary", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t ret = 
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&esave[0]);
  H5Dclose(dataset);
  H5Sclose(dataspace);
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
    LOGMSG("Summary is found \n BranchEngine is initialized  E_T=" << E_T << " EavgSum="<<EavgSum << " WgtSum=" << WgtSum)
  } else {
    LOGMSG("Summary is not found. Starting from scratch")
  }
}

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

