//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_WALKER_CONTROL_BASE_H
#define QMCPLUSPLUS_WALKER_CONTROL_BASE_H

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Message/MPIObjectBase.h"
#include "Message/CommOperators.h"
// #include "QMCDrivers/ForwardWalking/ForwardWalkingStructure.h"

//#include <boost/archive/binary_oarchive.hpp>

namespace qmcplusplus
{

/** Base class to control the walkers for DMC simulations.
 *
 * The virtual functions are implemented for a serial execution with a usual birth/death
 * process. Inherited classes implement other WalkerControl algorithms by implementing
 * branch function.
 * - 2006/10/18
 * -- curData and accumData are added.
 * -- branch function handles averaging of the local energies to apply the weight over the
 *   current walker sample correctly.
 * -- removes an extra global reduction.
 *   --- curData should be used for global reductions
 * - 2007/11/07
 * -- added functions to dump dmc histograms
 */
class WalkerControlBase: public MPIObjectBase
{

public:
  ///typedef of Walker_t
  typedef MCWalkerConfiguration::Walker_t  Walker_t;
  ///typedef of RealType
  typedef QMCTraits::EstimatorRealType     RealType;
  ///typedef of IndexType
  typedef QMCTraits::IndexType             IndexType;

  ///@enum enumeration to access curData and accumData for reduction
  enum {ENERGY_INDEX=0
        , ENERGY_SQ_INDEX
        , WALKERSIZE_INDEX
        , WEIGHT_INDEX
        , EREF_INDEX
        , R2ACCEPTED_INDEX
        , R2PROPOSED_INDEX
        , FNSIZE_INDEX
        , RNONESIZE_INDEX
        , RNSIZE_INDEX
        , B_ENERGY_INDEX
        , B_WGT_INDEX
        , SENTWALKERS_INDEX
        , LE_MAX
       };

  ///id for the method
  IndexType MyMethod;
  ///context id
  IndexType MyContext;
  ///number of contexts
  IndexType NumContexts;
  ///0 is default
  IndexType SwapMode;
  ///minimum number of walkers
  IndexType Nmin;
  ///maximum number of walkers
  IndexType Nmax;
  ///maximum copy per walker
  IndexType MaxCopy;
  ///current number of walkers per processor
  IndexType NumWalkers;
  ///Number of walkers created by this node
  IndexType NumWalkersCreated;
  ///Number of walkers sent during the exchange
  IndexType NumWalkersSent;
  ///trial energy energy
  RealType trialEnergy;
  ///target average energy
  RealType targetAvg;
  ///target average variance
  RealType targetVar;
  ///target sigma to limit fluctuations of the trial energy
  RealType targetSigma;
  ///bound of the energy window
  RealType targetEnergyBound;
  ///current variance
  RealType curVar;
  ///number of particle per node
  std::vector<int> NumPerNode;
  ///offset of the particle index
  std::vector<int> OffSet;
  ///offset of the particle index for a fair distribution
  std::vector<int> FairOffSet;

  ///ensenble properties
  MCDataType<RealType> EnsembleProperty;

  ///filename for dmc.dat
  std::string dmcFname;
  ///file to save energy histogram
  std::ofstream* dmcStream;
  ///archive
  //boost::archive::binary_oarchive *oa;

  ///any accumulated data over a block
  std::vector<RealType> accumData;
  ///any temporary data
  std::vector<RealType> curData;
  ///temporary storage for good and bad walkers
  std::vector<Walker_t*> good_w, bad_w;
  ///temporary storage for copy counters
  std::vector<int> ncopy_w;
  ///Add released-node fields to .dmc.dat file
  bool WriteRN;
  ///Use non-blocking isend/irecv
  bool use_nonblocking;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControlBase(Communicate* c, bool rn=false);

  /** empty destructor to clean up the derived classes */
  virtual ~WalkerControlBase();

  /** start a block */
  void start();

  /** start controller  and initialize the IDs of walkers*/
  void setWalkerID(MCWalkerConfiguration& walkers);

  /** take averages and writes to a file */
  void measureProperties(int iter);

  /** set the trial energy
   */
  inline void setTrialEnergy(RealType et)
  {
    trialEnergy=et;
  }

  /** return a value accumulated during a block
   * @param i index of the data
   *
   * use enum for i, see DMCEnergyEstimator
   */
  inline RealType getValue(int i)
  {
    return accumData[i];
  }

  /** return a current value
   * @param i index of the data
   *
   * use enum for i, see DMCEnergyEstimator
   */
  inline RealType getCurrentValue(int i)
  {
    return curData[i];
  }

  /** set the target average and variance
   */
  inline void setEnergyAndVariance(RealType e, RealType v)
  {
    trialEnergy=e;
    targetAvg=e;
    targetVar=v;
  }

  /** update properties without branching */
  int doNotBranch(int iter, MCWalkerConfiguration& W);

  /** sort Walkers between good and bad and prepare branching
   */
  int sortWalkers(MCWalkerConfiguration& W);
  /** copy good walkers to W
   */
  int copyWalkers(MCWalkerConfiguration& W);

  void Write2XYZ(MCWalkerConfiguration& W);

  /** reset to accumulate data */
  virtual void reset();

  /** perform branch and swap walkers as required */
  virtual int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  virtual RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 1.0/(static_cast<RealType>(ngen)*tau);
  }

  bool put(xmlNodePtr cur);

  void setMinMax(int nw_in, int nmax_in);

//     struct ForwardWalkingData
//     {
//       typedef TinyVector<float,DIM>       StoredPosType;
//       typedef ParticleAttrib<StoredPosType>       StoredPosVector;
//       long ID;
//       long ParentID;
//       StoredPosVector Pos;
//
//       inline ForwardWalkingData()
//       {
//       }
//
//       inline ForwardWalkingData(const Walker_t& a)
//       {
//         Pos.resize(a.R.size());
//         Pos = a.R;
//         ID = a.ID;
//         ParentID = a.ParentID;
//       }
//
//       inline ForwardWalkingData(const ForwardWalkingData& a)
//       {
//         Pos.resize(a.Pos.size());
//         Pos = a.Pos;
//         ID = a.ID;
//         ParentID = a.ParentID;
//       }
//
//       inline ForwardWalkingData(const int a)
//       {
//         Pos.resize(a);
//       }
//
//       inline ForwardWalkingData& operator=(const Walker_t& a) {
//         Pos.resize(a.R.size());
//         Pos = a.R;
//         ID = a.ID;
//         ParentID = a.ParentID;
//         return *this;
//       }
//
//       inline ForwardWalkingData& operator=(const ForwardWalkingData& a)
//       {
//         Pos.resize(a.Pos.size());
//         Pos = a.Pos;
//         ID = a.ID;
//         ParentID = a.ParentID;
//         return *this;
//       }
//
//       inline int SizeOf()
//       {
//         return sizeof(long)*2 + Pos.size()*DIM*sizeof(float);
//       }
//
//     };



//     typedef std::vector<ForwardWalkingData> ForwardWalkingConfiguration;
//     std::vector<ForwardWalkingConfiguration> ForwardWalkingHistory;
//     inline void storeConfigsForForwardWalking(MCWalkerConfiguration& W)
//     {
//       std::vector<ForwardWalkingData> ForwardWalkingHere;
//
//       for(std::vector<Walker_t*>::iterator Wit(W.begin()); Wit != W.end(); Wit++ )
//       {
//         ForwardWalkingData fwhere( *(*Wit) );
//         ForwardWalkingHere.push_back(fwhere);
//       }
//
//       ForwardWalkingHistory.push_back(ForwardWalkingHere);
//     }
//
//     inline void clearConfigsForForwardWalking()
//     {
//       ForwardWalkingHistory.clear();
//     }
//
//     inline int sizeOfConfigsForForwardWalking()
//     {
//       int szeFW(0);
//       int singleSize = ForwardWalkingHistory[0][0].SizeOf();
//       for(int i=0;i<ForwardWalkingHistory.size();i++) szeFW += ForwardWalkingHistory[i].size() * singleSize;
//       return szeFW;
//     }
//
//     inline void layoutOfConfigsForForwardWalking(std::vector<int>& returnVal)
//     {
//       returnVal.resize(ForwardWalkingHistory.size()+1,0);
//       returnVal[0]=0;
//       for(int i=0;i<ForwardWalkingHistory.size();i++) returnVal[i+1]=ForwardWalkingHistory[i].size();
//     }



};


;


}
#endif

