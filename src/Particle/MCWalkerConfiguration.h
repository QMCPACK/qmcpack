//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file MCWalkerConfiguration.h
 * @brief Declaration of a MCWalkerConfiguration
 */
#ifndef QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#define QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#include "Particle/ParticleSet.h"
#include "Particle/Walker.h"
#include "Utilities/IteratorUtility.h"
//#include "Particle/Reptile.h"

namespace qmcplusplus
{

//Forward declaration
class MultiChain;
class MCSample;
class HDFWalkerOutput;
class Reptile;

/** A set of walkers that are to be advanced by Metropolis Monte Carlo.
 *
 *As a derived class from ParticleSet, MCWalkerConfiguration interacts with
 *QMCHamiltonian and TrialWaveFunction as a ParticleSet, while QMCDrivers
 *use it as multiple walkers whose configurations are advanced according
 to MC algorithms.
 *
 Each walker is represented by Walker<PosVector_t> and
 *MCWalkerConfiguration contains a list of
 *the walkers.  This class enables two possible moves:
 *<ul>
 *<li> move the entire active walkers, similarly to molecu. Suitable for
 *small and big systems with a small time step.
 *<li> move a particle for each walker. Suitable for large systems.

 *</ul>
 */
class MCWalkerConfiguration: public ParticleSet
{

public:

  /**enumeration for update*/
  enum {Update_All = 0, ///move all the active walkers
        Update_Walker,  ///move a walker by walker
        Update_Particle ///move a particle by particle
       };

  ///container type of the Properties of a Walker
  typedef Walker_t::PropertyContainer_t  PropertyContainer_t;
  ///container type of Walkers
  typedef std::vector<Walker_t*>         WalkerList_t;
  ///iterator of Walker container
  typedef WalkerList_t::iterator         iterator;
  ///const_iterator of Walker container
  typedef WalkerList_t::const_iterator   const_iterator;
 
  typedef std::vector<Reptile*>          ReptileList_t;
   /** starting index of the walkers in a processor group
   *
   * WalkerOffsets[0]=0 and WalkerOffsets[WalkerOffsets.size()-1]=total number of walkers in a group
   * WalkerOffsets[processorid+1]-WalkerOffsets[processorid] is equal to the number of walkers on a processor,
   * i.e., W.getActiveWalkers().
   * WalkerOffsets is added to handle parallel I/O with hdf5
   */
  std::vector<int> WalkerOffsets;

  // Data for GPU-acceleration via CUDA
  // These hold a list of pointers to the positions, gradients, and
  // laplacians for each walker.  These vectors .data() is often
  // passed to GPU kernels.
#ifdef QMC_CUDA

  gpu::device_vector<CudaRealType*>  RList_GPU;
  gpu::device_vector<CudaValueType*> GradList_GPU, LapList_GPU;
  // First index is the species.  The second index is the walker
  std::vector<gpu::device_vector<CUDA_PRECISION_FULL*> > RhokLists_GPU;
  gpu::device_vector<CudaValueType*>  DataList_GPU;
  gpu::device_vector<CudaPosType>    Rnew_GPU;
  gpu::host_vector<CudaPosType>      Rnew_host;
  std::vector<PosType>                    Rnew;
  gpu::device_vector<CudaRealType*>  NLlist_GPU;
  gpu::host_vector<CudaRealType*>    NLlist_host;
  gpu::host_vector<CudaRealType*>    hostlist;
  gpu::host_vector<CudaValueType*>   hostlist_valueType;
  gpu::host_vector<CUDA_PRECISION_FULL*> hostlist_AA; 
  gpu::host_vector<CudaPosType>      R_host;
  gpu::host_vector<CudaGradType>     Grad_host;
  gpu::device_vector<int> iatList_GPU;
  gpu::host_vector<int> iatList_host;
  gpu::device_vector<int> AcceptList_GPU;
  gpu::host_vector<int> AcceptList_host;

  void allocateGPU(size_t buffersize);
  void copyWalkersToGPU(bool copyGrad=false);
  void copyWalkerGradToGPU();
  void updateLists_GPU();
  int CurrentParticle;
  void proposeMove_GPU
  (std::vector<PosType> &newPos, int iat);
  void acceptMove_GPU(std::vector<bool> &toAccept);
  void NLMove_GPU (std::vector<Walker_t*> &walkers,
                   std::vector<PosType> &Rnew,
                   std::vector<int> &iat);
#endif

  ///default constructor
  MCWalkerConfiguration();

  ///default constructor: copy only ParticleSet
  MCWalkerConfiguration(const MCWalkerConfiguration& mcw);

  ///default destructor
  ~MCWalkerConfiguration();

  /** create numWalkers Walkers
   *
   * Append Walkers to WalkerList.
   */
  void createWalkers(int numWalkers);
  /** create walkers
   * @param first walker iterator
   * @param last walker iterator
   */
  void createWalkers(iterator first, iterator last);
  /** copy walkers
   * @param first input walker iterator
   * @param last input walker iterator
   * @param start first target iterator
   *
   * No memory allocation is allowed.
   */
  void copyWalkers(iterator first, iterator last, iterator start);

  /** destroy Walkers from itstart to itend
   *@param first starting iterator of the walkers
   *@param last ending iterator of the walkers
   */
  iterator destroyWalkers(iterator first, iterator last);

  /** destroy Walkers
   *@param nw number of walkers to be destroyed
   */
  void destroyWalkers(int nw);

  /** copy the pointers to the Walkers to WalkerList
   * @param head pointer to the head walker
   * @param tail pointer to the tail walker
   *
   * Special function introduced to work with Reptation method.
   * Clear the current WalkerList and add two walkers, head and tail.
   * OwnWalkers are set to false.
   */
  void copyWalkerRefs(Walker_t* head, Walker_t* tail);

  ///clean up the walker list and make a new list
  void resize(int numWalkers, int numPtcls);

  ///make random moves for all the walkers
  //void sample(iterator first, iterator last, value_type tauinv);
  ///make a random move for a walker
  void sample(iterator it, RealType tauinv);

  ///return the number of active walkers
  inline int getActiveWalkers() const
  {
    return WalkerList.size();
  }
  ///return the total number of active walkers among a MPI group
  inline int getGlobalNumWalkers() const
  {
    return GlobalNumWalkers;
  }
  ///return the total number of active walkers among a MPI group
  inline void setGlobalNumWalkers(int nw)
  {
    GlobalNumWalkers=nw;
    EnsembleProperty.NumSamples=nw;
    EnsembleProperty.Weight=nw;
  }

  inline void setWalkerOffsets(const std::vector<int>& o)
  {
    WalkerOffsets=o;
  }

  ///return the number of particles per walker
  inline int getParticleNum() const
  {
    return R.size();
  }

  /// return the first iterator
  inline iterator begin()
  {
    return WalkerList.begin();
  }
  /// return the last iterator, [begin(), end())
  inline iterator end()
  {
    return WalkerList.end();
  }

  /// return the first const_iterator
  inline const_iterator begin() const
  {
    return WalkerList.begin();
  }

  /// return the last const_iterator  [begin(), end())
  inline const_iterator end() const
  {
    return WalkerList.end();
  }
  /**@}*/

  /** clear the WalkerList without destroying them
   *
   * Provide std::vector::clear interface
   */
  inline void clear()
  {
    WalkerList.clear();
  }

  /** insert elements
   * @param it locator where the inserting begins
   * @param first starting iterator
   * @param last ending iterator
   *
   * Provide std::vector::insert interface
   */
  template<class INPUT_ITER>
  inline void insert(iterator it, INPUT_ITER first, INPUT_ITER last)
  {
    WalkerList.insert(it,first,last);
  }

  /** add Walker_t* at the end
   * @param awalker pointer to a walker
   *
   * Provide std::vector::push_back interface
   */
  inline void push_back(Walker_t* awalker)
  {
    WalkerList.push_back(awalker);
  }

  ///resize Walker::PropertyHistory and Walker::PHindex:
  void resizeWalkerHistories();

  /** delete the last Walker_t*
   *
   * Provide std::vector::pop_back interface
   */
  inline void pop_back()
  {
    delete WalkerList.back();
    WalkerList.pop_back();
  }

  inline Walker_t* operator[](int i)
  {
    return WalkerList[i];
  }

  inline const Walker_t* operator[](int i) const
  {
    return WalkerList[i];
  }

  /** set LocalEnergy
   * @param e current average Local Energy
   */
  inline void setLocalEnergy(RealType e)
  {
    LocalEnergy = e;
  }

  /** return LocalEnergy
   */
  inline RealType getLocalEnergy() const
  {
    return LocalEnergy;
  }

  inline MultiChain* getPolymer()
  {
    return Polymer;
  }

  inline void setPolymer(MultiChain *chain)
  {
    Polymer=chain;
  }

  /** reset the Walkers
   */
  void reset();

  void resetWalkerProperty(int ncopy=1);

  inline bool updatePbyP() const
  {
    return ReadyForPbyP;
  }

  //@{save/load/clear function for optimization
  inline int numSamples() const
  {
    return CurSampleCount;
  }
  ///set the number of max samples
  void setNumSamples(int n);
  ///save the position of current walkers to SampleStack
  void saveEnsemble();
  ///save the position of current walkers
  void saveEnsemble(iterator first, iterator last);
  /// load a single sample from SampleStack
  void loadSample(ParticleSet::ParticlePos_t &Pos, size_t iw) const;
  /** load SampleStack data to current walkers
   */
  void loadEnsemble();
  //void loadEnsemble(const Walker_t& wcopy);
  /** load SampleStack from others
    */
  void loadEnsemble(std::vector<MCWalkerConfiguration*>& others, bool doclean=true);
  /** dump Samples to a file
   * @param others MCWalkerConfigurations whose samples will be collected
   * @param out engine to write the samples to state_0/walkers
   * @param np number of processors
   * @return true with non-zero samples
   */
  bool dumpEnsemble(std::vector<MCWalkerConfiguration*>& others, HDFWalkerOutput* out, int np, int nBlock);
  ///clear the ensemble
  void clearEnsemble();
  //@}

  template<typename ForwardIter>
  inline void putConfigurations(ForwardIter target)
  {
    int ds=OHMMS_DIM*TotalNum;
    for(iterator it=WalkerList.begin(); it!= WalkerList.end(); ++it,target+=ds)
    {
      copy(get_first_address((*it)->R),get_last_address((*it)->R),target);
    }
  }

  inline void resetWalkerParents( )
  {
    for(iterator it=WalkerList.begin(); it!= WalkerList.end(); ++it )
    {
      (*it)->ParentID = (*it)->ID;
    }
  }

  ///overwrite make_clones function
  void make_clones(int n);

protected:

  ///boolean for cleanup
  bool OwnWalkers;
  ///true if the buffer is ready for particle-by-particle updates
  bool ReadyForPbyP;
  ///number of walkers on a node
  int LocalNumWalkers;
  ///number of walkers shared by a MPI group
  int GlobalNumWalkers;
  ///update-mode index
  int UpdateMode;

  RealType LocalEnergy;

public:
  ///a collection of walkers
  WalkerList_t WalkerList;
  ///a collection of reptiles contained in MCWalkerConfiguration.
  ReptileList_t ReptileList;
  Reptile* reptile;

private:

  MultiChain *Polymer;

  int MaxSamples;
  int CurSampleCount;
  //add samples
  std::vector<MCSample*> SampleStack;

  /** initialize the PropertyList
   *
   * Add the named values of the default properties
  void initPropertyList();
   */
};
}
#endif
