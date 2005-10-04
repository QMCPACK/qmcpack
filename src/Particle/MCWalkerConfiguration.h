//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_QMC_MCWALKERCONFIGURATION_H
#define OHMMS_QMC_MCWALKERCONFIGURATION_H
#include "Particle/ParticleSet.h"
#include "Particle/Walker.h"

namespace ohmmsqmc {

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
  class MCWalkerConfiguration: public ParticleSet {

  public:

    /**enumeration for update*/
    enum {Update_All = 0, ///move all the active walkers
	  Update_Walker,  ///move a walker by walker
	  Update_Particle ///move a particle by particle
    };
    
    ///type of Walker class
    typedef Walker<RealType,ParticlePos_t> Walker_t;
    ///container type of the Properties of a Walker
    typedef Walker_t::PropertyContainer_t  PropertyContainer_t;
    ///container type of Walkers
    typedef std::vector<Walker_t*>         WalkerList_t;
    ///iterator of Walker container
    typedef WalkerList_t::iterator         iterator;
    ///const_iterator of Walker container
    typedef WalkerList_t::const_iterator   const_iterator;

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

    /** destroy Walkers from itstart to itend
     *@param first starting iterator of the walkers
     *@param last ending iterator of the walkers
     */
    iterator destroyWalkers(iterator first, iterator last);

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
    inline int getActiveWalkers() const { return WalkerList.size();}

    ///return the number of particles per walker
    inline int getParticleNum() const { return R.size();}

    /// return the first iterator
    inline iterator begin() { return WalkerList.begin();}
    /// return the last iterator, [begin(), end())
    inline iterator end() { return WalkerList.end();}

    /// return the first const_iterator
    inline const_iterator begin() const { return WalkerList.begin();}

    /// return the last const_iterator  [begin(), end())
    inline const_iterator end() const { return WalkerList.end();}
    /**@}*/

    /** set LocalEnergy
     * @param e current average Local Energy
     */
    inline void setLocalEnergy(RealType e) { LocalEnergy = e;}

    /** return LocalEnergy
     */
    inline RealType getLocalEnergy() const {return LocalEnergy;}

    /** destroy/create walkers  using Multiplicity of the Walkers
     * @param maxcopy maximum number of copies per Walker
     * @param Nmax maximum number of Walkers 
     * @param Nmin minimum number of Walkers
     * @return the current number of Walkers after branching
     *
     * Branching algorithm determines Walker_t::Multiplicity and branch
     * function destroy and copy walkers based on the Multiplicity
     */
    int branch(int maxcopy, int Nmax, int Nmin);

    /** reset the Walkers 
     */
    void reset();

    /**load a Walker_t to the current ParticleSet
     *@param awalker the reference to the walker to be loaded
     */
    void loadWalker(Walker_t& awalker);
    bool createAuxDataSet(int nfield=256);
    void registerData(Walker_t& awalker, PooledData<RealType>& buf);
    void updateBuffer(Walker_t& awalker, PooledData<RealType>& buf);
    void copyToBuffer(PooledData<RealType>& buf);
    void copyFromBuffer(PooledData<RealType>& buf);

    void resetWalkerProperty(int ncopy=1);

    /** swap walkers among MPI nodes
     * @return total number of walkers
     */
#if defined(HAVE_MPI)
    int swapWalkers();
#else
    inline int swapWalkers() { return WalkerList.size();}
#endif

  protected:

    bool OwnWalkers;

    bool ReadyForPbyP;

    int UpdateMode;

    RealType LocalEnergy;

    ///a collection of walkers
    WalkerList_t WalkerList;

   private:


    /** initialize the PropertyList
     *
     * Add the named values of the default properties
     */
    void initPropertyList();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
