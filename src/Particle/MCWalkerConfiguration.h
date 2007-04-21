//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file MCWalkerConfiguration.h
 * @brief Declaration of a MCWalkerConfiguration
 */
#ifndef QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#define QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#include "Particle/ParticleSet.h"
#include "Particle/Walker.h"

namespace qmcplusplus {

  //Forward declaration
  struct MultiChain;

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

    /** clear the WalkerList without destroying them
     *
     * Provide std::vector::clear interface
     */
    inline void clear() { WalkerList.clear(); }

    /** insert elements
     * @param it locator where the inserting begins
     * @param first starting iterator
     * @param last ending iterator
     *
     * Provide std::vector::insert interface
     */
    template<class INPUT_ITER>
    inline void insert(iterator it, INPUT_ITER first, INPUT_ITER last) {
      WalkerList.insert(it,first,last);
    }

    /** add Walker_t* at the end
     * @param awalker pointer to a walker
     *
     * Provide std::vector::push_back interface
     */
    inline void push_back(Walker_t* awalker) {
      WalkerList.push_back(awalker);
    }

    /** delete the last Walker_t*
     *
     * Provide std::vector::pop_back interface
     */
    inline void pop_back() {
      delete WalkerList.back();
      WalkerList.pop_back();
    }

    inline Walker_t* operator[](int i) { return WalkerList[i];}
    inline const Walker_t* operator[](int i) const  { return WalkerList[i];}

    /** set LocalEnergy
     * @param e current average Local Energy
     */
    inline void setLocalEnergy(RealType e) { LocalEnergy = e;}

    /** return LocalEnergy
     */
    inline RealType getLocalEnergy() const {return LocalEnergy;}

    inline MultiChain* getPolymer(){
      return Polymer;
    }

    inline void setPolymer(MultiChain *chain){
      Polymer=chain;
    }

    /** reset the Walkers 
     */
    void reset();

    /**load a Walker_t to the current ParticleSet
     *@param awalker the reference to the walker to be loaded
     */
    void loadWalker(Walker_t& awalker);
    void clearAuxDataSet();
    bool createAuxDataSet(int nfield=256);
    //void registerData(Walker_t& awalker, PooledData<RealType>& buf);
    //void updateBuffer(Walker_t& awalker, PooledData<RealType>& buf);
    //void copyToBuffer(PooledData<RealType>& buf);
    //void copyFromBuffer(PooledData<RealType>& buf);

    void resetWalkerProperty(int ncopy=1);

  protected:

    bool OwnWalkers;

    bool ReadyForPbyP;

    int UpdateMode;

    RealType LocalEnergy;

    ///a collection of walkers
    WalkerList_t WalkerList;

   private:

     MultiChain *Polymer;

    /** initialize the PropertyList
     *
     * Add the named values of the default properties
    void initPropertyList();
     */
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
