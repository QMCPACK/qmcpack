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
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/PooledData.h"
#include "Particle/Walker.h"
#include <list>

namespace ohmmsqmc {

  /** A set of walkers that are to be advanced by Metropolis Monte Carlo.  
   *
   *Each walker is represented by Walker<PosVector_t> and 
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
    
    typedef Walker<RealType,ParticlePos_t> Walker_t;
    typedef Walker_t::PropertyContainer_t  PropertyContainer_t;
    typedef vector<Walker_t*>              WalkerList_t;
    typedef WalkerList_t::iterator         iterator;
    typedef WalkerList_t::const_iterator   const_iterator;
    typedef PooledData<RealType>           WalkerData_t;

    std::vector<WalkerData_t*> DataSet;
    Matrix<ValueType> Energy;

    ///default constructor
    MCWalkerConfiguration();

    ///default constructor: create 2 null walkers
    MCWalkerConfiguration(const MCWalkerConfiguration& mcw, int nw=2);

    ///default destructor
    ~MCWalkerConfiguration();

    void createWalkers(int n);
    iterator destroyWalker(iterator walkerit);
    iterator destroyWalker(iterator itstart, iterator itend);
    void copyWalker(iterator walkerit, int n);

    /** copy the pointers to the Walkers to WalkerList */
    void copyEndWalkers(Walker_t* first, Walker_t* last) {
      WalkerList[0]=first;
      WalkerList[1]=last;
    }

    void addWalkers(vector<int>& Copies,
		    vector<Walker_t*>& InactiveList);

    ///clean up the walker list and make a new list
    void resize(int nw, int np);

    ///make random moves for all the walkers
    //void sample(iterator first, iterator last, value_type tauinv);
    ///make a random move for a walker
    void sample(iterator it, RealType tauinv);

    ///return the number of active walkers
    inline int getActiveWalkers() const { return WalkerList.size();}

    ///return the number of particles per walker
    inline int getParticleNum() const { return R.size();}

    /**@defgroup iterators to access WalkerList
     *@{*/
    /// return the first iterator
    inline iterator begin() { return WalkerList.begin();}
    /// return the last iterator, [begin(), end())
    inline iterator end() { return WalkerList.end();}

    /// return the first const_iterator
    inline const_iterator begin() const { return WalkerList.begin();}

    /// return the last const_iterator  [begin(), end())
    inline const_iterator end() const { return WalkerList.end();}

    /// return the active iterator
    inline iterator theWalker() { return WorkingWalker;}

    /// return the active const_iterator
    inline const_iterator theWalker() const { return WorkingWalker;}
    /**@}*/

    void setUpdateMode(int updatemode);

    inline void setLocalEnergy(RealType e) { LocalEnergy = e;}
    inline RealType getLocalEnergy() const {return LocalEnergy;}

    int branch(int maxcopy, int Nmax, int Nmin);
    int branch2(int maxcopy, int Nmax, int Nmin);

    void clear();

    void copy(iterator first, iterator last);

    void reset();


    /**load a Walker_t to the current ParticleSet
     *@param awalker the reference to the walker to be loaded
     */
    void loadWalker(Walker_t& awalker);

    /**move a particle
     *@param iat the index of the particle to be moved
     *@param newpos new position of the iat-th particle
     */
    SingleParticlePos_t makeMove(int iat, const SingleParticlePos_t& displ);

    /**accept the move
     *@param iat the index of the particle whose position and other attributes to be updated
     */
    void acceptMove(int iat);

    bool createAuxDataSet(int nfield=256);
    void registerData(Walker_t& awalker, PooledData<RealType>& buf);
    void copyToBuffer(PooledData<RealType>& buf);
    void copyFromBuffer(PooledData<RealType>& buf);

  protected:

    RealType LocalEnergy;

    ///the position of the active particle for particle-by-particle moves
    SingleParticlePos_t activePos;

    ///the indexp of the active particle for particle-by-particle moves
    Index_t             activePtcl;

    int UpdateMode;

    ///list of Walker
    WalkerList_t WalkerList;

    ///save a WorkingWalker
    iterator WorkingWalker;

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
