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
#include "OhmmsPETE/OhmmsMatrix.h"
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
    typedef list<Walker_t*>              WalkerList_t;
    typedef WalkerList_t::iterator         iterator;
    typedef WalkerList_t::const_iterator   const_iterator;

    Matrix<ValueType> Energy;

    ///default constructor
    MCWalkerConfiguration();

    ///default destructor
    ~MCWalkerConfiguration();

    void createWalkers(int n);

    iterator destroyWalker(iterator walkerit);
    iterator destroyWalker(iterator itstart, iterator itend);
    void copyWalker(iterator walkerit, int n);

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

    void setUpdateMode(int updatemode) { UpdateMode = updatemode;}


    inline void setLocalEnergy(RealType e) { LocalEnergy = e;}
    inline RealType getLocalEnergy() const {return LocalEnergy;}

    int branch(int maxcopy, int Nmax, int Nmin);

    void clear();

    void copy(iterator first, iterator last);

    void reset();

  protected:

    RealType LocalEnergy;

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
