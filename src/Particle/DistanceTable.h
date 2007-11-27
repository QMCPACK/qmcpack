//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_DISTANCETABLE_H
#define QMCPLUSPLUS_DISTANCETABLE_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus {

  /**@ingroup nnlist
   * @brief Class to manage multiple DistanceTableData objects.
   * 
   * There is only one instance of the data memebers of 
   * DistanceTable in an application and the data are shared by many objects.
   * Note that static data members and functions are used 
   * (based on singleton and factory patterns).
   *\todo DistanceTable should work as a factory, as well, to instantiate DistanceTableData
   * subject to different boundary conditions. 
   * Lattice/CrystalLattice.h and Lattice/CrystalLattice.cpp can be owned by DistanceTable 
   * to generically control the crystalline structure.
   */
  class DistanceTable {
    
  public:
    
    enum {DIM=ParticleSet::DIM};

    typedef ParticleSet::ParticleLayout_t    ParticleLayout_t;
    typedef ParticleSet::Scalar_t            RealType;
    typedef ParticleSet::SingleParticlePos_t PosType;

    ///add a named DistanceTableData_t of Symmectric type
    static DistanceTableData* add(ParticleSet& s, const char* aname = NULL);
    
    ///add a named DistanceTableData_t of Asymmectric type
    static DistanceTableData* add(const ParticleSet& s, ParticleSet& t, const char* aname = NULL);

    /** returns the pointer to SimulationCell
     */
    static const ParticleLayout_t* getSimulationCell() {
      return SimulationCell;
    }

    /** create a global SimulationCell referenced by ParticleSet objects
     * @param cur xml node
     */
    static void createSimulationCell(xmlNodePtr cur);


    ///reset the internal values, mainly Updated flags to prepare new series
    static void reset();

    /** resize the containers
     * @param walkers number of walkers
     */
    static void create(int walkers);
    
    ///remove the distance table with the name
    static void removeTable(const string& tname);

  private:
    
    ///Center_map[name] returns DistanceTableData*
    static map<string,DistanceTableData*> TableMap;
    
    ///Global object to define a simulation cell
    static ParticleLayout_t* SimulationCell;

    /// Default constructor.
    DistanceTable(){ }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
