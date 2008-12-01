//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_DISTANCETABLE_H
#define QMCPLUSPLUS_DISTANCETABLE_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus {

  /** Class to manage multiple DistanceTableData objects.
   * 
   * \date  2008-09-19
   * static data members are removed. DistanceTable::add functions
   * are kept for compatibility only. New codes should use a member function
   * of ParticleSet to add a distance table
   * int ParticleSet::addTable(const ParticleSet& source)
   *
   * \deprecated There is only one instance of the data memebers of 
   * DistanceTable in an application and the data are shared by many objects.
   * Note that static data members and functions are used 
   * (based on singleton and factory patterns).
   *\todo DistanceTable should work as a factory, as well, to instantiate DistanceTableData
   * subject to different boundary conditions. 
   * Lattice/CrystalLattice.h and Lattice/CrystalLattice.cpp can be owned by DistanceTable 
   * to generically control the crystalline structure.
   */
  struct DistanceTable {
    
    ///add a named DistanceTableData_t of Symmectric type
    static DistanceTableData* add(ParticleSet& s);//, const char* aname = NULL);
    
    ///add a named DistanceTableData_t of Asymmectric type
    static DistanceTableData* add(const ParticleSet& s, ParticleSet& t);//, const char* aname = NULL);

  };

  ///free function to create a distable table of s-s
  DistanceTableData* createDistanceTable(ParticleSet& s);

  ///free function create a distable table of s-t
  DistanceTableData* createDistanceTable(const ParticleSet& s, ParticleSet& t);
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
