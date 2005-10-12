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
#ifndef OHMMS_QMC_DISTANCETABLE_H
#define OHMMS_QMC_DISTANCETABLE_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace ohmmsqmc {

  class DistanceTableData;

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
    
    typedef ParticleSet::ParticleLayout_t    ParticleLayout_t;
    typedef ParticleSet::Scalar_t            RealType;
    typedef ParticleSet::SingleParticlePos_t PosType;

    enum {SUPERCELL_OPEN=0, SUPERCELL_WIRE=1, 
      SUPERCELL_SLAB=3, SUPERCELL_BULK=7};

    ///add a named DistanceTableData_t of Symmectric type
    static int add(ParticleSet& s, const char* aname = NULL);
    
    ///add a named DistanceTableData_t of Asymmectric type
    static int add(const ParticleSet& s, ParticleSet& t, const char* aname = NULL);
    
    /** returns a pointer to a DistanceTableData_t
     * @param i index to access TableList
     */
    static DistanceTableData* getTable(int i){
      return TableList[i];
    }

    /** returns a pointer to a DistanceTableData of the pair
     *@param atable name of the distance table
     */
    static DistanceTableData* getTable(const char* atable) {
      map<string,int>::iterator it = TableMap.find(atable);
      if(it == TableMap.end()) 
	return NULL;
      else 
        return TableList[(*it).second];
    }

    /** returns the pointer to SimulationCell
     */
    static const ParticleLayout_t* getSimulationCell() {
      return SimulationCell;
    }

    /** create a global SimulationCell referenced by ParticleSet objects
     * @param cur xml node
     */
    static void createSimulationCell(xmlNodePtr cur);

    /** select DistanceTableData objects whose visitor tag matches ptag
     *@param ptag the tag of a ParticleSet
     *@param tables The objects related to the particle set
     */
    static void getTables(int ptag, vector<DistanceTableData*>& tables);

    ///reset the internal values, mainly Updated flags to prepare new series
    static void reset();

    static void create(int walkers);
    
    ///return true if ith table has been updated
    static bool updated(int i) { return Updated[i];}
    
    //static void update(ParticleSet& t);

    //static void update(WalkerSetRef& t);

    static void registerData(PooledData<RealType>& buf);

    static void copyFromBuffer(PooledData<RealType>& buf);

    static void copyToBuffer(PooledData<RealType>& buf);

  private:
    
    ///a list of update flags
    static vector<bool> Updated;
    
    ///a list of DistanceTableData_t
    static vector<DistanceTableData*> TableList;

    static vector<int> VisitorID;

    ///Center_map[name] returns the Table index
    static map<string,int>  TableMap;
    
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
