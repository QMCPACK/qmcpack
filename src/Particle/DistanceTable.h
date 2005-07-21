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

  class ParticleSet;
  class WalkerSetRef;
  class DistanceTableData;

  /** Class to store tabulated data for particle distances, cosine vectors and inverses.
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
    
    ///derive the real type from ParticleSet::Scalar_t
    typedef ParticleSet::Scalar_t RealType;
    typedef ParticleSet::SingleParticlePos_t PosType;

    ///add a named DistanceTableData_t of Symmectric type
    static int add(const ParticleSet& s, const char* aname = NULL);
    
    ///add a named DistanceTableData_t of Asymmectric type
    static int add(const ParticleSet& s, const ParticleSet& t, 
		   const char* aname = NULL);
    
    /*!\fn DistanceTableData_t* getTable(int i)
     * \param i index to access TableList
     * \brief returns a pointer to a DistanceTableData_t
     */
    static DistanceTableData* getTable(int i){
      return TableList[i];
    }

    /*!\fn DistanceTableData_t* getTable(const char* atable)
     * \param atable name of the distance table
     * \brief returns a pointer to a DistanceTableData_t of the pair
     */
    static DistanceTableData* getTable(const char* atable) {
      map<string,int>::iterator it = TableMap.find(atable);
      if(it == TableMap.end()) 
	return NULL;
      else 
        return TableList[(*it).second];
    }

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
    
    static void update(ParticleSet& t);

    static void update(WalkerSetRef& t);

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
