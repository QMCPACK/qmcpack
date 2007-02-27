//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Jeongnim Kim
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
/** @file SparseSparseLocalizsedBasisSet.h
 * @author Jeongnim Kim
 * @brief A derived class from BasisSetBase which implements a sparse vector method to handle locazed basis sets.
 */
#ifndef QMCPLUSPLUS_SPARSELOCALIZEDBASISSET_H
#define QMCPLUSPLUS_SPARSELOCALIZEDBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus {

  /** A localized basis set derived from BasisSetBase<typename COT::value_type>
   *
   * This class performs the evaluation of the basis functions and their
   * derivatives for each of the N-particles in a configuration. 
   * The template parameter COT denotes Centered-Orbital-Type which provides
   * a set of localized orbitals associated with a center. Unlike LocalizedBasisSet<COT>
   * this class does not assume that all the centers have basis functions, as the name
   * SparseLocalizedBasisSet implies.
   */
  template<class COT>
  struct SparseLocalizedBasisSet: public BasisSetBase<typename COT::value_type> 
  {
    typedef BasisSetBase<typename COT::value_type> BasisSetType;
    typedef typename BasisSetType::RealType      RealType;
    typedef typename BasisSetType::ValueType     ValueType;
    typedef typename BasisSetType::IndexType     IndexType;
    typedef typename BasisSetType::IndexVector_t IndexVector_t;
    typedef typename BasisSetType::ValueVector_t ValueVector_t;
    typedef typename BasisSetType::ValueMatrix_t ValueMatrix_t;
    typedef typename BasisSetType::GradVector_t  GradVector_t;
    typedef typename BasisSetType::GradMatrix_t  GradMatrix_t;
    typedef typename BasisSetType::OptimizableSetType OptimizableSetType;

    using BasisSetType::BasisSetSize;
    using BasisSetType::Phi;
    using BasisSetType::dPhi;
    using BasisSetType::d2Phi;
    using BasisSetType::Y;
    using BasisSetType::dY;
    using BasisSetType::d2Y;
    ///Reference to the center
    const ParticleSet& CenterSys;
    ///number of centers, e.g., ions
    int NumCenters;
    ///number of quantum particles
    int NumTargets;
    ///number of groups
    int NumGroups;
    ///BasisOffset
    vector<int> BasisOffset;

    /** class to handle a COT and particle indices which share the COT
     */
    struct LOForCenter
    {
      ///Centered Orbital Type for this center
      COT* myO;
      ///IDs for the particle which belong to the same center
      vector<int> myP;
      inline LOForCenter(COT* o=0):myO(o) {}
    };

    typedef vector<LOForCenter*> ContainerType;
    ContainerType LOBasisSet;

    /** distance table, e.g., ion-electron
     *
     * Localized basis sets require a pair relationship between CenterSys 
     * and the quantum particle set. 
     */
    const DistanceTableData* myTable;

    /** constructor
     * @param ions ionic system
     * @param els electronic system
     */
    SparseLocalizedBasisSet(ParticleSet& ions, ParticleSet& els): 
      CenterSys(ions), myTable(0)
    { 
      myTable = DistanceTable::add(ions,els);
      NumCenters=CenterSys.getTotalNum();
      NumTargets=els.getTotalNum();
      NumGroups=CenterSys.getSpeciesSet().getTotalNum();
      BasisOffset.resize(NumCenters+1);
      LOBasisSet.resize(NumGroups,0);
    }

    /**
     @param atable the distance table (ion-electron)
     @brief Assign the distance table (ion-electron) and
     *determine the total number of basis states.
    */
    void setBasisSetSize(int nbs) { 
      if(nbs == BasisSetSize) return;

      if(myTable ==0) {
        app_error() << "SparseLocalizsedBasisSet cannot function without a distance table. Abort" << endl;
      }

      for(int ig=0; ig<NumGroups; ig++)
        if(LOBasisSet[ig]) LOBasisSet[ig]->myO->setTable(myTable);

      BasisOffset[0]=0;
      for(int c=0; c<NumCenters; c++)
      {
        int ig=CenterSys.GroupID[c];
        if(LOBasisSet[ig])
        {
          LOBasisSet[ig]->myP.push_back(c);//add this center to the basis set group
          BasisOffset[c+1] = BasisOffset[c]+LOBasisSet[ig]->myO->getBasisSetSize();
        }
        else
        {
          BasisOffset[c+1] = BasisOffset[c];
        }
      }
      BasisSetSize = BasisOffset[NumCenters];

      for(int ig=0; ig<NumGroups; ig++)
      {
        if(LOBasisSet[ig])
        {
          cout << "### SparseLOBasis " << ig << " : ";
          const vector<int>& plist(LOBasisSet[ig]->myP);
          for(int i=0; i<plist.size(); i++)
          {
            cout << plist[i] << "[" << BasisOffset[plist[i]] << "] ";
          }
          cout << endl;
        }
      }

      this->resize(NumTargets);
    }

    void resetParameters(OptimizableSetType& optVariables) 
    {
      for(int ig=0; ig<NumGroups; ig++)
        if(LOBasisSet[ig]) LOBasisSet[ig]->myO->resetParameters(optVariables);
    }
    
    /** reset the distance table with a new target P
     */
    void resetTargetParticleSet(ParticleSet& P) {
      LOGMSG("SparseLocalizsedBasisSet::resetTargetParticleSet")
      myTable = DistanceTable::add(CenterSys,P);
      for(int ig=0; ig<NumGroups; ig++)
        if(LOBasisSet[ig]) LOBasisSet[ig]->myO->setTable(myTable);
    }

    inline void 
    evaluateForWalkerMove(const ParticleSet& P) 
    {
      for(int ig=0; ig<NumGroups; ig++)
      {
        if(LOBasisSet[ig])
        {
          const vector<int>& plist(LOBasisSet[ig]->myP);
          COT& o(*(LOBasisSet[ig]->myO));
          for(int i=0; i<plist.size(); i++)
          {
            o.evaluateForWalkerMove(plist[i],0,NumTargets,BasisOffset[plist[i]],Y,dY,d2Y);
          }
        }
      }
    }

    inline void 
    evaluateForWalkerMove(const ParticleSet& P, int iat) 
    {
      for(int ig=0; ig<NumGroups; ig++)
      {
        if(LOBasisSet[ig])
        {
          const vector<int>& plist(LOBasisSet[ig]->myP);
          COT& o(*(LOBasisSet[ig]->myO));
          for(int i=0; i<plist.size(); i++)
          {
            o.evaluateForWalkerMove(plist[i],iat,BasisOffset[plist[i]],Phi,dPhi,d2Phi);
          }
        }
      }
    }

    inline void 
    evaluateForPtclMove(const ParticleSet& P, int iat)  
    {
      for(int ig=0; ig<NumGroups; ig++)
      {
        if(LOBasisSet[ig])
        {
          const vector<int>& plist(LOBasisSet[ig]->myP);
          COT& o(*(LOBasisSet[ig]->myO));
          for(int i=0; i<plist.size(); i++)
          {
            o.evaluateForPtclMove(plist[i],iat,BasisOffset[plist[i]],Phi);
          }
        }
      }
    }

    inline void 
    evaluateAllForPtclMove(const ParticleSet& P, int iat)  
    {
      for(int ig=0; ig<NumGroups; ig++)
      {
        if(LOBasisSet[ig])
        {
          const vector<int>& plist(LOBasisSet[ig]->myP);
          COT& o(*(LOBasisSet[ig]->myO));
          for(int i=0; i<plist.size(); i++)
          {
            o.evaluateAllForPtclMove(plist[i],iat,BasisOffset[plist[i]],Phi,dPhi,d2Phi);
          }
        }
      }
    }
    

    /** add a new set of Centered Atomic Orbitals
     * @param icenter the index of the center
     * @param aos a set of Centered Atomic Orbitals
     */
    void add(int icenter, COT* aos) {
      aos->setTable(myTable);
      LOBasisSet[icenter]=new LOForCenter(aos);
    }
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1817 $   $Date: 2007-02-26 12:37:58 -0600 (Mon, 26 Feb 2007) $
 * $Id: SparseLocalizsedBasisSet.h 1817 2007-02-26 18:37:58Z jnkim $ 
 ***************************************************************************/

