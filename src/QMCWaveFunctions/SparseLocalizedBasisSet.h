//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SparseSparseLocalizsedBasisSet.h
 * @author Jeongnim Kim
 * @brief A derived class from BasisSetBase which implements a sparse vector method to handle locazed basis sets.
 */
#ifndef QMCPLUSPLUS_SPARSELOCALIZEDBASISSET_H
#define QMCPLUSPLUS_SPARSELOCALIZEDBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

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

  using BasisSetType::BasisSetSize;
  using BasisSetType::Phi;
  using BasisSetType::dPhi;
  using BasisSetType::d2Phi;
  using BasisSetType::grad_grad_Phi;
  using BasisSetType::Y;
  using BasisSetType::dY;
  using BasisSetType::d2Y;
  ///Reference to the center
  const ParticleSet& CenterRef;
  ///number of centers, e.g., ions
  int NumCenters;
  ///number of quantum particles
  int NumTargets;
  ///number of groups
  int NumGroups;
  ///BasisOffset
  std::vector<int> BasisOffset;

  /** class to handle a COT and particle indices which share the COT
   */
  struct LOForCenter
  {
    ///Centered Orbital Type for this center
    COT* myO;
    ///IDs for the particle which belong to the same center
    std::vector<int> myP;
    inline LOForCenter(COT* o=0):myO(o) {}
    LOForCenter(const LOForCenter& old):myP(old.myP)
    {
      myO=old.myO->makeClone();
    }
  };

  typedef std::vector<LOForCenter*> ContainerType;
  ContainerType LOBasisSet;

  /** distance table, e.g., ion-electron
   *
   * Localized basis sets require a pair relationship between CenterRef
   * and the quantum particle set.
   */
  const DistanceTableData* myTable;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SparseLocalizedBasisSet(const ParticleSet& ions, ParticleSet& els):
    CenterRef(ions), myTable(0)
  {
    myTable = DistanceTable::add(ions,els,DT_AOS);
    NumCenters=CenterRef.getTotalNum();
    NumTargets=els.getTotalNum();
    NumGroups=CenterRef.getSpeciesSet().getTotalNum();
    BasisOffset.resize(NumCenters+1);
    LOBasisSet.resize(NumGroups,0);
  }

  BasisSetBase<typename COT::value_type>* makeClone() const
  {
    SparseLocalizedBasisSet<COT>* myclone=new SparseLocalizedBasisSet<COT>(*this);
    for(int i=0; i<LOBasisSet.size(); ++i)
      myclone->LOBasisSet[i]=new LOForCenter(*LOBasisSet[i]);
    return myclone;
  }

  void checkInVariables(opt_variables_type& active)
  {
    for(int i=0; i<LOBasisSet.size(); ++i)
      LOBasisSet[i]->myO->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    for(int i=0; i<LOBasisSet.size(); ++i)
      LOBasisSet[i]->myO->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<LOBasisSet.size(); ++i)
      LOBasisSet[i]->myO->resetParameters(active);
  }

  /**
   @param atable the distance table (ion-electron)
   @brief Assign the distance table (ion-electron) and
   *determine the total number of basis states.
  */
  void setBasisSetSize(int nbs)
  {
    if(nbs == BasisSetSize)
      return;
    if(myTable ==0)
    {
      APP_ABORT("SparseLocalizsedBasisSet cannot function without a distance table. Abort");
    }
    for(int ig=0; ig<NumGroups; ig++)
      if(LOBasisSet[ig])
        LOBasisSet[ig]->myO->setTable(myTable);
    BasisOffset[0]=0;
    for(int c=0; c<NumCenters; c++)
    {
      int ig=CenterRef.GroupID[c];
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
    //for(int ig=0; ig<NumGroups; ig++)
    //{
    //  if(LOBasisSet[ig])
    //  {
    //    std::cout << "### SparseLOBasis " << ig << " : ";
    //    const std::vector<int>& plist(LOBasisSet[ig]->myP);
    //    for(int i=0; i<plist.size(); i++)
    //    {
    //      std::cout << plist[i] << "[" << BasisOffset[plist[i]] << "] ";
    //    }
    //    std::cout << std::endl;
    //  }
    //}
    this->resize(NumTargets);
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  {
    LOGMSG("SparseLocalizsedBasisSet::resetTargetParticleSet")
    myTable = DistanceTable::add(CenterRef,P,DT_AOS);
    for(int ig=0; ig<NumGroups; ig++)
      if(LOBasisSet[ig])
        LOBasisSet[ig]->myO->setTable(myTable);
  }

  inline void
  evaluateWithHessian(const ParticleSet& P, int iat)
  {
    APP_ABORT("SparseLocalizsedBasisSet::evaluateWithHessian has not been implemented. \n");
  }

  void evaluateWithThirdDeriv(const ParticleSet& P, int iat)
  {
    APP_ABORT("SparseLocalizsedBasisSet::evaluateWithThirdDeriv has not been implemented. \n");
  }

  void evaluateThirdDerivOnly(const ParticleSet& P, int iat)
  {
    APP_ABORT("SparseLocalizsedBasisSet::evaluateThirdDerivOnly has not been implemented. \n");
  }

  inline void
  evaluateForWalkerMove(const ParticleSet& P)
  {
    for(int ig=0; ig<NumGroups; ig++)
    {
      if(LOBasisSet[ig])
      {
        const std::vector<int>& plist(LOBasisSet[ig]->myP);
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
        const std::vector<int>& plist(LOBasisSet[ig]->myP);
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
        const std::vector<int>& plist(LOBasisSet[ig]->myP);
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
        const std::vector<int>& plist(LOBasisSet[ig]->myP);
        COT& o(*(LOBasisSet[ig]->myO));
        for(int i=0; i<plist.size(); i++)
        {
          o.evaluateAllForPtclMove(plist[i],iat,BasisOffset[plist[i]],Phi,dPhi,d2Phi);
        }
      }
    }
  }

  inline void
  evaluateForPtclMoveWithHessian(const ParticleSet& P, int iat)
  {
    APP_ABORT("SparseLocalizsedBasisSet::evaluateForPtclMoveWithHessian has not been implemented. \n");
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    aos->setTable(myTable);
    LOBasisSet[icenter]=new LOForCenter(aos);
  }
};

}
#endif


