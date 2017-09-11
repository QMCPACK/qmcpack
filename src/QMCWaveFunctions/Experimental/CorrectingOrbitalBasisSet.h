//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file CorrectingOrbitalBasisSet.h
 * @brief A derived class from BasisSetBase
 *
 * This class calculates correcting orbitals (s-type NGO functions)
 * which are meant to go on top of a modified LCOrbitalSet
 * Used in cusp correction algorithms
 */
#ifndef QMCPLUSPLUS_CORRECTINGORBITALBASISSET_H
#define QMCPLUSPLUS_CORRECTINGORBITALBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/** A localized basis set derived from BasisSetBase<typename COT::value_type>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
template<class COT>
struct CorrectingOrbitalBasisSet: public BasisSetBase<typename COT::value_type>
{
  typedef BasisSetBase<typename COT::value_type> BasisSetType;
  typedef typename BasisSetType::RealType      RealType;
  typedef typename BasisSetType::ValueType     ValueType;
  typedef typename BasisSetType::IndexType     IndexType;
  typedef typename BasisSetType::HessType      HessType;
  typedef typename BasisSetType::IndexVector_t IndexVector_t;
  typedef typename BasisSetType::ValueVector_t ValueVector_t;
  typedef typename BasisSetType::ValueMatrix_t ValueMatrix_t;
  typedef typename BasisSetType::GradVector_t  GradVector_t;
  typedef typename BasisSetType::GradMatrix_t  GradMatrix_t;
  typedef typename BasisSetType::HessVector_t  HessVector_t;
  typedef typename BasisSetType::HessMatrix_t  HessMatrix_t;


  using BasisSetType::ActivePtcl;
  using BasisSetType::Counter;
  using BasisSetType::BasisSetSize;
  using BasisSetType::Phi;
  using BasisSetType::dPhi;
  using BasisSetType::d2Phi;
  using BasisSetType::grad_grad_Phi;
  using BasisSetType::grad_grad_grad_Phi;
  using BasisSetType::Y;
  using BasisSetType::dY;
  using BasisSetType::d2Y;
  ///Reference to the center
  const ParticleSet& CenterSys;
  ///number of centers, e.g., ions
  int NumCenters;
  ///number of quantum particles
  int NumTargets;

  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  std::vector<int>  BasisOffset;

  /** container of the pointers to the Atomic Orbitals
   *
   * size of LOBasis =  number  of centers (e.g., ions)
   * AO[i] returns a Centered Orbital for an ion i
   */
  std::vector<COT*> LOBasis;

  /** distance table, e.g., ion-electron
   *
   * CorrectingOrbitalBasisSet basis sets require a pair relationship between CenterSys
   * and the quantum particle set.
   */
  const DistanceTableData* myTable;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  CorrectingOrbitalBasisSet(const DistanceTableData* tbl, ParticleSet &ions, int nEls): myTable(tbl), CenterSys(ions), NumTargets(nEls)
  {
    NumCenters=CenterSys.getTotalNum();
    LOBasis.resize(NumCenters,0);
    BasisOffset.resize(NumCenters+1);
  }

  CorrectingOrbitalBasisSet<COT>* makeClone() const
  {
    CorrectingOrbitalBasisSet<COT>* myclone = new CorrectingOrbitalBasisSet<COT>(*this);
    for(int i=0; i<LOBasis.size(); ++i)
    {
      COT* cc=LOBasis[i]->makeClone();
      myclone->LOBasis[i]=cc;
    }
    myclone->setBasisSetSize(BasisSetSize);
    return myclone;
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
      app_error() << "CorrectingOrbitalBasisSet cannot function without a distance table. Abort" << std::endl;
    }
    //reset the distance table for the atomic orbitals
    for(int i=0; i<LOBasis.size(); i++)
      LOBasis[i]->setTable(myTable);
    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for(int c=0; c<NumCenters; c++)
      BasisOffset[c+1] = BasisOffset[c]+LOBasis[c]->getBasisSetSize();
    BasisSetSize = BasisOffset[NumCenters];
    this->resize(NumTargets);
  }

  void resetParameters(const opt_variables_type& active)
  {
    //reset each unique basis functions
    for(int i=0; i<LOBasis.size(); i++)
      LOBasis[i]->resetParameters(active);
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  {
    myTable = DistanceTable::add(CenterSys,P,DT_AOS);
    for(int i=0; i<LOBasis.size(); i++)
      LOBasis[i]->setTable(myTable);
  }

  inline void
  evaluateForWalkerMove(const ParticleSet& P)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,0,NumTargets,BasisOffset[c],Y,dY,d2Y);
    Counter++; // increment a conter
  }

  inline void
  evaluateForWalkerMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
    Counter++;
  }

  inline void
  evaluateForPtclMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForPtclMove(c,iat,BasisOffset[c],Phi);
    Counter++;
    ActivePtcl=iat;
  }

  inline void
  evaluateAllForPtclMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateAllForPtclMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
    Counter++;
    ActivePtcl=iat;
  }

  inline void
  evaluateForPtclMoveWithHessian(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateAllForPtclMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi);
    Counter++;
    ActivePtcl=iat;
  }

  inline void
  evaluateWithHessian(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi);
    Counter++; // increment a conter
  }

  inline void
  evaluateWithThirdDeriv(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi,grad_grad_grad_Phi);
    Counter++; // increment a conter
  }

  inline void
  evaluateThirdDerivOnly(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateThirdDerivOnly(c,iat,BasisOffset[c],grad_grad_grad_Phi);
    Counter++; // increment a conter
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    aos->setTable(myTable);
    LOBasis[icenter]=aos;
  }
};
}
#endif


