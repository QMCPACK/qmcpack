//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SoaLocalizedBasisSet.h
 * @brief A derived class from BasisSetBase
 *
 * This is intended as a replacement for MolecularOrbitalBase and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H
#define QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H

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
struct SoaLocalizedBasisSet: public BasisSetBase
{
  typedef COT                            ThisCOT_t;
  typedef typename COT::RadialOrbital_t  ThisRadialOrbital_t;
  typedef typename COT::value_type       value_type;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///number of quantum particles
  size_t myTableIndex;
  ///size of the basis set
  size_t BasisSetSize;
  ///Reference to the center
  const ParticleSet::ParticleIndex_t& IonID;

  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  aligned_vector<size_t>  BasisOffset;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  aligned_vector<COT*> LOBasisSet;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els): IonID(ions.GroupID), myTable(0)
  {
    myTableIndex=els.addTable(ions,DT_SOA);
    NumCenters=ions.getTotalNum();
    NumTargets=els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum(),0);
    BasisOffset.resize(NumCenters+1);
    BasisSetSize=-1;
  }

  /** copy constructor */
  SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a)=default;

  /** makeClone */
  SoaLocalizedBasisSet<COT>* makeClone() const
  {
    return new LocalizedBasisSet<COT>(*this);
  }


  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(size_t nbs)
  {
    if(nbs == BasisSetSize) return;

    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for(int c=0; c<NumCenters; c++)
      BasisOffset[c+1] = BasisOffset[c]+LOBasisSet[IonID[c]]->getBasisSetSize();
    BasisSetSize = BasisOffset[NumCenters];
  }

  void resetParameters(const opt_variables_type& active)
  {
    //reset each unique basis functions
    for(int i=0; i<LOBasisSet.size(); i++)
      LOBasisSet[i]->resetParameters(active);
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  { }

  inline void
  evaluateWithHessian(const ParticleSet& P, int iat) { }

  inline void
  evaluateWithThirdDeriv(const ParticleSet& P, int iat) { }

  inline void
  evaluateThirdDerivOnly(const ParticleSet& P, int iat) { }

  /** compute VGL 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, Matrix<value_type>& vgl, bool trialMove)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const value_type* restrict  dist = (trialMove)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto&  displ= (trialMove)? d_table->Temp_dr: d_table->Displacements[iat];
    for(int c=0; c<NumCenters; c++)
      LOBasisSet[IonID[c]]->evaluateVGL(dist[c],displ[c],BasisOffset[c],vgl);
  }

  /** compute values for the iat-paricle move
   *
   * Always uses Temp_r and Temp_dr
   */
  inline void evaluateV(const ParticleSet& P, int iat, T* restrict vals)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const value_type* restrict  dist=d_table->Temp_r.data();
    const auto&  displ=d_table->Temp_dr;
    for(int c=0; c<NumCenters; c++)
      LOBasisSet[IonID[c]]->evaluateV(dist[c],displ[c],vals+BasisOffset[c]);
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    LOBasisSet[icenter]=aos;
  }
};
}
#endif
