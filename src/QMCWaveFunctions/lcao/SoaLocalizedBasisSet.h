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
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

/** A localized basis set derived from RealBasisSetBase<typename COT::value_type>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
template<class COT>
struct SoaLocalizedBasisSet: public RealBasisSetBase<typename COT::value_type>
{
  typedef typename COT::value_type value_type;
  typedef typename RealBasisSetBase<value_type>::vgl_type vgl_type;

  using RealBasisSetBase<value_type>::BasisSetSize;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///number of quantum particles
  int myTableIndex;
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
  SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els): IonID(ions.GroupID)
  {
    myTableIndex=els.addTable(ions,DT_SOA);
    NumCenters=ions.getTotalNum();
    NumTargets=els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum(),0);
    BasisOffset.resize(NumCenters+1);
    BasisSetSize=0;
  }

  /** copy constructor */
  SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a)=default;

  /** makeClone */
  //SoaLocalizedBasisSet<COT>* makeClone() const
  RealBasisSetBase<value_type>* makeClone() const
  {
    SoaLocalizedBasisSet<COT>* myclone=new SoaLocalizedBasisSet<COT>(*this);
    for(int i=0; i<LOBasisSet.size(); ++i)
      myclone->LOBasisSet[i]=LOBasisSet[i]->makeClone();
    return myclone;
  }


  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    if(BasisSetSize>0 && nbs == BasisSetSize) return;

    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for(int c=0; c<NumCenters; c++)
    {
      BasisOffset[c+1] = BasisOffset[c]+LOBasisSet[IonID[c]]->getBasisSetSize();
    }
    BasisSetSize = BasisOffset[NumCenters];
  }

#if 0
  inline int getBasisSetSize()
  {
    return BasisSetSize;
  }

  void resetParameters(const opt_variables_type& active)
  {
    //reset each unique basis functions
    for(int i=0; i<LOBasisSet.size(); i++)
      LOBasisSet[i]->resetParameters(active);
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P) { }
#endif

  /** compute VGL 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, vgl_type& vgl)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const value_type* restrict  dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto& displ= (P.activePtcl==iat)? d_table->Temp_dr: d_table->Displacements[iat];
    for(int c=0; c<NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGL(dist[c],displ[c],BasisOffset[c],vgl);
    }
  }

  /** compute values for the iat-paricle move
   *
   * Always uses Temp_r and Temp_dr
   */
  inline void evaluateV(const ParticleSet& P, int iat, value_type* restrict vals)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const value_type* restrict  dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto& displ= (P.activePtcl==iat)? d_table->Temp_dr: d_table->Displacements[iat];
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
