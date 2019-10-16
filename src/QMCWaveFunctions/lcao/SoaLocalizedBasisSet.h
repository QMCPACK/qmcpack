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
 * This is intended as a replacement for MolecularWaveFunctionComponent and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H
#define QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
/** A localized basis set derived from SoaBasisSetBase<ORBT>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 * The template parameter ORBT denotes the orbital value return type
 */
template<class COT, typename ORBT>
struct SoaLocalizedBasisSet : public SoaBasisSetBase<ORBT>
{
  typedef typename COT::RealType RealType;
  typedef SoaBasisSetBase<ORBT> BaseType;
  typedef typename BaseType::vgl_type vgl_type;
  typedef typename BaseType::vgh_type vgh_type;
  typedef typename BaseType::vghgh_type vghgh_type;

  using BaseType::BasisSetSize;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///Reference to the center
  const ParticleSet::ParticleIndex_t& IonID;
  ///number of quantum particles
  const int myTableIndex;

  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  aligned_vector<size_t> BasisOffset;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  aligned_vector<COT*> LOBasisSet;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els)
    : IonID(ions.GroupID), myTableIndex(els.addTable(ions, DT_SOA))
  {
    NumCenters   = ions.getTotalNum();
    NumTargets   = els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum(), 0);
    BasisOffset.resize(NumCenters + 1);
    BasisSetSize = 0;
  }

  /** copy constructor */
  SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a) = default;

  /** makeClone */
  //SoaLocalizedBasisSet<COT>* makeClone() const
  BaseType* makeClone() const
  {
    SoaLocalizedBasisSet<COT, ORBT>* myclone = new SoaLocalizedBasisSet<COT, ORBT>(*this);
    for (int i = 0; i < LOBasisSet.size(); ++i)
      myclone->LOBasisSet[i] = LOBasisSet[i]->makeClone();
    return myclone;
  }
  /** set Number of periodic Images to evaluate the orbitals. 
      Set to 0 for non-PBC, and set manually in the input.
      Passes the pre-computed phase factor for evaluation of complex wavefunction. If WF is real Phase_factor is real and equals 1 if gamma or -1 if non-Gamma.  
  */
  void setPBCParams(const TinyVector<int, 3>& PBCImages,const std::vector<QMCTraits::ValueType>& phase_factor)
  {
    for (int i = 0; i < LOBasisSet.size(); ++i)
      LOBasisSet[i]->setPBCParams(PBCImages,phase_factor);
  }
  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    if (BasisSetSize > 0 && nbs == BasisSetSize)
      return;

    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for (int c = 0; c < NumCenters; c++)
    {
      BasisOffset[c + 1] = BasisOffset[c] + LOBasisSet[IonID[c]]->getBasisSetSize();
    }
    BasisSetSize = BasisOffset[NumCenters];
  }

  /**  Determine which orbitals are S-type.  Used by cusp correction.
    */
  void queryOrbitalsForSType(const std::vector<bool>& corrCenter, std::vector<bool>& is_s_orbital) const
  {
    int idx = 0;
    for (int c = 0; c < NumCenters; c++)
    {
      int bss = LOBasisSet[IonID[c]]->BasisSetSize;
      std::vector<bool> local_is_s_orbital(bss);
      LOBasisSet[IonID[c]]->queryOrbitalsForSType(local_is_s_orbital);
      for (int k = 0; k < bss; k++)
      {
        if (corrCenter[c])
        {
          is_s_orbital[idx++] = local_is_s_orbital[k];
        }
        else
        {
          is_s_orbital[idx++] = false;
        }
      }
    }
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
    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGL(P.Lattice, dist[c], displ[c], BasisOffset[c], vgl);
    }
    /*std::vector<double> K {0.333,0.333,0.333};
    RealType s,c;
    RealType vec_scalar;
    vec_scalar=(((P.activePtcl == iat) ? P.activePos : P.R[iat])[0]*K[0]+((P.activePtcl == iat) ? P.activePos : P.R[iat])[1]*K[1]+((P.activePtcl == iat) ? P.activePos : P.R[iat])[2]*K[2]);  
    sincos(-2*M_PI*vec_scalar, &s,&c);
    QMCTraits::ValueType PhaseFactor(c,s);
    for (int i =0; i<BasisSetSize;i++)
    {
      vgl.data(0)[i]*=PhaseFactor;
      vgl.data(1)[i]*=PhaseFactor;
      vgl.data(2)[i]*=PhaseFactor;
      vgl.data(3)[i]*=PhaseFactor;
      vgl.data(4)[i]*=PhaseFactor;
    }*/
  }


  /** compute VGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(10,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGH(const ParticleSet& P, int iat, vgh_type& vgh)
{
    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGH(P.Lattice, dist[c], displ[c], BasisOffset[c], vgh);
    }
    
  }

  /** compute VGHGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vghgh Matrix(20,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGHGH(const ParticleSet& P, int iat, vghgh_type& vghgh)
  {
   // APP_ABORT("SoaLocalizedBasisSet::evaluateVGH() not implemented\n");
    
    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGHGH(P.Lattice, dist[c], displ[c], BasisOffset[c], vghgh);
    }
    
  }

  /** compute values for the iat-paricle move
   *
   * Always uses Temp_r and Temp_dr
   */
  inline void evaluateV(const ParticleSet& P, int iat, ORBT* restrict vals)
  {
    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateV(P.Lattice, dist[c], displ[c], vals + BasisOffset[c]);
    }
    /*std::vector<double> K {0.333,0.333,0.333};
    RealType s,c;
    RealType vec_scalar;
    vec_scalar=(((P.activePtcl == iat) ? P.activePos : P.R[iat])[0]*K[0]+((P.activePtcl == iat) ? P.activePos : P.R[iat])[1]*K[1]+((P.activePtcl == iat) ? P.activePos : P.R[iat])[2]*K[2]);  
    sincos(-2*M_PI*vec_scalar, &s,&c);
    QMCTraits::ValueType PhaseFactor(c,s);
    for (int i =0; i<BasisSetSize;i++)
      vals[i]*=PhaseFactor;
    */
  }
  inline void evaluateGradSourceV(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vgl_type& vgl)
  {
    //We need to zero out the temporary array vgl.  
    auto* restrict gx  = vgl.data(1);
    auto* restrict gy  = vgl.data(2);
    auto* restrict gz  = vgl.data(3);

    for(int ib=0; ib<BasisSetSize; ib++)
    {
      gx[ib]=0;
      gy[ib]=0;
      gz[ib]=0;
    }

    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
  
    //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
    //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
    //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.
    LOBasisSet[IonID[jion]]->evaluateVGL(P.Lattice, dist[jion], displ[jion], BasisOffset[jion], vgl);

  }

  inline void evaluateGradSourceVGL(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vghgh_type& vghgh)
  {
    //We need to zero out the temporary array vghgh.
    auto* restrict gx  = vghgh.data(1);
    auto* restrict gy  = vghgh.data(2);
    auto* restrict gz  = vghgh.data(3);

    auto* restrict hxx = vghgh.data(4);
    auto* restrict hxy = vghgh.data(5);
    auto* restrict hxz = vghgh.data(6);
    auto* restrict hyy = vghgh.data(7);
    auto* restrict hyz = vghgh.data(8);
    auto* restrict hzz = vghgh.data(9);
    
    auto* restrict gxxx = vghgh.data(10);
    auto* restrict gxxy = vghgh.data(11);
    auto* restrict gxxz = vghgh.data(12);
    auto* restrict gxyy = vghgh.data(13);
    auto* restrict gxyz = vghgh.data(14);
    auto* restrict gxzz = vghgh.data(15);
    auto* restrict gyyy = vghgh.data(16);
    auto* restrict gyyz = vghgh.data(17);
    auto* restrict gyzz = vghgh.data(18);
    auto* restrict gzzz = vghgh.data(19);


    for(int ib=0; ib<BasisSetSize; ib++)
    {
      gx[ib]=0;
      gy[ib]=0;
      gz[ib]=0;

      hxx[ib]=0;
      hxy[ib]=0;
      hxz[ib]=0;
      hyy[ib]=0;
      hyz[ib]=0;
      hzz[ib]=0;

      gxxx[ib]=0;
      gxxy[ib]=0;
      gxxz[ib]=0;
      gxyy[ib]=0;
      gxyz[ib]=0;
      gxzz[ib]=0;
      gyyy[ib]=0;
      gyyz[ib]=0;
      gyzz[ib]=0;
      gzzz[ib]=0;
    }

    const auto& d_table = P.getDistTable(myTableIndex);
    const RealType* restrict dist    = (P.activePtcl == iat) ? d_table.Temp_r.data() : d_table.Distances[iat];
    const auto& displ                = (P.activePtcl == iat) ? d_table.Temp_dr : d_table.Displacements[iat];
    
    //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
    //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
    //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.
    LOBasisSet[IonID[jion]]->evaluateVGHGH(P.Lattice, dist[jion], displ[jion], BasisOffset[jion], vghgh);
    
  }
  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos) { LOBasisSet[icenter] = aos; }
};
} // namespace qmcplusplus
#endif
