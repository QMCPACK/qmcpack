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

#include <memory>

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
  typedef typename ParticleSet::PosType PosType;

  using BaseType::BasisSetSize;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///ion particle set
  const ParticleSet& ions_;
  ///number of quantum particles
  const int myTableIndex;
  ///Global Coordinate of Supertwist read from HDF5
  PosType SuperTwist;


  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  aligned_vector<size_t> BasisOffset;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  aligned_vector<std::unique_ptr<COT>> LOBasisSet;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els)
      : ions_(ions), myTableIndex(els.addTable(ions, true)), SuperTwist(0.0)
  {
    NumCenters = ions.getTotalNum();
    NumTargets = els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum());
    BasisOffset.resize(NumCenters + 1);
    BasisSetSize = 0;
  }

  /** copy constructor */
  SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a)
      : SoaBasisSetBase<ORBT>(a),
        NumCenters(a.NumCenters),
        NumTargets(a.NumTargets),
        ions_(a.ions_),
        myTableIndex(a.myTableIndex),
        SuperTwist(a.SuperTwist),
        BasisOffset(a.BasisOffset)
  {
    LOBasisSet.reserve(a.LOBasisSet.size());
    for (auto& elem : a.LOBasisSet)
      LOBasisSet.push_back(std::make_unique<COT>(*elem));
  }

  /** makeClone */
  //SoaLocalizedBasisSet<COT>* makeClone() const
  BaseType* makeClone() const
  {
    return new SoaLocalizedBasisSet<COT, ORBT>(*this);
  }
  /** set Number of periodic Images to evaluate the orbitals. 
      Set to 0 for non-PBC, and set manually in the input.
      Passes the pre-computed phase factor for evaluation of complex wavefunction. If WF is real Phase_factor is real and equals 1 if gamma or -1 if non-Gamma.  
  */
  void setPBCParams(const TinyVector<int, 3>& PBCImages,
                    const TinyVector<double, 3> Sup_Twist,
                    const std::vector<QMCTraits::ValueType>& phase_factor)
  {
    for (int i = 0; i < LOBasisSet.size(); ++i)
      LOBasisSet[i]->setPBCParams(PBCImages, Sup_Twist, phase_factor);

    SuperTwist = Sup_Twist;
  }
  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    const auto& IonID(ions_.GroupID);
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
    const auto& IonID(ions_.GroupID);
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

  /** compute VGL 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, vgl_type& vgl)
  {
    const auto& IonID(ions_.GroupID);
    const auto& coordR  = P.activeR(iat);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++)
    {
      Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
      Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
      Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
      LOBasisSet[IonID[c]]->evaluateVGL(P.Lattice, dist[c], displ[c], BasisOffset[c], vgl, Tv);
    }
  }


  /** compute VGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(10,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGH(const ParticleSet& P, int iat, vgh_type& vgh)
  {
    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGH(P.Lattice, dist[c], displ[c], BasisOffset[c], vgh);
    }
  }

  /** compute VGHGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vghgh Matrix(20,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGHGH(const ParticleSet& P, int iat, vghgh_type& vghgh)
  {
    // APP_ABORT("SoaLocalizedBasisSet::evaluateVGH() not implemented\n");

    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluateVGHGH(P.Lattice, dist[c], displ[c], BasisOffset[c], vghgh);
    }
  }

  /** compute values for the iat-paricle move
   *
   * Always uses getTempDists() and getTempDispls()
   * Tv is a translation vector; In PBC, in order to reduce the number
   * of images that need to be summed over when generating the AO the 
   * nearest image displacement, dr, is used. Tv corresponds to the 
   * translation that takes the 'general displacement' (displacement
   * between ion position and electron position) to the nearest image 
   * displacement. We need to keep track of Tv because it must be add
   * as a phase factor, i.e., exp(i*k*Tv).
   */
  inline void evaluateV(const ParticleSet& P, int iat, ORBT* restrict vals)
  {
    const auto& IonID(ions_.GroupID);
    const auto& coordR  = P.activeR(iat);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++)
    {
      Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
      Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
      Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
      LOBasisSet[IonID[c]]->evaluateV(P.Lattice, dist[c], displ[c], vals + BasisOffset[c], Tv);
    }
  }

  inline void evaluateGradSourceV(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vgl_type& vgl)
  {
    //We need to zero out the temporary array vgl.
    auto* restrict gx = vgl.data(1);
    auto* restrict gy = vgl.data(2);
    auto* restrict gz = vgl.data(3);

    for (int ib = 0; ib < BasisSetSize; ib++)
    {
      gx[ib] = 0;
      gy[ib] = 0;
      gz[ib] = 0;
    }

    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);


    PosType Tv;
    Tv[0] = Tv[1] = Tv[2] = 0;
    //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
    //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
    //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.
    LOBasisSet[IonID[jion]]->evaluateVGL(P.Lattice, dist[jion], displ[jion], BasisOffset[jion], vgl, Tv);
  }

  inline void evaluateGradSourceVGL(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vghgh_type& vghgh)
  {
    //We need to zero out the temporary array vghgh.
    auto* restrict gx = vghgh.data(1);
    auto* restrict gy = vghgh.data(2);
    auto* restrict gz = vghgh.data(3);

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


    for (int ib = 0; ib < BasisSetSize; ib++)
    {
      gx[ib] = 0;
      gy[ib] = 0;
      gz[ib] = 0;

      hxx[ib] = 0;
      hxy[ib] = 0;
      hxz[ib] = 0;
      hyy[ib] = 0;
      hyz[ib] = 0;
      hzz[ib] = 0;

      gxxx[ib] = 0;
      gxxy[ib] = 0;
      gxxz[ib] = 0;
      gxyy[ib] = 0;
      gxyz[ib] = 0;
      gxzz[ib] = 0;
      gyyy[ib] = 0;
      gyyz[ib] = 0;
      gyzz[ib] = 0;
      gzzz[ib] = 0;
    }

    // Since jion is indexed on the source ions not the ions_ the distinction between
    // ions_ and ions is extremely important.
    const auto& IonID(ions.GroupID);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
    //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
    //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.

    LOBasisSet[IonID[jion]]->evaluateVGHGH(P.Lattice, dist[jion], displ[jion], BasisOffset[jion], vghgh);
  }
  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, std::unique_ptr<COT> aos) { LOBasisSet[icenter] = std::move(aos); }
};
} // namespace qmcplusplus
#endif
