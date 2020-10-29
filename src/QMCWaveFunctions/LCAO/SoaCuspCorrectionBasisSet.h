//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file SoaCuspCorrectionBasisSet.h
 *
 * Convert CorrectingOrbitalBasisSet using MultiQuinticSpline1D<T>
 */
#ifndef QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H
#define QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/LCAO/MultiQuinticSpline1D.h"

namespace qmcplusplus
{
// /** Handles a set of correction orbitals per atom
//  *
//  * Reduction over the orbitals - beware of the reduction problem
//  */
template<typename T>
class CuspCorrectionAtomicBasis
{
  using QMCT = QMCTraits;
  typedef MultiQuinticSpline1D<T> RadialSetType;
  typedef ParticleSet::PosType PosType;

  QMCT::RealType r_max_ = 100;
  RadialSetType AOs;
  aligned_vector<size_t> ID;

public:
  CuspCorrectionAtomicBasis(){};

  /** copy constructor */
  CuspCorrectionAtomicBasis(const CuspCorrectionAtomicBasis& a) = default;

  inline void initializeRadialSet(LogGrid<T>& radial_grid, QMCT::IndexType orbital_set_size)
  {
    r_max_ = radial_grid.rmax();
    AOs.initialize(radial_grid, orbital_set_size);
  }

  template<class T1>
  inline void addSpline(int mo_idx, OneDimQuinticSpline<T1>& radial_spline)
  {
    AOs.add_spline(mo_idx, radial_spline);
  }

  inline void evaluate(const T r, T* restrict vals)
  {
    //assume output vars are zero'd
    if (r >= r_max_)
      return;

    size_t nr = AOs.getNumSplines();
    std::vector<T> phi(nr);

    AOs.evaluate(r, phi.data());
    for (size_t i = 0; i < nr; ++i)
      //vals[ID[i]]+=phi[i];
      vals[i] += phi[i];
  }

  inline void evaluate_vgl(const T r,
                           const PosType& dr,
                           T* restrict u,
                           T* restrict du_x,
                           T* restrict du_y,
                           T* restrict du_z,
                           T* restrict d2u)
  {
    //assume output vars are zero'd
    if (r >= r_max_)
      return;

    size_t nr = AOs.getNumSplines();
    std::vector<T> phi(nr);
    std::vector<T> dphi(nr);
    std::vector<T> d2phi(nr);

    AOs.evaluate(r, phi.data(), dphi.data(), d2phi.data());
    for (size_t i = 0; i < nr; ++i)
    {
      const size_t j = i; //ID[i];
      u[j] += phi[i];
      du_x[j] -= dphi[i] * dr[0] / r; // Displacements have opposite sign (relative to AOS)
      du_y[j] -= dphi[i] * dr[1] / r;
      du_z[j] -= dphi[i] * dr[2] / r;
      d2u[j] += d2phi[i] + 2 * dphi[i] / r;
    }
  }
};

/** A localized basis set derived from BasisSetBase<typename COT::ValueType>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
struct SoaCuspCorrection
{
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::RealType RealType;
  typedef VectorSoaContainer<ValueType, 5> VGLVector_t;
  typedef SPOSet::ValueMatrix_t ValueMatrix_t;
  typedef SPOSet::GradMatrix_t GradMatrix_t;
  typedef SPOSet::GradVector_t GradVector_t;
  typedef SPOSet::ValueVector_t ValueVector_t;
  typedef ParticleSet::PosType PosType;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///number of quantum particles
  const int myTableIndex;
  ///size of the basis set
  int BasisSetSize;

  ///COMPLEX WON'T WORK
  typedef CuspCorrectionAtomicBasis<RealType> COT;

  int unused = 1;
  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number of centers (atoms)
   */
  aligned_vector<COT*> LOBasisSet;

  Matrix<RealType> myVGL;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaCuspCorrection(ParticleSet& ions, ParticleSet& els) : myTableIndex(els.addTable(ions))
  {
    NumCenters = ions.getTotalNum();
    NumTargets = els.getTotalNum();
    LOBasisSet.resize(NumCenters, 0);
  }

  /** copy constructor */
  SoaCuspCorrection(const SoaCuspCorrection& a) = default;


  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    BasisSetSize = nbs;
    //THIS NEEDS TO BE FIXE for OpenMP
    myVGL.resize(5, BasisSetSize);
  }

  /** compute VGL
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
  {
    myVGL = 0.0;

    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++)
    {
      if (LOBasisSet[c])
      {
        LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
      }
    }

    {
      const auto v_in  = myVGL[0];
      const auto gx_in = myVGL[1];
      const auto gy_in = myVGL[2];
      const auto gz_in = myVGL[3];
      const auto l_in  = myVGL[4];
      auto v_out       = vgl.data(0);
      auto gx_out      = vgl.data(1);
      auto gy_out      = vgl.data(2);
      auto gz_out      = vgl.data(3);
      auto l_out       = vgl.data(4);
      for (size_t i = 0; i < BasisSetSize; ++i)
      {
        v_out[i] += v_in[i];
        gx_out[i] += gx_in[i];
        gy_out[i] += gy_in[i];
        gz_out[i] += gz_in[i];
        l_out[i] += l_in[i];
      }
    }
  }

  inline void evaluate_vgl(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    myVGL = 0.0;

    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++)
    {
      if (LOBasisSet[c])
      {
        LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
      }
    }

    const auto v_in  = myVGL[0];
    const auto gx_in = myVGL[1];
    const auto gy_in = myVGL[2];
    const auto gz_in = myVGL[3];
    const auto l_in  = myVGL[4];
    for (size_t i = 0; i < BasisSetSize; ++i)
    {
      psi[i] += v_in[i];
      dpsi[i][0] += gx_in[i];
      dpsi[i][1] += gy_in[i];
      dpsi[i][2] += gz_in[i];
      d2psi[i] += l_in[i];
    }
  }

  inline void evaluate_vgl(const ParticleSet& P,
                           int iat,
                           int idx,
                           ValueMatrix_t& psi,
                           GradMatrix_t& dpsi,
                           ValueMatrix_t& d2psi)
  {
    myVGL = 0.0;

    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++)
    {
      if (LOBasisSet[c])
      {
        LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
      }
    }

    const auto v_in  = myVGL[0];
    const auto gx_in = myVGL[1];
    const auto gy_in = myVGL[2];
    const auto gz_in = myVGL[3];
    const auto l_in  = myVGL[4];
    for (size_t i = 0; i < BasisSetSize; ++i)
    {
      psi[idx][i] += v_in[i];
      dpsi[idx][i][0] += gx_in[i];
      dpsi[idx][i][1] += gy_in[i];
      dpsi[idx][i][2] += gz_in[i];
      d2psi[idx][i] += l_in[i];
    }
  }

  /** compute values for the iat-paricle move
   *
   * Always uses getTempDists() and getTempDispls()
   */
  inline void evaluateV(const ParticleSet& P, int iat, ValueType* restrict vals)
  {
    ValueType* tmp_vals = myVGL[0];

    std::fill_n(tmp_vals, myVGL.size(), 0.0);

    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);

    //THIS IS SERIAL, only way to avoid this is to use myVGL
    for (int c = 0; c < NumCenters; c++)
    {
      if (LOBasisSet[c])
      {
        LOBasisSet[c]->evaluate(dist[c], tmp_vals);
      }
    }

    { //collect
      const auto v_in = myVGL[0];
      for (size_t i = 0; i < BasisSetSize; ++i)
      {
        vals[i] += v_in[i];
      }
    }
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos) { LOBasisSet[icenter] = aos; }
  void addVGL(const ParticleSet& P, int iat, VGLVector_t& vgl) { evaluateVGL(P, iat, vgl); }
  void addV(const ParticleSet& P, int iat, ValueType* restrict vals) { evaluateV(P, iat, vals); }
  void add_vgl(const ParticleSet& P, int iat, int idx, ValueMatrix_t& vals, GradMatrix_t& dpsi, ValueMatrix_t& d2psi)
  {
    evaluate_vgl(P, iat, idx, vals, dpsi, d2psi);
  }
  void add_vector_vgl(const ParticleSet& P, int iat, ValueVector_t& vals, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    evaluate_vgl(P, iat, vals, dpsi, d2psi);
  }
};
} // namespace qmcplusplus
#endif
