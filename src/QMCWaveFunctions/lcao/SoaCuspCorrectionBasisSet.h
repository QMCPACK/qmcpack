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

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include <QMCWaveFunctions/lcao/MultiQuinticSpline1D.h>

namespace qmcplusplus
{

  /** Handles a set of correction orbitals per atom
   *
   * Reduction over the orbitals - beware of the reduction problem
   */
  template<typename T>
  class CuspCorrectionAtomicBasis
  {
    typedef MultiQuinticSpline1D<T> RadialSetType;
    typedef ParticleSet::PosType PosType;


  public:
    RadialSetType AOs;
    aligned_vector<size_t> ID;

    CuspCorrectionAtomicBasis() {};

    /** copy constructor */
    CuspCorrectionAtomicBasis(const CuspCorrectionAtomicBasis& a)=default;

    inline void evaluate(const T r, T* restrict vals)
    {
      size_t nr=AOs.num_splines;
      T phi[nr];
      AOs.evaluate(r,phi);
      for(size_t i=0; i<nr; ++i)
        vals[ID[i]]+=phi[i];
    }

    inline void evaluate_vgl(const T r, const PosType& dr, T* restrict u, T* restrict du_x,
                             T* restrict du_y, T* restrict du_z, T* restrict d2u)
    {
      size_t nr=AOs.num_splines;
      T phi[nr];
      T dphi[nr];
      T d2phi[nr];
      AOs.evaluate(r,phi,dphi,d2phi);
      for(size_t i=0; i<nr; ++i)
      {
        const size_t j=ID[i];
        u[j]  += phi[i];
        du_x[j] -= dphi[i]*dr[0]/r; // Displacements have opposite sign (relative to AOS)
        du_y[j] -= dphi[i]*dr[1]/r;
        du_z[j] -= dphi[i]*dr[2]/r;
        d2u[j]+= d2phi[i] + 2*dphi[i]/r;
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
  int myTableIndex;
  ///size of the basis set
  int BasisSetSize;

  ///COMPLEX WON'T WORK
  typedef CuspCorrectionAtomicBasis<RealType> COT;

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
  SoaCuspCorrection(ParticleSet& ions, ParticleSet& els)
  {
    myTableIndex=els.addTable(ions,DT_SOA);
    NumCenters=ions.getTotalNum();
    NumTargets=els.getTotalNum();
    LOBasisSet.resize(NumCenters, 0);
  }

  /** copy constructor */
  SoaCuspCorrection(const SoaCuspCorrection& a)=default;


  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    BasisSetSize=nbs;
    //THIS NEEDS TO BE FIXE for OpenMP
    myVGL.resize(5,BasisSetSize);
  }

  /** compute VGL
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
  {

    constexpr RealType czero(0);

    myVGL=czero;

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto displ = (P.activePtcl==iat)? d_table->Temp_dr: d_table->Displacements[iat];
    for(int c=0; c<NumCenters; c++)
    {
      if (LOBasisSet[c]) {
        LOBasisSet[c]->evaluate_vgl(dist[c],displ[c],myVGL[0],myVGL[1],myVGL[2],myVGL[3],myVGL[4]);
      }
    }

    {
      const auto v_in=myVGL[0];
      const auto gx_in=myVGL[1];
      const auto gy_in=myVGL[2];
      const auto gz_in=myVGL[3];
      const auto l_in=myVGL[4];
      auto v_out =vgl.data(0);
      auto gx_out=vgl.data(1);
      auto gy_out=vgl.data(2);
      auto gz_out=vgl.data(3);
      auto l_out =vgl.data(4);
      for(size_t i=0; i<BasisSetSize; ++i)
      {
        v_out[i] +=v_in[i];
        gx_out[i]+=gx_in[i];
        gy_out[i]+=gy_in[i];
        gz_out[i]+=gz_in[i];
        l_out[i] +=l_in[i];
      }
    }
  }

  inline void evaluate_vgl(const ParticleSet& P, int iat, ValueVector_t &psi, GradVector_t &dpsi, ValueVector_t &d2psi)
  {

    constexpr RealType czero(0);

    myVGL=czero;

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto displ = (P.activePtcl==iat)? d_table->Temp_dr: d_table->Displacements[iat];
    for(int c=0; c<NumCenters; c++)
    {
      if (LOBasisSet[c]) {
        LOBasisSet[c]->evaluate_vgl(dist[c],displ[c],myVGL[0],myVGL[1],myVGL[2],myVGL[3],myVGL[4]);
      }
    }

    const auto v_in=myVGL[0];
    const auto gx_in=myVGL[1];
    const auto gy_in=myVGL[2];
    const auto gz_in=myVGL[3];
    const auto l_in=myVGL[4];
    for(size_t i=0; i<BasisSetSize; ++i)
    {
      psi[i] += v_in[i];
      dpsi[i][0] += gx_in[i];
      dpsi[i][1] += gy_in[i];
      dpsi[i][2] += gz_in[i];
      d2psi[i] += l_in[i];
    }
  }

  inline void evaluate_vgl(const ParticleSet& P, int iat, int idx, ValueMatrix_t &psi, GradMatrix_t &dpsi, ValueMatrix_t &d2psi)
  {

    constexpr RealType czero(0);

    myVGL=czero;

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    const auto displ = (P.activePtcl==iat)? d_table->Temp_dr: d_table->Displacements[iat];
    for(int c=0; c<NumCenters; c++)
    {
      if (LOBasisSet[c]) {
        LOBasisSet[c]->evaluate_vgl(dist[c],displ[c],myVGL[0],myVGL[1],myVGL[2],myVGL[3],myVGL[4]);
      }
    }

    const auto v_in=myVGL[0];
    const auto gx_in=myVGL[1];
    const auto gy_in=myVGL[2];
    const auto gz_in=myVGL[3];
    const auto l_in=myVGL[4];
    for(size_t i=0; i<BasisSetSize; ++i)
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
   * Always uses Temp_r and Temp_dr
   */
  inline void evaluateV(const ParticleSet& P, int iat, ValueType* restrict vals)
  {
    constexpr RealType czero(0);
    ValueType* tmp_vals=myVGL[0];

    std::fill_n(tmp_vals, myVGL.size(), czero);

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];

    //THIS IS SERIAL, only way to avoid this is to use myVGL
    for(int c=0; c<NumCenters; c++)
    {
      if (LOBasisSet[c]) {
        LOBasisSet[c]->evaluate(dist[c],tmp_vals);
      }
    }

    {//collect
      const auto v_in=myVGL[0];
      for(size_t i=0; i<BasisSetSize; ++i)
      {
        vals[i] +=v_in[i];
      }
    }

  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    LOBasisSet[icenter]=aos;
  }
  void addVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
  {
    evaluateVGL(P, iat, vgl);
  }
  void addV(const ParticleSet& P, int iat, ValueType* restrict vals)
  {
    evaluateV(P, iat, vals);
  }
  void add_vgl(const ParticleSet& P, int iat, int idx, ValueMatrix_t &vals, GradMatrix_t &dpsi, ValueMatrix_t &d2psi)
  {
    evaluate_vgl(P, iat, idx, vals, dpsi, d2psi);
  }
  void add_vector_vgl(const ParticleSet& P, int iat, ValueVector_t &vals, GradVector_t &dpsi, ValueVector_t &d2psi)
  {
    evaluate_vgl(P, iat, vals, dpsi, d2psi);
  }
};
}
#endif
