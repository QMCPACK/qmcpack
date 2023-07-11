//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file SoaCuspCorrection.h
 */
#ifndef QMCPLUSPLUS_SOA_CUSPCORRECTION_H
#define QMCPLUSPLUS_SOA_CUSPCORRECTION_H

#include "Configuration.h"
#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
template<typename T>
class CuspCorrectionAtomicBasis;

/** A localized basis set derived from BasisSetBase<typename COT::ValueType>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
class SoaCuspCorrection
{
  using ValueType   = QMCTraits::ValueType;
  using RealType    = QMCTraits::RealType;
  using VGLVector   = VectorSoaContainer<ValueType, 5>;
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradMatrix  = SPOSet::GradMatrix;
  using GradVector  = SPOSet::GradVector;
  using ValueVector = SPOSet::ValueVector;
  using PosType     = ParticleSet::PosType;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///number of quantum particles
  const int myTableIndex;
  /** Maximal number of supported MOs
   * this is not the AO basis because cusp correction is applied on the MO directly.
   */
  int MaxOrbSize = 0;

  ///COMPLEX WON'T WORK
  using COT = CuspCorrectionAtomicBasis<RealType>;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number of centers (atoms)
   * should use unique_ptr once COT is fixed for better performance
   */
  std::vector<std::shared_ptr<const COT>> LOBasisSet;

  Matrix<RealType> myVGL;

public:
  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaCuspCorrection(ParticleSet& ions, ParticleSet& els);

  /** copy constructor */
  SoaCuspCorrection(const SoaCuspCorrection& a);

  /** set the number of orbitals this cusp correction may serve. call this before adding any correction centers.
   */
  void setOrbitalSetSize(int norbs);

  /** compute VGL
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  void evaluateVGL(const ParticleSet& P, int iat, VGLVector& vgl);

  void evaluate_vgl(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  void evaluate_vgl(const ParticleSet& P, int iat, int idx, ValueMatrix& psi, GradMatrix& dpsi, ValueMatrix& d2psi);

  /** compute values for the iat-paricle move
   *
   * Always uses getTempDists() and getTempDispls()
   */
  void evaluateV(const ParticleSet& P, int iat, ValueVector& psi);

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, std::unique_ptr<COT> aos);

  void addVGL(const ParticleSet& P, int iat, VGLVector& vgl) { evaluateVGL(P, iat, vgl); }
  void addV(const ParticleSet& P, int iat, ValueVector& psi) { evaluateV(P, iat, psi); }
  void add_vgl(const ParticleSet& P, int iat, int idx, ValueMatrix& vals, GradMatrix& dpsi, ValueMatrix& d2psi)
  {
    evaluate_vgl(P, iat, idx, vals, dpsi, d2psi);
  }
  void add_vector_vgl(const ParticleSet& P, int iat, ValueVector& vals, GradVector& dpsi, ValueVector& d2psi)
  {
    evaluate_vgl(P, iat, vals, dpsi, d2psi);
  }
};
} // namespace qmcplusplus
#endif
