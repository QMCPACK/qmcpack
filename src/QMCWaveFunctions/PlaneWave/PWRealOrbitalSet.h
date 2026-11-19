//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWRealOrbitalSet.h
 * @brief Real and imaginary projections of complex plane-wave SPOs for real builds.
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_REALORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_REALORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
/** Real-valued SPOs formed by projecting a complex primitive-cell expansion.
 *
 * Each orbital selects its real or imaginary projection when coefficients are
 * loaded. The same projection is applied to values and every derivative, allowing
 * a retained time-reversal representative to supply two distinct real orbitals.
 * Call resize() and load coefficients before evaluation.
 *
 * Output vectors hold all orbitals. Matrix outputs use (particle-first, orbital)
 * indexing; legacy logdet names denote orbital values, not logarithms or determinants.
 */
class PWRealOrbitalSet : public SPOSet
{
public:
  using BasisSet_t = PWBasis;
  using PWBasisPtr = PWBasis*;

  /** inherit the enum of BasisSet_t */
  enum
  {
    PW_VALUE    = BasisSet_t::PW_VALUE,
    PW_LAP      = BasisSet_t::PW_LAP,
    PW_GRADX    = BasisSet_t::PW_GRADX,
    PW_GRADY    = BasisSet_t::PW_GRADY,
    PW_GRADZ    = BasisSet_t::PW_GRADZ,
    PW_MAXINDEX = BasisSet_t::PW_MAXINDEX
  };

  /// Create an SPO of the requested size without attaching a basis or coefficients.
  PWRealOrbitalSet(const std::string& my_name, size_t size)
      : SPOSet(my_name, size), OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0)
  {}

  std::string getClassName() const override { return "PWRealOrbitalSet"; }

  /// Delete myBasisSet only when OwnBasisSet is true.
  ~PWRealOrbitalSet() override;

  /// Copy the SPO and give the clone its own PWBasis evaluation state.
  std::unique_ptr<SPOSet> makeClone() const override;

  /** Attach an initialized basis and allocate coefficient and derivative buffers.
   * Projection flags are reset to select real parts.
   * @param bset basis with loaded G vectors and prepared evaluation arrays
   * @param cleanup whether this SPO should delete bset on destruction
   */
  void resize(PWBasisPtr bset, bool cleanup = false);

  /** Load real coefficients for an orbital in input G-vector order.
   * @param coefs one coefficient per entry of the basis inputmap
   * @param jorb valid zero-based orbital index
   * @param use_imaginary_part select imaginary rather than real values and derivatives
   */
  void addVector(const std::vector<RealType>& coefs, int jorb, bool use_imaginary_part = false);

  /** Load complex coefficients for an orbital in input G-vector order.
   * @param coefs one coefficient per entry of the basis inputmap
   * @param jorb valid zero-based orbital index
   * @param use_imaginary_part select imaginary rather than real values and derivatives
   */
  void addVector(const std::vector<ComplexType>& coefs, int jorb, bool use_imaginary_part = false);


  /// Evaluate the selected projection of orbital ib at a Cartesian position.
  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return project(BLAS::dot(BasisSetSize, CC[ib], myBasisSet->Zv.data()), ib);
  }

  /// Evaluate all projected orbital values at P.activeR(iat).
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  /// Evaluate projected values, Cartesian gradients, and Laplacians.
  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  /** Evaluate projected VGL for particles in [first, last).
   * Output matrices require last-first rows and OrbitalSetSize columns.
   */
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

  /** Evaluate projected values, gradients, and Cartesian Hessians for [first, last).
   * Hessian elements use the same per-orbital projection as values and gradients.
   */
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override;


  /// Whether the destructor deletes myBasisSet.
  bool OwnBasisSet;
  /// Legacy twist metadata; evaluation uses the twist stored in myBasisSet.
  PosType TwistAngle;
  /// Basis supplying primitive-cell phases and Cartesian derivatives.
  PWBasisPtr myBasisSet;
  /// Number of retained plane waves per orbital.
  IndexType BasisSetSize;
  /// Complex coefficients indexed by (orbital, G-vector), before projection.
  Matrix<ComplexType> CC;
  /// Contracted complex VGL data indexed by (orbital, PW_* column).
  Matrix<ComplexType> Temp;
  /// Contracted complex Hessians indexed by (orbital, row*OHMMS_DIM+column).
  Matrix<ComplexType> TempHess;
  /// Contracted complex values before per-orbital projection.
  Vector<ComplexType> tempPsi;
  ///select the imaginary rather than real projection for each orbital
  std::vector<bool> useImaginaryPart;

private:
  /// Apply the projection chosen when this orbital's coefficients were loaded.
  RealType project(const ComplexType& value, int orbital_index) const
  {
    return useImaginaryPart[orbital_index] ? value.imag() : value.real();
  }
};
} // namespace qmcplusplus
#endif
