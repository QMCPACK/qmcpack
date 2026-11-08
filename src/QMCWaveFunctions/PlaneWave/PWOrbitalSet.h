//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWOrbitalSet.h
 * @brief Complex plane-wave SPOs with values, VGL, and Cartesian Hessians.
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
/** Contract complex orbital coefficients with a single primitive-cell PWBasis.
 *
 * Used by complex builds; tiled SPOs combine these components in a CompositeSPOSet.
 * Call resize() and load coefficients before evaluation. Output vectors hold all
 * orbitals, and matrix outputs use (particle-first, orbital) indexing; the legacy
 * logdet parameter names denote orbital values, not logarithms or determinants.
 */
class PWOrbitalSet : public SPOSet
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
  PWOrbitalSet(const std::string& my_name, size_t size)
      : SPOSet(my_name, size), OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0), C(nullptr), IsCloned(false)
  {}

  std::string getClassName() const override { return "PWOrbitalSet"; }


  /// Release an owned basis and, for a non-clone, the coefficient matrix.
  ~PWOrbitalSet() override;

  /** Copy the basis evaluation state while sharing orbital coefficients.
   * The source coefficient storage must outlive the clone.
   */
  std::unique_ptr<SPOSet> makeClone() const override;
  /** Attach an initialized basis and allocate coefficient and derivative buffers.
   * @param bset basis with loaded G vectors and prepared evaluation arrays
   * @param cleanup whether this SPO should delete bset on destruction
   */
  void resize(PWBasisPtr bset, bool cleanup = false);

  /** Load one orbital in input G-vector order, applying the basis inputmap.
   * @pre jorb is a valid orbital index and coefs.size() equals inputmap.size().
   */
  void addVector(const std::vector<ComplexType>& coefs, int jorb);
  /// Load real coefficients using the same input ordering as the complex overload.
  void addVector(const std::vector<RealType>& coefs, int jorb);


  /// Evaluate orbital ib at a Cartesian position, updating the basis value buffer.
  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return BLAS::dot(BasisSetSize, (*C)[ib], myBasisSet->Zv.data());
  }

  /// Evaluate all orbital values at P.activeR(iat).
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  /// Evaluate all orbital values, Cartesian gradients, and Laplacians.
  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  /** Evaluate VGL for particles in [first, last).
   * Output matrices require last-first rows and OrbitalSetSize columns.
   */
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

  /** Evaluate values, Cartesian gradients, and full Hessians for [first, last).
   * Output matrices use the same particle/ orbital layout as the VGL overload.
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
  /// Coefficients indexed by (orbital, G-vector); makeClone() shares this matrix.
  ValueMatrix* C;
  /// Whether C is shared with a source SPO and must not be deleted here.
  bool IsCloned;
  /////Plane-wave coefficients: (iband,g-vector)
  //Matrix<ValueType> Coefs;
  /// Contracted VGL data indexed by (orbital, PW_* column).
  Matrix<ValueType> Temp;
  /// Contracted Hessians indexed by (orbital, row*OHMMS_DIM+column).
  Matrix<ValueType> TempHess;
};
} // namespace qmcplusplus
#endif
