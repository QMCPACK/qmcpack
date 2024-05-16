//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file LinearMethod.h
 * @brief LinearMethod related functions.
 */
#ifndef QMCPLUSPLUS_LINEARMETHOD_H
#define QMCPLUSPLUS_LINEARMETHOD_H

#include <Configuration.h>
#include <NewTimer.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{
///forward declaration of a cost function
class QMCCostFunctionBase;

class LinearMethod
{
  using Real = QMCTraits::RealType;

  // obtain the range of non-linear parameters
  void getNonLinearRange(int& first, int& last, const QMCCostFunctionBase& optTarget) const;

public:
  /** Use generalized eigenvalue solver
   *  @param[in]  A Hamiltonian matrix
   *  @param[in]  B Overlap matrix
   *  @param[out] eigenvals Eigenvalues
   *  @param[out] eigenvectors Eigenvectors corresponding to the eigenvalues
   *
   */
  static void solveGeneralizedEigenvalues(Matrix<Real>& A,
                                          Matrix<Real>& B,
                                          std::vector<Real>& eigenvals,
                                          Matrix<Real>& eigenvectors);

  /** Solve by explicitly inverting the overlap matrix
   *  @param[in]  A Hamiltonian matrix
   *  @param[in]  B Overlap matrix
   *  @param[out] eigenvals Eigenvalues
   *  @param[out] eigenvectors Eigenvectors corresponding to the eigenvalues
   *
   */
  static void solveGeneralizedEigenvalues_Inv(Matrix<Real>& A,
                                              Matrix<Real>& B,
                                              std::vector<Real>& eigenvals,
                                              Matrix<Real>& eigenvectors);

  /** Select eigenvalue and return corresponding scaled eigenvector
   *  @param[in] eigenvals Eigenvalues
   *  @param[in] eigenvectors Eigenvectors
   *  @param[in] zerozero The H(0,0) element, used to guide eigenvalue selection
   *  @param[out] ev The eigenvector scaled by the reciprocal of the first element
   *  @return The selected eigenvalue
   */
  static Real selectEigenvalue(std::vector<Real>& eigenvals,
                               Matrix<Real>& eigenvectors,
                               Real zerozero,
                               std::vector<Real>& ev);

  /** Solve the generalized eigenvalue problem and return a scaled eigenvector corresponding to the selected eigenvalue
   * @param[in] A Hamiltonian matrix
   * @param[in] B Overlap matrix
   * @param[out] ev Scaled eigenvector
   *
   * This uses a regular eigenvalue solver for B^-1 A (LAPACK geev).
   * In theory using the generalized eigenvalue solver is more numerically stable, but
   * in practice this solver is sufficiently stable, and is faster than the generalized solver.
   */
  static Real getLowestEigenvector_Inv(Matrix<Real>& A, Matrix<Real>& B, std::vector<Real>& ev);

  /** Solve the generalized eigenvalue problem and return a scaled eigenvector corresponding to the selected eigenvalue
   * @param[in] A Hamiltonian matrix
   * @param[in] B Overlap matrix
   * @param[out] ev Scaled eigenvector
   *
   * This uses a generalized eigenvalue solver (LAPACK ggev).
   */
  static Real getLowestEigenvector_Gen(Matrix<Real>& A, Matrix<Real>& B, std::vector<Real>& ev);

  //asymmetric EV
  Real getLowestEigenvector(Matrix<Real>& A, std::vector<Real>& ev) const;
  // compute a rescale factor. Ye: Where is the method from?
  Real getNonLinearRescale(std::vector<Real>& dP, Matrix<Real>& S, const QMCCostFunctionBase& optTarget) const;
};
} // namespace qmcplusplus
#endif
