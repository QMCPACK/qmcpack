//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com
//////////////////////////////////////////////////////////////////////////////////////


/** @file EigenSolver.h
 * @brief EigenSolver related functions.
 */
#ifndef QMCPLUSPLUS_EIGENSOLVER_H
#define QMCPLUSPLUS_EIGENSOLVER_H

#include <Configuration.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{

class Eigensolver
{
  using Real = QMCTraits::RealType;

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
};
} // namespace qmcplusplus
#endif
