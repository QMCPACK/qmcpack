//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Eigensolver.h"
#include <vector>
#include <CPU/BLAS.hpp>
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{


// A - Hamiltonian
// B - Overlap
void Eigensolver::solveGeneralizedEigenvalues(Matrix<Real>& A,
                                              Matrix<Real>& B,
                                              std::vector<Real>& eigenvals,
                                              Matrix<Real>& eigenvectors)
{
  int Nl = A.rows();
  assert(A.rows() == A.cols());
  assert(Nl == eigenvals.size());
  assert(Nl == eigenvectors.rows());
  assert(Nl == eigenvectors.cols());

  // transpose the A and B (row-major vs. column-major)
  for (int i = 0; i < Nl; i++)
    for (int j = i + 1; j < Nl; j++)
    {
      std::swap(A(i, j), A(j, i));
      std::swap(B(i, j), B(j, i));
    }

  //   Getting the optimal worksize
  char jl('N');
  char jr('V');
  std::vector<Real> alphar(Nl), alphai(Nl), beta(Nl);
  //Matrix<Real> eigenT(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<Real> work(1);
  Real tt(0);
  int t(1);
  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, alphar.data(), alphai.data(), beta.data(), &tt, &t,
               eigenvectors.data(), &Nl, work.data(), &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, alphar.data(), alphai.data(), beta.data(), &tt, &t,
               eigenvectors.data(), &Nl, work.data(), &lwork, &info);
  if (info != 0)
    throw std::runtime_error("Invalid Matrix Diagonalization Function, ggev info = " + std::to_string(info));

  for (int i = 0; i < Nl; i++)
    eigenvals[i] = alphar[i] / beta[i];
}


// A - Hamiltonian
// B - Overlap
void Eigensolver::solveGeneralizedEigenvalues_Inv(Matrix<Real>& A,
                                                  Matrix<Real>& B,
                                                  std::vector<Real>& eigenvals,
                                                  Matrix<Real>& eigenvectors)
{
  int Nl = A.rows();
  assert(A.rows() == A.cols());
  assert(Nl == eigenvals.size());
  assert(Nl == eigenvectors.rows());
  assert(Nl == eigenvectors.cols());

  invert_matrix(B, false);

  Matrix<Real> prdMat(Nl, Nl);

  MatrixOperators::product(B, A, prdMat);

  // transpose the result (why?)
  for (int i = 0; i < Nl; i++)
    for (int j = i + 1; j < Nl; j++)
      std::swap(prdMat(i, j), prdMat(j, i));

  //   Getting the optimal worksize
  char jl('N');
  char jr('V');
  std::vector<Real> alphai(Nl);
  Matrix<Real> eigenD(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<Real> work(1);
  LAPACK::geev(&jl, &jr, &Nl, prdMat.data(), &Nl, eigenvals.data(), alphai.data(), eigenD.data(), &Nl,
               eigenvectors.data(), &Nl, work.data(), &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::geev(&jl, &jr, &Nl, prdMat.data(), &Nl, eigenvals.data(), alphai.data(), eigenD.data(), &Nl,
               eigenvectors.data(), &Nl, work.data(), &lwork, &info);
  if (info != 0)
    throw std::runtime_error("Invalid Matrix Diagonalization Function, geev info = " + std::to_string(info));
}

} // namespace qmcplusplus
