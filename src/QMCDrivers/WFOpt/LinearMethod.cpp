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


#include "LinearMethod.h"
#include <vector>
#include "QMCCostFunctionBase.h"
#include <CPU/BLAS.hpp>
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{

LinearMethod::Real LinearMethod::getLowestEigenvector_Inv(Matrix<Real>& A, Matrix<Real>& B, std::vector<Real>& ev)
{
  int Nl(ev.size());
  std::vector<Real> alphar(Nl);
  Matrix<Real> eigenT(Nl, Nl);
  Real zerozero = A(0, 0);
  solveGeneralizedEigenvalues_Inv(A, B, alphar, eigenT);
  return selectEigenvalue(alphar, eigenT, zerozero, ev);
}

LinearMethod::Real LinearMethod::getLowestEigenvector_Gen(Matrix<Real>& A, Matrix<Real>& B, std::vector<Real>& ev)
{
  int Nl(ev.size());
  std::vector<Real> alphar(Nl);
  Matrix<Real> eigenT(Nl, Nl);
  Real zerozero = A(0, 0);
  solveGeneralizedEigenvalues(A, B, alphar, eigenT);
  return selectEigenvalue(alphar, eigenT, zerozero, ev);
}

// Input
//  -eigenvals
//  -eigenvectors
//  Returns selected eigenvalue
//   - ev - scaled eigenvector corresponding to selected eigenvalue

LinearMethod::Real LinearMethod::selectEigenvalue(std::vector<Real>& eigenvals,
                                                  Matrix<Real>& eigenvectors,
                                                  Real zerozero,
                                                  std::vector<Real>& ev)
{
  // Filter and sort to find desired eigenvalue.
  // Filter accepts eigenvalues between E_0 and E_0 - 100.0,
  // where E_0 is H(0,0), the current estimate for the VMC energy.
  // Sort searches for eigenvalue closest to E_0 - 2.0
  int Nl = eigenvals.size();

  bool found_any_eigenvalue = false;
  std::vector<std::pair<Real, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    Real evi(eigenvals[i]);
    if ((evi < zerozero) && (evi > (zerozero - 1e2)))
    {
      mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
      mappedEigenvalues[i].second = i;
      found_any_eigenvalue        = true;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<Real>::max();
      mappedEigenvalues[i].second = i;
    }
  }

  // Sometimes there is no eigenvalue less than E_0, but there is one slightly higher.
  // Run filter and sort again, except this time accept eigenvalues between E_0 + 100.0 and E_0 - 100.0.
  // Since it's already been determined there are no eigenvalues below E_0, the sort
  // finds the eigenvalue closest to E_0.
  if (!found_any_eigenvalue)
  {
    app_log() << "No eigenvalues passed initial filter. Trying broader 100 a.u. filter" << std::endl;

    bool found_higher_eigenvalue;
    for (int i = 0; i < Nl; i++)
    {
      Real evi(eigenvals[i]);
      if ((evi < zerozero + 1e2) && (evi > (zerozero - 1e2)))
      {
        mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
        mappedEigenvalues[i].second = i;
        found_higher_eigenvalue     = true;
      }
      else
      {
        mappedEigenvalues[i].first  = std::numeric_limits<Real>::max();
        mappedEigenvalues[i].second = i;
      }
    }
    if (!found_higher_eigenvalue)
    {
      app_log() << "No eigenvalues passed second filter. Optimization is likely to fail." << std::endl;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  //         for (int i=0; i<4; i++) app_log()<<i<<": "<<alphar[mappedEigenvalues[i].second]<< std::endl;
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenvectors(mappedEigenvalues[0].second, i) / eigenvectors(mappedEigenvalues[0].second, 0);
  return eigenvals[mappedEigenvalues[0].second];
}

// A - Hamiltonian
// B - Overlap
void LinearMethod::solveGeneralizedEigenvalues(Matrix<Real>& A,
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
  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t,
               eigenvectors.data(), &Nl, &work[0], &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t,
               eigenvectors.data(), &Nl, &work[0], &lwork, &info);
  if (info != 0)
    throw std::runtime_error("Invalid Matrix Diagonalization Function, ggev info = " + std::to_string(info));

  for (int i = 0; i < Nl; i++)
  {
    eigenvals[i] = alphar[i] / beta[i];
  }
}


// A - Hamiltonian
// B - Overlap
void LinearMethod::solveGeneralizedEigenvalues_Inv(Matrix<Real>& A,
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
  std::vector<Real> alphar(Nl), alphai(Nl);
  //Matrix<Real> eigenT(Nl, Nl);
  Matrix<Real> eigenD(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<Real> work(1);
  LAPACK::geev(&jl, &jr, &Nl, prdMat.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenvectors.data(), &Nl,
               &work[0], &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::geev(&jl, &jr, &Nl, prdMat.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenvectors.data(), &Nl,
               &work[0], &lwork, &info);
  if (info != 0)
    throw std::runtime_error("Invalid Matrix Diagonalization Function, geev info = " + std::to_string(info));

  for (int i = 0; i < Nl; i++)
  {
    eigenvals[i] = alphar[i];
  }
}


LinearMethod::Real LinearMethod::getLowestEigenvector(Matrix<Real>& A, std::vector<Real>& ev) const
{
  int Nl(ev.size());
  //   Getting the optimal worksize
  Real zerozero = A(0, 0);
  char jl('N');
  char jr('V');
  std::vector<Real> alphar(Nl), alphai(Nl), beta(Nl);
  Matrix<Real> eigenT(Nl, Nl);
  Matrix<Real> eigenD(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<Real> work(1);
  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  if (info != 0)
  {
    APP_ABORT("Invalid Matrix Diagonalization Function!");
  }

  // Filter and sort to find desired eigenvalue.
  // Filter accepts eigenvalues between E_0 and E_0 - 100.0,
  // where E_0 is H(0,0), the current estimate for the VMC energy.
  // Sort searches for eigenvalue closest to E_0 - 2.0

  bool found_any_eigenvalue = false;
  std::vector<std::pair<Real, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    Real evi(alphar[i]);
    if ((evi < zerozero) && (evi > (zerozero - 1e2)))
    {
      mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
      mappedEigenvalues[i].second = i;
      found_any_eigenvalue        = true;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<Real>::max();
      mappedEigenvalues[i].second = i;
    }
  }

  // Sometimes there is no eigenvalue less than E_0, but there is one slightly higher.
  // Run filter and sort again, except this time accept eigenvalues between E_0 + 100.0 and E_0 - 100.0.
  // Since it's already been determined there are no eigenvalues below E_0, the sort
  // finds the eigenvalue closest to E_0.
  if (!found_any_eigenvalue)
  {
    app_log() << "No eigenvalues passed initial filter. Trying broader 100 a.u. filter" << std::endl;

    bool found_higher_eigenvalue;
    for (int i = 0; i < Nl; i++)
    {
      Real evi(alphar[i]);
      if ((evi < zerozero + 1e2) && (evi > (zerozero - 1e2)))
      {
        mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
        mappedEigenvalues[i].second = i;
        found_higher_eigenvalue     = true;
      }
      else
      {
        mappedEigenvalues[i].first  = std::numeric_limits<Real>::max();
        mappedEigenvalues[i].second = i;
      }
    }
    if (!found_higher_eigenvalue)
    {
      app_log() << "No eigenvalues passed second filter. Optimization is likely to fail." << std::endl;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  //         for (int i=0; i<4; i++) app_log()<<i<<": "<<alphar[mappedEigenvalues[i].second]<< std::endl;
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenT(mappedEigenvalues[0].second, i) / eigenT(mappedEigenvalues[0].second, 0);
  return alphar[mappedEigenvalues[0].second];
  //     }
}

void LinearMethod::getNonLinearRange(int& first, int& last, const QMCCostFunctionBase& optTarget) const
{
  std::vector<int> types;
  optTarget.getParameterTypes(types);
  first = 0;
  last  = types.size();
  //assume all non-linear coeffs are together.
  if (types[0] == optimize::LINEAR_P)
  {
    int i(0);
    while (i < types.size())
    {
      if (types[i] == optimize::LINEAR_P)
        first = i;
      i++;
    }
    first++;
  }
  else
  {
    int i(types.size() - 1);
    while (i >= 0)
    {
      if (types[i] == optimize::LINEAR_P)
        last = i;
      i--;
    }
  }
  //     returns the number of non-linear parameters.
  //    app_log()<<"line params: "<<first<<" "<<last<< std::endl;
}

LinearMethod::Real LinearMethod::getNonLinearRescale(std::vector<Real>& dP,
                                                     Matrix<Real>& S,
                                                     const QMCCostFunctionBase& optTarget) const
{
  int first(0), last(0);
  getNonLinearRange(first, last, optTarget);
  if (first == last)
    return 1.0;
  Real rescale(1.0);
  Real xi(0.5);
  Real D(0.0);
  for (int i = first; i < last; i++)
    for (int j = first; j < last; j++)
      D += S(i + 1, j + 1) * dP[i + 1] * dP[j + 1];
  rescale = (1 - xi) * D / ((1 - xi) + xi * std::sqrt(1 + D));
  rescale = 1.0 / (1.0 - rescale);
  //     app_log()<<"rescale: "<<rescale<< std::endl;
  return rescale;
}


} // namespace qmcplusplus
