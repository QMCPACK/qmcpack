//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SPARSE_MATRIX_HELPER_H
#define AFQMC_SPARSE_MATRIX_HELPER_H

#include "AFQMC/Matrix/SparseMatrix.h"

#include <stdio.h>
#include <string>
#include <complex>

using std::string;
using std::complex;

namespace qmcplusplus
{

template<typename T> void output_data(T *data, int size)
{
  std::cout << "[ ";
  for (int i= 0; i < size; i++)
  {
    std::cout << data[i] << " ";
  }
  std::cout << "]" << std::endl;
}

template<typename T> void output_matrix(SparseMatrix<T> &M)
{
  std::cout << "Colms = ";
  output_data(M.column_data(), M.size());
  std::cout << "Row index = ";
  output_data(M.row_index(), M.rows()+1);
  std::cout << "Values = ";
  output_data(M.values(), M.size());
}

template <typename T> double realPart(T &a) { REQUIRE(false); return 0.0; }

template<> double realPart(const double &a)
{
  return a;
}

template<> double realPart(double &a)
{
  return a;
}

template<> double realPart(const complex<double> &a)
{
  return a.real();
}

template<> double realPart(complex<double> &a)
{
  return a.real();
}

}

#endif
