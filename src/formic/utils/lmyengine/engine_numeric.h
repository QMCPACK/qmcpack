///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/numeric/numeric.h
///
/// \brief   header file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CQMC_NUMERIC_HEADER
#define CQMC_NUMERIC_HEADER

#include<string>
#include<cassert>
#include<cmath>
#include<vector>
#include<set>

namespace cqmc {

  void unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, const bool correct, double & r, double & v);
  void mpi_unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, const bool correct, double & r, double & v);
  int my_round(const double d);


  
}

#endif
