///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SPARSE_HIP_CATCH_ALL_HPP
#define AFQMC_SPARSE_HIP_CATCH_ALL_HPP

#include <cassert>
#include <typeinfo>
#include "AFQMC/config.0.h"

namespace device
{
template<typename T, typename Q, class ptrA, class ptrI, class ptrI2, class ptrB, class ptrC>
void csrmm(const char transa,
           const int M,
           const int N,
           const int K,
           T alpha,
           const char* matdescra,
           ptrA A,
           ptrI indx,
           ptrI2 pntrb,
           ptrI2 pntre,
           ptrB B,
           const int ldb,
           Q beta,
           ptrC C,
           const int ldc)
{
  std::cout << " types: "
            << "  T: " << typeid(alpha).name() << "\n"
            << "  Q: " << typeid(beta).name() << "\n"
            << "  ptrA: " << typeid(A).name() << "\n"
            << "  ptrI: " << typeid(indx).name() << "\n"
            << "  ptrI2: " << typeid(pntrb).name() << "\n"
            << "  ptrB: " << typeid(B).name() << "\n"
            << "  ptrC: " << typeid(C).name() << std::endl;
  print_stacktrace throw std::runtime_error("Error: Calling qmc_hip::csrmm catch all.");
}

template<typename T, typename Q, class ptrA, class ptrI, class ptrI2, class ptrB, class ptrC>
void csrmv(const char transa,
           const int M,
           const int K,
           T alpha,
           const char* matdescra,
           ptrA A,
           ptrI indx,
           ptrI2 pntrb,
           ptrI2 pntre,
           ptrB x,
           Q beta,
           ptrC y)
{
  std::cout << " types: "
            << "  T: " << typeid(alpha).name() << "\n"
            << "  Q: " << typeid(beta).name() << "\n"
            << "  ptrA: " << typeid(A).name() << "\n"
            << "  ptrI: " << typeid(indx).name() << "\n"
            << "  ptrI2: " << typeid(pntrb).name() << "\n"
            << "  ptrB: " << typeid(x).name() << "\n"
            << "  ptrC: " << typeid(y).name() << std::endl;
  print_stacktrace throw std::runtime_error("Error: Calling qmc_hip::csrmv catch all.");
}


} // namespace device

#endif
