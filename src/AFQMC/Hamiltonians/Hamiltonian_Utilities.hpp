//////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIAN_UTILITIES_H
#define QMCPLUSPLUS_AFQMC_HAMILTONIAN_UTILITIES_H

#include <cstdlib>
#include <complex>
#include <iostream>
#include <vector>

#include "Configuration.h"
#include "AFQMC/config.h"

namespace qmcplusplus
{
namespace afqmc
{
inline void print_tuple(s1D<ValueType>& t)
{
  std::cout << "  -  " << std::get<0>(t) << " " << std::get<1>(t) << std::endl;
}

inline void print_tuple(s4D<ValueType>& v)
{
  std::cout << std::get<4>(v) << " " << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << " "
            << std::get<3>(v) << std::endl;
}

inline void print_Vs4D(std::vector<s4D<ValueType>>& v)
{
  for (int i = 0; i < v.size(); i++)
    std::cout << std::get<4>(v[i]) << " " << std::get<0>(v[i]) << " " << std::get<1>(v[i]) << " " << std::get<2>(v[i])
              << " " << std::get<3>(v[i]) << std::endl;
  std::cout << std::endl;
}

inline long mapUT(long i, long j, long N)
{
  if (j >= i)
    return N * i + j - (i * (i + 1)) / 2;
  else
    return N * j + i - (j * (j + 1)) / 2;
}

inline long mapUT_woD(long i, long j, long N)
{
  if (j == i)
  {
    APP_ABORT(" Error in mapUT_woD: This should not happen. \n");
  }
  else if (j > i)
    return N * i + j - (i * (i + 1)) / 2 - i - 1;
  return N * j + i - (j * (j + 1)) / 2 - j - 1;
}

inline int mapUT(int i, int j, int N)
{
  if (j >= i)
    return N * i + j - (i * (i + 1)) / 2;
  return N * j + i - (j * (j + 1)) / 2;
}

inline int mapUT_woD(int i, int j, int N)
{
  if (j == i)
  {
    APP_ABORT(" Error in mapUT_woD: This should not happen. \n");
  }
  else if (j > i)
    return N * i + j - (i * (i + 1)) / 2 - i - 1;
  return N * j + i - (j * (j + 1)) / 2 - j - 1;
}

inline IndexType Index2Mat(int NMO, IndexType i, IndexType j, bool GHF = false)
{
  if (GHF)
    return i * 2 * NMO + j;
  else
    return (i < NMO) ? (i * NMO + j) : (NMO * NMO + (i - NMO) * NMO + j - NMO);
}

inline int getSpinSector(const int NMO, const IndexType& i, const IndexType& j, const IndexType& k, const IndexType& l)
{
  if (i < NMO)
  {
    if (j < NMO)
      return 0; // <alpha,alpha | alpha,alpha>
    else
      return 1; // <alpha,beta  | alpha,beta >
  }
  else
  {
    if (j < NMO)
      return 2; // <beta,alpha | beta,alpha>
    else
      return 3; // <beta,beta  | beta,beta >
  }
}

inline bool goodSpinSector(const IndexType& i, const IndexType& j, const IndexType& k, const IndexType& l, int NT)
{
  if (i < NT)
  {
    if (j < NT) // <alpha,alpha | alpha,alpha>
      return (k < NT && l < NT);
    else // <alpha,beta  | alpha,beta >
      return (k < NT && l >= NT);
  }
  else
  {
    if (j < NT) // <beta,alpha | beta,alpha>
      return (k >= NT && l < NT);
    else // <beta,beta  | beta,beta >
      return (k >= NT && l >= NT);
  }
}

inline int getSpinSector(const int NMO, const IndexType& i, const IndexType& j)
{
  if (i < NMO)
    return 0;
  return 1;
}

inline IndexType Index2Col(const int NMO, const IndexType i)
{
#if defined(AFQMC_DEBUG)
// assert( ((i<NMO)&&(j<NMO)) || ((i>NMO)&&(j>NMO))   )
#endif
  return (i < NMO) ? (i) : (i - NMO);
}

inline bool find_smallest_permutation(s4D<ValueType>& val)
{
#if defined(QMC_COMPLEX)
  // jl < ik
  if (std::forward_as_tuple(std::get<1>(val), std::get<3>(val)) <
      std::forward_as_tuple(std::get<0>(val), std::get<2>(val)))
  {
    std::swap(std::get<0>(val), std::get<1>(val));
    std::swap(std::get<2>(val), std::get<3>(val));
  }
  // kl < ij
  if (std::forward_as_tuple(std::get<2>(val), std::get<3>(val)) <
      std::forward_as_tuple(std::get<0>(val), std::get<1>(val)))
  {
    std::swap(std::get<0>(val), std::get<2>(val));
    std::swap(std::get<1>(val), std::get<3>(val));
    std::get<4>(val) = ma::conj(std::get<4>(val));
    // jl < ik again since ij<->kl swap occurred
    if (std::forward_as_tuple(std::get<1>(val), std::get<3>(val)) <
        std::forward_as_tuple(std::get<0>(val), std::get<2>(val)))
    {
      std::swap(std::get<0>(val), std::get<1>(val));
      std::swap(std::get<2>(val), std::get<3>(val));
    }
    return true;
  }
  else
  {
    // only possibility is that l < i, since I know that the current i is smaller than j and k
    //
    if (std::forward_as_tuple(std::get<3>(val), std::get<2>(val)) <
        std::forward_as_tuple(std::get<0>(val), std::get<1>(val)))
    {
      std::swap(std::get<0>(val), std::get<3>(val));
      std::swap(std::get<2>(val), std::get<1>(val));
      std::get<4>(val) = ma::conj(std::get<4>(val));
      return true;
    }
    return false;
  }
#else
  // i < k
  if (std::get<2>(val) < std::get<0>(val))
    std::swap(std::get<0>(val), std::get<2>(val));
  // j < l
  if (std::get<3>(val) < std::get<1>(val))
    std::swap(std::get<1>(val), std::get<3>(val));
  // ik < jl
  if (std::get<1>(val) < std::get<0>(val))
  {
    std::swap(std::get<0>(val), std::get<1>(val));
    std::swap(std::get<2>(val), std::get<3>(val));
  }
  else if ((std::get<1>(val) == std::get<0>(val)) && (std::get<3>(val) < std::get<2>(val)))
    std::swap(std::get<2>(val), std::get<3>(val));
  return false;
#endif
}

} // namespace afqmc
} // namespace qmcplusplus

#endif
