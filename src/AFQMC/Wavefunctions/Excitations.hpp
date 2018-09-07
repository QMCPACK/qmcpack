//////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_EXCITATIONS_HPP
#define QMCPLUSPLUS_AFQMC_EXCITATIONS_HPP

#include "AFQMC/config.h"
#include "AFQMC/Matrix/csr_matrix.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// try putting this in shared memory later on
template<class I = int,
         class Alloc = boost::mpi3::intranode::allocator<T>,
         class is_root = boost::mpi3::intranode::is_root>
struct ph_excitations
{
  public:
  using integer_type = I;

  private:
  using std::size_t;
  using IAllocator = Alloc;
  using CAllocator = typename Alloc::template rebind<ComplexType>::other;
  using coeff_ucsr = ma::sparse::ucsr_matrix<ComplexType,integer_type,int,
                                CAllocator,is_root>;
  using index_ucsr = ma::sparse::ucsr_matrix<integer_type,integer_type,int,
                                IAllocator,is_root>;

  IAllocator i_allocator_;
  CAllocator c_allocator_;

  public:

  ph_excitations() = delete;

  ph_excitations(size_t max_, std::vector<size_t>& terms_per_excitation, Alloc alloc_ = Alloc{}):
                i_allocator_(alloc_),
                c_allocator_(alloc_),
                coefficients({max_,size_t(MAX_PH_EXCITATION_PER_CHANNEL)}, {0,0}, 
                             terms_per_excitation,c_allocator_),
                excitations({max_,size_t(MAX_PH_EXCITATION_PER_CHANNEL)}, {0,0}, 
                             terms_per_excitation,i_allocator_)
  {}

  ph_excitations(ph_excitations const& other) = default;
  ph_excitations(ph_excitations && other) = default;

  ph_excitations& operator=(ph_excitations const& other) = default;
  ph_excitations& operator=(ph_excitations && other) = default;

  int number_of_excitations(int n) { return coeffs[n].size(); }

  private:

  coeff_ucsr coefficients;
  index_ucsr excitations; 

}

}

}

#endif
