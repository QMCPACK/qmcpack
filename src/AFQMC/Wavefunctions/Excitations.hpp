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

#include <boost/iterator/iterator_facade.hpp>

#include "AFQMC/config.h"
#include "AFQMC/Matrix/array_of_sequences.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// try putting this in shared memory later on
template<class I = int,
	 class VType = std::complex<double>,	
         class Alloc = boost::mpi3::intranode::allocator<T>,
         class is_root = boost::mpi3::intranode::is_root>
struct ph_excitations
{
  public:
  using integer_type = I;

  private:
  using std::size_t;
  using IAllocator = Alloc;
  using CAllocator = typename Alloc::template rebind<VType>::other;
  using coeff_aos = ma::sparse::array_of_sequences<VType,int,CAllocator,is_root>;
  using index_aos = ma::sparse::array_of_sequences<integer_type,int,IAllocator,is_root>;

  IAllocator i_allocator_;
  CAllocator c_allocator_;

  template< typename Integer>
  class Iterator: public boost::iterator_facade< Iterator<Integer>,
                                               Integer*,
                                               std::random_access_iterator_tag,
                                               Integer*,
                                               std::ptrdiff_t
                                             >
  { 
    public: 
      using difference_type = std::ptrdiff_t;
      using reference = Integer*;
      using value_tupe = Integer*;
    
      Iterator(Integer* index, size_t d_) : p_index(index),D(d_) {}
        
        // What we implement is determined by the boost::forward_traversal_tag
    private: 
      friend class boost::iterator_core_access;
    
      void increment() { p_index+=2*D; }
    
      bool equal(Iterator const& other) const
      {   
        return this->p_index == other.p_index;
      }
    
      reference dereference() const
      {   
        return reference(p_index);
      }
    
      void decrement() { p_index-=2*D; }
    
      void advance(int n) { p_index+=2*D*n; }
    
      difference_type distance_to(Iterator const& z) const { return ((z.p_index-p_index)/2/D); }
  
    private:
      Integer* p_index;
      size_t D;
  };

  public:

  using Excitation_Iterator = Iterator<integer_type>;

  ph_excitations() = delete;

  ph_excitations(size_t max_, std::vector<size_t>& terms_per_excitation, Alloc alloc_ = Alloc{}):
                i_allocator_(alloc_),
                c_allocator_(alloc_),
                coefficients(max_,terms_per_excitation,c_allocator_),
                excitations(max_,terms_per_excitation,i_allocator_)
  {}

  ph_excitations(ph_excitations const& other) = delete;
  ph_excitations(ph_excitations && other) = default;

  ph_excitations& operator=(ph_excitations const& other) = delete;
  ph_excitations& operator=(ph_excitations && other) = default;

  size_t maximum_excitation_number() { return coefficients.size(); } 
  size_t number_of_excitations(int n=0) { return coefficients.num_elements(n); } 

  template<typename integer1, typename integer2, typename value> 
  void emplace_back(integer n, integer2* indx, value ci) {
    assert(n < maximum_excitation_number());
    coefficients.emplace_back(n,static_cast<VType>(ci));
    for(int i=0; i<2*n; i++, ++indx)	
      excitations.emplace_back(i,static_cast<integer_type>(*indx));	
  }

  VType* coefficients_begin(int n) {
    assert(n < maximum_excitation_number());
    return coefficients.values(n);
  }

  VType* coefficients_end(int n) {
    assert(n < maximum_excitation_number());
    return coefficients.values() + (*coefficients.pointers_end(n));
  }

  Excitation_Iterator excitations_begin(int n) {
    assert(n < maximum_excitation_number());
    return Excitation_Iterator(excitations.values(n),n); 
  }

  Excitation_Iterator excitations_end(int n) {
    assert(n < maximum_excitation_number());
    return Excitation_Iterator(excitations.values()+(*excitations.pointers_end(n)),n);
  }

  private:

  coeff_aos coefficients;
  index_aos excitations; 

}

}

}

#endif
