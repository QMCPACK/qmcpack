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

template<class Vector, class T>
int get_excitation_number(bool getIndx, Vector& refc, Vector& confg, Vector* exct, T& ci, Vector& loc, Vector& Iwork) {
  int NE = refc.size();
  auto ref = refc.data();
  exct.clear();
  int cnt=0;
  if(getIndx) std::copy(refc.begin(),refc.end(),Iwork.begin());
  auto it = ref;
  for(int i=0; i<NE; i++,it++)
    if(!std::binary_search(confg.begin(),confg.end(),*it)) {
      if(getIndx) {
        exct.emplace_back(*it);
        loc[cnt]=i;
      }
      cnt++;
    }
  if(!getIndx)
    return cnt;
  it = confg.data();
  int cnt2=0;
  for(int i=0; i<NE; i++,it++)
    if(!std::binary_search(refc.begin(),refc.end(),*it)) {
      exct.emplace_back(*it);
      Iwork[loc[cnt2]]=*it;
      cnt2++;
    }
  assert(cnt==cnt2);
  // sort Iwork and count number of exchanges to determine permutation sign
  // sooo slow but sooo simple too
  for(int i=0; i<NE; i++)
    for(int j=i+1; j<NE; j++)
    {
      if(Iwork[j] < Iwork[i])
      {
        ci *= T(-1.0);
        std::swap(Iwork[i],Iwork[j]);
      }
    }
  return cnt;
}

// try putting this in shared memory later on
template<class I = int,
	 class VType = std::complex<double>,	
         class Alloc = boost::mpi3::intranode::allocator<I>,
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
      using const_reference = Integer const*;
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

  template< typename Integer>
  class Iterator_const: public boost::iterator_facade< Iterator<Integer>,
                                               Integer const*,
                                               std::random_access_iterator_tag,
                                               Integer const*,
                                               std::ptrdiff_t
                                             >
  {
    public:
      using difference_type = std::ptrdiff_t;
      using reference = Integer const*;
      using const_reference = Integer const*;
      using value_tupe = Integer*;

      Iterator(Integer * index, size_t d_) : p_index(index),D(d_) {}
      Iterator(Integer const* index, size_t d_) : p_index(index),D(d_) {}

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
  using Excitation_const_Iterator = Iterator_const<integer_type>;

  ph_excitations() = delete;

  // Note: terms_per_excitation[0] has special meaning, the number of electrons in the calculation.
  // coefficients[0] will store the reference configuration itself.
  ph_excitations(std::vector<size_t>& terms_per_excitation, Alloc alloc_ = Alloc{}):
                i_allocator_(alloc_),
                c_allocator_(alloc_),
                excitations(terms_per_excitation.size(),terms_per_excitation,i_allocator_),
                coefficients(terms_per_excitation.size(),0,c_allocator_)
  { 
    // to accomodate a different treatment of [0]
    terms_per_excitation[0]=1;
    coefficients.reserve(terms_per_excitation);
  }

  ph_excitations(ph_excitations const& other) = delete;
  ph_excitations(ph_excitations && other) = default;

  ph_excitations& operator=(ph_excitations const& other) = delete;
  ph_excitations& operator=(ph_excitations && other) = default;

  size_t maximum_excitation_number() const{ return coefficients.size(); } 
  size_t number_of_excitations(int n) const{ return coefficients.num_elements(n); } 
  size_t number_of_excitations() const{ return coefficients.num_elements(); } 

  template<typename integer1, typename integer2, typename value> 
  void emplace_back(integer n, integer2* indx, value ci) {
    assert(n < maximum_excitation_number());
    assert(n>0);
    coefficients.emplace_back(n,static_cast<VType>(ci));
    for(int i=0; i<2*n; i++, ++indx)	
      excitations.emplace_back(n,static_cast<integer_type>(*indx));	
  }

  template<typename integer1, typename integer2, typename value>
  void emplace_reference(integer1 N, integer2* indx, value ci) {
    coefficients.emplace_back(0,static_cast<VType>(ci));
    for(int i=0; i<N; i++, ++indx)  
      excitations.emplace_back(0,static_cast<integer_type>(*indx));  
  }

  Excitation_Iterator::const_reference reference_configuration() const{
    return excitations.values(0);
  }

  VType const* coefficients_begin(int n) const {
    assert(n < maximum_excitation_number());
    return coefficients.values(n);
  }

  VType const* coefficients_end(int n) const {
    assert(n < maximum_excitation_number());
    return coefficients.values() + (*coefficients.pointers_end(n));
  }

  Excitation_Iterator excitations_begin(int n) {
    assert(n < maximum_excitation_number());
    assert(n>0);
    return Excitation_Iterator(excitations.values(n),n); 
  }

  Excitation_Iterator excitations_end(int n) {
    assert(n < maximum_excitation_number());
    assert(n>0);
    return Excitation_Iterator(excitations.values()+(*excitations.pointers_end(n)),n);
  }

  Excitation_const_Iterator excitations_begin(int n) {
    assert(n < maximum_excitation_number());
    assert(n>0);
    return Excitation_const_Iterator(excitations.values(n),n);
  }

  Excitation_const_Iterator excitations_end(int n) {
    assert(n < maximum_excitation_number());
    assert(n>0);
    return Excitation_const_Iterator(excitations.values()+(*excitations.pointers_end(n)),n);
  }

  private:

  coeff_aos coefficients;
  index_aos excitations; 

}

}

}

#endif
