//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERS_HPP
#define QMCPLUSPLUS_AFQMC_WALKERS_HPP

#include <random>
#include <type_traits>
#include <memory>

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
  template<class Ptr>
  struct const_walker {

    public:

      using pointer = const typename std::decay<Ptr>::type; 
      using element = typename std::pointer_traits<typename std::decay<Ptr>::type>::element_type;
//      using element = __element;
      using SMType = boost::multi::array_ref<element,2,pointer>;
    
      template<class ma>
      const_walker(ma const& a, const wlk_indices& i_, const wlk_descriptor& d_): 
        w_(boost::multi::array_ref<element,1,pointer>(a.origin(),iextensions<1u>{a.size()})),indx(i_),desc(d_) 
      {
	static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
	assert(stride(w_)==1);
      }
      
      ~const_walker() {}

      const_walker(const_walker&& other) = default;  
      const_walker(const_walker const& other) = default;  
      const_walker& operator=(const_walker&& other) = delete;  
      const_walker& operator=(const_walker const& other) = delete; 

      pointer base() const {return w_.origin(); }
      int size() const {return w_.size(0); }
      SMType SlaterMatrix(SpinTypes s) const{ 
	if(desc[2] <= 0 && s!=Alpha)
	  APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
	return (s==Alpha)?(SMType(get_(SM),{desc[0],desc[1]})):
			  (SMType(get_(SM)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      SMType SlaterMatrixN(SpinTypes s) const {
        if(indx[SMN] < 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType(get_(SMN),{desc[0],desc[1]})):
                          (SMType(get_(SMN)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      pointer weight() const { return get_(WEIGHT); } 
      pointer phase() const { return get_(PHASE); } 
      pointer pseudo_energy() const { return get_(PSEUDO_ELOC_); } 
      pointer onebody_energy() const { return get_(E1_); } 
      pointer exchange_energy() const { return get_(EXX_); } 
      pointer coulomb_energy() const { return get_(EJ_); } 
      pointer E1() const { return get_(E1_); } 
      pointer EXX() const { return get_(EXX_); } 
      pointer EJ() const { return get_(EJ_); } 
      element energy() const { return *get_(E1_)+*get_(EXX_)+*get_(EJ_); } 
      pointer overlap() const { return get_(OVLP); } 
      // propagators can not be spin dependent
      SMType BMatrix(int ip) const {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
        return SMType( get_(PROPAGATORS) + desc[0]*desc[0]*ip ,
                                                        {desc[0],desc[0]});
      }
      bool isBMatrixBufferFull() const {
        return getHead()==0;
      }
      int NumBackProp() const {
        return desc[3];
      }
      pointer BPWeightFactor() const {
        if(indx[WEIGHT_FAC] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        return get_(WEIGHT_FAC);
      }
      void copy_to_buffer(Ptr data) const {
        copy(base(),base()+size(),data);
      }	

    private:

      // needs new strategy for gpu!!!
      int getHead() const { return static_cast<int>(*get_(HEAD).real()); }

      boost::multi::array_ref<element,1,pointer> w_;
      const wlk_indices& indx;
      const wlk_descriptor& desc;	 

      pointer get_(int P) const { return w_.origin() + indx[P]; }   
  };
*/

  template<class Ptr>
  struct walker {

    public:

      using pointer = Ptr;
      using element = typename std::pointer_traits<pointer>::element_type;
      using SMType = boost::multi::array_ref<element,2,pointer>;
    
      template<class ma>
      walker(ma&& a, const wlk_indices& i_, const wlk_descriptor& d_): 
        w_(boost::multi::array_ref<element,1,Ptr>(a.origin(),iextensions<1u>{a.size()})),indx(i_),desc(d_) 
      {
	static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
	assert(stride(w_)==1);
      }
      
      ~walker() {}

      walker(walker&& other) = default;  
      walker(walker const& other) = default;  
      walker& operator=(walker&& other) = delete;  
      walker& operator=(walker const& other) = delete; 

      pointer base() {return w_.origin(); }
      int size() const {return w_.size(0); }
      SMType SlaterMatrix(SpinTypes s) {
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType(get_(SM),{desc[0],desc[1]})):
              (SMType(get_(SM)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      SMType SlaterMatrixN(SpinTypes s) {
        if(indx[SMN] < 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType(get_(SMN),{desc[0],desc[1]})):
                          (SMType(get_(SMN)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      pointer weight() { return get_(WEIGHT); } 
      pointer phase() { return get_(PHASE); } 
      pointer pseudo_energy() { return get_(PSEUDO_ELOC_); } 
      pointer onebody_energy() { return get_(E1_); } 
      pointer exchange_energy() { return get_(EXX_); } 
      pointer coulomb_energy() { return get_(EJ_); } 
      pointer E1() { return get_(E1_); } 
      pointer EXX() { return get_(EXX_); } 
      pointer EJ() { return get_(EJ_); } 
      element energy() { return *get_(E1_)+*get_(EXX_)+*get_(EJ_); } 
      pointer overlap() { return get_(OVLP); } 
      // propagators can not be spin dependent
      SMType BMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3]) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        return SMType(get_(PROPAGATORS)+desc[0]*desc[0]*ip,
                      {desc[0],desc[0]});
      }
      SMType BMatrix(int ip) const {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
        return SMType( get_(PROPAGATORS) + desc[0]*desc[0]*ip ,
                                                        {desc[0],desc[0]});
      }
      int NumBackProp() const {
        return desc[3];
      }
      void incrementBMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
// problems on GPU
        *get_(HEAD) = element((ip+1)%desc[3],0);
      }
      void decrementBMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
// problems on GPU
        *get_(HEAD) = element((ip-1+desc[3])%desc[3],0);
      }
      bool isBMatrixBufferFull() const {
        return getHead() == 0;
      }
      // Reset back propagation information. B = I, weight factors = 1.0.
      void resetForBackPropagation() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        int nbp = desc[3];
        for(int ip = 0; ip < nbp; ip++) {
          SMType B = SMType(get_(PROPAGATORS)+desc[0]*desc[0]*ip,
                            {desc[0],desc[0]});
          for(int i = 0; i < desc[0]; i++) {
            for(int j = 0; j < desc[0]; j++) {
              B[i][j] = ((i==j)?element(1.0):element(0.0));
            }
          }
        }
        *BPWeightFactor() = element(1.0);
        setSlaterMatrixN();
      }
      // Weight factors for partial path restoration approximation.
// problems on GPU
      pointer BPWeightFactor() {
        if(indx[WEIGHT_FAC] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        return get_(WEIGHT_FAC);
      }
      void copy_to_buffer(Ptr data) {
        copy(base(),base()+size(),data);
      }	
      void copy_from_buffer(Ptr data) {
        copy(data,data+size(),base());
      }	
      // replaces Slater Matrix at timestep M+N to timestep N for back propagation.
      void setSlaterMatrixN() {
        SlaterMatrixN(Alpha) = SlaterMatrix(Alpha);
        if(desc[2] > 0) {
          SlaterMatrixN(Beta) = SlaterMatrix(Beta);
        }
      }

    private:

      int getHead() const { return static_cast<int>(get_(HEAD)->real()); }

      boost::multi::array_ref<element,1,Ptr> w_;
      const wlk_indices& indx;
      const wlk_descriptor& desc;	 

      pointer get_(int P) const { return w_.origin() + indx[P]; }   

  };

  template<class Ptr>
  struct walker_iterator: 
	public boost::iterator_facade<
	  walker_iterator<Ptr>,
	  void,
	  std::random_access_iterator_tag,
	  walker<Ptr>,
	  std::ptrdiff_t 
	>
  {
    public:
    template<class WBuff>
    walker_iterator(int k, WBuff&& w_, const wlk_indices& i_, const wlk_descriptor& d_): 
	pos(k),W(w_.origin(),w_.extensions()),indx(&i_),desc(&d_)  {}

    using pointer = Ptr;
    using element = typename std::pointer_traits<pointer>::element_type;
    using Wlk_Buff = boost::multi::array_ref<element,2,Ptr>;
    using difference_type = std::ptrdiff_t;
    using reference = walker<Ptr>;

    private:

    int pos;
    Wlk_Buff W;
    wlk_indices const* indx;
    wlk_descriptor const* desc; 

    friend class boost::iterator_core_access;

    void increment(){++pos;}
    void decrement(){--pos;}
    bool equal(walker_iterator const& other) const{ return pos == other.pos; }
    reference dereference() const { return reference(W[pos],*indx,*desc);}
    void advance(difference_type n){pos += n;}
    difference_type distance_to(walker_iterator other) const{ return other.pos - pos; }
  };

/*
  template<class Ptr>
  struct const_walker_iterator:
        public boost::iterator_facade<
          const_walker_iterator<Ptr>,
          void,
          std::random_access_iterator_tag,
          const_walker<Ptr>,
          std::ptrdiff_t
        >
  {
    public:
    template<class WBuff>
    const_walker_iterator(int k, WBuff&& w_, const wlk_indices& i_, const wlk_descriptor& d_):
        pos(k),W(to_address(w_.origin()),w_.extensions()),indx(&i_),desc(&d_)  {}

    using pointer = const Ptr;
    using element = typename std::pointer_traits<pointer>::element_type;
    using Wlk_Buff = boost::multi::array_ref<element,2,pointer>;
    using difference_type = std::ptrdiff_t;
    using reference = const_walker<Ptr>;

    private:

    int pos;
    Wlk_Buff W;
    wlk_indices const* indx;
    wlk_descriptor const* desc;

    friend class boost::iterator_core_access;

    void increment(){++pos;}
    void decrement(){--pos;}
    bool equal(const_walker_iterator const& other) const{  return pos == other.pos;  }
    reference dereference() const { return reference(W[pos],*indx,*desc);}
    void advance(difference_type n){pos += n;}
    difference_type distance_to(const_walker_iterator other) const{ return other.pos - pos; }
  };
*/
}

}

#endif
