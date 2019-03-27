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
#include "AFQMC/Numerics/ma_blas_extensions.hpp"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{

namespace afqmc
{

  template<class Ptr, class BPPtr>
  struct walker {

    public:

      using pointer = Ptr;
      using element = typename std::pointer_traits<pointer>::element_type;
      using SMType = boost::multi::array_ref<element,2,pointer>;

      using bp_pointer = BPPtr;
      using bp_element = typename std::pointer_traits<bp_pointer>::element_type;
      using BPSMType = boost::multi::array_ref<bp_element,2,bp_pointer>;
      using BPVector = boost::multi::array_ref<bp_element,1,bp_pointer>;
    
      template<class ma, class bp_ma>
      walker(ma&& a, bp_ma&& bp_a, const wlk_indices& i_, const wlk_descriptor& d_): 
        w_(boost::multi::array_ref<element,1,pointer>(a.origin(),iextensions<1u>{a.size()})),
        bp_(boost::multi::array_ref<bp_element,1,bp_pointer>(bp_a.origin(),
                                                             iextensions<1u>{bp_a.size()})),
        indx(i_),desc(d_) 
      {
	static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
	static_assert(std::decay<bp_ma>::type::dimensionality == 1, "Wrong dimensionality");
	assert(stride(w_)==1);
	assert(stride(bp_)==1);
      }
      
      ~walker() {}

      walker(walker&& other) = default;  
      walker(walker const& other) = default;  
      walker& operator=(walker&& other) = delete;  
      walker& operator=(walker const& other) = delete; 

      pointer base() {return w_.origin(); }
      int size() const {return w_.size(0); }
      int bp_size() const {return bp_.size(0); }
      SMType SlaterMatrix(SpinTypes s) {
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType(getw_(SM),{desc[0],desc[1]})):
              (SMType(getw_(SM)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      SMType SlaterMatrixN(SpinTypes s) {
        if(indx[SMN] < 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType(getw_(SMN),{desc[0],desc[1]})):
                          (SMType(getw_(SMN)+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      pointer weight() { return getw_(WEIGHT); } 
      pointer phase() { return getw_(PHASE); } 
      pointer pseudo_energy() { return getw_(PSEUDO_ELOC_); } 
      pointer onebody_energy() { return getw_(E1_); } 
      pointer exchange_energy() { return getw_(EXX_); } 
      pointer coulomb_energy() { return getw_(EJ_); } 
      pointer E1() { return getw_(E1_); } 
      pointer EXX() { return getw_(EXX_); } 
      pointer EJ() { return getw_(EJ_); } 
      element energy() { return *getw_(E1_)+*getw_(EXX_)+*getw_(EJ_); } 
      pointer overlap() { return getw_(OVLP); } 
      int NumBackProp() const {
        return desc[3];
      }
      int NumCholVecs() const {
        return desc[4];
      }
      int NumReferences() const {
        return desc[5];
      }
      BPVector getField(int ip) {
        if(indx[FIELDS] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        if(ip < 0 || ip >= desc[3]) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        return BPVector(getbp_(FIELDS)+desc[4]*ip,iextensions<1u>{desc[4]});
      }
      template<class CVec>
      void setField(CVec&& X, int ip) {
	static_assert(std::decay<CVec>::type::dimensionality == 1, "Wrong dimensionality");
        if(indx[FIELDS] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        if(ip < 0) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        if( ip >= desc[3]) return;  // don't store if outside buffer length
        BPVector V(getbp_(FIELDS)+desc[4]*ip,iextensions<1u>{desc[4]});
        V = X;
      }
      // Weight factors for partial path restoration approximation.
      bp_pointer BPWeightFactor() {
        if(indx[WEIGHT_FAC] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        return get_(WEIGHT_FAC);
      }
      // replaces Slater Matrix at timestep M+N to timestep N for back propagation.
      void setSlaterMatrixN() {
        SlaterMatrixN(Alpha) = SlaterMatrix(Alpha);
        if(desc[2] > 0) {
          SlaterMatrixN(Beta) = SlaterMatrix(Beta);
        }
      }
      BPSMType getReference(int ip, SpinTypes s) {
        if(indx[SM_REFS] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        if(ip < 0 || ip >= desc[5]) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in getReference(ip,SpinType).\n");
        int skip(ip*desc[0]*(desc[1]+desc[2]);
        return (s==Alpha)?(BPSMType(getbp_(SM_REFS)+skip,{desc[0],desc[1]})):
              (BPSMType(getbp_(SM_REFS)+skip+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      template<class CMat>
      void setReference(int ip, SpinTypes s, CMat&& A) {
        if(indx[SM_REFS] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        if(ip < 0 || ip >= desc[5]) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in getReference(ip,SpinType).\n");
        int skip(ip*desc[0]*(desc[1]+desc[2]);
        return (s==Alpha)?(BPSMType(getbp_(SM_REFS)+skip,{desc[0],desc[1]})):
              (BPSMType(getbp_(SM_REFS)+skip+desc[0]*desc[1],{desc[0],desc[2]}));
      }


      template<class CVec>
      void setField(CVec&& X, int ip) {
        static_assert(std::decay<CVec>::type::dimensionality == 1, "Wrong dimensionality");
        if(indx[FIELDS] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        if(ip < 0) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        if( ip >= desc[3]) return;  // don't store if outside buffer length
        BPVector V(getbp_(FIELDS)+desc[4]*ip,iextensions<1u>{desc[4]});
        V = X;
      }

    private:

      boost::multi::array_ref<element,1,pointer> w_;
      boost::multi::array_ref<bp_element,1,bp_pointer> bp_;
      const wlk_indices& indx;
      const wlk_descriptor& desc;	 

      pointer getw_(int P) const { return w_.origin() + indx[P]; }   
      bp_pointer getbp_(int P) const { return bp_.origin() + indx[P]; }   

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

}

}

#endif
