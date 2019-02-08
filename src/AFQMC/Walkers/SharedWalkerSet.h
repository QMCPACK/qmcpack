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

#ifndef QMCPLUSPLUS_AFQMC_SHAREDWALKERSET_H
#define QMCPLUSPLUS_AFQMC_SHAREDWALKERSET_H

#include <random>
#include <type_traits>
#include <memory>

#include "Configuration.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "Utilities/NewTimer.h"
#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/Walkers/WalkerControl.hpp"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class SharedWalkerSet: public AFQMCInfo 
{
  enum SharedWalkerSetTimers { LoadBalance, PopControl };

  enum walker_data { SM, WEIGHT, PHASE, PSEUDO_ELOC_, E1_, EXX_, EJ_, OVLP, PROPAGATORS, HEAD, TAIL, SMN, COS_FAC, WEIGHT_FAC };

//#define MA_TEST

  using Wlk_Buff = boost::multi::array_ref<ComplexType,2>; 
  using const_Wlk_Buff = boost::multi::const_array_ref<ComplexType,2>; 

  // wlk_descriptor: {nmo, naea, naeb, nback_prop} 
  using wlk_descriptor = std::array<int,4>;
  using wlk_indices = std::array<int,14>;
  using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using CMatrix = boost::multi::array<ComplexType,2>;

  public:

  // contiguous_walker = true means that all the data of a walker is continguous in memory
  static const bool contiguous_walker = true;
  // contiguous_storage = true means that the data of all walker is continguous in memory
  static const bool contiguous_storage = true;
  static const bool fixed_population = true;

  using SMType = boost::multi::array_ref<ComplexType,2>;
  using const_SMType = boost::multi::array_ref<const ComplexType,2>;

  struct const_walker {

    public:
    
      template<class ma>
      const_walker(ma const& a, const wlk_indices& i_, const wlk_descriptor& d_): 
        w_(boost::multi::const_array_ref<ComplexType,1>(std::addressof(*a.origin()),extensions<1u>{a.size()})),indx(i_),desc(d_) 
      {

	static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
	assert(w_.strides()[0]==1);
      }
      
      ~const_walker() {}

      const_walker(const_walker&& other) = default;  
      const_walker(const_walker const& other) = default;  
      const_walker& operator=(const_walker&& other) = delete;  
      const_walker& operator=(const_walker const& other) = delete; 

      ComplexType const* base() const {return &w_[0]; }
      int size() const {return w_.shape()[0]; }
      const_SMType SlaterMatrix(SpinTypes s) const{ 
	if(desc[2] <= 0 && s!=Alpha)
	  APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
	return (s==Alpha)?(const_SMType((&w_[indx[SM]]),{desc[0],desc[1]})):
			  (const_SMType((&w_[indx[SM]])+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      const_SMType SlaterMatrixN(SpinTypes s) const {
        if(indx[SMN] < 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(const_SMType((&w_[indx[SMN]]),{desc[0],desc[1]})):
                          (const_SMType((&w_[indx[SMN]])+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      ComplexType weight() const { return w_[indx[WEIGHT]]; } 
      ComplexType phase() const { return w_[indx[PHASE]]; } 
      ComplexType pseudo_energy() const { return w_[indx[PSEUDO_ELOC_]]; } 
      ComplexType onebody_energy() const { return w_[indx[E1_]]; } 
      ComplexType exchange_energy() const { return w_[indx[EXX_]]; } 
      ComplexType coulomb_energy() const { return w_[indx[EJ_]]; } 
      ComplexType E1() const { return w_[indx[E1_]]; } 
      ComplexType EXX() const { return w_[indx[EXX_]]; } 
      ComplexType EJ() const { return w_[indx[EJ_]]; } 
      ComplexType energy() const { return w_[indx[E1_]]+w_[indx[EXX_]]+w_[indx[EJ_]]; } 
      ComplexType overlap() const { return w_[indx[OVLP]]; } 
      // propagators can not be spin dependent
      const_SMType BMatrix(int ip) const {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
        return const_SMType( &(w_[indx[PROPAGATORS] + desc[0]*desc[0]*ip]) ,
                                                        {desc[0],desc[0]});
      }
      //void incrementBMatrix() const {
        //if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          //APP_ABORT("error: access to uninitialized BP sector. \n");
        //}
        //auto ip = getHead();
        //if(ip < 0 || ip >= desc[3]) {
          //APP_ABORT("error: Index out of bounds.\n");
        //}
        //w_[indx[HEAD]] = ComplexType((ip+1)%desc[3],0);	
      //}
      //void decrementBMatrix() const {
        //if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          //APP_ABORT("error: access to uninitialized BP sector. \n");
        //}
        //auto ip = getHead();
        //if(ip < 0 || ip >= desc[3]) {
          //APP_ABORT("error: Index out of bounds.\n");
        //}
        //w_[indx[HEAD]] = ComplexType((ip-1+desc[3])%desc[3],0);  
      //}
      bool isBMatrixBufferFull() const {
        return getHead()==0;
      }
      int NumBackProp() const {
        return desc[3];
      }
      ComplexType BPWeightFactor() const {
        if(indx[WEIGHT_FAC] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        return w_[indx[WEIGHT_FAC]];
      }
      //ComplexType weightFactor() const {
        //if(indx[COS_FAC]) {
          //APP_ABORT("error: access to uninitialized BP sector. \n");
        //}
        //return w_[indx[WEIGHT_FAC] + getHead()];
      //}
      void copy_to_buffer(ComplexType* data) const {
        std::copy(base(),base()+size(),data);
      }	

    private:

      int getHead() const { return static_cast<int>(w_[indx[HEAD]].real()); }

      boost::multi::const_array_ref<ComplexType,1> w_;
      const wlk_indices& indx;
      const wlk_descriptor& desc;	 
  };

  struct walker {

    public:
    
      template<class ma>
      walker(ma&& a, const wlk_indices& i_, const wlk_descriptor& d_): 
        w_(boost::multi::array_ref<ComplexType,1>(std::addressof(*a.origin()),extensions<1u>{a.size()})),indx(i_),desc(d_) 
      {

	static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
	assert(w_.strides()[0]==1);
      }
      
      ~walker() {}

      walker(walker&& other) = default;  
      walker(walker const& other) = default;  
      walker& operator=(walker&& other) = delete;  
      walker& operator=(walker const& other) = delete; 

      ComplexType* base() {return &w_[0]; }
      int size() const {return w_.shape()[0]; }
      SMType SlaterMatrix(SpinTypes s) {
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType((&w_[indx[SM]]),{desc[0],desc[1]})):
              (SMType((&w_[indx[SM]])+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      SMType SlaterMatrixN(SpinTypes s) {
        if(indx[SMN] < 0)
          APP_ABORT("error: access to uninitialized BP sector. \n");
        if(desc[2] <= 0 && s!=Alpha)
          APP_ABORT("error:walker spin out of range in SM(SpinType).\n");
        return (s==Alpha)?(SMType((&w_[indx[SMN]]),{desc[0],desc[1]})):
                          (SMType((&w_[indx[SMN]])+desc[0]*desc[1],{desc[0],desc[2]}));
      }
      ComplexType& weight() { return w_[indx[WEIGHT]]; } 
      ComplexType& phase() { return w_[indx[PHASE]]; } 
      ComplexType& pseudo_energy() { return w_[indx[PSEUDO_ELOC_]]; } 
      ComplexType& onebody_energy() { return w_[indx[E1_]]; } 
      ComplexType& exchange_energy() { return w_[indx[EXX_]]; } 
      ComplexType& coulomb_energy() { return w_[indx[EJ_]]; } 
      ComplexType& E1() { return w_[indx[E1_]]; } 
      ComplexType& EXX() { return w_[indx[EXX_]]; } 
      ComplexType& EJ() { return w_[indx[EJ_]]; } 
      ComplexType energy() { return w_[indx[E1_]]+w_[indx[EXX_]]+w_[indx[EJ_]]; } 
      ComplexType& overlap() { return w_[indx[OVLP]]; } 
      // propagators can not be spin dependent
      SMType BMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3]) {
          APP_ABORT("error: Index out of bounds.\n");
        }
        return SMType(&(w_[indx[PROPAGATORS]+desc[0]*desc[0]*ip]),
                      {desc[0],desc[0]});
      }
      int NumBackProp() {
        return desc[3];
      }
      void incrementBMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
        w_[indx[HEAD]] = ComplexType((ip+1)%desc[3],0);
      }
      void decrementBMatrix() {
        if(indx[PROPAGATORS] < 0 || indx[HEAD] < 0 || desc[3] <= 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        auto ip = getHead();
        if(ip < 0 || ip >= desc[3])
          APP_ABORT("error: Index out of bounds.\n");
        w_[indx[HEAD]] = ComplexType((ip-1+desc[3])%desc[3],0);
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
          SMType B = SMType(&(w_[indx[PROPAGATORS]+desc[0]*desc[0]*ip]),
                            {desc[0],desc[0]});
          for(int i = 0; i < desc[0]; i++) {
            for(int j = 0; j < desc[0]; j++) {
              B[i][j] = ((i==j)?ComplexType(1.0,0.0):ComplexType(0.0,0.0));
            }
          }
        }
        BPWeightFactor() = ComplexType(1.0,0.0);
        setSlaterMatrixN();
      }
      // Weight factors for partial path restoration approximation.
      ComplexType& BPWeightFactor() {
        if(indx[WEIGHT_FAC] < 0) {
          APP_ABORT("error: access to uninitialized BP sector. \n");
        }
        return w_[indx[WEIGHT_FAC]];
      }
      void copy_to_buffer(ComplexType* data) {
        std::copy(base(),base()+size(),data);
      }	
      void copy_from_buffer(ComplexType* data) {
        std::copy(data,data+size(),base());
      }	
      // replaces Slater Matrix at timestep M+N to timestep N for back propagation.
      void setSlaterMatrixN() {
        SlaterMatrixN(Alpha) = SlaterMatrix(Alpha);
        if(desc[2] > 0) {
          SlaterMatrixN(Beta) = SlaterMatrix(Beta);
        }
      }

    private:

      int getHead() const { return static_cast<int>(w_[indx[HEAD]].real()); }

      boost::multi::array_ref<ComplexType,1> w_;
      const wlk_indices& indx;
      const wlk_descriptor& desc;	 
      //wlk_indices indx;
      //wlk_descriptor desc;	 

  };

  struct walker_iterator: 
	public boost::iterator_facade<
	  walker_iterator,
	  void,
	  std::random_access_iterator_tag,
	  walker,
	  std::ptrdiff_t 
	>
  {
    public:
    template<class WBuff>
    walker_iterator(int k, WBuff&& w_, const wlk_indices& i_, const wlk_descriptor& d_): 
	pos(k),W(std::addressof(*w_.origin()),w_.extensions()),indx(&i_),desc(&d_)  {}

    using difference_type = std::ptrdiff_t;
    using reference = walker;

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

  struct const_walker_iterator:
        public boost::iterator_facade<
          const_walker_iterator,
          void,
          std::random_access_iterator_tag,
          const_walker,
          std::ptrdiff_t
        >
  {
    public:
    template<class WBuff>
    const_walker_iterator(int k, WBuff&& w_, const wlk_indices& i_, const wlk_descriptor& d_):
        pos(k),W(std::addressof(*w_.origin()),w_.extensions()),indx(&i_),desc(&d_)  {}

    using difference_type = std::ptrdiff_t;
    using reference = const_walker;

    private:

    int pos;
    const_Wlk_Buff W;
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

  using reference = walker; 
  using const_reference = const_walker; 
  using iterator = walker_iterator; 
  using const_iterator = const_walker_iterator; 
//  using reverse_iterator = walker_reverse_iterator; 

  /// constructor
  SharedWalkerSet(afqmc::TaskGroup_& tg_, xmlNodePtr cur, AFQMCInfo& info, 
        RandomGenerator_t* r):
                TG(tg_),AFQMCInfo(info),rng(r),
                walker_memory_usage(0),tot_num_walkers(0),
		walker_buffer({1,1}),
		//walker_buffer({1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                load_balance(UNDEFINED_LOAD_BALANCE),
                pop_control(UNDEFINED_BRANCHING),min_weight(0.05),max_weight(4.0),
                walkerType(UNDEFINED_WALKER_TYPE),nback_prop(0)
  {
    parse(cur);
    setup(); 
  }

  /// destructor
  ~SharedWalkerSet() {}

  SharedWalkerSet(SharedWalkerSet const& other) = delete;
  SharedWalkerSet(SharedWalkerSet&& other) = default;
  SharedWalkerSet& operator=(SharedWalkerSet const& other) = delete;
  SharedWalkerSet& operator=(SharedWalkerSet&& other) = default;

  /*
   * Returns the current number of walkers in the set.
   */	
  int size() const { 
    return tot_num_walkers; 
  } 

  /*
   * Returns the maximum number of walkers in the set that can be stored without reallocation.
   */	
  int capacity() const{ 
    return int(walker_buffer.shape()[0]); 
  } 

  /*
   * Returns iterator to the first walker in the set
   */
  iterator begin() {
    assert(walker_buffer.shape()[1] == walker_size);
    return iterator(0,walker_buffer,data_displ,wlk_desc);
  }

  /*
   * Returns iterator to the past-the-end walker in the set
   */
  iterator end() {
    assert(walker_buffer.shape()[1] == walker_size);
    return iterator(tot_num_walkers,walker_buffer,data_displ,wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  reference operator[](int i) {
    if(i<0 || i>tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(walker_buffer.shape()[1] == walker_size);
    return walker(walker_buffer[i],data_displ,wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  const_reference operator[](int i) const {
    if(i<0 || i>tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(walker_buffer.shape()[1] == walker_size);
    return const_walker(walker_buffer[i],data_displ,wlk_desc);
  }

  // cleans state of object. 
  //   -erases allocated memory 
  bool clean(); 

  /*
   * Increases the capacity of the containers to n.
   */
  void reserve(int n);

  /*
   * Adds/removes the number of walkers in the set to match the requested value.
   * Walkers are removed from the end of the set 
   *     and buffer capacity remains unchanged in this case.
   * New walkers are initialized from already existing walkers in a round-robin fashion. 
   * If the set is empty, calling this routine will abort. 
   * Capacity is increased if necessary.
   * Target Populations are set to n.
   */  
  void resize(int n);

  /*
   * Adds/removes the number of walkers in the set to match the requested value.
   * Walkers are removed from the end of the set 
   *     and buffer capacity remains unchanged in this case.
   * New walkers are initialized from the supplied matrix. 
   * Capacity is increased if necessary.
   * Target Populations are set to n.
   */ 
  template<class MatA, class MatB>
  void resize(int n, MatA&& A, MatB&& B);

  // perform and report tests/timings
  void benchmark(std::string& blist,int maxnW,int delnW,int repeat);

  int get_TG_target_population() const{ return targetN_per_TG; }
  int get_global_target_population() const{ return targetN; }

  int getNBackProp() const { return nback_prop; }
  std::pair<int,int> walker_dims() const {
    return std::pair<int,int> {wlk_desc[0], wlk_desc[1]};
  }

  int GlobalPopulation() const{
    int res=0;
    assert(walker_buffer.shape()[1] == walker_size);
    if(TG.TG_local().root())
      res += tot_num_walkers;
    return (TG.Global() += res);
  }

  RealType GlobalWeight() const {
    RealType res=0;
    assert(walker_buffer.shape()[1] == walker_size);
    if(TG.TG_local().root()) {
      for(int i=0; i<tot_num_walkers; i++) 
        res += std::abs(walker_buffer[i][data_displ[WEIGHT]]);
    }
    return (TG.Global() += res);
  }

  // population control algorithm
  void popControl(std::vector<ComplexType>& curData); 

  template<class Mat>
  void push_walkers(Mat&& M);
  template<class Mat>
  void pop_walkers(Mat&& M);

  // given a list of new weights and counts, add/remove walkers and reassign weight accordingly
  template<class Mat>
  void branch(std::vector<std::pair<double,int>>::iterator itb,
              std::vector<std::pair<double,int>>::iterator ite,
              Mat& M  
              );

  template<class T>
  void scaleWeight(const T& w0) {
    if(!TG.TG_local().root()) return;
    assert(walker_buffer.shape()[1] == walker_size);
    for(int i=0; i<tot_num_walkers; i++) 
      walker_buffer[i][data_displ[WEIGHT]]*=w0; 
  } 

  void scaleWeightsByOverlap() {
    if(!TG.TG_local().root()) return;
    assert(walker_buffer.shape()[1] == walker_size);
    for(int i=0; i<tot_num_walkers; i++) {
      walker_buffer[i][data_displ[WEIGHT]] *= ComplexType(1.0/std::abs(walker_buffer[i][data_displ[OVLP]]),0.0);
      walker_buffer[i][data_displ[PHASE]] *= std::exp(ComplexType(0.0, -std::arg(walker_buffer[i][data_displ[OVLP]]))); 
    }
  }

  afqmc::TaskGroup_& getTG() { return TG; } 

  int single_walker_memory_usage() const{ return walker_memory_usage; } 
  int single_walker_size() const{ return walker_size; }  
 
  WALKER_TYPES getWalkerType() { return walkerType; } 

  int walkerSizeIO() {
    return wlk_desc[0]*(wlk_desc[1]+wlk_desc[2])+7;
  }

  template<class Vec>
  void copyToIO(Vec&& x, int n) {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(walker_buffer.shape()[1] == walker_size);
    auto ptr = walker_buffer[n].origin(); // pointer to walker data
    auto xd = std::addressof(x[0]);
    int sm_size = wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]);
    std::copy_n(ptr,sm_size,xd);
    xd[data_displ[WEIGHT]] = ptr[data_displ[WEIGHT]];
    xd[data_displ[PHASE]] = ptr[data_displ[PHASE]];
    xd[data_displ[PSEUDO_ELOC_]] = ptr[data_displ[PSEUDO_ELOC_]];
    xd[data_displ[E1_]] = ptr[data_displ[E1_]];
    xd[data_displ[EXX_]] = ptr[data_displ[EXX_]];
    xd[data_displ[EJ_]] = ptr[data_displ[EJ_]];
    xd[data_displ[OVLP]] = ptr[data_displ[OVLP]];
  }

  template<class Vec>
  void copyFromIO(Vec&& x, int n) {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(walker_buffer.shape()[1] == walker_size);
    auto ptr = walker_buffer[n].origin(); // pointer to walker data
    auto xd = std::addressof(x[0]);
    int sm_size = wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]);
    std::copy_n(xd,sm_size,ptr);
    ptr[data_displ[WEIGHT]] = xd[data_displ[WEIGHT]];
    ptr[data_displ[PHASE]] = xd[data_displ[PHASE]];
    ptr[data_displ[PSEUDO_ELOC_]] = xd[data_displ[PSEUDO_ELOC_]];
    ptr[data_displ[E1_]] = xd[data_displ[E1_]];
    ptr[data_displ[EXX_]] = xd[data_displ[EXX_]];
    ptr[data_displ[EJ_]] = xd[data_displ[EJ_]];
    ptr[data_displ[OVLP]] = xd[data_displ[OVLP]];
  }

  private:

  RandomGenerator_t* rng;

  int nback_prop;
  int walker_size, walker_memory_usage;

  // wlk_descriptor: {nmo, naea, naeb, nback_prop} 
  wlk_descriptor wlk_desc; 
  wlk_indices data_displ; 

  WALKER_TYPES walkerType; 

  int targetN_per_TG;
  int targetN;
  int tot_num_walkers;

  TimerList_t Timers;

  afqmc::TaskGroup_& TG;  

  CMatrix walker_buffer;
  //shmCMatrix walker_buffer;

  // reads xml and performs setup
  void parse(xmlNodePtr cur); 

  // performs setup
  void setup();

  // load balance algorithm
  LOAD_BALANCE_ALGORITHM load_balance; 

  // load balancing algorithm
  template<class Mat>
  void loadBalance(Mat&& M); 

  // branching algorithm
  BRANCHING_ALGORITHM pop_control; 
  double min_weight, max_weight;

  std::vector<int> nwalk_counts_new, nwalk_counts_old;

};

}

}

#include "AFQMC/Walkers/SharedWalkerSet.icc"

#endif
