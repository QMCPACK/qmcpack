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

#include "AFQMC/Walkers/Walkers.hpp"
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

//  // wlk_descriptor: {nmo, naea, naeb, nback_prop} 
//  using wlk_descriptor = std::array<int,4>;
//  using wlk_indices = std::array<int,14>;
  using element = ComplexType;
  using pointer = element *; 
  using const_element = const ComplexType;
  using const_pointer = ComplexType const*; 

  using CMatrix = boost::multi::array<element,2>;
  using shmCMatrix = boost::multi::array<element,2,shared_allocator<element>>;

  public:

  // contiguous_walker = true means that all the data of a walker is continguous in memory
  static const bool contiguous_walker = true;
  // contiguous_storage = true means that the data of all walker is continguous in memory
  static const bool contiguous_storage = true;
  static const bool fixed_population = true;

  using reference = walker<pointer>; 
  using iterator = walker_iterator<pointer>; 
  //using const_reference = const_walker<const_pointer>; 
  //using const_iterator = const_walker_iterator<const_pointer>; 
  using const_reference = walker<pointer>; 
  using const_iterator = walker_iterator<pointer>; 

  /// constructor
  SharedWalkerSet(afqmc::TaskGroup_& tg_, xmlNodePtr cur, AFQMCInfo& info, 
        RandomGenerator_t* r):
                TG(tg_),AFQMCInfo(info),rng(r),
                walker_memory_usage(0),tot_num_walkers(0),
		walker_buffer({1,1},shared_allocator<element>{TG.TG_local()}),
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
    return int(walker_buffer.size(0)); 
  } 

  /*
   * Returns iterator to the first walker in the set
   */
  iterator begin() {
    assert(walker_buffer.size(1) == walker_size);
    return iterator(0,boost::multi::static_array_cast<element, pointer>(walker_buffer),data_displ,wlk_desc);
  }

  /*
   * Returns iterator to the past-the-end walker in the set
   */
  iterator end() {
    assert(walker_buffer.size(1) == walker_size);
    return iterator(tot_num_walkers,boost::multi::static_array_cast<element, pointer>(walker_buffer),data_displ,wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  reference operator[](int i) {
    if(i<0 || i>tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(walker_buffer.size(1) == walker_size);
    return reference(boost::multi::static_array_cast<element, pointer>(walker_buffer)[i],data_displ,wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  const_reference operator[](int i) const {
    if(i<0 || i>tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(walker_buffer.size(1) == walker_size);
    return const_reference(boost::multi::static_array_cast<element, pointer>(walker_buffer)[i],data_displ,wlk_desc);
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
    assert(walker_buffer.size(1) == walker_size);
    if(TG.TG_local().root())
      res += tot_num_walkers;
    return (TG.Global() += res);
  }

  RealType GlobalWeight() const {
    RealType res=0;
    assert(walker_buffer.size(1) == walker_size);
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
    assert(walker_buffer.size(1) == walker_size);
    for(int i=0; i<tot_num_walkers; i++) 
      walker_buffer[i][data_displ[WEIGHT]]*=w0; 
  } 

  void scaleWeightsByOverlap() {
    if(!TG.TG_local().root()) return;
    assert(walker_buffer.size(1) == walker_size);
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
    if(walkerType==CLOSED)
      return wlk_desc[0]*wlk_desc[1]+7;
    else if(walkerType==COLLINEAR)
      return wlk_desc[0]*(wlk_desc[1]+wlk_desc[2])+7;
    else if(walkerType==NONCOLLINEAR)
      return 2*wlk_desc[0]*(wlk_desc[1]+wlk_desc[2])+7;
    return 0; 
  }

  template<class Vec>
  void copyToIO(Vec&& x, int n) {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(walker_buffer.size(1) == walker_size);
    auto ptr = walker_buffer[n].origin(); // pointer to walker data
    auto xd = std::addressof(x[0]);
    xd[0] = ptr[WEIGHT];
    xd[1] = ptr[PHASE];
    xd[2] = ptr[PSEUDO_ELOC_];
    xd[3] = ptr[E1_];
    xd[4] = ptr[EXX_];
    xd[5] = ptr[EJ_];
    xd[6] = ptr[OVLP];
    if(walkerType==CLOSED) {
      std::copy_n(ptr+SM,wlk_desc[0]*wlk_desc[1],xd+7);
    } else if(walkerType==COLLINEAR) {
      std::copy_n(ptr+SM,wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]),xd+7);
    } else if(walkerType==NONCOLLINEAR) {
      std::copy_n(ptr+SM,2*wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]),xd+7);
    } else {
      APP_ABORT(" Error: Unknown walkerType.\n");
    }
  }

  template<class Vec>
  void copyFromIO(Vec&& x, int n) {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(walker_buffer.size(1) == walker_size);
    auto ptr = walker_buffer[n].origin(); // pointer to walker data
    auto xd = std::addressof(x[0]);
    ptr[WEIGHT] = xd[0];
    ptr[PHASE] = xd[1];
    ptr[PSEUDO_ELOC_] = xd[2];
    ptr[E1_] = xd[3];
    ptr[EXX_] = xd[4];
    ptr[EJ_] = xd[5];
    ptr[OVLP] = xd[6];
    if(walkerType==CLOSED) {
      std::copy_n(xd+7,wlk_desc[0]*wlk_desc[1],ptr+SM);
    } else if(walkerType==COLLINEAR) {
      std::copy_n(xd+7,wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]),ptr+SM);
    } else if(walkerType==NONCOLLINEAR) {
      std::copy_n(xd+7,2*wlk_desc[0]*(wlk_desc[1]+wlk_desc[2]),ptr+SM);
    } else {
      APP_ABORT(" Error: Unknown walkerType.\n");
    }
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

  shmCMatrix walker_buffer;

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
