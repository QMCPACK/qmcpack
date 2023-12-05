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

#ifndef QMCPLUSPLUS_AFQMC_WALKERSETBASE_H
#define QMCPLUSPLUS_AFQMC_WALKERSETBASE_H

#include <random>
#include <type_traits>
#include <memory>

#include "Configuration.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/TimerManager.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Numerics/ma_blas.hpp"

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
template<class Alloc, class Ptr> //, class BPAlloc, class BPPtr>
class WalkerSetBase : public AFQMCInfo
{
protected:
  enum WalkerSetBaseTimers
  {
    LoadBalance_t,
    Branching_t
  };

  inline static const TimerNameList_t<WalkerSetBaseTimers> WalkerSetBaseTimerNames =
      {{LoadBalance_t, "WalkerSetBase::loadBalance"}, {Branching_t, "WalkerSetBase::branching"}};

  using element       = typename std::pointer_traits<Ptr>::element_type;
  using pointer       = Ptr;
  using const_element = const element;
  using const_pointer = const Ptr;
  using Allocator     = Alloc;

  using bp_element       = SPComplexType;  //typename std::pointer_traits<BPPtr>::element_type;
  using bp_pointer       = SPComplexType*; //BPPtr;
  using const_bp_element = const bp_element;
  using const_bp_pointer = const bp_pointer;
  using BPAllocator      = shared_allocator<bp_element>; //BPAlloc;

  using CMatrix       = boost::multi::array<element, 2, Allocator>;
  using BPCMatrix     = boost::multi::array<bp_element, 2, BPAllocator>;
  using BPCVector_ref = boost::multi::array_ref<bp_element, 1, bp_pointer>;
  using BPCMatrix_ref = boost::multi::array_ref<bp_element, 2, bp_pointer>;
  using BPCTensor_ref = boost::multi::array_ref<bp_element, 3, bp_pointer>;

  using stdCMatrix_ptr = boost::multi::array_ptr<bp_element, 2>;
  using stdCTensor_ptr = boost::multi::array_ptr<bp_element, 3>;

public:
  // contiguous_walker = true means that all the data of a walker is continguous in memory
  static const bool contiguous_walker = true;
  // contiguous_storage = true means that the data of all walker is continguous in memory
  static const bool contiguous_storage = true;
  static const bool fixed_population   = true;

  using reference = walker<pointer>;
  using iterator  = walker_iterator<pointer>;
  //using const_reference = const_walker<const_pointer>;
  //using const_iterator = const_walker_iterator<const_pointer>;
  using const_reference = walker<pointer>;
  using const_iterator  = walker_iterator<pointer>;

  /// constructor
  WalkerSetBase(afqmc::TaskGroup_& tg_,
                xmlNodePtr cur,
                AFQMCInfo& info,
                RandomBase<RealType>& r,
                Allocator alloc_,
                BPAllocator bpalloc_)
      : AFQMCInfo(info),
        TG(tg_),
        rng(r),
        walker_size(1),
        walker_memory_usage(0),
        bp_walker_size(0),
        bp_walker_memory_usage(0),
        bp_pos(-1),
        history_pos(0),
        walkerType(UNDEFINED_WALKER_TYPE),
        tot_num_walkers(0),
        Timers(getGlobalTimerManager(), WalkerSetBaseTimerNames, timer_level_coarse),
        walker_buffer({0, 1}, alloc_),
        bp_buffer({0, 0}, bpalloc_),
        load_balance(UNDEFINED_LOAD_BALANCE),
        pop_control(UNDEFINED_BRANCHING),
        min_weight(0.05),
        max_weight(4.0)
  {
    parse(cur);
    setup();
  }

  /// destructor
  ~WalkerSetBase() {}

  WalkerSetBase(WalkerSetBase const& other)            = delete;
  WalkerSetBase(WalkerSetBase&& other)                 = default;
  WalkerSetBase& operator=(WalkerSetBase const& other) = delete;
  WalkerSetBase& operator=(WalkerSetBase&& other)      = delete;

  /*
   * Returns the current number of walkers in the set.
   */
  int size() const { return tot_num_walkers; }

  /*
   * Returns the maximum number of walkers in the set that can be stored without reallocation.
   */
  int capacity() const { return int(std::get<0>(walker_buffer.sizes())); }

  /*
   * Returns the maximum number of fields in the set that can be stored without reallocation. 
   */
  int NumBackProp() const { return wlk_desc[3]; }
  /*
   * Returns the maximum number of cholesky vectors in the set that can be stored without reallocation. 
   */
  int NumCholVecs() const { return wlk_desc[4]; }
  /*
   * Returns the length of the history buffers. 
   */
  int HistoryBufferLength() const { return wlk_desc[6]; }

  /*
   * Returns the position of the insertion point in the BP stack. 
   */
  int getBPPos() const { return bp_pos; }
  void setBPPos(int p) { bp_pos = p; }
  void advanceBPPos() { bp_pos++; }

  /*
   * Returns, sets and advances the position of the insertion point in the History circular buffers. 
   */
  int getHistoryPos() const { return history_pos; }
  void setHistoryPos(int p) { history_pos = p % wlk_desc[6]; }
  void advanceHistoryPos() { history_pos = (history_pos + 1) % wlk_desc[6]; }


  /*
   * Returns iterator to the first walker in the set
   */
  iterator begin()
  {
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    return iterator(0, boost::multi::static_array_cast<element, pointer>(walker_buffer), data_displ, wlk_desc);
  }

  /*
   * Returns iterator to the first walker in the set
   */
  const_iterator begin() const
  {
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    return const_iterator(0, boost::multi::static_array_cast<element, pointer>(walker_buffer), data_displ, wlk_desc);
  }


  /*
   * Returns iterator to the past-the-end walker in the set
   */
  iterator end()
  {
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    return iterator(tot_num_walkers, boost::multi::static_array_cast<element, pointer>(walker_buffer), data_displ,
                    wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  reference operator[](int i)
  {
    if (i < 0 || i > tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    return reference(boost::multi::static_array_cast<element, pointer>(walker_buffer)[i], data_displ, wlk_desc);
  }

  /*
   * Returns a reference to a walker
   */
  const_reference operator[](int i) const
  {
    if (i < 0 || i > tot_num_walkers)
      APP_ABORT("error: index out of bounds.\n");
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    return const_reference(boost::multi::static_array_cast<element, pointer>(walker_buffer)[i], data_displ, wlk_desc);
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
  void resize(int n, MatA&& A, MatB&& B)
  {
    assert(std::get<0>(A.sizes()) == wlk_desc[0]);
    assert(std::get<1>(A.sizes()) == wlk_desc[1]);
    if (walkerType == COLLINEAR)
    {
      assert(std::get<0>(B.sizes()) == wlk_desc[0]);
      assert(std::get<1>(B.sizes()) == wlk_desc[2]);
    }
    reserve(n);
    if (n > tot_num_walkers)
    {
      if (TG.TG_local().root())
      {
        auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
        auto pos = tot_num_walkers;
        // careful here!!!
        while (pos < n)
        {
          using std::fill_n;
          fill_n(W[pos].origin(), W[pos].size(), ComplexType(0, 0));
          reference w0(W[pos], data_displ, wlk_desc);
          //w0.SlaterMatrix(Alpha) = A;
          auto&& SM_(*w0.SlaterMatrix(Alpha));
          SM_ = A;
          if (walkerType == COLLINEAR)
          {
            //w0.SlaterMatrix(Beta) = B;
            auto&& SMB_(*w0.SlaterMatrix(Beta));
            SMB_ = B;
          }
          pos++;
        }
        // use operator= or assign when ready!!!
        boost::multi::array<ComplexType, 1> buff(iextensions<1u>{n - tot_num_walkers}, ComplexType(1.0));
        ma::copy(buff, W({tot_num_walkers, n}, data_displ[WEIGHT]));
        ma::copy(buff, W({tot_num_walkers, n}, data_displ[OVLP]));
        ma::copy(buff, W({tot_num_walkers, n}, data_displ[PHASE]));
      }
    }
    tot_num_walkers = n;
    targetN_per_TG  = tot_num_walkers;
    targetN         = GlobalPopulation();
    if (targetN != targetN_per_TG * TG.getNumberOfTGs())
    {
      std::cerr << " targetN, targetN_per_TG, # of TGs: " << targetN << " " << targetN_per_TG << " "
                << TG.getNumberOfTGs() << std::endl;
      std::cout << " targetN, targetN_per_TG, # of TGs: " << targetN << " " << targetN_per_TG << " "
                << TG.getNumberOfTGs() << std::endl;
      APP_ABORT("Error in WalkerSetBase::resize(n,A,B).\n");
    }
  }

  void resize_bp(int nbp, int nCV, int nref)
  {
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    assert(bp_buffer.size() == bp_walker_size);
    assert(walker_buffer.size() == std::get<1>(bp_buffer.sizes()));
    // wlk_descriptor: {nmo, naea, naeb, nback_prop, nCV, nRefs, nHist}
    wlk_desc[3] = nbp;
    wlk_desc[4] = nCV;
    wlk_desc[5] = nref;
    wlk_desc[6] = 3 * nbp;
    int ncol    = NAEA;
    int nrow    = NMO;
    if (walkerType != CLOSED)
    {
      if (walkerType == COLLINEAR)
      {
        ncol += NAEB;
      }
      else if (walkerType == NONCOLLINEAR)
      {
        nrow += NMO;
        ncol += NAEB;
      }
      else
      {
        app_error() << " Error: Incorrect walker_type on WalkerSetBase::setup \n";
        APP_ABORT("");
      }
    }
    // store nbpx3 history of weights and factors in circular buffer
    int cnt            = 0;
    data_displ[FIELDS] = cnt;
    cnt += nbp * nCV;
    data_displ[WEIGHT_FAC] = cnt;
    cnt += wlk_desc[6];
    data_displ[WEIGHT_HISTORY] = cnt;
    cnt += wlk_desc[6];
    bp_walker_size = cnt;
    if (std::get<0>(bp_buffer.sizes()) != bp_walker_size)
    {
      bp_buffer.reextent({bp_walker_size, std::get<0>(walker_buffer.sizes())});
      using std::fill_n;
      fill_n(bp_buffer.origin() + data_displ[WEIGHT_FAC] * std::get<1>(bp_buffer.sizes()),
             wlk_desc[6] * std::get<1>(bp_buffer.sizes()), bp_element(1.0));
    }
    if (nbp > 0 && (data_displ[SMN] < 0 || data_displ[SM_AUX] < 0))
    {
      auto sz(walker_size);
      data_displ[SMN] = walker_size;
      walker_size += nrow * ncol;
      data_displ[SM_AUX] = walker_size;
      walker_size += nrow * ncol;
      CMatrix wb({std::get<0>(walker_buffer.sizes()), walker_size}, walker_buffer.get_allocator());
      ma::copy(walker_buffer, wb(wb.extension(0), {0, sz}));
      walker_buffer = std::move(wb);
    }
  }

  // perform and report tests/timings
  void benchmark(std::string& blist, int maxnW, int delnW, int repeat);

  int get_TG_target_population() const { return targetN_per_TG; }
  int get_global_target_population() const { return targetN; }

  std::pair<int, int> walker_dims() const { return std::pair<int, int>{wlk_desc[0], wlk_desc[1]}; }

  int GlobalPopulation() const
  {
    int res = 0;
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    if (TG.TG_local().root())
      res += tot_num_walkers;
    return (TG.Global() += res);
  }

  RealType GlobalWeight() const
  {
    RealType res = 0;
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    if (TG.TG_local().root())
    {
      boost::multi::array<ComplexType, 1> buff(iextensions<1u>{tot_num_walkers});
      getProperty(WEIGHT, buff);
      for (int i = 0; i < tot_num_walkers; i++)
        res += std::abs(buff[i]);
    }
    return (TG.Global() += res);
  }

  // population control algorithm
  void popControl(std::vector<ComplexType>& curData);

  template<class Mat>
  void push_walkers(Mat&& M)
  {
    static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
    if (tot_num_walkers + M.size() > capacity())
      APP_ABORT("Insufficient capacity");
    if (single_walker_size() + single_walker_bp_size() != std::get<1>(M.sizes()))
      APP_ABORT("Incorrect dimensions.");
    if (M.stride(1) != 1)
      APP_ABORT("Incorrect strides.");
    if (!TG.TG_local().root())
    {
      tot_num_walkers += M.size();
      return;
    }
    auto&& W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    auto&& BPW(boost::multi::static_array_cast<bp_element, bp_pointer>(bp_buffer));
    for (int i = 0; i < M.size(); i++)
    {
      W[tot_num_walkers] = M[i].sliced(0, walker_size);
      if (wlk_desc[3] > 0)
        BPW(BPW.extension(0), tot_num_walkers) = M[i].sliced(walker_size, walker_size + bp_walker_size);
      tot_num_walkers++;
    }
  }

  template<class Mat>
  void pop_walkers(Mat&& M)
  {
    static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
    if (tot_num_walkers < int(M.size()))
      APP_ABORT("Insufficient walkers");
    if (wlk_desc[3] > 0)
    {
      if (walker_size + bp_walker_size != int(std::get<1>(M.sizes())))
        APP_ABORT("Incorrect dimensions.");
    }
    else
    {
      if (walker_size != int(std::get<1>(M.sizes())))
        APP_ABORT("Incorrect dimensions.");
    }
    if (M.stride(1) != 1)
      APP_ABORT("Incorrect strides.");

    if (!TG.TG_local().root())
    {
      tot_num_walkers -= int(M.size());
      return;
    }
    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    auto BPW(boost::multi::static_array_cast<bp_element, bp_pointer>(bp_buffer));
    for (int i = 0; i < M.size(); i++)
    {
      M[i].sliced(0, walker_size) = W[tot_num_walkers - 1];
      if (wlk_desc[3] > 0)
        M[i].sliced(walker_size, walker_size + bp_walker_size) = BPW(BPW.extension(0), tot_num_walkers - 1);
      tot_num_walkers--;
    }
  }

  // given a list of new weights and counts, add/remove walkers and reassign weight accordingly
  template<class Mat>
  void branch(std::vector<std::pair<double, int>>::iterator itbegin,
              std::vector<std::pair<double, int>>::iterator itend,
              Mat& M)
  {
    if (std::distance(itbegin, itend) != tot_num_walkers)
      APP_ABORT("Error in WalkerSetBase::branch(): ptr_range != # walkers. \n");

    // checking purposes
    int nW = 0;
    for (auto it = itbegin; it != itend; ++it)
      nW += it->second;
    if (int(std::get<0>(M.sizes())) < std::max(0, nW - targetN_per_TG))
    {
      std::cout << " Error in WalkerSetBase::branch(): Not enough space in excess matrix. \n"
                << std::get<0>(M.sizes()) << " " << nW << " " << targetN_per_TG << std::endl;
      APP_ABORT("Error in WalkerSetBase::branch(): Not enough space in excess matrix.\n");
    }
    if (int(std::get<1>(M.sizes())) < walker_size + ((wlk_desc[3] > 0) ? bp_walker_size : 0))
      APP_ABORT("Error in WalkerSetBase::branch(): Wrong dimensions in excess matrix.\n");

    // if all walkers are dead, don't bother with routine, reset tot_num_walkers and return
    if (nW == 0)
    {
      tot_num_walkers = 0;
      return;
    }

    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    auto BPW(boost::multi::static_array_cast<bp_element, bp_pointer>(bp_buffer));

    //1. push/swap all dead walkers to the end and adjust tot_num_walkers
    {
      auto kill = itbegin;
      auto keep = itend - 1;

      while (keep > kill)
      {
        // 1. look for next keep
        while (keep->second == 0 && keep > kill)
        {
          tot_num_walkers--;
          --keep;
        }
        if (keep == kill)
          break;

        // 2. look for next kill
        while (kill->second != 0 && kill < keep)
          ++kill;
        if (keep == kill)
          break;

        // 3. swap
        std::swap(*kill, *keep);
        W[std::distance(itbegin, kill)] = W[tot_num_walkers - 1];
        if (wlk_desc[3] > 0)
          BPW(BPW.extension(0), std::distance(itbegin, kill)) = BPW(BPW.extension(0), tot_num_walkers - 1);
        --tot_num_walkers;
        --keep;
      }

      // check
      int n = 0;
      for (auto it = itbegin; it != itbegin + tot_num_walkers; ++it)
        n += it->second;
      if (n != nW)
        APP_ABORT("Error in WalkerSetBase::branch(): Problems with walker counts after sort.\n");
      for (auto it = itbegin + tot_num_walkers; it != itend; ++it)
        if (it->second != 0)
          APP_ABORT("Error in WalkerSetBase::branch(): Problems after sort.\n");
      // you can also check the energy if things look incorrect
    }

    //2. Adjust weights and replicate walkers. Those beyond targetN_per_TG go in M
    itend   = itbegin + tot_num_walkers;
    int pos = 0;
    int cnt = 0;
    // circular buffer
    int his_pos = ((history_pos == 0) ? wlk_desc[6] - 1 : history_pos - 1);
    for (; itbegin != itend; ++itbegin, ++pos)
    {
      if (itbegin->second <= 0)
      { // just checking
        APP_ABORT("Error in WalkerSetBase::branch(): Problems during branch.\n");
      }
      else if (itbegin->second == 1)
      {
        //walker_buffer[pos][data_displ[WEIGHT]] = ComplexType(itbegin->first,0.0);
        // need synthetic references to make this easier!!!
        using std::fill_n;
        fill_n(W[pos].origin() + data_displ[WEIGHT], 1, ComplexType(itbegin->first, 0.0));
        if (wlk_desc[6] > 0 && his_pos >= 0 && his_pos < wlk_desc[6])
          fill_n(BPW[data_displ[WEIGHT_HISTORY] + his_pos].origin() + pos, 1, ComplexType(itbegin->first, 0.0));
      }
      else
      {
        // if there is space, branch within walker set
        // otherwise send excess to M
        int n = std::min(targetN_per_TG - tot_num_walkers, itbegin->second - 1);
        //walker_buffer[pos][data_displ[WEIGHT]] = ComplexType(itbegin->first,0.0);
        // need synthetic references to make this easier!!!
        using std::fill_n;
        fill_n(W[pos].origin() + data_displ[WEIGHT], 1, ComplexType(itbegin->first, 0.0));
        if (wlk_desc[6] > 0 && his_pos >= 0 && his_pos < wlk_desc[6])
          fill_n(BPW[data_displ[WEIGHT_HISTORY] + his_pos].origin() + pos, 1, ComplexType(itbegin->first, 0.0));
        for (int i = 0; i < n; i++)
        {
          W[tot_num_walkers] = W[pos];
          if (wlk_desc[3] > 0)
            BPW(BPW.extension(0), tot_num_walkers) = BPW(BPW.extension(0), pos);
          tot_num_walkers++;
        }
        for (int i = 0, in = itbegin->second - 1 - n; i < in; i++, cnt++)
        {
          M[cnt].sliced(0, walker_size) = W[pos];
          if (wlk_desc[3] > 0)
            M[cnt].sliced(walker_size, walker_size + bp_walker_size) = BPW(BPW.extension(0), pos);
        }
      }
    }
  }

  template<class T>
  void scaleWeight(const T& w0, bool scale_last_history = false)
  {
    if (!TG.TG_local().root())
      return;
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    ma::scal(ComplexType(w0), W({0, tot_num_walkers}, data_displ[WEIGHT]));
    if (scale_last_history)
    {
      int his_pos = ((history_pos == 0) ? wlk_desc[6] - 1 : history_pos - 1);
      if (wlk_desc[6] > 0 && his_pos >= 0 && his_pos < wlk_desc[6])
      {
        auto BPW(boost::multi::static_array_cast<bp_element, bp_pointer>(bp_buffer));
        ma::scal(bp_element(w0), BPW[data_displ[WEIGHT_HISTORY] + his_pos]);
      }
    }
  }

  void scaleWeightsByOverlap()
  {
    if (!TG.TG_local().root())
      return;
    assert(walker_buffer.size(1) == walker_size);
    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    boost::multi::array<ComplexType, 1> ov(iextensions<1u>{tot_num_walkers});
    boost::multi::array<ComplexType, 1> buff(iextensions<1u>{tot_num_walkers});
    getProperty(OVLP, ov);
    for (int i = 0; i < tot_num_walkers; i++)
      buff[i] = ComplexType(1.0 / std::abs(ov[i]), 0.0);
    ma::axty(ComplexType(1.0), buff, W({0, tot_num_walkers}, data_displ[WEIGHT]));
    for (int i = 0; i < tot_num_walkers; i++)
      buff[i] = std::exp(ComplexType(0.0, -std::arg(ov[i])));
    ma::axty(ComplexType(1.0), buff, W({0, tot_num_walkers}, data_displ[PHASE]));
  }

  afqmc::TaskGroup_& getTG() { return TG; }

  int single_walker_memory_usage() const { return walker_memory_usage; }
  int single_walker_size() const { return walker_size; }
  int single_walker_bp_memory_usage() const { return (wlk_desc[3] > 0) ? bp_walker_memory_usage : 0; }
  int single_walker_bp_size() const { return (wlk_desc[3] > 0) ? bp_walker_size : 0; }

  WALKER_TYPES getWalkerType() { return walkerType; }

  int walkerSizeIO()
  {
    if (walkerType == COLLINEAR)
      return wlk_desc[0] * (wlk_desc[1] + wlk_desc[2]) + 7;
    else // since NAEB = 0 in both CLOSED and NONCOLLINEAR cases
      return wlk_desc[0] * wlk_desc[1] + 7;
    return 0;
  }

  // I am going to assume that the relevant data to be copied is continuous,
  // careful not to break this in the future
  template<class Vec>
  void copyToIO(Vec&& x, int n)
  {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    using std::copy_n;
    copy_n(W[n].origin(), walkerSizeIO(), x.origin());
  }

  template<class Vec>
  void copyFromIO(Vec&& x, int n)
  {
    assert(n < tot_num_walkers);
    assert(x.size() >= walkerSizeIO());
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    using std::copy_n;
    copy_n(x.origin(), walkerSizeIO(), W[n].origin());
  }

  template<class TVec>
  void getProperty(walker_data id, TVec&& v) const
  {
    static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
    if (v.num_elements() < tot_num_walkers)
      APP_ABORT("Error: getProperty(v):: v.size < tot_num_walkers.\n");
    auto W_(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    ma::copy(W_({0, tot_num_walkers}, data_displ[id]), v.sliced(0, tot_num_walkers));
  }

  template<class TVec>
  void setProperty(walker_data id, TVec&& v)
  {
    static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
    if (v.num_elements() < tot_num_walkers)
      APP_ABORT("Error: setProperty(v):: v.size < tot_num_walkers.\n");
    auto W_(boost::multi::static_array_cast<element, pointer>(walker_buffer));
    ma::copy(v.sliced(0, tot_num_walkers), W_({0, tot_num_walkers}, data_displ[id]));
  }

  void resetWeights()
  {
    TG.TG_local().barrier();
    if (TG.TG_local().root())
    {
      boost::multi::array<element, 1> w_(iextensions<1u>{tot_num_walkers}, ComplexType(1.0));
      setProperty(WEIGHT, w_);
    }
    TG.TG_local().barrier();
  }

  // Careful!!! This matrix returns an array_ref, NOT a copy!!!
  stdCMatrix_ptr getFields(int ip)
  {
    if (ip < 0 || ip > wlk_desc[3])
      APP_ABORT(" Error: index out of bounds in getFields. \n");
    int skip = (data_displ[FIELDS] + ip * wlk_desc[4]) * std::get<1>(bp_buffer.sizes());
    return stdCMatrix_ptr(to_address(bp_buffer.origin()) + skip, {wlk_desc[4], std::get<1>(bp_buffer.sizes())});
  }

  stdCTensor_ptr getFields()
  {
    return stdCTensor_ptr(to_address(bp_buffer.origin()) + data_displ[FIELDS] * std::get<1>(bp_buffer.sizes()),
                          {wlk_desc[3], wlk_desc[4], std::get<1>(bp_buffer.sizes())});
  }

  template<class Mat>
  void storeFields(int ip, Mat&& V)
  {
    static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
    auto&& F(*getFields(ip));
    if (V.stride(0) == std::get<1>(V.sizes()))
    {
      using std::copy_n;
      copy_n(V.origin(), F.num_elements(), F.origin());
    }
    else
      F = V;
  }

  stdCMatrix_ptr getWeightFactors()
  {
    return stdCMatrix_ptr(to_address(bp_buffer.origin()) + data_displ[WEIGHT_FAC] * std::get<1>(bp_buffer.sizes()),
                          {wlk_desc[6], std::get<1>(bp_buffer.sizes())});
  }

  stdCMatrix_ptr getWeightHistory()
  {
    return stdCMatrix_ptr(to_address(bp_buffer.origin()) + data_displ[WEIGHT_HISTORY] * std::get<1>(bp_buffer.sizes()),
                          {wlk_desc[6], std::get<1>(bp_buffer.sizes())});
  }

  double getLogOverlapFactor() const { return LogOverlapFactor; }
  // nx= {2:CLOSED&&COLLINEAR, 1:NONCOLLINEAR }
  // before: OV_full = exp( nx*LogOverlapFactor ) * OVold
  // after: OV_full = exp( nx*LogOverlapFactor+f ) * OVnew
  // OVnew = OVold * exp( -f )
  // LogOverlapFactor_new = LogOverlapFactor + f/nx
  void adjustLogOverlapFactor(const double f)
  {
    assert(std::get<1>(walker_buffer.sizes()) == walker_size);
    double nx = (walkerType == NONCOLLINEAR ? 1.0 : 2.0);
    if (TG.TG_local().root())
    {
      auto W(boost::multi::static_array_cast<element, pointer>(walker_buffer));
      ma::scal(ComplexType(std::exp(-f)), W({0, tot_num_walkers}, data_displ[OVLP]));
    }
    LogOverlapFactor += f / nx;
    TG.TG_local().barrier();
  }

protected:
  afqmc::TaskGroup_& TG;

  RandomBase<RealType>& rng;

  int walker_size, walker_memory_usage;
  int bp_walker_size, bp_walker_memory_usage;
  int bp_pos;
  int history_pos;

  // wlk_descriptor: {nmo, naea, naeb, nback_prop, nCV, nRefs, nHist}
  wlk_descriptor wlk_desc;
  wlk_indices data_displ;

  WALKER_TYPES walkerType;

  int targetN_per_TG;
  int targetN;
  int tot_num_walkers;

  // Stores an overall scaling factor for walker weights (assumed to be
  // consistent over all walker groups).
  // The actual overlap of a walker is exp(OverlapFactor)*wset[i].weight()
  // Notice that overlap ratios (which are always what matters) are
  // independent of this value.
  // If this value is changed, the overlaps of the walkers must be adjusted
  // This is needed for stability reasons in large systems
  // Note that this is stored on a "per spin" capacity
  double LogOverlapFactor = 0.0;

  TimerList_t Timers;

  // Contains main walker data needed for propagation
  CMatrix walker_buffer;

  // Contains stack of fields and slater matrix references for back propagation
  BPCMatrix bp_buffer;

  // reads xml and performs setup
  void parse(xmlNodePtr cur);

  // performs setup
  void setup();

  // load balance algorithm
  LOAD_BALANCE_ALGORITHM load_balance;

  // load balancing algorithm
  template<class Mat>
  void loadBalance(Mat&& M)
  {
    if (load_balance == SIMPLE)
    {
      if (TG.TG_local().root())
        afqmc::swapWalkersSimple(*this, std::forward<Mat>(M), nwalk_counts_old, nwalk_counts_new, TG.TG_heads());
    }
    else if (load_balance == ASYNC)
    {
      if (TG.TG_local().root())
        afqmc::swapWalkersAsync(*this, std::forward<Mat>(M), nwalk_counts_old, nwalk_counts_new, TG.TG_heads());
    }
    TG.local_barrier();
    // since tot_num_walkers is local, you need to sync it
    if (TG.TG_local().size() > 1)
      TG.TG_local().broadcast_n(&tot_num_walkers, 1, 0);
  }

  // branching algorithm
  BRANCHING_ALGORITHM pop_control;
  double min_weight, max_weight;

  std::vector<int> nwalk_counts_new, nwalk_counts_old;
};

} // namespace afqmc

} // namespace qmcplusplus

#include "AFQMC/Walkers/WalkerSetBase.icc"

#endif
