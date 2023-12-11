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
template<class Ptr>
struct walker
{
public:
  using pointer = Ptr;
  using element = typename std::pointer_traits<pointer>::element_type;
  using SMType  = boost::multi::array_ptr<element, 2, pointer>;

  template<class ma>
  walker(ma&& a, const wlk_indices& i_, const wlk_descriptor& d_)
      : w_(a.origin(), iextensions<1u>{a.size()}), indx(i_), desc(d_)
  {
    static_assert(std::decay<ma>::type::dimensionality == 1, "Wrong dimensionality");
  }

  ~walker() {}

  /*
      walker(walker&& other): w_(other.w_.origin(), iextensions<1u>{other.w_.size()}), 
                              indx(other.indx),desc(other.desc)  {} 
      walker(walker const& other): w_(other.w_.origin(),iextensions<1u>{other.w_.size()}), 
                              indx(other.indx),desc(other.desc)  {} 
*/
  // no copy/move assignment
  walker(walker&& other)      = default;
  walker(walker const& other) = default;
  walker& operator=(walker&& other) = delete;
  walker& operator=(walker const& other) = delete;

  pointer base() { return (*w_).origin(); }
  int size() const { return (*w_).size(0); }
  SMType SlaterMatrix(SpinTypes s)
  {
    if (desc[2] <= 0 && s != Alpha)
      APP_ABORT("error:walker spin out of range in SlaterMatrix(SpinType).\n");
    return (s == Alpha) ? (SMType(getw_(SM), {desc[0], desc[1]}))
                        : (SMType(getw_(SM) + desc[0] * desc[1], {desc[0], desc[2]}));
  }
  SMType SlaterMatrixN(SpinTypes s)
  {
    if (indx[SMN] < 0)
      APP_ABORT("error: access to uninitialized BP sector. \n");
    if (desc[2] <= 0 && s != Alpha)
      APP_ABORT("error:walker spin out of range in SlaterMatrixN(SpinType).\n");
    return (s == Alpha) ? (SMType(getw_(SMN), {desc[0], desc[1]}))
                        : (SMType(getw_(SMN) + desc[0] * desc[1], {desc[0], desc[2]}));
  }
  SMType SlaterMatrixAux(SpinTypes s)
  {
    if (indx[SM_AUX] < 0)
      APP_ABORT("error: access to uninitialized BP sector. \n");
    if (desc[2] <= 0 && s != Alpha)
      APP_ABORT("error:walker spin out of range in SlaterMatrixAux(SpinType).\n");
    return (s == Alpha) ? (SMType(getw_(SM_AUX), {desc[0], desc[1]}))
                        : (SMType(getw_(SM_AUX) + desc[0] * desc[1], {desc[0], desc[2]}));
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
  element energy() { return *getw_(E1_) + *getw_(EXX_) + *getw_(EJ_); }
  pointer overlap() { return getw_(OVLP); }
  // replaces Slater Matrix at timestep M+N to timestep N for back propagation.
  void setSlaterMatrixN()
  {
    *SlaterMatrixN(Alpha) = *SlaterMatrix(Alpha);
    if (desc[2] > 0)
      *SlaterMatrixN(Beta) = *SlaterMatrix(Beta);
  }

private:
  boost::multi::array_ptr<element, 1, pointer> w_;
  const wlk_indices& indx;
  const wlk_descriptor& desc;

  pointer getw_(int P) const { return (*w_).origin() + indx[P]; }
};

template<class Ptr>
struct walker_iterator
    : public boost::
          iterator_facade<walker_iterator<Ptr>, void, std::random_access_iterator_tag, walker<Ptr>, std::ptrdiff_t>
{
public:
  template<class WBuff>
  walker_iterator(int k, WBuff&& w_, const wlk_indices& i_, const wlk_descriptor& d_)
      : pos(k), W(w_.origin(), w_.extensions()), indx(&i_), desc(&d_)
  {}

  using pointer         = Ptr;
  using element         = typename std::pointer_traits<pointer>::element_type;
  using Wlk_Buff        = boost::multi::array_ptr<element, 2, Ptr>;
  using difference_type = std::ptrdiff_t;
  using reference       = walker<Ptr>;

  /*
    walker_iterator(walker_iterator const& it):
        pos(it.pos),W(it.W.origin(),it.W.extensions()),indx(it.indx),desc(it.desc)
    {}

    walker_iterator(walker_iterator && it):
        pos(it.pos),W(it.W.origin(),it.W.extensions()),indx(it.indx),desc(it.desc)
    {}
*/

private:
  int pos;
  Wlk_Buff W;
  wlk_indices const* indx;
  wlk_descriptor const* desc;

  friend class boost::iterator_core_access;

  void increment() { ++pos; }
  void decrement() { --pos; }
  bool equal(walker_iterator const& other) const { return pos == other.pos; }
  reference dereference() const { return reference((*W)[pos], *indx, *desc); }
  void advance(difference_type n) { pos += n; }
  difference_type distance_to(walker_iterator other) const { return other.pos - pos; }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
