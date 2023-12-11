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
#include <map>
#include "AFQMC/config.h"
#include "AFQMC/Matrix/array_of_sequences.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// CHEAT!!!
template<class TP, class integer>
size_t get_index(TP const& tp_, integer loc)
{
  if (loc == 0)
    return size_t(std::get<0>(tp_));
  else if (loc == 1)
    return size_t(std::get<1>(tp_));
  else
    throw std::runtime_error(" Error in qmcplusplus::afqmc::get_index<TP,integer>(). \n");
}

template<class Vector>
void push_excitation(Vector const& abij, Vector& v)
{
  if (abij.size() == 0)
    return;
  assert(v.size() % abij.size() == 0);
  size_t n = abij.size();
  for (typename Vector::iterator it = v.begin(); it < v.end(); it += n)
    if (std::equal(abij.begin(), abij.end(), it))
      return;
  v.insert(v.end(), abij.begin(), abij.end());
}

template<class Vector>
size_t find_excitation(Vector const& abij, Vector& v)
{
  if (abij.size() == 0)
    return 0; // this assumes that the reference is 0
  assert(v.size() % abij.size() == 0);
  size_t n   = abij.size();
  size_t loc = 0;
  for (typename Vector::iterator it = v.begin(); it < v.end(); it += n, loc++)
    if (std::equal(abij.begin(), abij.end(), it))
      return loc;
  APP_ABORT("Error: Sequence not found in find_excitation.\n");
  return 0;
}

template<class excitations>
std::map<int, int> find_active_space(bool single_list, excitations const& abij, int NMO, int NAEA, int NAEB)
{
  std::map<int, int> mo2active;
  for (int i = 0; i < 2 * NMO; i++)
    mo2active[i] = -1;
  std::vector<size_t> count(2 * NMO);
  // reference first
  auto refc = abij.reference_configuration();
  for (int i = 0; i < NAEA + NAEB; i++, ++refc)
    ++count[*refc];
  auto nex = abij.maximum_excitation_number();
  for (int n = 1; n < nex[0]; ++n)
  {
    auto it    = abij.alpha_begin(n);
    auto itend = abij.alpha_end(n);
    for (; it < itend; ++it)
    {
      auto exct = *it + n; // skip locations
      for (int ak = 0; ak < n; ++ak, ++exct)
        ++count[*exct];
    }
  }
  for (int n = 1; n < nex[1]; ++n)
  {
    auto it    = abij.beta_begin(n);
    auto itend = abij.beta_end(n);
    for (; it < itend; ++it)
    {
      auto exct = *it + n; // skip locations
      for (int ak = 0; ak < n; ++ak, ++exct)
        ++count[*exct];
    }
  }
  if (not single_list)
  {
    int ik = 0;
    for (int i = 0; i < NMO; ++i)
      if (count[i] > 0 || count[i + NMO] > 0)
      {
        if (count[i] > 0)
          mo2active[i] = ik;
        if (count[i + NMO] > 0)
          mo2active[i + NMO] = ik;
        ++ik;
      }
  }
  else
  {
    int ik = 0;
    for (int i = 0; i < NMO; ++i)
      if (count[i] > 0)
        mo2active[i] = ik++;
    ik = 0;
    for (int i = NMO; i < 2 * NMO; ++i)
      if (count[i] > 0)
        mo2active[i] = ik++;
  }
  return mo2active;
}

/*
 * - exct stores, for each electron excitaion, the location of the orbital in the reference 
 *    being excited and the index of the excited orbital.
 */
template<class Vector, class T>
int get_excitation_number(bool getIndx, Vector& refc, Vector& confg, Vector& exct, T& ci, Vector& Iwork)
{
  int NE = refc.size();
  exct.clear();
  int cnt = 0;
  if (getIndx)
    std::copy(refc.begin(), refc.end(), Iwork.begin());
  auto it = refc.data();
  assert(Iwork.size() >= refc.size());
  for (int i = 0; i < NE; i++, it++)
    if (!std::binary_search(confg.begin(), confg.end(), *it))
    {
      if (getIndx)
      {
        // store the location, NOT the index!!!
        //exct.emplace_back(*it);
        exct.emplace_back(i);
      }
      cnt++;
    }
  if (!getIndx)
    return cnt;
  it       = confg.data();
  int cnt2 = 0;
  for (int i = 0; i < NE; i++, it++)
    if (!std::binary_search(refc.begin(), refc.end(), *it))
    {
      exct.emplace_back(*it);
      Iwork[exct[cnt2]] = *it;
      cnt2++;
    }
  assert(cnt == cnt2);
  // sort Iwork and count number of exchanges to determine permutation sign
  // sooo slow but sooo simple too
  for (int i = 0; i < NE; i++)
    for (int j = i + 1; j < NE; j++)
    {
      if (Iwork[j] < Iwork[i])
      {
        ci *= T(-1.0);
        std::swap(Iwork[i], Iwork[j]);
      }
    }
  return cnt;
}

template<typename intT, class csr>
inline std::vector<size_t> get_nnz(csr const& PsiT_MO, intT* refc, size_t N, size_t shift)
{
  std::vector<size_t> res(N);
  for (size_t i = 0; i < N; i++)
    res[i] = PsiT_MO.num_non_zero_elements(*(refc + i) - shift);
  return res;
}


// try putting this in shared memory later on
template<class I       = int,
         class VType   = std::complex<double>,
         class Alloc   = shared_allocator<I>,
         class is_root = ma::sparse::is_root>
struct ph_excitations
{
public:
  using integer_type       = I;
  using configuration_type = std::tuple<int, int, VType>;

private:
  using IAllocator = Alloc;
  using CAllocator = typename Alloc::template rebind<configuration_type>::other;
  using confg_aos  = ma::sparse::array_of_sequences<configuration_type, int, CAllocator, is_root>;
  using index_aos  = ma::sparse::array_of_sequences<integer_type, int, IAllocator, is_root>;

  IAllocator i_allocator_;
  CAllocator c_allocator_;

  template<typename Integer>
  class Iterator
      : public boost::
            iterator_facade<Iterator<Integer>, Integer*, std::random_access_iterator_tag, Integer*, std::ptrdiff_t>
  {
  public:
    using difference_type = std::ptrdiff_t;
    using reference       = Integer*;
    using const_reference = Integer const*;
    using value_tupe      = Integer*;

    Iterator(Integer* index, size_t d_) : p_index(index), D(d_) {}

    // What we implement is determined by the boost::forward_traversal_tag
  private:
    friend class boost::iterator_core_access;

    void increment() { p_index += 2 * D; }

    bool equal(Iterator const& other) const { return this->p_index == other.p_index; }

    reference dereference() const { return reference(p_index); }

    void decrement() { p_index -= 2 * D; }

    void advance(int n) { p_index += 2 * D * n; }

    difference_type distance_to(Iterator const& z) const { return ((z.p_index - p_index) / 2 / D); }

  private:
    Integer* p_index;
    size_t D;
  };

  template<typename Integer>
  class Iterator_const : public boost::iterator_facade<Iterator_const<Integer>,
                                                       Integer const*,
                                                       std::random_access_iterator_tag,
                                                       Integer const*,
                                                       std::ptrdiff_t>
  {
  public:
    using difference_type = std::ptrdiff_t;
    using reference       = Integer const*;
    using const_reference = Integer const*;
    using value_tupe      = Integer*;

    Iterator_const(Integer* index, size_t d_) : p_index(index), D(d_) {}
    Iterator_const(Integer const* index, size_t d_) : p_index(index), D(d_) {}

    // What we implement is determined by the boost::forward_traversal_tag
  private:
    friend class boost::iterator_core_access;

    void increment() { p_index += 2 * D; }

    bool equal(Iterator_const const& other) const { return this->p_index == other.p_index; }

    reference dereference() const { return reference(p_index); }

    void decrement() { p_index -= 2 * D; }

    void advance(int n) { p_index += 2 * D * n; }

    difference_type distance_to(Iterator_const const& z) const { return ((z.p_index - p_index) / 2 / D); }

  private:
    Integer* p_index;
    size_t D;
  };

public:
  using Excitation_Iterator       = Iterator<integer_type>;
  using Excitation_const_Iterator = Iterator_const<integer_type>;

  ph_excitations() = delete;

  // Note: terms_per_excitation[0] has special meaning, the number of electrons in the calculation.
  // coefficients[0] will store the reference configuration itself.
  ph_excitations(size_t number_of_configurations,
                 int na_,
                 int nb_,
                 std::vector<size_t>& unique_alpha_counts,
                 std::vector<size_t>& unique_beta_counts,
                 Alloc alloc_ = Alloc{})
      : i_allocator_(alloc_),
        c_allocator_(alloc_),
        NAEA(na_),
        NAEB(nb_),
        configurations(1, number_of_configurations, c_allocator_),
        reference(1, NAEA + NAEB, i_allocator_),
        unique_alpha(unique_alpha_counts.size(), unique_alpha_counts, i_allocator_),
        unique_beta(unique_beta_counts.size(), unique_beta_counts, i_allocator_)
  {
    size_t emax = std::max(unique_alpha_counts.size(), unique_beta_counts.size());
    sum_of_exct.resize(emax + 1);
    sum_of_exct[0] = {0, 0};
    sum_of_exct[1] = {1, 1};
    for (size_t n = 1; n < unique_alpha.size(); ++n)
      sum_of_exct[n + 1][0] = sum_of_exct[n][0] + unique_alpha_counts[n] / size_t(2) / n;
    for (size_t n = unique_alpha.size() + 1; n <= emax; n++)
      sum_of_exct[n][0] = sum_of_exct[n - 1][0];
    for (size_t n = 1; n < unique_beta.size(); ++n)
      sum_of_exct[n + 1][1] = sum_of_exct[n][1] + unique_beta_counts[n] / size_t(2) / n;
    for (size_t n = unique_beta.size() + 1; n <= emax; n++)
      sum_of_exct[n][1] = sum_of_exct[n - 1][1];
  }

  ph_excitations(ph_excitations const& other) = delete;
  //ph_excitations(ph_excitations && other) = default;
  ph_excitations(ph_excitations&& other)
      : i_allocator_(other.i_allocator_),
        c_allocator_(other.c_allocator_),
        NAEA(other.NAEA),
        NAEB(other.NAEB),
        configurations(std::move(other.configurations)),
        reference(std::move(other.reference)),
        unique_alpha(std::move(other.unique_alpha)),
        unique_beta(std::move(other.unique_beta)),
        sum_of_exct(std::move(other.sum_of_exct))
  {}

  ph_excitations& operator=(ph_excitations const& other) = delete;
  ph_excitations& operator=(ph_excitations&& other) = default;

  std::array<size_t, 2> maximum_excitation_number() const { return {unique_alpha.size(), unique_beta.size()}; }
  size_t number_of_unique_alpha_excitations(int n) const
  {
    if (n == 0)
      return 1;
    return unique_alpha.num_elements(n) / 2 / n;
  }
  size_t number_of_unique_beta_excitations(int n) const
  {
    if (n == 0)
      return 1;
    return unique_beta.num_elements(n) / 2 / n;
  }
  std::array<size_t, 2> number_of_unique_excitations(int n) const
  {
    if (n == 0)
      return {1, 1};
    std::array<size_t, 2> res{0, 0};
    if (n < unique_alpha.size())
      res[0] = unique_alpha.num_elements(n) / 2 / n;
    if (n < unique_beta.size())
      res[1] = unique_beta.num_elements(n) / 2 / n;
    return res;
  }
  std::array<size_t, 2> number_of_unique_excitations() const { return sum_of_exct.back(); }

  size_t number_of_configurations() const { return configurations.num_elements(0); }

  // returns the number of unique excitations with particle number less than n
  std::array<size_t, 2> number_of_unique_smaller_than(int n) const { return sum_of_exct[n]; }

  template<class intIt>
  void add_alpha(size_t n, intIt indx)
  {
    assert(n < unique_alpha.size());
    assert(n > 0);
    for (int i = 0; i < 2 * n; i++, ++indx)
      unique_alpha.emplace_back(n, static_cast<integer_type>(*indx));
  }

  template<class intIt>
  void add_beta(size_t n, intIt indx)
  {
    assert(n < unique_beta.size());
    assert(n > 0);
    for (int i = 0; i < 2 * n; i++, ++indx)
      unique_beta.emplace_back(n, static_cast<integer_type>(*indx));
  }

  template<class Vector>
  void add_reference(Vector& refa, Vector& refb)
  {
    for (auto k : refa)
      reference.emplace_back(0, static_cast<integer_type>(k));
    for (auto k : refb)
      reference.emplace_back(0, static_cast<integer_type>(k));
  }

  // index=0 is reserved for the reference!!!
  template<typename integer, typename value>
  void add_configuration(integer alpha_indx, integer beta_index, value ci)
  {
    configurations.emplace_back(0, configuration_type{alpha_indx, beta_index, ci});
  }

  typename Excitation_Iterator::const_reference reference_configuration(int spin = 0) const
  {
    return to_address(reference.values(0)) + (spin == 0 ? 0 : NAEA);
  }

  typename Excitation_Iterator::reference reference_configuration(int spin = 0)
  {
    return to_address(reference.values(0)) + (spin == 0 ? 0 : NAEA);
  }

  configuration_type const* configurations_begin() const { return to_address(configurations.values(0)); }

  configuration_type const* configurations_end() const
  {
    return to_address(configurations.values()) + (*configurations.pointers_end(0));
  }

  configuration_type const* configuration(int i) const { return to_address(configurations.values()) + i; }

  Excitation_Iterator alpha_begin(int n)
  {
    assert(n > 0);
    if (n < unique_alpha.size())
    {
      return Excitation_Iterator(to_address(unique_alpha.values(n)), n);
    }
    else
      return alpha_end(n);
  }

  Excitation_Iterator alpha_end(int n)
  {
    assert(n > 0);
    if (n < unique_alpha.size())
      return Excitation_Iterator(to_address(unique_alpha.values()) + (*unique_alpha.pointers_end(n)), n);
    else
      return Excitation_Iterator(to_address(unique_alpha.values()) +
                                     (*unique_alpha.pointers_end(unique_alpha.size() - 1)),
                                 1);
  }

  Excitation_const_Iterator alpha_begin(int n) const
  {
    assert(n > 0);
    if (n < unique_alpha.size())
    {
      return Excitation_const_Iterator(to_address(unique_alpha.values(n)), n);
    }
    else
      return alpha_end(n);
  }

  Excitation_const_Iterator alpha_end(int n) const
  {
    assert(n > 0);
    if (n < unique_alpha.size())
      return Excitation_const_Iterator(to_address(unique_alpha.values()) + (*unique_alpha.pointers_end(n)), n);
    else
      return Excitation_const_Iterator(to_address(unique_alpha.values()) +
                                           (*unique_alpha.pointers_end(unique_alpha.size() - 1)),
                                       1);
  }

  Excitation_Iterator beta_begin(int n)
  {
    assert(n > 0);
    if (n < unique_beta.size())
      return Excitation_Iterator(to_address(unique_beta.values(n)), n);
    else
      return beta_end(n);
  }

  Excitation_Iterator beta_end(int n)
  {
    assert(n > 0);
    if (n < unique_beta.size())
      return Excitation_Iterator(to_address(unique_beta.values()) + (*unique_beta.pointers_end(n)), n);
    else
      return Excitation_Iterator(to_address(unique_beta.values()) + (*unique_beta.pointers_end(unique_beta.size() - 1)),
                                 1);
  }

  Excitation_const_Iterator beta_begin(int n) const
  {
    assert(n > 0);
    if (n < unique_beta.size())
      return Excitation_const_Iterator(to_address(unique_beta.values(n)), n);
    else
      return beta_end(n);
  }

  Excitation_const_Iterator beta_end(int n) const
  {
    assert(n > 0);
    if (n < unique_beta.size())
      return Excitation_const_Iterator(to_address(unique_beta.values()) + (*unique_beta.pointers_end(n)), n);
    else
      return Excitation_const_Iterator(to_address(unique_beta.values()) +
                                           (*unique_beta.pointers_end(unique_beta.size() - 1)),
                                       1);
  }

  // for generic access
  std::array<Excitation_Iterator, 2> unique_begin(int n)
  {
    return std::array<Excitation_Iterator, 2>{alpha_begin(n), beta_begin(n)};
  }
  std::array<Excitation_Iterator, 2> unique_end(int n)
  {
    return std::array<Excitation_Iterator, 2>{alpha_end(n), beta_end(n)};
  }
  std::array<Excitation_const_Iterator, 2> unique_begin(int n) const
  {
    return std::array<Excitation_const_Iterator, 2>{alpha_begin(n), beta_begin(n)};
  }
  std::array<Excitation_const_Iterator, 2> unique_end(int n) const
  {
    return std::array<Excitation_const_Iterator, 2>{alpha_end(n), beta_end(n)};
  }

  template<class Vector>
  void get_alpha_configuration(size_t index, Vector& confg) const
  {
    assert(confg.size() >= NAEA);
    std::copy_n(to_address(reference.values(0)), NAEA, confg.data());
    if (index == 0)
      return;
    // could use lower bound
    for (int i = 1; i < unique_alpha.size(); i++)
    {
      if (index >= sum_of_exct[i][0] && index < sum_of_exct[i + 1][0])
      {
        size_t dn = index - sum_of_exct[i][0];
        auto exct = unique_alpha.values(i) + 2 * i * dn;
        for (int n = 0; n < i; n++)
          confg[exct[n]] = exct[n + i];
        return;
      }
    }
    APP_ABORT(" Error in ph_excitations::get_alpha_configuration() \n");
  }

  template<class Vector>
  void get_beta_configuration(size_t index, Vector& confg) const
  {
    assert(confg.size() >= NAEB);
    std::copy_n(to_address(reference.values(0)) + NAEA, NAEB, confg.data());
    if (index == 0)
      return;
    // could use lower bound
    for (int i = 1; i < unique_beta.size(); i++)
    {
      if (index >= sum_of_exct[i][1] && index < sum_of_exct[i + 1][1])
      {
        size_t dn = index - sum_of_exct[i][1];
        auto exct = unique_beta.values(i) + 2 * i * dn;
        for (int n = 0; n < i; n++)
          confg[exct[n]] = exct[n + i];
        return;
      }
    }
    APP_ABORT(" Error in ph_excitations::get_beta_configuration() \n");
  }

  template<class Vector>
  void get_configuration(int spin, size_t index, Vector& confg) const
  {
    if (spin == 0)
      get_alpha_configuration(index, confg);
    else
      get_beta_configuration(index, confg);
  }

private:
  // using array_of_seq until I switch to Boost.Multi to be able to use shared_allocator
  int NAEA, NAEB;
  confg_aos configurations;
  index_aos reference;
  index_aos unique_alpha;
  index_aos unique_beta;
  std::vector<std::array<size_t, 2>> sum_of_exct;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
