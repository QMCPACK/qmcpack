#if COMPILATION_INSTRUCTIONS
(echo "#include<" $0 ">" > $0x.cpp) && time clang++ - O3 - std = c++ 1z - Wall `# - Wfatal
        - errors` - I..- D_TEST_SPARSE_COO_MATRIX $0x.cpp - lstdc++ fs - lboost_system - lboost_timer - o $0x.x
    && time $0x.x $ @ && rm - f $0x.cpp;
exit
#endif

////////////////////////////////////////////////////////////////////////////////
// File developed by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef SPARSE_COO_MATRIX_HPP
#define SPARSE_COO_MATRIX_HPP

#include <array>
#include <cassert>
#include <cstddef> // ptrdiff_t
#include <vector>
#include <tuple>

    namespace ma
{
  namespace sparse
  {
  using size_type = std::size_t;

  template<class T, class index, class Alloc = std::allocator<T>>
  class coo_matrix
  {
  protected:
    using alloc_ts        = std::allocator_traits<Alloc>;
    using index_allocator = typename alloc_ts::template rebind_alloc<index>;
    std::vector<index, index_allocator> is_;
    std::vector<index, index_allocator> js_;
    std::vector<T, Alloc> vs_;
    size_type cols_;
    size_type rows_;

  public:
    using element                            = T;
    coo_matrix(coo_matrix const&)            = default;
    coo_matrix(coo_matrix&&)                 = default;
    bool operator==(coo_matrix const&) const = delete;
    coo_matrix(std::tuple<index, index> const& arr                                    = std::tuple<index, index>{0, 0},
               std::initializer_list<std::pair<std::tuple<index, index>, element>> il = {})
        : cols_(std::get<0>(arr)), rows_(std::get<1>(arr))
    {
      for (auto e : il)
        emplace(e.first, e.second);
    }
    void reserve(size_type s)
    {
      is_.reserve(s);
      js_.reserve(s);
      vs_.reserve(s);
    }
    auto size() const { return rows_; }
    auto num_elements() const { return size() * cols_; }
    std::array<index, 2> shape() const { return {{size(), cols_}}; }
    auto num_non_zero_elements() const { return vs_.size(); }
    T* non_zero_values_data() { return vs_.data(); }
    index* non_zero_indices1_data() { return (index*)is_.data(); }
    index* non_zero_indices2_data() const { return (index*)js_.data(); }
    decltype(auto) move_non_zero_values() && { return std::move(vs_); }
    decltype(auto) move_non_zero_indices1() && { return std::move(is_); }
    decltype(auto) move_non_zero_indices2() && { return std::move(js_); }
    void clear()
    {
      is_.clear();
      js_.clear();
      vs_.clear();
      cols_ = 0;
      rows_ = 0;
    }
    template<class Pair = std::array<index, 2>, class TT>
    void emplace(Pair&& indices, TT&& tt)
    {
      using std::get;
      is_.emplace_back(get<0>(std::forward<Pair>(indices)));
      js_.emplace_back(get<1>(std::forward<Pair>(indices)));
      vs_.emplace_back(std::forward<TT>(tt));
    }

  protected:
    struct row_reference
    {
      coo_matrix& self_;
      index i_;
      struct element_reference
      {
        row_reference& self_;
        index j_;
        template<class TT>
        element_reference&& operator=(TT&& tt) &&
        {
          self_.self_.emplace({{self_.i_, j_}}, std::forward<TT>(tt));
          return std::move(*this);
        }
      };
      using reference = element_reference;
      reference operator[](index i) && { return reference{*this, i}; }
    };

  public:
    using reference = row_reference;
    reference operator[](index i) { return reference{*this, i}; }
    friend decltype(auto) size(coo_matrix const& s) { return s.size(); }
    friend decltype(auto) shape(coo_matrix const& s) { return s.shape(); }
    friend decltype(auto) clear(coo_matrix& s) { s.clear(); }
  };
  /*
template<class... Ts>
std::array<index, 2> index_bases(coo_matrix<Ts...> const&){return {{0,0}};}
template<class... Ts>
auto num_non_zero_elements(coo_matrix<Ts...> const& s){
	return s.num_non_zero_elements();
}
template<class... Ts>
auto non_zero_values_data(coo_matrix<Ts...>& s){return s.non_zero_values_data();}
template<class... Ts>
auto non_zero_indices1_data(coo_matrix<Ts...>& s){
	return s.non_zero_indices1_data();
}
template<class... Ts>
auto non_zero_indices2_data(coo_matrix<Ts...>& s){
	return s.non_zero_indices2_data();
}
*/
  } // namespace sparse
}

#ifdef _TEST_SPARSE_COO_MATRIX

#include "iterator/zipper.hpp"
#include "timer/timed.hpp"

#include <boost/timer/timer.hpp>

#include <algorithm> // std::sort
#include <cassert>
#include <iostream>
#include <random>

using std::cerr;
using std::cout;
using std::get;

int main()
{
  using ma::sparse::coo_matrix;

  auto const M          = 40000;
  auto const N          = 40000;
  double const sparsity = 0.01;

  // generate tuples for reference
  auto const csource = [&]() {
    std::vector<std::tuple<int, int, double>> source;
    source.reserve(M * N * sparsity);
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist;
    std::uniform_real_distribution<double> dist2(0, 10);
    for (auto i = M; i != 0; --i)
      for (auto j = N; j != 0; --j)
        if (dist(gen) < sparsity)
          source.emplace_back(i - 1, j - 1, dist2(gen));
    return source;
  }();

  cerr << "value size " << csource.size() * sizeof(double) / 1000000 << " MB\n";
  {
    coo_matrix<double> coom({M, N});
    {
      boost::timer::auto_cpu_timer t("ugly syntax: %t seconds\n");
      for (auto& s : csource)
        coom.emplace({{get<0>(s), get<1>(s)}}, get<2>(s));
    }
  }
  {
    coo_matrix<double> coom({M, N});
    {
      boost::timer::auto_cpu_timer t("generic syntax: %t seconds\n");
      for (auto& s : csource)
        coom[get<0>(s)][get<1>(s)] = get<2>(s);
    }
  }

  coo_matrix<double> small({4, 4}, {{{3, 3}, 1}, {{2, 1}, 3}, {{0, 1}, 9}});
}

#endif
#endif
