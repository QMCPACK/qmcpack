////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_TTI_HPP 
#define AFQMC_TTI_HPP 

#include "boost/hana.hpp"

namespace qmcplusplus
{
namespace afqmc
{

/*
// checks if clas has a member function called reserve that accepts a vector of size_t
template<class T, typename = decltype(std::declval<T>().shape())>
std::true_type  has_shape_aux(T);
std::false_type has_shape_aux(...);
template<class V> struct has_shape : decltype(has_shape_aux(std::declval<V&>())){};

template<class Container,
    class IArr = std::array<int,2>, 
    typename = std::enable_if_t<has_shape<Container>{}>
>
inline bool check_shape(Container const& C, IArr const& sz) 
{
  if(C.shape()[0] != sz[0] || C.shape()[1] != sz[1])
    return false;
  return true; 
}

template<class Container,
    class IArr = std::array<int,2>, 
    typename = std::enable_if_t<not has_shape<Container>{}>,
    typename = void
>
inline bool check_shape(Container const& C, IArr const& sz)
{
  return true;
}
*/

template<class Container,
         class IArr = std::array<int,2>> 
bool check_shape(Container const& C, IArr const& sz)
{
  auto has_ = boost::hana::is_valid([](auto const& p) -> decltype((void)p.shape()) { });
  return boost::hana::if_(has_(C),
    [] (auto const& C, auto const& sz) {
        if(C.shape()[0] != sz[0] || C.shape()[1] != sz[1])
            return false;
        return true; 
    },
    [] (auto const& C, auto const& a) { return true; } 
  )(C,sz);
}

/*
// checks if clas has a member function called reserve that accepts a vector of size_t
template<class T, typename = decltype(std::declval<T>().reserve(std::vector<std::size_t>{}))>
std::true_type  has_reserve_with_vector_aux(T);
std::false_type has_reserve_with_vector_aux(...);
template<class V> struct has_reserve_with_vector : decltype(has_reserve_with_vector_aux(std::declval<V&>())){};

// reserve with either vector or size_t
template<
    class Container, typename integer,
    typename = std::enable_if_t<has_reserve_with_vector<Container>{}>
>
void reserve_to_fit(Container& C, std::vector<integer> const& v){
    C.reserve(v);
}

template<
    class Container, typename integer,
    typename = std::enable_if_t<not has_reserve_with_vector<Container>{}>
>
void reserve_to_fit(Container& C, std::vector<integer> const& v, double = 0){
    C.reserve(std::accumulate(v.begin(), v.end(), std::size_t(0)));
}
*/

template<class Container,
         class IVec = std::vector<int>>
void reserve_to_fit(Container& C, IVec const& sz)
{
  auto has_reserve = boost::hana::is_valid([](auto& p, auto const& v) -> decltype((void)p.reserve(v)) { });
  boost::hana::if_(has_reserve(C,sz),
    [] (auto& C, auto const& v) {C.reserve(v);},
    [] (auto& C, auto const& v) {
        C.reserve( std::size_t(std::accumulate(v.begin(), v.end(), std::size_t(0)) ) ); }
  )(C,sz);
}

/*
// checks for emplace_back(tuple<int,int,SPComplexType>)  
using tp = std::tuple<int,int,SPComplexType>;
template<class T, typename = decltype(std::declval<T>().emplace_back(tp{}))>
std::true_type  has_emplace_back_tp_aux(T);
std::false_type has_emplace_back_tp_aux(...);
template<class V> struct has_emplace_back_tp : decltype(has_emplace_back_tp_aux(std::declval<V&>())){};

// checks for emplace(tuple<int,int,SPComplexType>)  
template<class T, typename = decltype(std::declval<T>().emplace(tp{}))>
std::true_type  has_emplace_tp_aux(T);
std::false_type has_emplace_tp_aux(...);
template<class V> struct has_emplace_tp : decltype(has_emplace_aux(std::declval<V&>())){};

// dispatch to emplace_back preferentially
template<
    class Container,
    typename = std::enable_if_t<has_emplace_back_tp<Container>{}>
>
void emplace(Container& C, tp const& a){
    C.emplace_back(a);
}

// dispatch to emplace if exists (and emplace_back doesn't) 
template<
    class Container,
    typename = std::enable_if_t<not has_emplace_back_tp<Container>{}>,
    typename = std::enable_if_t<has_emplace_tp<Container>{}>
>
void emplace(Container& C, tp const& a){
    C.emplace(a);
}
*/

template<class Container,
          class tp>
void emplace(Container& C, tp const& a)
{
  auto has_ = boost::hana::is_valid([](auto& p, auto const& a) -> decltype((void)p.emplace_back(a)) { });
  boost::hana::if_(has_(C,a),
    [] (auto& C, auto const& a) {C.emplace_back(a);},
    [] (auto& C, auto const& a) {
        C.emplace(a); }
  )(C,a);
}

}

}

#endif

