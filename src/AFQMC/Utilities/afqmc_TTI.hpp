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

namespace qmcplusplus
{
namespace afqmc
{
// checks if class has a member function called reserve that accepts a vector of size_t
template<class T, typename = decltype(std::declval<T>().reserve(std::vector<std::size_t>{}))>
std::true_type has_reserve_with_vector_aux(T);
std::false_type has_reserve_with_vector_aux(...);
template<class V>
struct has_reserve_with_vector : decltype(has_reserve_with_vector_aux(std::declval<V>()))
{};

// reserve with either vector or size_t
template<class Container, typename integer, typename = std::enable_if_t<has_reserve_with_vector<Container>{}>>
void reserve_to_fit(Container& C, std::vector<integer> const& v)
{
  C.reserve(v);
}

template<class Container, typename integer, typename = std::enable_if_t<not has_reserve_with_vector<Container>{}>>
void reserve_to_fit(Container& C, std::vector<integer> const& v, double = 0)
{
  C.reserve(std::accumulate(v.begin(), v.end(), std::size_t(0)));
}


// checks for emplace_back(tuple<int,int,SPComplexType>)
using tp = std::tuple<int, int, SPComplexType>;
template<class T, typename = decltype(std::declval<T>().emplace_back(tp{}))>
std::true_type has_emplace_back_tp_aux(T);
std::false_type has_emplace_back_tp_aux(...);
template<class V>
struct has_emplace_back_tp : decltype(has_emplace_back_tp_aux(std::declval<V>()))
{};

// checks for emplace(tuple<int,int,SPComplexType>)
template<class T, typename = decltype(std::declval<T>().emplace(tp{}))>
std::true_type has_emplace_tp_aux(T);
std::false_type has_emplace_tp_aux(...);
template<class V>
struct has_emplace_tp : decltype(has_emplace_tp_aux(std::declval<V>()))
{};

// dispatch to emplace_back preferentially
template<class Container, class tp_, typename = std::enable_if_t<has_emplace_back_tp<Container>{}>>
void emplace(Container& C, tp_ const& a)
{
  C.emplace_back(a);
}

// dispatch to emplace if exists (and emplace_back doesn't)
template<class Container,
         class tp_,
         typename = std::enable_if_t<not has_emplace_back_tp<Container>{}>,
         typename = std::enable_if_t<has_emplace_tp<Container>{}>>
void emplace(Container& C, tp_ const& a)
{
  C.emplace(a);
}

/*
// checks for emplace_back(tuple<int,int,SPComplexType>)  
template<typename T> using tp = std::tuple<int,int,T>;
template<class Q, class T, typename = decltype(std::declval<T>().emplace_back(tp<Q>{}))>
std::true_type  has_emplace_back_tp_aux(T);
template<class Q>
std::false_type has_emplace_back_tp_aux(...);
template<class V, class Q> struct has_emplace_back_tp : decltype(has_emplace_back_tp_aux<Q>(std::declval<V>())){};

// checks for emplace(tuple<int,int,SPComplexType>)  
template<class Q, class T, typename = decltype(std::declval<T>().emplace(tp<Q>{}))>
std::true_type  has_emplace_tp_aux(T);
template<class Q>
std::false_type has_emplace_tp_aux(...);
template<class V, class Q> struct has_emplace_tp : decltype(has_emplace_tp_aux<Q>(std::declval<V>())){};

// dispatch to emplace_back preferentially
template<
    class Container,
    typename T,
    typename = std::enable_if_t<has_emplace_back_tp<Container,T>{}>
>
void emplace(Container& C, tp<T> const& a){
    C.emplace_back(a);
}

// dispatch to emplace if exists (and emplace_back doesn't) 
template<
    class Container,
    typename T,
    typename = std::enable_if_t<not has_emplace_back_tp<Container,T>{}>,
    typename = std::enable_if_t<has_emplace_tp<Container,T>{}>
>
void emplace(Container& C, tp<T> const& a){
    C.emplace(a);
}
*/
} // namespace afqmc

} // namespace qmcplusplus

#endif
