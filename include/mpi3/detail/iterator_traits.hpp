//  -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef MPI3_DETAIL_ITERATOR_TRAITS_HPP
#define MPI3_DETAIL_ITERATOR_TRAITS_HPP
#pragma once

#include "./iterator.hpp"

#include "./just.hpp"

#include<iterator>
#include<type_traits>

namespace boost {
namespace mpi3 {
namespace detail {

template<class Iterator>
struct iterator_traits : std::iterator_traits<Iterator> {};

struct unspecified{};
struct output_iterator_tag {using base = unspecified;};
struct input_iterator_tag {using base = unspecified;};
struct forward_iterator_tag : input_iterator_tag {using base = input_iterator_tag;};
struct random_access_iterator_tag : forward_iterator_tag {
	using base = forward_iterator_tag;
};
struct contiguous_iterator_tag : random_access_iterator_tag{
	using base = random_access_iterator_tag;
};
struct strided_contiguous_iterator_tag : std::random_access_iterator_tag{};

template<class T> struct std_translate;

template<> struct std_translate<std::output_iterator_tag>{using type = output_iterator_tag;};
template<> struct std_translate<std::input_iterator_tag>{using type = input_iterator_tag;};
template<> struct std_translate<std::forward_iterator_tag>{using type = forward_iterator_tag;};
template<> struct std_translate<std::bidirectional_iterator_tag>{using type = forward_iterator_tag;};
template<> struct std_translate<std::random_access_iterator_tag>{using type = random_access_iterator_tag;};

//template<class T> struct is_declared_contiguous : std::false_type{};
template<class T> struct is_contiguous;

template<class It, class = decltype(detail::data(std::declval<It>()))>
std::true_type is_contiguous_aux(It);
std::false_type is_contiguous_aux(...);

template<class T> struct is_contiguous 
: decltype(is_contiguous_aux(std::declval<T>())){};

template<class T, typename = decltype(detail::data(T{}))> 
std::true_type  has_data_aux(T  );
std::false_type has_data_aux(...);

template<class T> struct has_data : decltype(has_data_aux(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class It, typename = std::enable_if_t<not has_data<It>::value> >
typename std_translate<typename std::iterator_traits<It>::iterator_category>::type iterator_category_aux(It);

//template<class It, typename = std::enable_if_t<std::is_same<typename std::iterator_traits<It>::iterator_category, std::output_iterator_tag>{}>>
//std_translate<std::output_iterator_tag>::type iterator_category_aux(It);

// intel compiler 17 needs this specialization
template<class T>
contiguous_iterator_tag iterator_category_aux(T*);
 
template<class Iter, typename = std::enable_if_t<has_data<Iter>{}>>
contiguous_iterator_tag iterator_category_aux(Iter);
template<class Iter, class = decltype(data(Iter{}.base())), class = decltype(Iter{}.stride()), class = std::enable_if_t<not has_data<Iter>{}>>
strided_contiguous_iterator_tag iterator_category_aux(Iter);

template<class Iter>
struct iterator_category {
	using type = decltype(iterator_category_aux(std::declval<std::decay_t<Iter>>()));
};
template<class Iter>
using iterator_category_t = typename iterator_category<Iter>::type;

template<class Iter> struct forward_iterator       : just_t<Iter> {
	using just_t<Iter>::just_t;
};
template<class Iter> struct random_access_iterator : forward_iterator<Iter> {
	using forward_iterator<Iter>::forward_iterator;
};

template<class Iter> struct contiguous_iterator    : random_access_iterator<Iter> {
	using random_access_iterator<Iter>::random_access_iterator;
};

template<class Iter>
forward_iterator<Iter&&> category_iterator_aux(Iter&& it, forward_iterator_tag /*forward*/) {
	return {std::forward<Iter>(it)};
}

template<class Iter>
random_access_iterator<Iter&&> category_iterator_aux(Iter&& it, random_access_iterator_tag /*random_access*/) {
	return {std::forward<Iter>(it)};
}

template<class Iter>
contiguous_iterator<Iter&&> category_iterator_aux(Iter&& it, contiguous_iterator_tag /*contiguous*/) {
	return {std::forward<Iter>(it)};
}

template<class Iter>
auto category_iterator(Iter&& it){
	return category_iterator_aux(std::forward<Iter>(it), typename iterator_category<Iter>::type{});
}

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_DETAIL_ITERATOR_TRAITS

//#include "../../mpi3/main.hpp"
//#include "../../mpi3/vector.hpp"

//#include<deque>
//#include<list>
//#include<vector>

//#include<iostream>

//namespace mpi3 = boost::mpi3;
//using std::cout;

//template<class It>
//std::string f(It it, mpi3::detail::forward_iterator_tag const&){
//	return "forward" + std::to_string(*it);
//};

//template<class It>
//std::string f(It it, mpi3::detail::random_access_iterator_tag const&){
//	return "random" + std::to_string(*it);
//};

//template<class It>
//std::string f(It it, mpi3::detail::contiguous_iterator_tag const&){
//	return "cont" + std::to_string(*it);
//};

//template<class It>
//std::string f(It&& it){
//	return f(
//		std::forward<It>(it),
//		typename boost::mpi3::detail::iterator_category<It>::type{}
//	);
//}

//#if 1
//namespace dispatch{
//template<class It>
//std::string g(mpi3::detail::forward_iterator<It>&& it){
//	return "forward" + std::to_string(**&it);
//}
//template<class It>
//std::string g(mpi3::detail::random_access_iterator<It>&& it){
//	return "random" + std::to_string(**&it);
//}
//template<class It>
//std::string g(mpi3::detail::contiguous_iterator<It>&& it){
//	return "contiguous" + std::to_string(**&it);
//}
//}

//template<class It>
//std::string g(It&& it){
//	return dispatch::g(mpi3::detail::category_iterator(std::forward<It>(it)));
//}
//#endif

//int mpi3::main(int, char*[], mpi3::communicator){
//	{
//		std::list<int> l = {11, 22};
//		assert( f(l.begin()) == "forward11" );
//		std::vector<int> v = {44, 33};
//		assert( f(v.begin()) == "cont44" );
//		std::deque<int> q{55,66};// q.push(5.);//,6.};
//		assert( f(q.begin()) == "random55" );
//		assert( g(l.begin()) == "forward11" );

//		std::istringstream iss("1 2 3");
//		std::istream_iterator<int> it(iss);
//		assert( *std::addressof(*it) == 1 );
//		++it;
//		assert( *std::addressof(*it) == 2 );		
//		cout << typeid(detail::iterator_category_t<std::vector<double>::iterator>).name() <<'\n';
//	}
//	return 0;
//	{
//		mpi3::uvector<int> v = {444, 333};
//	//	assert( f(v.begin()) == "cont444" );
//		cout << *detail::data(v.begin()) << '\n';
//	//	static_assert( detail::has_data<mpi3::uvector<int>::iterator>{} );
//	}
//	return 0;
//}
//#endif
#endif

