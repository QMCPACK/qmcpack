#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_POINTER_TRAITS_HPP
#define MULTI_POINTER_TRAITS_HPP

#include<memory>
#include<type_traits> // std::conditional_t

namespace boost{
namespace multi{

template<std::size_t I> struct Priority : std::conditional_t<I==0, std::true_type, struct Priority<I-1>>{}; 

template<class Pointer> std::allocator<typename std::iterator_traits<Pointer>::value_type> dat_aux(Priority<0>, Pointer);
template<class T> std::allocator<typename std::iterator_traits<T*>::value_type> dat_aux(Priority<1>, T*);
template<class FancyPtr> typename FancyPtr::default_allocator_type dat_aux(Priority<2>, FancyPtr);

template<class Pointer>
struct pointer_traits/*, typename Pointer::default_allocator_type>*/ : std::pointer_traits<Pointer>{
	using default_allocator_type = decltype(dat_aux(Priority<2>{}, std::declval<Pointer>()));
};

}}

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi memory pointer traits"
#include<boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_pointer_traits){
	static_assert(std::is_same<multi::pointer_traits<double*>::default_allocator_type, std::allocator<double>>{}, "!");
}

#endif
#endif

