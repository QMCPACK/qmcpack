#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MULTI_DETAIL_OPERATORS $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_OPERATORS_HPP
#define BOOST_MULTI_DETAIL_OPERATORS_HPP

//#include "layout.hpp"
#include<type_traits>

namespace boost{
namespace multi{

struct empty_base{};

template<class T, class V, class B = empty_base> struct equality_comparable2;

template <class T, class B>
struct equality_comparable2<T, void, B> : B{
//	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
//	friend bool operator==(const U& y, const T& x){return x == y;}
//	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
//	friend bool operator!=(const U& y, const T& x){return not (x == y);}
	template<class U>
	friend bool operator!=(const T& y, const U& x){return not (y == x);}
};

template<class T, class V, class B = empty_base> struct partially_ordered2;

template <class T, class B>
struct partially_ordered2<T, void, B> : B{
	template<class U>
	friend bool operator<=(const T& x, const U& y){return (x < y) or (x == y);}
	template<class U>
	friend bool operator>=(const T& x, const U& y){return (x > y) or (x == y);}
	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
	friend bool operator>(const U& x, const T& y){return y < x;}
	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
	friend bool operator<(const U& x, const T& y){return y > x;}
	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
	friend bool operator<=(const U& x, const T& y){return (y > x) or (y == x);}
	template<class U, typename = std::enable_if_t<not std::is_same<U, T>{}> >
	friend bool operator>=(const U& x, const T& y){return (y < x) or (y == x);}
};

template<class T, class V, class B = empty_base> struct totally_ordered2;

template<class T, class B>
struct totally_ordered2<T, void, B> : B{
	template<class U>
	friend auto operator<=(const T& x, const U& y){return (x < y) or (x == y);}
	template<class U>
	friend auto operator>=(const T& x, const U& y){return (y < x) or (x == y);}
	template<class U>
	friend auto operator>(const T& x, const U& y){return y < x;}
};

template<class T, class B = empty_base>
struct random_iterable : B{
//	using iterator = Iterator;
///	template<typename U = T>
//	typename U::const_iterator
//	 cbegin() const{return typename T::const_iterator{static_cast<T const&>(*this).begin()};}
//	template<typename U = T>
//	typename U::const_iterator
//	 cend() const{return typename T::const_iterator{static_cast<T const&>(*this).end()};}
//	template<typename U = T>
//	friend typename U::const_iterator cbegin(U const& s, T* = 0, ...){
//		return static_cast<random_iterable const&>(s).cbegin();}
//	template<typename U = T>
//	friend typename U::const_iterator cend  (U const& s, T* = 0, ...){
//		return static_cast<random_iterable const&>(s).cend  ();}
//	auto cend()   const{return typename T::const_iterator{static_cast<T const*>(this)->end()};}
//	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
//	friend auto cbegin(TT const& s)->typename TT::const_iterator
//	{return typename TT::const_iterator{static_cast<T const&>(s).begin()};}
//	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
//	friend auto cend  (TT const& s)->typename TT::const_iterator
//	{return typename TT::const_iterator{static_cast<T const&>(s)->end()};}

//	auto crbegin() const{return typename T::const_reverse_iterator{cend  ()};}
//	auto crend  () const{return typename T::const_reverse_iterator{cbegin()};}
//	friend auto crbegin(T const& s){return static_cast<random_iterable const&>(s).cbegin();}
//	friend auto crend  (T const& s){return static_cast<random_iterable const&>(s).cend  ();}

	friend auto begin(T const& t){return t.begin();}
	friend auto end  (T const& t){return t.end();}

	friend auto begin(T& t){return t.begin();}
	friend auto end  (T& t){return t.end();}

	auto rbegin(){return typename T::reverse_iterator{static_cast<T&>(*this).end  ()};}
	auto rend  (){return typename T::reverse_iterator{static_cast<T&>(*this).begin()};}
	friend auto rbegin(T& s){return static_cast<random_iterable&>(s).rbegin();}
	friend auto rend  (T& s){return static_cast<random_iterable&>(s).rend  ();}

	decltype(auto) cfront() const{return static_cast<T const&>(*this).front();}
	decltype(auto) cback()  const{return static_cast<T const&>(*this).back() ;}
	friend auto cfront(T const& s){return s.cfront();}
	friend auto cback (T const& s){return s.cback() ;}
};

//template<class T, class B>
//typename T::const_iterator cbegin(random_iterable<T, B> const& c){return c.cbegin();}
//template<class T, class B>
//typename T::const_iterator cend(random_iterable<T, B> const& c){return c.cend();}

}}

#if _TEST_BOOST_MULTI_DETAIL_OPERATORS

#include<iostream>
using std::cout;
namespace multi = boost::multi;

int main(){}

#endif
#endif

