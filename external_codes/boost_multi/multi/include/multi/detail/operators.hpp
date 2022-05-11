// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MULTI_DETAIL_OPERATORS_HPP
#define BOOST_MULTI_DETAIL_OPERATORS_HPP

#include<type_traits>  // enable_if
#include<utility>  // forward

namespace boost {
namespace multi {

struct empty_base {};

template<class T, class V, class B = empty_base> struct equality_comparable2;

template <class T, class B>
struct equality_comparable2<T, void, B> : B {
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}> >
	friend constexpr auto operator!=(const T& y, const U& x) -> bool {return not(y==x);}
};

template<class T, class V> struct partially_ordered2;

template <class T>
struct partially_ordered2<T, void>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator>(const U& x, const T& y) -> bool{return y < x;}
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator<(const U& x, const T& y) -> bool{return y > x;}

	template<class U>
	friend constexpr auto operator<=(T&& x, U&& y) -> bool{return (std::forward<T>(x) < std::forward<T>(y)) or (std::forward<T>(x) == std::forward<T>(y));}
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator<=(const U& x, const T& y) -> bool{return (y > x) or (y == x);}
	template<class U>
	friend constexpr auto operator>=(const T& x, const U& y) -> bool{return (x > y) or (x == y);}
	template<class U, typename = std ::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator>=(const U& x, const T& y) -> bool{return (y < x) or (y == x);}
};

template<class T, class V, class B = empty_base> struct totally_ordered2;

template<class T, class B>
struct totally_ordered2<T, void, B> : B {
	template<class U>
	friend constexpr auto operator<=(const T& x, const U& y) {return (x < y) or (x == y);}
	template<class U>
	friend constexpr auto operator>=(const T& x, const U& y) {return (y < x) or (x == y);}
	template<class U>
	friend constexpr auto operator> (const T& x, const U& y) {return y < x;}
};

template<class T>
struct copy_constructible{};

template<class T>
struct weakly_incrementable{
//  friend T& operator++(weakly_incrementable& t){return ++static_cast<T&>(t);}
};

template<class T>
struct weakly_decrementable{
//  friend T& operator--(weakly_decrementable& t){return --static_cast<T&>(t);}
};

template<class T>
struct incrementable : weakly_incrementable<T>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator++(U& self, int) -> T{T tmp{self}; ++self; return tmp;}
};

template<class T>
struct decrementable : weakly_decrementable<T>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator--(U& self, int) -> T{T tmp{self}; --self; return tmp;}
};

template<class T>
struct steppable : incrementable<T>, decrementable<T>{};

template<class T, class Reference>
struct dereferenceable{
	using reference = Reference;
	friend constexpr auto operator*(dereferenceable const& t) -> reference{return *static_cast<T const&>(t);}
};

template<class T, class D>
struct addable2{
	using difference_type = D;
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr auto operator+(TT&& t, difference_type const& d) -> T{T tmp{std::forward<TT>(t)}; tmp+=d; return tmp;}
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr auto operator+(difference_type const& d, TT&& t) -> T {return std::forward<TT>(t) + d;}
};

template<class T, class D>
struct subtractable2{
	using difference_type = D;
	template<class TT, class = T>
	friend auto operator-(TT&& t, difference_type const& d) -> T {T tmp{std::forward<TT>(t)}; tmp-=d; return tmp;}
};

template<class T, class Difference>
struct affine : addable2<T, Difference>, subtractable2<T, Difference>{
	using difference_type = Difference;
};

template<class T>
struct random_iterable {
	       constexpr auto rbegin()&    {return typename T::reverse_iterator{static_cast<T&>(*this).end  ()};}
	       constexpr auto rend  ()&    {return typename T::reverse_iterator{static_cast<T&>(*this).begin()};}
	friend /*consx*/ auto rbegin(T& s) {return static_cast<random_iterable&>(s).rbegin();}
	friend /*consx*/ auto rend  (T& s) {return static_cast<random_iterable&>(s).rend  ();}

	       constexpr auto cfront() const&    -> decltype(auto) {return static_cast<T const&>(*this).front();}
	       constexpr auto cback () const&    -> decltype(auto) {return static_cast<T const&>(*this).back() ;}
	friend constexpr auto cfront(T const& s) -> decltype(auto) {return s.cfront();}
	friend constexpr auto cback (T const& s) -> decltype(auto) {return s.cback() ;}
};

}  // end namespace multi
}  // end namespace boost

#endif

