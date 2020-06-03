#if defined(COMPILATION)
$CXX $0 -o $0x&&$0x&&rm $0x $0.cpp;exit
#endif
#ifndef BOOST_MULTI_DETAIL_OPERATORS_HPP
#define BOOST_MULTI_DETAIL_OPERATORS_HPP

#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD 
#endif

#include<type_traits> // enable_if
#include<utility> // forward

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
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>> 
	friend bool operator!=(const T& y, const U& x){return not(y==x);}
};

template<class T, class V> struct partially_ordered2;

template <class T>
struct partially_ordered2<T, void>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend bool operator>(const U& x, const T& y){return y < x;}
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend bool operator<(const U& x, const T& y){return y > x;}

	template<class U>
	friend bool operator<=(T&& x, U&& y){return (std::forward<T>(x) < std::forward<T>(y)) or (std::forward<T>(x) == std::forward<T>(y));}
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend bool operator<=(const U& x, const T& y){return (y > x) or (y == x);}
	template<class U>
	friend bool operator>=(const T& x, const U& y){return (x > y) or (x == y);}
	template<class U, typename = std ::enable_if_t<not std::is_base_of<T, U>{}>>
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

template<class T>
struct copy_constructible{};

template<class T>
struct weakly_incrementable{
//	friend T& operator++(weakly_incrementable& t){return ++static_cast<T&>(t);}
};

template<class T>
struct weakly_decrementable{
//	friend T& operator--(weakly_decrementable& t){return --static_cast<T&>(t);}
};

template<class T>
struct incrementable : weakly_incrementable<T>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend T operator++(U& self, int){T tmp{self}; ++self; return tmp;}
};

template<class T>
struct decrementable : weakly_decrementable<T>{
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend T operator--(U& self, int){T tmp{self}; --self; return tmp;}
};

template<class T>
struct steppable : incrementable<T>, decrementable<T>{};

template<class T, class Reference>//, typename Reference = decltype(*std::declval<T const&>())>
struct dereferenceable{
	using reference = Reference;
	friend reference operator*(dereferenceable const& t){return *static_cast<T const&>(t);}
};

template<class T, class D>
struct addable2{
	using difference_type = D;
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr T operator+(TT&& t, difference_type const& d){T tmp{std::forward<TT>(t)}; tmp+=d; return tmp;}
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr T operator+(difference_type const& d, TT&& t) HD{return std::forward<TT>(t) + d;}
};

template<class T, class D>
struct subtractable2{
	using difference_type = D;
	template<class TT, class = T>
	friend T operator-(TT&& t, difference_type const& d){T tmp{std::forward<TT>(t)}; tmp-=d; return tmp;}
};

template<class T, class D = typename T::difference_type>
struct affine : addable2<T, D>, subtractable2<T, D>{
	using difference_type = D;
};

template<class T>
struct random_iterable{

//	friend constexpr decltype(auto) begin(T& t){return t.begin();}
//	friend constexpr decltype(auto) end  (T& t){return t.end();}

	auto rbegin(){return typename T::reverse_iterator{static_cast<T&>(*this).end  ()};}
	auto rend  (){return typename T::reverse_iterator{static_cast<T&>(*this).begin()};}
	friend auto rbegin(T& s){return static_cast<random_iterable&>(s).rbegin();}
	friend auto rend  (T& s){return static_cast<random_iterable&>(s).rend  ();}

	decltype(auto) cfront() const{return static_cast<T const&>(*this).front();}
	decltype(auto) cback()  const{return static_cast<T const&>(*this).back() ;}
	friend auto cfront(T const& s){return s.cfront();}
	friend auto cback (T const& s){return s.cback() ;}

};

#if 0
// TODO random_iterable_container ??
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

};
#endif

//template<class T, class B>
//typename T::const_iterator cbegin(random_iterable<T, B> const& c){return c.cbegin();}
//template<class T, class B>
//typename T::const_iterator cend(random_iterable<T, B> const& c){return c.cend();}

}}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MULTI_DETAIL_OPERATORS

#include<iostream>
using std::cout;
namespace multi = boost::multi;

int main(){}

#endif
#endif

