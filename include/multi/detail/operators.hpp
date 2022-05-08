// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_OPERATORS_HPP
#define MULTI_DETAIL_OPERATORS_HPP

#include<type_traits>  // enable_if
#include<utility>  // forward

namespace boost::multi {

struct empty_base {};

template<class Self> struct selfable {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto self()       -> self_type      & {return static_cast<self_type      &>(*this);}
	friend constexpr auto self(selfable const& s) -> self_type const& {return s.self();}
};

template<class Self, class U> struct equality_comparable2;

//template<class Self, class Other>
//struct equality_comparable2 : virtual selfable<Self> {
//	friend constexpr auto operator==(equality_comparable2 const& s, Other const& o) {return     s.self() == o ;}
//	friend constexpr auto operator!=(equality_comparable2 const& s, Other const& o) {return not(s.self() == o);}

//	friend constexpr auto operator!=(Other const& o, equality_comparable2 const& s) {return not(s.self() == o);}
//	friend constexpr auto operator==(Other const& o, equality_comparable2 const& s) {return     s.self() == o ;}
//};

template<class Self>
struct equality_comparable2<Self, Self> : selfable<Self> {
//	friend constexpr auto operator==(equality_comparable2 const& s, equality_comparable2 const& o) {return     s.self() == o.self() ;}
	friend constexpr auto operator!=(equality_comparable2 const& s, equality_comparable2 const& o) {return not(s.self() == o.self());}
};

template<class Self> struct equality_comparable : equality_comparable2<Self, Self> {};

//template<class T>
//struct equality_comparable2<T, void> {
//	template<class Other> friend constexpr auto operator==(equality_comparable2 const& s, Other const& o) {return     s.self() == o ;}
//	template<class Other> friend constexpr auto operator!=(equality_comparable2 const& s, Other const& o) {return not(s.self() == o);}

//	template<class Other> friend constexpr auto operator!=(Other const& o, equality_comparable2 const& s) {return not(s.self() == o);}
//	template<class Other> friend constexpr auto operator==(Other const& o, equality_comparable2 const& s) {return     s.self() == o ;}
//};

//template<class T, class V> struct partially_ordered2;

//template <class T>
//struct partially_ordered2<T, void>{
//	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
//	friend constexpr auto operator>(const U& x, const T& y) -> bool {return y < x;}
//	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
//	friend constexpr auto operator<(const U& x, const T& y) -> bool {return y > x;}

//	template<class U>
//	friend constexpr auto operator<=(T&& x, U&& y) -> bool{return (std::forward<T>(x) < std::forward<T>(y)) or (std::forward<T>(x) == std::forward<T>(y));}
//	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
//	friend constexpr auto operator<=(const U& x, const T& y) -> bool {return (y > x) or (y == x);}
//	template<class U>
//	friend constexpr auto operator>=(const T& x, const U& y) -> bool {return (x > y) or (x == y);}
//	template<class U, typename = std ::enable_if_t<not std::is_base_of<T, U>{}>>
//	friend constexpr auto operator>=(const U& x, const T& y) -> bool {return (y < x) or (y == x);}
//};

template<class T, class V> struct totally_ordered2;

template<class Self>
struct totally_ordered2<Self, Self> : equality_comparable2<totally_ordered2<Self, Self>, totally_ordered2<Self, Self>> {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}

//	friend auto operator< (totally_ordered2 const& s, totally_ordered2 const& o) -> bool {return     s.self() < o.self() ;}
	friend auto operator==(totally_ordered2 const& s, totally_ordered2 const& o) -> bool {return not(s.self() < o.self()) and not(o.self() < s.self());}
//	friend auto operator!=(totally_ordered2 const& s, totally_ordered2 const& o) {return    (s.self() < o.self()) or     (o.self() < s.self());}

	friend auto operator<=(totally_ordered2 const& s, totally_ordered2 const& o) -> bool {return not(o.self() < s.self());}

	friend auto operator> (totally_ordered2 const& s, totally_ordered2 const& o) -> bool {return not(s.self() < o.self()) and not(s.self() == o.self());}
	friend auto operator>=(totally_ordered2 const& s, totally_ordered2 const& o) -> bool {return not(s.self() < o.self());}
};

template<class Self> using totally_ordered = totally_ordered2<Self, Self>;

template<class T>
struct totally_ordered2<T, void> {
	template<class U>
	friend constexpr auto operator<=(const T& x, const U& y) {return (x < y) or (x == y);}
	template<class U>
	friend constexpr auto operator>=(const T& x, const U& y) {return (y < x) or (x == y);}
	template<class U>
	friend constexpr auto operator> (const T& x, const U& y) {return y < x;}
};

template<class T>
struct copy_constructible {};

template<class T>
struct weakly_incrementable {
//  friend T& operator++(weakly_incrementable& t){return ++static_cast<T&>(t);}
};

template<class T>
struct weakly_decrementable {
//  friend T& operator--(weakly_decrementable& t){return --static_cast<T&>(t);}
};

//template<class T>
//struct incrementable : weakly_incrementable<T> {
//	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
//	friend constexpr auto operator++(U& self, int) -> T {T tmp{self}; ++self; return tmp;}
//};

template<class Self> struct incrementable : totally_ordered<Self> {//, self_mutable<Self> {
	friend constexpr auto operator++(incrementable& s, int) -> Self {Self tmp{s.self()}; ++s.self(); assert(s.self() > tmp); return tmp;}
};

template<class T>
struct decrementable : weakly_decrementable<T> {
	template<class U, typename = std::enable_if_t<not std::is_base_of<T, U>{}>>
	friend constexpr auto operator--(U& self, int) -> T {T tmp{self}; --self; return tmp;}
};

template<class Self>
struct steppable : totally_ordered<Self> {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto self()       -> self_type      & {return static_cast<self_type      &>(*this);}

	friend constexpr auto operator++(steppable& s, int) -> Self {Self tmp{s.self()}; ++s.self(); return tmp;}
	friend constexpr auto operator--(steppable& s, int) -> Self {Self tmp{s.self()}; --s.self(); return tmp;}
};

template<class Self, typename Difference>
struct affine_with_unit : steppable<Self> {//affine_with_unit<Self, Difference> > {
	using self_type = Self;
	constexpr auto cself() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto  self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto  self()       -> self_type      & {return static_cast<self_type      &>(*this);}

	using difference_type = Difference;
	friend constexpr auto operator++(affine_with_unit& s) -> Self& {return s.self() += difference_type{1};}
	friend constexpr auto operator--(affine_with_unit& s) -> Self& {return s.self() -= difference_type{1};}

	friend constexpr auto operator-(affine_with_unit const& s, difference_type const& d) -> Self {
		auto ret{s.self()};
		ret += (-d);
		return ret;
	}
	constexpr auto operator+(difference_type const& d) const -> Self {
		auto ret{cself()};
		ret += d;
		return ret;
	}
//	friend
////  constexpr
//	auto operator+(affine_with_unit const& s, difference_type const& d) -> Self {
//		auto ret{s.self()};
//		ret += d;
//		return ret;
//	}
	friend constexpr auto operator+(difference_type const& d, affine_with_unit const& s) -> Self {
		auto ret{s.self()};
		ret += d;
		return ret;
	}
	friend constexpr auto operator<(affine_with_unit const& s, affine_with_unit const& o) -> bool {
		return difference_type{0} < o.self() - s.self();
	}
};

template<class Self, typename Reference>
struct dereferenceable {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto self()       -> self_type      & {return static_cast<self_type      &>(*this);}

	using reference = Reference;

	constexpr auto operator*() const -> reference {return *(self().operator->());}
};

template<class Self, typename Difference, typename Reference>
struct random_accessable  // NOLINT(fuchsia-multiple-inheritance)
: affine_with_unit<Self, Difference>
, dereferenceable<Self, Reference> {
	using difference_type = Difference;
	using reference = Reference;
	using iterator_category = std::random_access_iterator_tag;

	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto self()       -> self_type      & {return static_cast<self_type      &>(*this);}

	constexpr auto operator[](difference_type d) const -> reference {return *(self() + d);}
};

//template<class T, class Reference>
//struct dereferenceable {
//	using reference = Reference;
//	friend constexpr auto operator*(dereferenceable const& t) -> reference {return *static_cast<T const&>(t);}
//};

template<class T, class D>
struct addable2 {
	using difference_type = D;
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr auto operator+(TT&& t, difference_type const& d) -> T {T tmp{std::forward<TT>(t)}; tmp+=d; return tmp;}
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr auto operator+(difference_type const& d, TT&& t) -> T {return std::forward<TT>(t) + d;}
};

template<class T, class D>
struct subtractable2 {
	using difference_type = D;
	template<class TT, class = T>
	friend auto operator-(TT&& t, difference_type const& d) -> T {T tmp{std::forward<TT>(t)}; tmp-=d; return tmp;}
};

template<class T, class Difference>
struct affine : addable2<T, Difference>, subtractable2<T, Difference> {
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

template<class Self, class Value, class Reference = Value&, class Pointer = Value*, class Difference = std::ptrdiff_t>
struct random_access_iterator : equality_comparable2<Self, Self> {
	using difference_type = Difference;
	using value_type = Value;
	using pointer = Pointer;
	using reference = Reference;
	using iterator_category = std::random_access_iterator_tag;
	auto operator*() const -> Reference {return *static_cast<Self const&>(*this);}
};

}  // end namespace boost::multi
#endif
