// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_OPERATORS_HPP
#define MULTI_DETAIL_OPERATORS_HPP

#include<type_traits>  // for enable_if
#include<utility>  // for forward

namespace boost::multi {

struct empty_base {};

template<class Self> struct selfable {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto self()       -> self_type      & {return static_cast<self_type      &>(*this);}
	friend constexpr auto self(selfable const& self) -> self_type const& {return self.self();}
};

template<class Self, class U> struct equality_comparable2;

template<class Self>
struct equality_comparable2<Self, Self> : selfable<Self> {
//	friend constexpr auto operator==(equality_comparable2 const& self, equality_comparable2 const& other) {return     self.self() == other.self() ;}
	friend constexpr auto operator!=(equality_comparable2 const& self, equality_comparable2 const& other) {return not(self.self() == other.self());}
};

template<class Self> struct equality_comparable : equality_comparable2<Self, Self> {};

template<class T, class V> struct totally_ordered2;

template<class Self>
struct totally_ordered2<Self, Self> : equality_comparable2<totally_ordered2<Self, Self>, totally_ordered2<Self, Self>> {
	using self_type = Self;
	constexpr auto self() const -> self_type const& {return static_cast<self_type const&>(*this);}

//	friend auto operator< (totally_ordered2 const& self, totally_ordered2 const& other) -> bool {return     self.self() < other.self() ;}
	friend auto operator==(totally_ordered2 const& self, totally_ordered2 const& other) -> bool {return not(self.self() < other.self()) and not(other.self() < self.self());}
//	friend auto operator!=(totally_ordered2 const& self, totally_ordered2 const& other) {return    (s.self() < o.self()) or     (o.self() < s.self());}

	friend auto operator<=(totally_ordered2 const& self, totally_ordered2 const& other) -> bool {return not(other.self() < self.self());}

	friend auto operator> (totally_ordered2 const& self, totally_ordered2 const& other) -> bool {return not(self.self() < other.self()) and not(self.self() == other.self());}
	friend auto operator>=(totally_ordered2 const& self, totally_ordered2 const& other) -> bool {return not(self.self() < other.self());}
};

template<class Self> using totally_ordered = totally_ordered2<Self, Self>;

template<class T>
struct totally_ordered2<T, void> {
	template<class U>
	friend constexpr auto operator<=(const T& self, const U& other) {return (self < other) or (self == other);}
	template<class U>
	friend constexpr auto operator>=(const T& self, const U& other) {return (other < self) or (self == other);}
	template<class U>
	friend constexpr auto operator> (const T& self, const U& other) {return  other < self;}
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

template<class Self> struct incrementable : totally_ordered<Self> {//, self_mutable<Self> {
	friend constexpr auto operator++(incrementable& self, int) -> Self {Self tmp{self.self()}; ++self.self(); assert(self.self() > tmp); return tmp;}
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

	friend constexpr auto operator++(steppable& self, int) -> Self {Self tmp{self.self()}; ++self.self(); return tmp;}
	friend constexpr auto operator--(steppable& self, int) -> Self {Self tmp{self.self()}; --self.self(); return tmp;}
};

template<class Self, typename Difference>
struct affine_with_unit : steppable<Self> {//affine_with_unit<Self, Difference> > {
	using self_type = Self;
	constexpr auto cself() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto  self() const -> self_type const& {return static_cast<self_type const&>(*this);}
	constexpr auto  self()       -> self_type      & {return static_cast<self_type      &>(*this);}

	using difference_type = Difference;
	friend constexpr auto operator++(affine_with_unit& self) -> Self& {return self.self() += difference_type{1};}
	friend constexpr auto operator--(affine_with_unit& self) -> Self& {return self.self() -= difference_type{1};}

	friend constexpr auto operator-(affine_with_unit const& self, difference_type const& diff) -> Self {
		auto ret{self.self()};
		ret += (-diff);
		return ret;
	}
	constexpr auto operator+(difference_type const& diff) const -> Self {
		auto ret{cself()};
		ret += diff;
		return ret;
	}
	friend constexpr auto operator+(difference_type const& diff, affine_with_unit const& self) -> Self {
		auto ret{self.self()};
		ret += diff;
		return ret;
	}
	friend constexpr auto operator<(affine_with_unit const& self, affine_with_unit const& other) -> bool {
		return difference_type{0} < other.self() - self.self();
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

	constexpr auto operator[](difference_type idx) const -> reference {return *(self() + idx);}
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
	friend constexpr auto operator+(TT&& self, difference_type const& diff) -> T {T tmp{std::forward<TT>(self)}; tmp += diff; return tmp;}
	template<class TT, typename = std::enable_if_t<std::is_base_of<T, TT>{}> >
	friend constexpr auto operator+(difference_type const& diff, TT&& self) -> T {return std::forward<TT>(self) + diff;}
};

template<class T, class D>
struct subtractable2 {
	using difference_type = D;
	template<class TT, class = T>
	friend auto operator-(TT&& self, difference_type const& diff) -> T {T tmp{std::forward<TT>(self)}; tmp -= diff; return tmp;}
};

template<class T, class Difference>
struct affine : addable2<T, Difference>, subtractable2<T, Difference> {
	using difference_type = Difference;
};

template<class T>
struct random_iterable {
	       constexpr auto cfront() const&       -> decltype(auto) {return static_cast<T const&>(*this).front();}
	       constexpr auto cback () const&       -> decltype(auto) {return static_cast<T const&>(*this).back() ;}
	friend constexpr auto cfront(T const& self) -> decltype(auto) {return self.cfront();}
	friend constexpr auto cback (T const& self) -> decltype(auto) {return self.cback() ;}
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
