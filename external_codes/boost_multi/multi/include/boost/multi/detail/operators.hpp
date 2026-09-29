// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_OPERATORS_HPP
#define BOOST_MULTI_DETAIL_OPERATORS_HPP
// #pragma once

#include <cstddef>      // for ptrdiff_t
#include <iterator>     // for random_access_iterator_tag
#include <type_traits>  // for enable_if_t, is_base_of
#include <utility>      // for forward

#ifdef __NVCC__
	#define BOOST_MULTI_HD __host__ __device__
#else
	#define BOOST_MULTI_HD
#endif

namespace boost::multi {

struct empty_base {};

template<class Self> struct selfable {
 protected:
	selfable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using self_type = Self;
	constexpr auto        self() const -> self_type const& {
		static_assert(std::is_base_of_v<selfable<Self>, Self>);
		return static_cast<self_type const&>(*this);
	}
	constexpr auto        self() -> self_type& {
		static_assert(std::is_base_of_v<selfable<Self>, Self>);
		return static_cast<self_type&>(*this);
	}
	friend constexpr auto self(selfable const& self) -> self_type const& {
		static_assert(std::is_base_of_v<selfable<Self>, Self>);
		return self.self();
	}
};

template<class Self>
class ra_iterable : selfable<Self> {
	ra_iterable() = default;
	friend Self;

	template<class Self2 = Self>
	using difference_type_t = decltype(std::declval<Self2 const&>() - std::declval<Self2 const&>());

 public:
	using iterator_category = std::random_access_iterator_tag;

	template<class Self2 = Self, std::enable_if_t<std::is_same_v<Self2, Self>, int> =0>  // NOLINT(modernize-use-constraints) for C++20
	friend BOOST_MULTI_HD constexpr auto operator+(Self2 self, difference_type_t<Self2> const& n) { return self += n; }
	// template<class Self2 = Self>
	// friend auto operator+(difference_type<Self2> const& n, Self2 const& self) { return self + n; }
	BOOST_MULTI_HD constexpr auto operator++(int) { Self tmp{*this}; ++this->self(); return tmp; }
	BOOST_MULTI_HD constexpr auto operator--(int) { Self tmp{*this}; --this->self(); return tmp; }  // NOLINT(cert-dcl21-cpp)

	template<class Self2 = Self, std::enable_if_t<std::is_same_v<Self2, Self>, int> =0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr auto friend operator!=(Self2 const& self, Self2 const& other) { return !(self == other); }
};

template<class Self, class U> struct equality_comparable2;

template<class Self>
struct equality_comparable2<Self, Self> : selfable<Self> {
	// friend constexpr auto operator==(equality_comparable2 const& self, equality_comparable2 const& other) {return     self.self() == other.self() ;}
	friend constexpr auto operator!=(equality_comparable2 const& self, equality_comparable2 const& other) { return !(self.self() == other.self()); }
};

template<class Self> struct equality_comparable : equality_comparable2<Self, Self> {
 protected:
	equality_comparable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;
};

template<class T, class V> struct totally_ordered2;

template<class Self>
struct totally_ordered2<Self, Self> : equality_comparable2<totally_ordered2<Self, Self>, totally_ordered2<Self, Self>> {
	using self_type = Self;
	BOOST_MULTI_HD constexpr auto self() const -> self_type const& { return static_cast<self_type const&>(*this); }

	friend BOOST_MULTI_HD constexpr auto operator==(totally_ordered2 const& myself, totally_ordered2 const& other) -> bool { return !(myself.self() < other.self()) && !(other.self() < myself.self()); }

	friend BOOST_MULTI_HD constexpr auto operator<=(totally_ordered2 const& myself, totally_ordered2 const& other) -> bool { return !(other.self() < myself.self()); }

	friend BOOST_MULTI_HD constexpr auto operator>(totally_ordered2 const& myself, totally_ordered2 const& other) -> bool { return !(myself.self() < other.self()) && !(myself.self() == other.self()); }
	friend BOOST_MULTI_HD constexpr auto operator>=(totally_ordered2 const& myself, totally_ordered2 const& other) -> bool { return !(myself.self() < other.self()); }
};

template<class Self> using totally_ordered = totally_ordered2<Self, Self>;

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4820 )  // '3' bytes padding added after data member
#endif
template<class T>
struct totally_ordered2<T, void> {
	// template<class U>
	// friend constexpr auto operator<=(T const& self, U const& other) { return (self < other) || (self == other); }
	// template<class U>
	// friend constexpr auto operator>=(T const& self, U const& other) { return (other < self) || (self == other); }
	// template<class U>
	// friend constexpr auto operator>(T const& self, U const& other) { return other < self; }
};
#ifdef _MSC_VER
#pragma warning( pop )
#endif

template<class T>
struct copy_constructible {};

template<class T>
struct weakly_incrementable : selfable<T> {
 protected:
	weakly_incrementable() = default;

 public:
	constexpr auto operator++(int) -> T {
		auto ret{this->self()}; ++ this->self(); return ret;
	}
};

template<class T>
struct weakly_decrementable {
 protected:
	weakly_decrementable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend T;
	// friend T& operator--(weakly_decrementable& t){return --static_cast<T&>(t);}
};

template<class Self>
struct incrementable : totally_ordered<Self> {
 protected:
	incrementable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	friend BOOST_MULTI_HD constexpr auto operator++(incrementable& self, int) -> Self {
		static_assert(std::is_base_of_v<incrementable<Self>, Self>);
		Self tmp{self.self()};
		++self.self();
		assert(self.self() > tmp);
		return tmp;
	}
};

template<class T>
struct decrementable : weakly_decrementable<T> {
 protected:
	decrementable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend T;

 public:
	template<class U, typename = std::enable_if_t<!std::is_base_of_v<T, U>>>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend constexpr auto operator--(U& self, int) -> T {
		T tmp{self};
		--self;
		return tmp;
	}
};

template<class Self>
struct steppable : totally_ordered<Self> {
 protected:
	steppable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using self_type = Self;
	BOOST_MULTI_HD constexpr auto self() const -> self_type const& { return static_cast<self_type const&>(*this); }
	BOOST_MULTI_HD constexpr auto self() -> self_type& { return static_cast<self_type&>(*this); }

	friend BOOST_MULTI_HD constexpr auto operator++(steppable& myself, int) -> Self {
		Self tmp{myself.self()};
		++myself.self();
		return tmp;
	}
	friend BOOST_MULTI_HD constexpr auto operator--(steppable& myself, int) -> Self {
		Self tmp{myself.self()};
		--myself.self();
		return tmp;
	}
};

template<class Self, typename Difference>
struct affine_with_unit : steppable<Self> {
 protected:
	affine_with_unit() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using self_type = Self;
	BOOST_MULTI_HD constexpr auto cself() const -> self_type const& { return static_cast<self_type const&>(*this); }
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto self() const -> self_type const& { return static_cast<self_type const&>(*this); }
	BOOST_MULTI_HD constexpr auto self() -> self_type& { return static_cast<self_type&>(*this); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	using difference_type = Difference;
	friend BOOST_MULTI_HD constexpr auto operator++(affine_with_unit& myself) -> Self& { return myself.self() += difference_type{1}; }
	friend BOOST_MULTI_HD constexpr auto operator--(affine_with_unit& myself) -> Self& { return myself.self() -= difference_type{1}; }

	BOOST_MULTI_HD constexpr auto operator+(difference_type const& diff) const -> Self {
		auto ret{cself()};
		ret += diff;
		return ret;
	}
	friend BOOST_MULTI_HD constexpr auto operator+(difference_type const& diff, affine_with_unit const& myself) -> Self {
		auto ret{myself.self()};
		ret += diff;
		return ret;
	}
	friend constexpr auto operator<(affine_with_unit const& myself, affine_with_unit const& other) -> bool {
		return difference_type{0} < other.self() - myself.self();
	}
};

template<class Self, typename Reference>
struct dereferenceable {
 protected:
	dereferenceable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using self_type = Self;
	constexpr auto self() const -> self_type const& { return static_cast<self_type const&>(*this); }
	constexpr auto self() -> self_type& { return static_cast<self_type&>(*this); }

	using reference = Reference;

	BOOST_MULTI_HD constexpr auto operator*() const -> reference { return * self().operator->(); }
};

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4820 )  // '7' bytes padding added after base class
#endif

template<class Self, typename Difference, typename Reference>
struct random_accessable
: affine_with_unit<Self, Difference>
, dereferenceable<Self, Reference> {
 protected:
	random_accessable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using difference_type   = Difference;
	using reference         = Reference;
	using iterator_category = std::random_access_iterator_tag;

	using self_type = Self;
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto self() const -> self_type const& { return static_cast<self_type const&>(*this); }
	BOOST_MULTI_HD constexpr auto self() -> self_type& { return static_cast<self_type&>(*this); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	BOOST_MULTI_HD constexpr auto operator[](difference_type idx) const -> reference { return *(self() + idx); }
};

#ifdef _MSC_VER
#pragma warning( pop )
#endif

template<class Self, class D>
class addable2 {
 protected:
	addable2() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using difference_type = D;

	template<class TT, typename = std::enable_if_t<std::is_base_of<Self, TT>{}>>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend BOOST_MULTI_HD constexpr auto operator+(TT&& myself, difference_type const& diff) -> Self { return Self{std::forward<TT>(myself)} += diff; }

	template<class TT, typename = std::enable_if_t<std::is_base_of<Self, TT>{}>>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend BOOST_MULTI_HD constexpr auto operator+(difference_type const& diff, TT&& myself) -> Self { return std::forward<TT>(myself) + diff; }
};

template<class T, class D>
class subtractable2 {
 protected:
	subtractable2() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend T;

 public:
	using difference_type = D;
	// TODO(correaa) clang 16 picks up this and converts the difference_type to TT !!
	// template<class TT, class = T>
	// friend auto operator-(TT&& self, difference_type const& diff) -> T {T tmp{std::forward<TT>(self)}; tmp -= diff; return tmp;}
};

template<class T, class Difference>
struct affine : addable2<T, Difference>
, subtractable2<T, Difference> {
 protected:
	affine() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend T;

 public:
	using difference_type = Difference;
};

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4820 )  // '3' bytes padding added after data member
#endif
template<class T>
class random_iterable {
 protected:
	random_iterable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend T;

 public:
	constexpr auto        cfront() const& -> decltype(auto) { return static_cast<T const&>(*this).front(); }
	constexpr auto        cback() const& -> decltype(auto) { return static_cast<T const&>(*this).back(); }
	friend constexpr auto cfront(T const& myself) -> decltype(auto) { return myself.cfront(); }
	friend constexpr auto cback(T const& myself) -> decltype(auto) { return myself.cback(); }
};
#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace detail {
template<class Self, class Value, class Reference = Value&, class Pointer = Value*, class Difference = std::ptrdiff_t>
struct random_access_iterator : equality_comparable2<Self, Self> {
 protected:
	random_access_iterator() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 public:
	using difference_type   = Difference;
	using value_type        = Value;
	using pointer           = Pointer;
	using reference         = Reference;
	using iterator_category = std::random_access_iterator_tag;
	BOOST_MULTI_HD constexpr auto operator*() const -> Reference { return *static_cast<Self const&>(*this); }
};
}  // end namespace detail

}  // end namespace boost::multi

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_OPERATORS_HPP
