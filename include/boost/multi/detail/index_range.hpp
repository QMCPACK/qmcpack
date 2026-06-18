// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
#define BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
// #pragma once

#include "boost/multi/detail/implicit_cast.hpp"
#include "boost/multi/detail/serialization.hpp"
#include "boost/multi/detail/tuple_zip.hpp"
#include "boost/multi/detail/types.hpp"

#include "boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp"

#include <algorithm>    // for min, max
#include <cassert>
#include <cstddef>      // for ptrdiff_t
#include <functional>   // for minus, plus
#include <iterator>     // for reverse_iterator, random_access_iterator_tag
#include <limits>       // for numeric_limits
#include <memory>       // for pointer_traits
#include <type_traits>  // for declval, true_type, decay_t, enable_if_t
#include <utility>      // for forward

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost::multi {

using boost::multi::detail::tuple;

template<
	class Self,
	class ValueType, class AccessCategory,
	class Reference = ValueType&, class DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, class Pointer = ValueType*>
class iterator_facade {
 protected:
	iterator_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	friend Self;

 private:
	using self_type = Self;
	[[nodiscard]] constexpr auto self_() & { return static_cast<self_type&>(*this); }
	[[nodiscard]] constexpr auto self_() const& { return static_cast<self_type const&>(*this); }

 public:
	using value_type        = ValueType;
	using reference         = Reference;
	using pointer           = Pointer;  // NOSONAR(cpp:S5008) false positive
	using difference_type   = DifferenceType;
	using iterator_category = AccessCategory;

	// friend constexpr auto operator!=(self_type const& self, self_type const& other) { return !(self == other); }

	friend constexpr auto operator<=(self_type const& self, self_type const& other) { return self < other || self == other; }
	friend constexpr auto operator>(self_type const& self, self_type const& other) { return !(self <= other); }
	friend constexpr auto operator>=(self_type const& self, self_type const& other) { return !(self < other); }

	constexpr auto        operator-(difference_type n) const { return self_type{self_()} -= n; }
	constexpr auto        operator+(difference_type n) const { return self_type{self_()} += n; }

	template<class = void>  // nvcc workaround 
	friend constexpr auto operator+(difference_type n, self_type const& self) { return self + n; }

	friend constexpr auto operator++(self_type& self, int) -> self_type {
		self_type ret = self;
		++self;
		return ret;
	}
	friend constexpr auto operator--(self_type& self, int) -> self_type {
		self_type ret = self;
		--self;
		return ret;
	}

	constexpr auto operator[](difference_type n) const { return *(self_() + n); }
};

template<typename IndexType = std::true_type, typename IndexTypeLast = IndexType, class Plus = std::plus<>, class Minus = std::minus<>>
class range {
	#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4820)	// 'boost::multi::range<IndexType,IndexTypeLast,std::plus<void>,std::minus<void>>': '3' bytes padding added after data member 'boost::multi::range<IndexType,IndexTypeLast,std::plus<void>,std::minus<void>>::first_'
	#endif

	BOOST_MULTI_NO_UNIQUE_ADDRESS
	IndexType first_;  //  = {};

	#ifdef _MSC_VER
	#pragma warning(pop)
	#endif

	#ifdef __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wpadded"
	#endif

	IndexTypeLast last_;  // = first_;  // TODO(correaa) check how to do partially initialzed

	#ifdef __clang__
	#pragma clang diagnostic pop
	#endif

 public:
	template<class Archive>  // , class ArT = multi::archive_traits<Ar>>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		arxiv
			& multi::archive_traits<Archive>::make_nvp("first", first_)
			& multi::archive_traits<Archive>::make_nvp("last", last_)
		;

		// arxiv & multi::archive_traits<Archive>::make_nvp("first", first_);
		// arxiv &               BOOST_SERIALIZATION_NVP(         first_);
		// arxiv &                     cereal:: make_nvp("first", first_);
		// arxiv &                            CEREAL_NVP(         first_);
		// arxiv &                                                first_ ;

		// arxiv & multi::archive_traits<Archive>::make_nvp("last", last_);
		// arxiv &                  BOOST_SERIALIZATION_NVP(         last_ );
		// arxiv &                        cereal:: make_nvp("last" , last_ );
		// arxiv &                               CEREAL_NVP(         last_ );
		// arxiv &                                                   last_  ;
	}

	using value_type      = decltype(IndexTypeLast{} + IndexType{});
	using difference_type = decltype(IndexTypeLast{} - IndexType{});  // std::make_signed_t<value_type>;
	using size_type       = difference_type;
	using const_reference = value_type;
	using reference       = const_reference;
	using const_pointer   = value_type;
	using pointer         = value_type;

	range() = default;  // cppcheck-suppress uninitMemberVar ;

	// range(range const&) = default;

	template<class Range,
	         std::enable_if_t<!std::is_base_of_v<range, std::decay_t<Range>>, int> = 0,  // NOLINT(modernize-type-traits) for C++20
	         decltype(detail::implicit_cast<IndexType>(std::declval<Range&&>().first()),
	                  detail::implicit_cast<IndexTypeLast>(std::declval<Range&&>().last())
	         )*                                                                    = nullptr>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ range(Range&& other)  // NOLINT(bugprone-forwarding-reference-overload,google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) // NOSONAR(cpp:S1709) ranges are implicitly convertible if elements are implicitly convertible
	: first_{std::forward<Range>(other).first()}, last_{std::forward<Range>(other).last()} {}  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)

	template<
		class Range,
		std::enable_if_t<!std::is_base_of_v<range, std::decay_t<Range>>, unsigned> = 0,
		decltype(detail::explicit_cast<IndexType>(std::declval<Range&&>().first()),
		         detail::explicit_cast<IndexTypeLast>(std::declval<Range&&>().last())
		)*                                                                    = nullptr>
	constexpr explicit range(Range&& other)  // NOLINT(bugprone-forwarding-reference-overload)
	: first_{std::forward<Range>(other).first()}, last_{std::forward<Range>(other).last()} {}

	BOOST_MULTI_HD constexpr range(IndexType first, IndexTypeLast last) : first_{first}, last_{last} {}

	// TODO(correaa) make this iterator SCARY
	class const_iterator : public boost::multi::iterator_facade<const_iterator, value_type, std::random_access_iterator_tag, const_reference, difference_type> {
		typename const_iterator::value_type curr_;
		constexpr explicit const_iterator(value_type current) : curr_{current} {}
		friend class range;

	 public:
		template<class T> using rebind = typename range<std::decay_t<T>>::const_iterator;
		using pointer = const_iterator;
		using element_type = IndexTypeLast;
 
		const_iterator() = default;

		template<class OtherConstIterator, class = decltype(std::declval<typename const_iterator::value_type&>() = *OtherConstIterator{})>
		// cppcheck-suppress noExplicitConstructor ; see below
		const_iterator(OtherConstIterator const& other) : curr_{*other} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

		BOOST_MULTI_HD constexpr auto operator==(const_iterator const& other) const -> bool { return curr_ == other.curr_; }
		BOOST_MULTI_HD constexpr auto operator!=(const_iterator const& other) const -> bool { return curr_ != other.curr_; }

		BOOST_MULTI_HD constexpr auto operator<(const_iterator const& other) const -> bool { return curr_ < other.curr_; }  // mull-ignore: cxx_lt_to_le

		constexpr auto operator++() -> const_iterator& {
			++curr_;
			return *this;
		}
		constexpr auto operator--() noexcept(noexcept(--curr_)) -> const_iterator& {
			--curr_;
			return *this;
		}

		constexpr auto operator-=(typename const_iterator::difference_type n) -> const_iterator& {
			curr_ -= n;
			return *this;
		}
		constexpr auto operator+=(typename const_iterator::difference_type n) -> const_iterator& {
			curr_ += n;
			return *this;
		}

		constexpr auto operator-(typename const_iterator::difference_type n) const -> const_iterator {
			return const_iterator{*this} -= n;
		}

		constexpr auto operator+(typename const_iterator::difference_type n) const -> const_iterator {
			return const_iterator{*this} += n;
		}

		constexpr auto operator-(const_iterator const& other) const { return curr_ - other.curr_; }
		constexpr auto operator*() const noexcept -> typename const_iterator::reference { return curr_; }
	};

	using iterator               = const_iterator;
	using reverse_iterator       = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	[[nodiscard]] BOOST_MULTI_HD constexpr auto first() const { return first_; }
	[[nodiscard]] BOOST_MULTI_HD constexpr auto last() const { return last_; }

	constexpr auto operator[](difference_type n) const -> const_reference {
		assert(n >= 0);
		assert(n < size());
		return first() + n;
	}

	[[nodiscard]] BOOST_MULTI_HD constexpr auto front() const -> value_type { return first(); }  // cppcheck-suppress functionStatic ;  // bug in cppcheck 2.19.0
	[[nodiscard]] BOOST_MULTI_HD constexpr auto back() const -> value_type { return last() - 1; }  // cppcheck-suppress functionStatic ;  // bug in cppcheck 2.19.0

	[[nodiscard]] constexpr auto cbegin() const { return const_iterator{first_}; }
	[[nodiscard]] constexpr auto cend() const { return const_iterator{last_}; }

	[[nodiscard]] constexpr auto rbegin() const { return reverse_iterator{end()}; }
	[[nodiscard]] constexpr auto rend() const { return reverse_iterator{begin()}; }

	[[nodiscard]] constexpr auto begin() const noexcept -> const_iterator { return cbegin(); }
	[[nodiscard]] constexpr auto end() const noexcept -> const_iterator { return cend(); }

	BOOST_MULTI_HD constexpr auto        is_empty() const& noexcept { return first_ == last_; }

	[[nodiscard]] BOOST_MULTI_HD constexpr auto empty() const& noexcept { return is_empty(); }

	#ifdef __NVCC__
	#pragma nv_diagnostic push
	#pragma nv_diag_suppress = 20013  // calling a constexpr __host__ function("operator std::streamoff") from a __host__ __device__ function("size") is not allowed.  // TODO(correaa) implement HD integral_constant
	#endif

	BOOST_MULTI_HD constexpr auto  size() const noexcept -> size_type { return last_ - first_; }
	BOOST_MULTI_HD constexpr auto ssize() const noexcept { return size(); }
	BOOST_MULTI_HD constexpr auto usize() const noexcept { return static_cast<std::size_t>(size()); }

	#ifdef __NVCC__
	#pragma nv_diagnostic pop
	#endif

	friend BOOST_MULTI_HD constexpr auto operator!=(range const& self, range const& other) { return !(self == other); }  // NOLINT(readability-redundant-parentheses) bug in clang-tidy
	friend BOOST_MULTI_HD constexpr auto operator==(range const& self, range const& other) {
		return (self.empty() && other.empty()) || (self.first_ == other.first_ && self.last_ == other.last_);
	}

	[[nodiscard]]  // ("find returns an iterator to the sequence, that is the only effect")]] for C++20
	constexpr auto find(value_type const& value) const -> const_iterator {
		if(value >= last_ || value < first_) {
			return end();
		}
		return begin() + (value - front());
	}
	template<class Value> [[nodiscard]] BOOST_MULTI_HD constexpr auto contains(Value const& value) const -> bool { return (first_ <= value) && (value < last_); }
	template<class Value> [[nodiscard]] BOOST_MULTI_HD constexpr auto count(Value const& value) const -> size_type { return contains(value); }

	friend constexpr auto intersection(range const& self, range const& other) {
		using std::max;
		using std::min;
		auto new_first = max(self.first(), other.first());
		auto new_last  = min(self.last(), other.last());
		new_first      = min(new_first, new_last);
		return range<decltype(new_first), decltype(new_last)>(new_first, new_last);
	}
};

#if defined(__cpp_deduction_guides) && (__cpp_deduction_guides >= 201703)
template<typename IndexType, typename IndexTypeLast = IndexType>     // , class Plus = std::plus<>, class Minus = std::minus<> >
range(IndexType, IndexTypeLast) -> range<IndexType, IndexTypeLast>;  // #3
#endif

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr auto make_range(IndexType first, IndexTypeLast last) -> range<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range {
	range<IndexType> impl_;
	
	constexpr intersecting_range() noexcept :  // MSVC 19.07 needs constexpr to initialize ALL later
		impl_(
			(std::numeric_limits<IndexType>::min)(),  // NOLINT(readability-redundant-parentheses) for MSVC min macros
			(std::numeric_limits<IndexType>::max)()   // NOLINT(readability-redundant-parentheses) for MSVC max macros
		)
	{}

	static constexpr auto make_(IndexType first, IndexType last) -> intersecting_range {
		intersecting_range ret;
		ret.impl_ = range<IndexType>{first, last};
		return ret;
	}
	friend constexpr auto intersection(intersecting_range const& self, range<IndexType> const& other) {
		return intersection(self.impl_, other);
	}
	friend constexpr auto intersection(range<IndexType> const& other, intersecting_range const& self) {
		return intersection(other, self.impl_);
	}
	friend constexpr auto operator<(intersecting_range const& self, IndexType end) {
		return intersecting_range::make_(self.impl_.first(), end);
	}
	friend constexpr auto operator<=(IndexType first, intersecting_range const& self) {
		return intersecting_range::make_(first, self.impl_.last());
	}

 public:
	constexpr auto        operator*() const& -> intersecting_range const& { return *this; }
	static constexpr auto all() noexcept { return intersecting_range{}; }
};

[[maybe_unused]] constexpr intersecting_range<> ALL = intersecting_range<>::all();
[[maybe_unused]] constexpr intersecting_range<> _   = ALL;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> U   = ALL;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> ooo = ALL;

[[maybe_unused]] constexpr intersecting_range<> V = U;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> A = V;  // NOLINT(readability-identifier-length)

#if !defined(__clang__) && !defined(CPPCHECK)
// cppcheck-suppress preprocessorErrorDirective ;  // unicode
// [[maybe_unused]] constexpr intersecting_range<> ∀ = V;  // not valid in clang or g++-15
// [[maybe_unused]] constexpr intersecting_range<> https://www.compart.com/en/unicode/U+2200 = V;
// [[maybe_unused]] constexpr intersecting_range<> ┄ = ALL;  // not valid in g++-15
#endif

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + IndexType{1})>
struct extent_t : public range<IndexType, IndexTypeLast> {
	using range<IndexType, IndexTypeLast>::range;

	BOOST_MULTI_HD constexpr extent_t(IndexType first, IndexTypeLast last) noexcept
	: range<IndexType, IndexTypeLast>{first, last} {}

//	BOOST_MULTI_HD constexpr extent_t(extent_t::size_type size) : extents_t(IndexType{}, IndexType{} + size) {}

	// cppcheck-suppress noExplicitConstructor ; because syntax convenience // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extent_t(IndexTypeLast last) noexcept  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) // NOSONAR(cpp:S1709) allow terse syntax
	: range<IndexType, IndexTypeLast>(IndexType{}, IndexType{} + last) {}

	template<
		class OtherExtension,
		decltype(
			detail::implicit_cast<IndexType>(std::declval<OtherExtension>().first()),
			detail::implicit_cast<IndexTypeLast>(std::declval<OtherExtension>().last())
		)* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extent_t(OtherExtension const& other) noexcept  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: extent_t{other.first(), other.last()} {}

	using index = IndexTypeLast;

	template<class OtherExtension>
	BOOST_MULTI_HD constexpr auto operator=(OtherExtension const& other) -> extent_t& {
		(*this) = extent_t{other};
		return *this;
	}

	// BOOST_MULTI_HD constexpr extent_t() noexcept : range<IndexType, IndexTypeLast>() {}
	constexpr extent_t() = default;

	friend constexpr auto intersection(extent_t const& ex1, extent_t const& ex2) -> extent_t {
		using std::max;
		using std::min;

		auto       first = max(ex1.first(), ex2.first());
		auto const last  = min(ex1.last(), ex2.last());

		first = min(first, last);

		return extent_t{first, last};
	}
};

#if defined(__cpp_deduction_guides) && (__cpp_deduction_guides >= 201703)
template<class IndexType, class IndexTypeLast>
extent_t(IndexType, IndexTypeLast) -> extent_t<IndexType, IndexTypeLast>;

template<class IndexType = multi::index>
extent_t(IndexType) -> extent_t<std::integral_constant<IndexType, 0>, IndexType>;
#endif

// deprecated spelling kept so old code keeps compiling; use extent_t
template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + IndexType{1})>
using extension_t [[deprecated("use extent_t")]] = extent_t<IndexType, IndexTypeLast>;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
constexpr auto make_extension_t(IndexType first, IndexTypeLast last) {
	return extent_t<IndexType, IndexTypeLast>{first, last};
}

template<class IndexType = multi::index>
constexpr auto make_extension_t(IndexType last) { return make_extension_t(std::integral_constant<IndexType, 0>{}, last); }

using index_range     = range<index>;
using index_extension = extent_t<index>;
using iextension      = index_extension;
using irange          = index_range;

namespace detail {

template<typename, typename>
struct append_to_type_seq {};

template<typename T, typename... Ts, template<typename...> class TT>
struct append_to_type_seq<T, TT<Ts...>> {
	using type = TT<Ts..., T>;
};

template<typename T, dimensionality_type N, template<typename...> class TT>
struct repeat {
	using type = typename append_to_type_seq<
		T,
		typename repeat<T, N - 1, TT>::type>::type;
};

template<typename T, template<typename...> class TT>
struct repeat<T, 0, TT> {
	using type = TT<>;
};

}  // end namespace detail

template<dimensionality_type D> using index_extensions = typename detail::repeat<index_extension, D, tuple>::type;

template<dimensionality_type D, class Tuple>
constexpr auto contains(index_extensions<D> const& iex, Tuple const& tup) {
	//  using detail::head;
	//  using detail::tail;
	return contains(head(iex), head(tup)) && contains(tail(iex), tail(tup));
}

}  // end namespace boost::multi

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
