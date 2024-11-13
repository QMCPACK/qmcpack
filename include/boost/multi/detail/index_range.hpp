// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
#define BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
#pragma once

#include <boost/multi/detail/implicit_cast.hpp>
#include <boost/multi/detail/serialization.hpp>
#include <boost/multi/detail/tuple_zip.hpp>
#include <boost/multi/detail/types.hpp>

#include <algorithm>    // for min, max
#include <cstddef>      // for ptrdiff_t
#include <functional>   // for minus, plus
#include <iterator>     // for reverse_iterator, random_access_iterator_tag
#include <limits>       // for numeric_limits
#include <memory>       // for pointer_traits
#include <type_traits>  // for declval, true_type, decay_t, enable_if_t
#include <utility>      // for forward

namespace boost::multi {

using boost::multi::detail::tuple;

template<
	class Self,
	class ValueType, class AccessCategory,
	class Reference = ValueType&, class DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, class Pointer = ValueType*>
class iterator_facade {
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

	friend constexpr auto operator<=(self_type const& self, self_type const& other) { return (self < other) || (self == other); }
	friend constexpr auto operator>(self_type const& self, self_type const& other) { return !(self <= other); }
	friend constexpr auto operator>=(self_type const& self, self_type const& other) { return !(self < other); }

	constexpr auto        operator-(difference_type n) const { return self_type{self_()} -= n; }
	constexpr auto        operator+(difference_type n) const { return self_type{self_()} += n; }
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
	IndexType     first_ = {};
	IndexTypeLast last_  = first_;  // TODO(correaa) check how to do partially initialzed

 public:
	template<class Archive>  // , class ArT = multi::archive_traits<Ar>>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		arxiv& multi::archive_traits<Archive>::make_nvp("first", first_);
		// arxiv &                  BOOST_SERIALIZATION_NVP(         first_);
		// arxiv &                        cereal:: make_nvp("first", first_);
		// arxiv &                               CEREAL_NVP(         first_);
		// arxiv &                                                   first_ ;

		arxiv& multi::archive_traits<Archive>::make_nvp("last", last_);
		// arxiv &                  BOOST_SERIALIZATION_NVP(         last_ );
		// arxiv &                        cereal:: make_nvp("last" , last_ );
		// arxiv &                               CEREAL_NVP(         last_ );
		// arxiv &                                                   last_  ;
	}

	using value_type      = IndexType;
	using difference_type = decltype(IndexTypeLast{} - IndexType{});  // std::make_signed_t<value_type>;
	using size_type       = difference_type;
	using const_reference = value_type;
	using reference       = const_reference;
	using const_pointer   = value_type;
	using pointer         = value_type;

	range() = default;

	// range(range const&) = default;

	template<class Range,
	         std::enable_if_t<!std::is_base_of_v<range, std::decay_t<Range>>, int> = 0,
	         decltype(detail::implicit_cast<IndexType>(std::declval<Range&&>().first()),
	                  detail::implicit_cast<IndexTypeLast>(std::declval<Range&&>().last())
	         )*                                                                    = nullptr>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ range(Range&& other)  // NOLINT(bugprone-forwarding-reference-overload,google-explicit-constructor,hicpp-explicit-conversions) // NOSONAR(cpp:S1709) ranges are implicitly convertible if elements are implicitly convertible
	: first_{std::forward<Range>(other).first()}, last_{std::forward<Range>(other).last()} {}

	template<
		class Range,
		std::enable_if_t<!std::is_base_of_v<range, std::decay_t<Range>>, int> = 0,
		decltype(detail::explicit_cast<IndexType>(std::declval<Range&&>().first()),
		         detail::explicit_cast<IndexTypeLast>(std::declval<Range&&>().last())
		)*                                                                    = nullptr>
	constexpr explicit range(Range&& other)  // NOLINT(bugprone-forwarding-reference-overload)
	: first_{std::forward<Range>(other).first()}, last_{std::forward<Range>(other).last()} {}

	constexpr range(IndexType first, IndexTypeLast last) : first_{first}, last_{last} {}

	class const_iterator : public boost::multi::iterator_facade<const_iterator, value_type, std::random_access_iterator_tag, const_reference, difference_type> {
		typename const_iterator::value_type curr_;
		constexpr explicit const_iterator(value_type current) : curr_{current} {}
		friend class range;

	 public:
		const_iterator() = default;

		constexpr auto operator==(const_iterator const& other) const -> bool { return curr_ == other.curr_; }
		constexpr auto operator!=(const_iterator const& other) const -> bool { return curr_ != other.curr_; }

		constexpr auto operator<(const_iterator const& other) const -> bool { return curr_ < other.curr_; }

		constexpr auto operator++() -> const_iterator& {
			++curr_;
			return *this;
		}
		constexpr auto operator--() -> const_iterator& {
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

	[[nodiscard]] constexpr auto first() const -> const_reference { return first_; }
	[[nodiscard]] constexpr auto last() const -> const_reference { return last_; }

	constexpr auto operator[](difference_type n) const -> const_reference { return first() + n; }

	[[nodiscard]] constexpr auto front() const -> const_reference { return first(); }
	[[nodiscard]] constexpr auto back() const -> const_reference { return last() - 1; }

	[[nodiscard]] constexpr auto cbegin() const { return const_iterator{first_}; }
	[[nodiscard]] constexpr auto cend() const { return const_iterator{last_}; }

	[[nodiscard]] constexpr auto rbegin() const { return reverse_iterator{end()}; }
	[[nodiscard]] constexpr auto rend() const { return reverse_iterator{begin()}; }

	[[nodiscard]] constexpr auto begin() const -> const_iterator { return cbegin(); }
	[[nodiscard]] constexpr auto end() const -> const_iterator { return cend(); }

	constexpr auto        is_empty() const& noexcept { return first_ == last_; }
	friend constexpr auto is_empty(range const& self) noexcept { return self.is_empty(); }

	[[nodiscard]] constexpr auto empty() const& noexcept { return is_empty(); }
	friend constexpr auto        empty(range const& self) noexcept { return self.empty(); }

	constexpr auto        size() const& noexcept -> size_type { return last_ - first_; }
	friend constexpr auto size(range const& self) noexcept -> size_type { return self.size(); }

	friend constexpr auto begin(range const& self) { return self.begin(); }
	friend constexpr auto end(range const& self) { return self.end(); }

	friend constexpr auto operator==(range const& self, range const& other) {
		return (self.empty() && other.empty()) || (self.first_ == other.first_ && self.last_ == other.last_);
	}
	friend constexpr auto operator!=(range const& self, range const& other) { return !(self == other); }

	[[nodiscard]] constexpr auto find(value_type const& value) const -> range::const_iterator {
		if(value >= last_ || value < first_) {
			return end();
		}
		return begin() + (value - front());
	}
	template<class Value> [[nodiscard]] constexpr auto contains(Value const& value) const -> bool { return (value >= first_) && (value < last_); }
	template<class Value> [[nodiscard]] constexpr auto count(Value const& value) const -> value_type { return contains(value); }

	friend constexpr auto intersection(range const& self, range const& other) {
		using std::max;
		using std::min;
		auto new_first = max(self.first(), other.first());
		auto new_last  = min(self.last(), other.last());
		new_first      = min(new_first, new_last);
		return range<decltype(new_first), decltype(new_last)>(new_first, new_last);
	}
	[[nodiscard]] constexpr auto contains(value_type const& value) const { return value >= first_ && value < last_; }
};

template<typename IndexType, typename IndexTypeLast = IndexType>     // , class Plus = std::plus<>, class Minus = std::minus<> >
range(IndexType, IndexTypeLast) -> range<IndexType, IndexTypeLast>;  // #3

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr auto make_range(IndexType first, IndexTypeLast last) -> range<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range {
	range<IndexType> impl_{
		(std::numeric_limits<IndexType>::min)(),  // parent needed for MSVC min/max macros
		(std::numeric_limits<IndexType>::max)()
	};

	constexpr intersecting_range() = default;  // MSVC 19.07 needs constexpr to initialize ALL later
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

// [[maybe_unused]] constexpr intersecting_range<> ∀ = V;
// [[maybe_unused]] constexpr intersecting_range<> https://www.compart.com/en/unicode/U+2200 = V;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
struct extension_t : public range<IndexType, IndexTypeLast> {
	using range<IndexType, IndexTypeLast>::range;

	constexpr extension_t(IndexType first, IndexTypeLast last) noexcept
	: range<IndexType, IndexTypeLast>{first, last} {}

	// cppcheck-suppress noExplicitConstructor ; because syntax convenience // NOLINTNEXTLINE(runtime/explicit)
	constexpr extension_t(IndexType last) noexcept  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) // NOSONAR(cpp:S1709) allow terse syntax
	: range<IndexType, IndexTypeLast>(0, last) {}

	constexpr extension_t() noexcept : range<IndexType, IndexTypeLast>() {}

	friend constexpr auto size(extension_t const& self) -> typename extension_t::size_type { return self.size(); }

	// constexpr auto operator==(extension_t const& other) const {return static_cast<range<IndexType> const&>(*this) == static_cast<range<IndexType> const&>(other);}
	// constexpr auto operator!=(extension_t const& other) const {return static_cast<range<IndexType> const&>(*this) != static_cast<range<IndexType> const&>(other);}

	// constexpr friend auto operator==(extension_t const& self, extension_t const& other) { return static_cast<range<IndexType> const&>(self) == static_cast<range<IndexType> const&>(other); }
	// constexpr friend auto operator!=(extension_t const& self, extension_t const& other) { return static_cast<range<IndexType> const&>(self) != static_cast<range<IndexType> const&>(other); }

	friend constexpr auto intersection(extension_t const& ex1, extension_t const& ex2) -> extension_t {
		using std::max;
		using std::min;

		auto       first = max(ex1.first(), ex2.first());
		auto const last  = min(ex1.last(), ex2.last());

		first = min(first, last);

		return extension_t{first, last};
	}
};

template<class IndexType, class IndexTypeLast>
extension_t(IndexType, IndexTypeLast) -> extension_t<IndexType, IndexTypeLast>;

template<class IndexType>
extension_t(IndexType) -> extension_t<IndexType>;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
constexpr auto make_extension_t(IndexType first, IndexTypeLast last) -> extension_t<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexTypeLast = std::ptrdiff_t>
constexpr auto make_extension_t(IndexTypeLast last) { return make_extension_t(IndexTypeLast{0}, last); }

using index_range     = range<index>;
using index_extension = extension_t<index>;
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
#endif  // BOOST_MULTI_DETAIL_INDEX_RANGE_HPP
