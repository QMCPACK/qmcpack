// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_INDEX_RANGE_HPP
#define MULTI_DETAIL_INDEX_RANGE_HPP

#include "multi/detail/serialization.hpp"
#include "multi/detail/tuple_zip.hpp"
#include "multi/detail/types.hpp"

#include <algorithm>  // for min
#include <iterator>   // for std::random_iterator_tag // std::reverse_iterator
#include <limits>     // for numeric_limits
#include <utility>    // for forward

namespace boost::multi {

using boost::multi::detail::tuple;

template<
	class Self,
	class ValueType, class AccessCategory,
	class Reference = ValueType&, class DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, class Pointer = ValueType*
>
class iterator_facade {
	using self_type = Self;
	[[nodiscard]] constexpr auto self()      & {return static_cast<self_type      &>(*this);}
	[[nodiscard]] constexpr auto self() const& {return static_cast<self_type const&>(*this);}

 public:
	using value_type        = ValueType;
	using reference         = Reference;
	using pointer           = Pointer;
	using difference_type   = DifferenceType;
	using iterator_category = AccessCategory;

	friend constexpr auto operator!=(self_type const& self, self_type const& other) {return not(self == other);}

	friend constexpr auto operator<=(self_type const& self, self_type const& other) {return (self < other) or (self == other);}
	friend constexpr auto operator> (self_type const& self, self_type const& other) {return not(self <= other);}
	friend constexpr auto operator>=(self_type const& self, self_type const& other) {return not(self <  other);}

	       constexpr auto operator-(difference_type n) const {return self_type{self()} -= n;}
	       constexpr auto operator+(difference_type n) const {return self_type{self()} += n;}
	friend constexpr auto operator+(difference_type n, self_type const& self) {return self + n;}

	friend constexpr auto operator++(self_type& self, int) -> self_type {self_type ret = self; ++self; return ret;}
	friend constexpr auto operator--(self_type& self, int) -> self_type {self_type ret = self; --self; return ret;}

	constexpr auto operator[](difference_type n) const {return *(self() + n);}
};

template<typename IndexType = std::true_type, typename IndexTypeLast = IndexType>
class range {
	IndexType     first_ = {};
	IndexTypeLast last_  = first_;

 public:
	template<class Archive>  // , class ArT = multi::archive_traits<Ar>>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		arxiv & multi::archive_traits<Archive>::make_nvp("first", first_);
	//	arxiv &                  BOOST_SERIALIZATION_NVP(         first_);
	//	arxiv &                        cereal:: make_nvp("first", first_);
	//	arxiv &                               CEREAL_NVP(         first_);
	//	arxiv &                                                   first_ ;

		arxiv & multi::archive_traits<Archive>::make_nvp("last" , last_ );
	//	arxiv &                  BOOST_SERIALIZATION_NVP(         last_ );
	//	arxiv &                        cereal:: make_nvp("last" , last_ );
	//	arxiv &                               CEREAL_NVP(         last_ );
	//	arxiv &                                                   last_  ;
	}

	using value_type      = IndexType;
	using difference_type = decltype(IndexTypeLast{} - IndexType{});  // std::make_signed_t<value_type>;
	using size_type       = difference_type;
	using const_reference = value_type;
	using reference       = const_reference;
	using const_pointer   = value_type;
	using pointer         = value_type;

	range() = default;

	template<class Range, typename = std::enable_if_t<std::is_same_v<std::decay_t<Range>, value_type>> >
	// cxxcheck-suppress internalAstError ; because bug in cppcheck
	constexpr explicit range(Range&& other)
	: first_{std::forward<Range>(other).first()}, last_{std::forward<Range>(other).last()} {}

	constexpr range(IndexType first, IndexTypeLast last) noexcept : first_{first}, last_{last} {}
	[[deprecated]] constexpr explicit range(IndexType first) : range{first, first + 1} {}

	class const_iterator : public boost::multi::iterator_facade<
		const_iterator,
		value_type, std::random_access_iterator_tag,
		const_reference, difference_type
	> {
		typename const_iterator::value_type curr_;
		constexpr explicit const_iterator(value_type current) : curr_{current} {}
		friend class range;

	 public:
		const_iterator() = default;

		constexpr auto operator==(const_iterator const& other) const -> bool {return curr_ == other.curr_;}
		constexpr auto operator< (const_iterator const& other) const -> bool {return curr_ <  other.curr_;}

		constexpr auto operator++() -> const_iterator& {++curr_; return *this;}
		constexpr auto operator--() -> const_iterator& {--curr_; return *this;}

		constexpr auto operator-=(typename const_iterator::difference_type n) -> const_iterator& {curr_ -= n; return *this;}
		constexpr auto operator+=(typename const_iterator::difference_type n) -> const_iterator& {curr_ += n; return *this;}

		constexpr auto operator-(const_iterator const& other) const {return curr_ - other.curr_;}
		constexpr auto operator*() const -> typename const_iterator::reference {return curr_;}
	};

	using               iterator =                       const_iterator ;
	using       reverse_iterator = std::reverse_iterator<      iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	[[nodiscard]] constexpr auto first() const -> const_reference {return first_;}
	[[nodiscard]] constexpr auto last()  const -> const_reference {return last_ ;}

	constexpr auto operator[](difference_type n) const -> const_reference {return first() + n;}

	[[nodiscard]] constexpr auto front() const -> const_reference {return first()   ;}
	[[nodiscard]] constexpr auto back()  const -> const_reference {return last() - 1;}

	[[nodiscard]] constexpr auto cbegin() const {return const_iterator{first_};}
	[[nodiscard]] constexpr auto cend()   const {return const_iterator{last_ };}

	[[nodiscard]] constexpr auto rbegin() const {return reverse_iterator{end()  };}
	[[nodiscard]] constexpr auto rend()   const {return reverse_iterator{begin()};}

	[[nodiscard]] constexpr auto begin() const -> const_iterator {return cbegin();}
	[[nodiscard]] constexpr auto end()   const -> const_iterator {return cend()  ;}

	       constexpr auto is_empty()     const&       noexcept {return first_ == last_;}
	friend constexpr auto is_empty(range const& self) noexcept {return self.is_empty();}

	[[nodiscard]]
	       constexpr auto empty()     const&       noexcept {return is_empty();}
	friend constexpr auto empty(range const& self) noexcept {return self.empty();}

	       constexpr auto size()     const&       noexcept -> size_type {return last_ - first_;}
	friend constexpr auto size(range const& self) noexcept -> size_type {return self.size();}

	friend constexpr auto begin(range const& self) {return self.begin();}
	friend constexpr auto end  (range const& self) {return self.end()  ;}

	friend constexpr auto operator==(range const& self, range const& other) {
		return (self.empty() and other.empty()) or (self.first_ == other.first_ and self.last_ == other.last_);
	}
	friend constexpr auto operator!=(range const& self, range const& other) {return not(self == other);}

	[[nodiscard]] constexpr auto find(value_type const& value) const -> range::const_iterator {
		if(value >= last_ or value < first_) {
			return end();
		}
		return begin() + (value - front());
	}
	template<class Value> [[nodiscard]] constexpr auto contains(Value const& value) const {return (value >=first_) and (value < last_);}
	template<class Value> [[nodiscard]] constexpr auto count   (Value const& value) const -> value_type {return contains(value);}

	friend constexpr auto intersection(range const& self, range const& other) {
		using std::max; using std::min;
		auto new_first = max(self.first(), other.first());
		auto new_last  = min(self.last() , other.last() );
		new_first = min(new_first, new_last);
		return range<decltype(new_first), decltype(new_last)>{new_first, new_last};
	}
	[[nodiscard]] constexpr auto contains(value_type const& value) const {return value >= first_ and value < last_;}
};

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr auto make_range(IndexType first, IndexTypeLast last) -> range<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range {
	range<IndexType> impl_{std::numeric_limits<IndexType>::min(), std::numeric_limits<IndexType>::max()};

	constexpr intersecting_range() = default;  // MSVC 19.07 needs constexpr to initialize ALL later
	static constexpr auto make(IndexType first, IndexType last) -> intersecting_range {
		intersecting_range ret; ret.impl_ = range<IndexType>{first, last}; return ret;
	}
	friend constexpr auto intersection(intersecting_range const& self, range<IndexType> const& other) {
		return intersection(self.impl_, other);
	}
	friend constexpr auto intersection(range<IndexType> const& other, intersecting_range const& self) {
		return intersection(other, self.impl_);
	}
	friend constexpr auto operator<(intersecting_range const& self, IndexType end) {
		return intersecting_range::make(self.impl_.first(), end);
	}
	friend constexpr auto operator<=(IndexType first, intersecting_range const& self) {
		return intersecting_range::make(first, self.impl_.last());
	}

 public:
	constexpr auto operator*() const& -> intersecting_range const& {return *this;}
	static constexpr auto all() noexcept {return intersecting_range{};}
};

[[maybe_unused]] constexpr intersecting_range<> ALL   = intersecting_range<>::all();
[[maybe_unused]] constexpr intersecting_range<> _     = ALL;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> U     = ALL;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> ooo   = ALL;

[[maybe_unused]] constexpr intersecting_range<> V     = U;  // NOLINT(readability-identifier-length)
[[maybe_unused]] constexpr intersecting_range<> A     = V;  // NOLINT(readability-identifier-length)

//[[maybe_unused]] constexpr intersecting_range<> https://www.compart.com/en/unicode/U+2200 = V;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
struct extension_t : public range<IndexType, IndexTypeLast> {
	using range<IndexType, IndexTypeLast>::range;

	constexpr extension_t(IndexType first, IndexTypeLast last) noexcept
	: range<IndexType, IndexTypeLast>{first, last} {}

	// cppcheck-suppress noExplicitConstructor ; because syntax convenience // NOLINTNEXTLINE(runtime/explicit)
	constexpr extension_t(IndexType last) noexcept  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) because syntax convenience
	: range<IndexType, IndexTypeLast>(0, last) {}

	constexpr extension_t() noexcept : range<IndexType, IndexTypeLast>() {}

	friend constexpr auto size(extension_t const& self) -> typename extension_t::size_type {return self.size();}

//  template<class OStream>
//  friend auto operator<<(OStream& os, extension_t const& self) -> decltype(os<<"[]") {
//  	if(self.empty()) {
//  		return os << static_cast<range<IndexType> const&>(self);
//  	}
//  	if(self.first() == 0) {
//  		return os <<"["<< self.last() <<"]";
//  	}
//  	return os << static_cast<range<IndexType> const&>(self);
//  }

	[[nodiscard]] constexpr auto start () const -> IndexType {return this->first();}
	[[nodiscard]] constexpr auto finish() const -> IndexType {return this->last ();}

	friend constexpr auto operator==(extension_t const& self, extension_t const& other) {return static_cast<range<IndexType> const&>(self) == static_cast<range<IndexType> const&>(other);}
	friend constexpr auto operator!=(extension_t const& self, extension_t const& other) {return static_cast<range<IndexType> const&>(self) != static_cast<range<IndexType> const&>(other);}

	friend constexpr auto intersection(extension_t const& ex1, extension_t const& ex2) -> extension_t {
		using std::max; using std::min;
		auto       first = max(ex1.first(), ex2.first());
		auto const last  = min(ex1.last() , ex2.last() );
		first = min(first, last);
		return extension_t{first, last};
	}
};

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
constexpr auto make_extension_t(IndexType first, IndexTypeLast last) -> extension_t<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexTypeLast = std::ptrdiff_t>
constexpr auto make_extension_t(IndexTypeLast last) {return make_extension_t(IndexTypeLast{0}, last);}

using index_range     = range<index>;
using index_extension = extension_t<index>;
using iextension      = index_extension;
using irange          = index_range;

namespace detail {

template<typename, typename>
struct append_to_type_seq {};

template<typename T, typename... Ts, template<typename...> class TT>
struct append_to_type_seq<T, TT<Ts...> > {
    using type = TT<Ts..., T>;
};

template<typename T, dimensionality_type N, template<typename...> class TT>
struct repeat {
    using type = typename
        append_to_type_seq<
            T,
            typename repeat<T, N-1, TT>::type
        >::type;
};

template<typename T, template<typename...> class TT>
struct repeat<T, 0, TT> {
	using type = TT<>;
};

//template<class T, std::size_t N>
//constexpr auto array_size_impl(const std::array<T, N>&)
//  -> std::integral_constant<std::size_t, N>;

//template<class... T>
//constexpr auto array_size_impl(const std::tuple<T...>&)
//    -> std::integral_constant<std::size_t, std::tuple_size<std::tuple<T...>>{}>;

//template<class Array>
//using array_size = decltype(array_size_impl(std::declval<const Array&>()));

//template<class Array>
//constexpr auto static_size() -> std::decay_t<decltype(array_size<Array>::value)> {
//	return array_size<Array>::value;
//}
//template<class Array>
//constexpr auto static_size(Array const& /*unused*/) -> decltype(static_size<Array>()) {
//	return static_size<Array>();
//}

//// TODO(correaa) consolidate with tuple_tail defined somewhere else
//template<class Tuple>
//constexpr auto head(Tuple&& t)
//->decltype(std::get<0>(std::forward<Tuple>(t))) {
//	return std::get<0>(std::forward<Tuple>(t)); }

//template<typename Tuple, std::size_t... Ns>
//constexpr auto tail_impl(std::index_sequence<Ns...> /*012*/, [[maybe_unused]] Tuple&& t) {  // [[maybe_unused]] needed by icpc "error #869: parameter "t" was never referenced"
//	using boost::multi::detail::get;
//	return boost::multi::detail::tuple{get<Ns + 1U>(std::forward<Tuple>(t))...};
////  return make_tuple(std::get<Ns + 1U>(std::forward<Tuple>(t))...);
//}

//template<class Tuple>
//constexpr auto tail(Tuple const& t) {
//	return tail_impl(std::make_index_sequence<std::tuple_size_v<Tuple> - 1U>(), t);
//}

}  // end namespace detail

template<dimensionality_type D> using index_extensions = typename detail::repeat<index_extension, D, tuple>::type;

template<dimensionality_type D, class Tuple>
constexpr auto contains(index_extensions<D> const& iex, Tuple const& tup) {
//  using detail::head;
//  using detail::tail;
	return contains(head(iex), head(tup)) and contains(tail(iex), tail(tup));
}

}  // end namespace boost::multi
#endif
