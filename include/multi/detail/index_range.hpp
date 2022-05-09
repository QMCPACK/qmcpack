// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_INDEX_RANGE_HPP
#define MULTI_DETAIL_INDEX_RANGE_HPP

//#include "../config/MAYBE_UNUSED.hpp"
//#include "../config/NODISCARD.hpp"

#include "../detail/serialization.hpp"
#include "../detail/tuple_zip.hpp"
#include "../detail/types.hpp"

#include <algorithm>  // for min

#include <iostream>   // TODO(correaa) remove, add include in QMCP

#include <iterator>   // for std::random_iterator_tag // std::reverse_iterator
#include <limits>     // for numeric_limits
#include <utility>    // for forward

namespace boost::multi {

//template<class... Ts>
//using tuple = std::tuple<Ts...>;

using boost::multi::detail::tuple;
//using std::tuple;
//using std::make_tuple;
//using std::tuple_cat;

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

	friend constexpr auto operator!=(self_type const& s, self_type const& o) {return not(s == o);}

	friend constexpr auto operator<=(self_type const& s, self_type const& o) {return (s < o) or (s == o);}
	friend constexpr auto operator> (self_type const& s, self_type const& o) {return not(s <= o);}
	friend constexpr auto operator>=(self_type const& s, self_type const& o) {return not(s < o);}

	       constexpr auto operator-(difference_type n) const {return self_type{self()} -= n;}
	       constexpr auto operator+(difference_type n) const {return self_type{self()} += n;}
	friend constexpr auto operator+(difference_type n, self_type const& s) {return s + n;}

	friend constexpr auto operator++(self_type& s, int) -> self_type {self_type r = s; ++s; return r;}
	friend constexpr auto operator--(self_type& s, int) -> self_type {self_type r = s; --s; return r;}

	constexpr auto operator[](difference_type n) const {return *(self() + n);}
};

template<typename IndexType = std::true_type, typename IndexTypeLast = IndexType>
class range {
	IndexType first_ = {};
	IndexTypeLast last_ = first_;

 public:
	template<class Ar>//, class ArT = multi::archive_traits<Ar>>
	void serialize(Ar& ar, unsigned /*version*/) {
		ar & multi::archive_traits<Ar>::make_nvp("first", first_);
	//	ar &             BOOST_SERIALIZATION_NVP(         first_);
	//	ar &                   cereal:: make_nvp("first", first_);
	//	ar &                          CEREAL_NVP(         first_);
	//	ar &                                              first_ ;

		ar & multi::archive_traits<Ar>::make_nvp("last" , last_ );
	//	ar &             BOOST_SERIALIZATION_NVP(         last_ );
	//	ar &                   cereal:: make_nvp("last" , last_ );
	//	ar &                          CEREAL_NVP(         last_ );
	//	ar &                                              last_  ;
	}

	using value_type      = IndexType;
	using difference_type = decltype(IndexTypeLast{} - IndexType{});  // std::make_signed_t<value_type>;
	using size_type       = difference_type;
	using const_reference = value_type;
	using reference       = const_reference;
	using const_pointer   = value_type;
	using pointer         = value_type;

	range() = default;

	template<class Range, typename = std::enable_if_t<std::is_same<std::decay_t<Range>, value_type>{}> >
	// cxxcheck-suppress internalAstError ; because bug in cppcheck
	constexpr explicit range(Range&& o) : first_{std::forward<Range>(o).first()}, last_{std::forward<Range>(o).last()} {}
	constexpr range(IndexType f, IndexTypeLast l) noexcept : first_{f}, last_{l} {}
	constexpr explicit range(IndexType f) : range{f, f + 1} {}

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

		constexpr auto operator==(const_iterator const& y) const -> bool {return curr_ == y.curr_;}
		constexpr auto operator< (const_iterator const& y) const -> bool {return curr_ <  y.curr_;}

		constexpr auto operator++() -> const_iterator& {++curr_; return *this;}
		constexpr auto operator--() -> const_iterator& {--curr_; return *this;}

		constexpr auto operator-=(typename const_iterator::difference_type n) -> const_iterator& {curr_ -= n; return *this;}
		constexpr auto operator+=(typename const_iterator::difference_type n) -> const_iterator& {curr_ += n; return *this;}

		constexpr auto operator-(const_iterator const& y) const {return curr_ - y.curr_;}
		constexpr auto operator*() const -> typename const_iterator::reference {return curr_;}
	};

	using               iterator =                       const_iterator ;
	using       reverse_iterator = std::reverse_iterator<      iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	[[nodiscard]] constexpr auto first() const -> const_reference {return first_;}
	[[nodiscard]] constexpr auto last()  const -> const_reference {return last_ ;}

	constexpr auto operator[](difference_type p) const -> const_reference {return first() + p;}

	[[nodiscard]] constexpr auto front() const -> const_reference {return first()   ;}
	[[nodiscard]] constexpr auto back()  const -> const_reference {return last() - 1;}

	[[nodiscard]] constexpr auto cbegin() const {return const_iterator{first_};}
	[[nodiscard]] constexpr auto cend()   const {return const_iterator{last_ };}

	[[nodiscard]] constexpr auto rbegin() const {return reverse_iterator{end()  };}
	[[nodiscard]] constexpr auto rend()   const {return reverse_iterator{begin()};}

	[[nodiscard]] constexpr auto begin() const -> const_iterator {return cbegin();}
	[[nodiscard]] constexpr auto end()   const -> const_iterator {return cend()  ;}

	       constexpr auto is_empty()     const&    noexcept {return first_ == last_;}
	friend constexpr auto is_empty(range const& s) noexcept {return s.is_empty();}

	[[nodiscard]]
	       constexpr auto empty()     const&    noexcept {return is_empty();}
	friend constexpr auto empty(range const& s) noexcept {return s.empty();}

	       constexpr auto size()     const&    noexcept -> size_type {return last_ - first_;}
	friend constexpr auto size(range const& s) noexcept -> size_type {return s.size();}

	friend constexpr auto begin(range const& self) {return self.begin();}
	friend constexpr auto end  (range const& self) {return self.end()  ;}

	friend constexpr auto operator==(range const& a, range const& b) {
		return (a.empty() and b.empty()) or (a.first_==b.first_ and a.last_==b.last_);
	}
	friend constexpr auto operator!=(range const& a, range const& b) {return not(a == b);}

	[[nodiscard]] constexpr auto find(value_type const& value) const -> range::const_iterator {
		if(value >= last_ or value < first_) {
			return end();
		}
		return begin() + (value - front());
	}
	template<class K> [[nodiscard]] constexpr auto contains(K const& k) const {return (k>=first_) and (k<last_);}
	template<class K>               constexpr auto count   (K const& k) const -> value_type {return contains(k);}

	friend constexpr auto intersection(range const& a, range const& b) {
		using std::max; using std::min;
		auto new_first = max(a.first(), b.first());
		auto new_last  = min(a.last() , b.last() );
		new_first = min(new_first, new_last);
		return range<decltype(new_first), decltype(new_last)>{new_first, new_last};
	}
	[[nodiscard]] constexpr auto contains(value_type const& v) const {return v >= first_ and v < last_;}
};

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr auto make_range(IndexType first, IndexTypeLast last) -> range<IndexType, IndexTypeLast> {
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range {
	range<IndexType> impl_{std::numeric_limits<IndexType>::min(), std::numeric_limits<IndexType>::max()};
	intersecting_range() = default;
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

[[maybe_unused]] constexpr intersecting_range<> const ALL   = intersecting_range<>::all();
[[maybe_unused]] constexpr intersecting_range<> const _     = ALL;
[[maybe_unused]] constexpr intersecting_range<> const U     = ALL;
[[maybe_unused]] constexpr intersecting_range<> const ooo   = ALL;

[[maybe_unused]] constexpr intersecting_range<> const V     = U;
[[maybe_unused]] constexpr intersecting_range<> const A     = V;
//  [[maybe_unused]] constexpr intersecting_range<> const ∀      = V;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
struct extension_t : public range<IndexType, IndexTypeLast> {
	using range<IndexType, IndexTypeLast>::range;

	constexpr extension_t(IndexType f, IndexTypeLast l) noexcept : range<IndexType, IndexTypeLast>{f, l} {}

	// cppcheck-suppress noExplicitConstructor ; because syntax convenience // NOLINTNEXTLINE(runtime/explicit)
	constexpr extension_t(IndexType last) noexcept : range<IndexType, IndexTypeLast>(0, last) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) because syntax convenience
	constexpr extension_t() noexcept : range<IndexType, IndexTypeLast>() {}

	friend constexpr auto size(extension_t const& s) -> typename extension_t::size_type {return s.size();}

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

	friend constexpr auto operator==(extension_t const& a, extension_t const& b) {return static_cast<range<IndexType> const&>(a) == static_cast<range<IndexType> const&>(b);}
	friend constexpr auto operator!=(extension_t const& a, extension_t const& b) {return static_cast<range<IndexType> const&>(a) != static_cast<range<IndexType> const&>(b);}

	friend constexpr auto intersection(extension_t const& r1, extension_t const& r2) -> extension_t {
		using std::max; using std::min;
		auto       first = max(r1.first(), r2.first());
		auto const last  = min(r1.last() , r2.last() );
		first = min(first, last);
		return extension_t{first, last};
	}
};

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
constexpr auto make_extension_t(IndexType f, IndexTypeLast l) -> extension_t<IndexType, IndexTypeLast> {return {f, l};}

template<class IndexTypeLast = std::ptrdiff_t>
constexpr auto make_extension_t(IndexTypeLast l) {return make_extension_t(IndexTypeLast{0}, l);}

using index_range     = range<index>;
using index_extension = extension_t<index>;
using iextension      = index_extension;
using irange          = index_range;

namespace detail {

template<typename, typename>
struct append_to_type_seq{};

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
constexpr auto contains(index_extensions<D> const& ie, Tuple const& tp) {
//  using detail::head;
//  using detail::tail;
	return contains(head(ie), head(tp)) and contains(tail(ie), tail(tp));
}

}  // end namespace boost::multi
#endif
