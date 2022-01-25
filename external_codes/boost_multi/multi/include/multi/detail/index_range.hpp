// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_DETAIL_INDEX_RANGE_HPP
#define MULTI_DETAIL_INDEX_RANGE_HPP

#include "../config/MAYBE_UNUSED.hpp"
#include "../config/NODISCARD.hpp"

#include "../detail/serialization.hpp"

#include <algorithm>  // for min

#include <iostream>   // TODO(correaa) remove, add include in QMCP

#include <iterator>   // for std::random_iterator_tag // std::reverse_iterator
#include <limits>     // for numeric_limits
#include <utility>    // for forward

namespace boost {
namespace multi {

template<
	class Self,
	class ValueType, class AccessCategory, class Reference = ValueType&, class DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, class Pointer = ValueType*
>
class iterator_facade {
	using self_type = Self;
	NODISCARD("") constexpr auto self()      & -> self_type      & {return static_cast<self_type      &>(*this);}
	NODISCARD("") constexpr auto self() const& -> self_type const& {return static_cast<self_type const&>(*this);}

 public:
	using value_type        = ValueType;
	using reference         = Reference;
	using pointer           = Pointer;
	using difference_type   = DifferenceType;
	using iterator_category = AccessCategory;

	constexpr auto operator==(self_type const& o) const {return o==self();}
	constexpr auto operator!=(self_type const& o) const {return not(o==self());}

	       constexpr auto operator+(difference_type n) const -> self_type {self_type r = self(); r += n; return r;}
	       constexpr auto operator-(difference_type n) const -> self_type {self_type r = self(); r -= n; return r;}

	friend constexpr auto operator+(difference_type n, self_type const& s) -> self_type {return s + n;}

	friend constexpr auto operator++(self_type& s, int) -> self_type {self_type r = s; ++s; return r;}
	friend constexpr auto operator--(self_type& s, int) -> self_type {self_type r = s; --s; return r;}
};

template<typename IndexType = std::true_type, typename IndexTypeLast = IndexType>
class range {
	IndexType first_ = {};
	IndexTypeLast last_ = first_;

 public:
	template<class Ar>//, class ArT = multi::archive_traits<Ar>>
	void serialize(Ar& ar, unsigned /*version*/) {
		{
		ar & multi::archive_traits<Ar>::make_nvp("first", first_);
	//	ar & BOOST_SERIALIZATION_NVP(first_);
	//	ar &       cereal:: make_nvp("first", first_);
	//	ar &              CEREAL_NVP(first_);   ///ArT::make_nvp("first", first_);
	//	ar &            first_ ;   ///ArT::make_nvp("first", first_);
		}
		{
		ar & multi::archive_traits<Ar>::make_nvp("last", last_);
	//	ar &             BOOST_SERIALIZATION_NVP        (last_);
	//	ar &                   cereal:: make_nvp("last", last_);
	//	ar &                          CEREAL_NVP        (last_);
	//	ar &                                             last_ ;   // ArT::make_nvp("last" , last_ );
		}
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
	constexpr range(IndexType f, IndexTypeLast l) : first_{f}, last_{l} {}
	constexpr explicit range(IndexType f) : range{f, f + 1} {}

	class const_iterator
	: public boost::multi::iterator_facade<const_iterator,
		value_type, std::random_access_iterator_tag,
		const_reference, difference_type
	> {
		typename const_iterator::value_type curr_;
		constexpr auto dereference() const -> typename const_iterator::reference {return curr_;}
		constexpr void increment() {++curr_;}
		constexpr void decrement() {--curr_;}
		constexpr void advance(typename const_iterator::difference_type n) {curr_+=n;}
		constexpr auto equal(const_iterator const& y) const -> bool {return curr_ == y.curr_;}
		constexpr auto distance_to(const_iterator const& z) const -> difference_type {return z.curr_-curr_;}
		constexpr explicit const_iterator(value_type current) : curr_(current) {}
		friend class range;

	public:
		using difference_type = std::ptrdiff_t;
		const_iterator() = default;
		constexpr auto operator==(const_iterator const& y) const -> bool {return curr_ == y.curr_;}
		constexpr auto operator< (const_iterator const& y) const -> bool {return curr_ < y.curr_;}

		constexpr auto operator++() -> const_iterator& {++curr_; return *this;}
		constexpr auto operator--() -> const_iterator& {--curr_; return *this;}

		constexpr auto operator-=(typename const_iterator::difference_type n) -> const_iterator& {curr_-=n; return *this;}
		constexpr auto operator+=(typename const_iterator::difference_type n) -> const_iterator& {curr_+=n; return *this;}

		constexpr auto operator-(const_iterator const& y) const {return curr_ - y.curr_;}
		constexpr auto operator-(typename const_iterator::difference_type n) const -> const_iterator {return curr_ - n;}
		constexpr auto operator*() const -> typename const_iterator::reference {return curr_;}
		constexpr auto operator[](typename const_iterator::difference_type n) const{return *((*this)+n);}
	};

	using iterator = const_iterator;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	NODISCARD("") constexpr auto first() const -> const_reference {return first_;}
	NODISCARD("") constexpr auto last()  const -> const_reference {return last_;}

	constexpr auto operator[](difference_type p) const -> const_reference {return first() + p;}

	NODISCARD("") constexpr auto front() const -> const_reference {return first();}
	NODISCARD("") constexpr auto back()  const -> const_reference {return last() - 1;}

	NODISCARD("") constexpr auto cbegin() const {return const_iterator{first_};}
	NODISCARD("") constexpr auto cend()   const {return const_iterator{last_};}

	NODISCARD("") constexpr auto rbegin() const {return reverse_iterator{end()};}
	NODISCARD("") constexpr auto rend()   const {return reverse_iterator{begin()};}

	NODISCARD("") constexpr auto begin() const -> const_iterator {return cbegin();}
	NODISCARD("") constexpr auto end()   const -> const_iterator {return cend();}

	       NODISCARD("") constexpr auto is_empty()     const&    noexcept -> bool {return first_ == last_;}
	friend               constexpr auto is_empty(range const& s) noexcept -> bool {return s.is_empty();}

	       NODISCARD("") constexpr auto empty()     const&    noexcept -> bool{return is_empty();}
	friend               constexpr auto empty(range const& s) noexcept -> bool{return s.empty();}

	       NODISCARD("") constexpr auto size()     const&    noexcept -> size_type {return last_ - first_;}
	friend               constexpr auto size(range const& s) noexcept -> size_type {return s.size();}

//  friend auto operator<<(std::ostream& os, range const& s) -> std::ostream& {
//  	return s.empty()?os<<"[)":os <<"["<< s.first() <<", "<< s.last() <<")";
//  }
	friend constexpr auto begin(range const& self) -> const_iterator {return self.begin();}
	friend constexpr auto end  (range const& self) -> const_iterator {return self.end()  ;}

	friend constexpr auto operator==(range const& a, range const& b) -> bool{
		return (a.empty() and b.empty()) or (a.first_==b.first_ and a.last_==b.last_);
	}
	friend constexpr auto operator!=(range const& r1, range const& r2) -> bool{return not(r1 == r2);}

	NODISCARD("") constexpr auto find(value_type const& value) const -> range::const_iterator{
		if(value >= last_ or value < first_) {
			return end();
		}
		return begin() + (value - front());
	}
	template<class K> NODISCARD("") constexpr auto contains(K const& k) const -> bool {return (k>=first_) and (k<last_);}
	template<class K>               constexpr auto count   (K const& k) const -> value_type {return contains(k);}
	friend constexpr auto intersection(range const& r1, range const& r2) {
		using std::max; using std::min;
		auto new_first = max(r1.first(), r2.first());
		auto new_last  = min(r1.last() , r2.last() );
		new_first = min(new_first, new_last);
		return range<decltype(new_first), decltype(new_last)>{new_first, new_last};
	}
	NODISCARD("") constexpr auto contains(value_type const& v) const {return v>=first_ and v<last_;}
};

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr auto make_range(IndexType first, IndexTypeLast last) -> range<IndexType, IndexTypeLast>{
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range{
	range<IndexType> impl_{std::numeric_limits<IndexType>::min(), std::numeric_limits<IndexType>::max()};
	intersecting_range() = default;
	static constexpr auto make(IndexType first, IndexType last) -> intersecting_range{
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
	static constexpr auto all() -> intersecting_range {return {};}
};

MAYBE_UNUSED constexpr intersecting_range<> const ALL = intersecting_range<>::all();
MAYBE_UNUSED constexpr intersecting_range<> const _   = ALL;
MAYBE_UNUSED constexpr intersecting_range<> const U   = ALL;

[[deprecated]] constexpr intersecting_range<> const all = ALL;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
struct extension_t : public range<IndexType, IndexTypeLast>{
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

	NODISCARD("") constexpr auto start () const -> IndexType {return this->first();}
	NODISCARD("") constexpr auto finish() const -> IndexType {return this->last ();}
	friend constexpr auto operator==(extension_t const& a, extension_t const& b) {return static_cast<range<IndexType> const&>(a)==static_cast<range<IndexType> const&>(b);}
	friend constexpr auto operator!=(extension_t const& a, extension_t const& b) {return not(a==b);}
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

}  // end namespace multi
}  // end namespace boost
#endif
