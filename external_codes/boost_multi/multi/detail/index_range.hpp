#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x &&$0x&&rm $0x;exit
#endif

#ifndef MULTI_INDEX_RANGE_HPP
#define MULTI_INDEX_RANGE_HPP

#include "../config/MAYBE_UNUSED.hpp"

#include<limits> // numeric_limits
#include<iterator> // std::random_iterator_tag // std::reverse_iterator
#include<iostream>

#if 0
//#include<boost/serialization/nvp.hpp>
namespace boost{
namespace serialization{
	template<class> struct nvp;
	template<class T> const nvp<T> make_nvp(char const* name, T& t) 
#if defined(BOOST_VERSION) and (BOOST_VERSION > 107100)
noexcept
#endif
	;
	template<class T> class array_wrapper;
	template<class T, class S> const array_wrapper<T> make_array(T* t, S s);

	template<class Archive>
	struct archive_traits{
		template<class T>
		static decltype(auto) make_nvp(char const* name, T&& t){
			return boost::serialization::make_nvp(name, std::forward<T>(t));
		}
		template<class P1, class P2>
		static decltype(auto) make_array(P1&& p1, P2&& p2){
			return boost::serialization::make_array(std::forward<P1>(p1), std::forward<P2>(p2));
		}
	};
}}
#endif

namespace boost{

namespace multi{

template<class Self, typename ValueType, class AccessCategory, typename Reference = ValueType&,  typename DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, typename Pointer = ValueType*>
class iterator_facade{
	using self_type = Self;
	constexpr self_type&       self()      {return static_cast<self_type&      >(*this);}
	constexpr self_type const& self() const{return static_cast<self_type const&>(*this);}
public:
	using value_type        = ValueType;
	using reference         = Reference;
	using pointer           = Pointer;
	using difference_type   = DifferenceType;
	using iterator_category = AccessCategory;
	constexpr auto operator==(self_type const& o) const{return o==self();}
	constexpr auto operator!=(self_type const& o) const{return not(o==self());}
	       constexpr self_type operator+(difference_type n) const{self_type r = self(); r += n; return r;}
	       constexpr self_type operator-(difference_type n) const{self_type r = self(); r -= n; return r;}
	friend constexpr self_type operator+(difference_type n, self_type const& s){return s + n;}
	friend constexpr self_type operator++(self_type& s, int){self_type r = s; ++s; return r;}
	friend constexpr self_type operator--(self_type& s, int){self_type r = s; --s; return r;}
};

//class iterator_core_access{};
}

namespace multi{

template<class T> struct archive_traits;

template<typename IndexType = std::true_type, typename IndexTypeLast = IndexType>
class range{
	IndexType first_ = {};
	IndexTypeLast last_ = first_;
public:
	template<class Archive>
	void serialize(Archive& ar, unsigned){
//		ar & boost::serialization::make_nvp("first", first_);//BOOST_SERIALIZATION_NVP(first);
//		ar & boost::serialization::make_nvp("last", last_);//BOOST_SERIALIZATION_NVP(last);
		ar & multi::archive_traits<std::decay_t<Archive>>::make_nvp("first", first_);//BOOST_SERIALIZATION_NVP(first); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if you get an error here you need to include the adators/serialization/xml_archive.hpp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
		ar & multi::archive_traits<std::decay_t<Archive>>::make_nvp("last", last_);//BOOST_SERIALIZATION_NVP(last);
	}
	using value_type      = IndexType;
	using difference_type = decltype(IndexTypeLast{} - IndexType{});// std::make_signed_t<value_type>;
	using size_type       = difference_type;
	using const_reference = value_type const;
	using reference       = const_reference;
	using const_pointer   = value_type;
	using pointer         = value_type;
	range() = default;
	template<class Range, typename = std::enable_if_t<std::is_same<std::decay_t<Range>, value_type>{}> >
	// cppcheck-suppress internalAstError ; because bug in cppcheck
	constexpr explicit range(Range&& o) : first_{std::forward<Range>(o).first()}, last_{std::forward<Range>(o).last()}{}
//	constexpr range(value_type const& fl) : first_{fl}, last_{fl + 1}{}
//	constexpr range(value_type f, value_type l) : first_{f}, last_{l}{}
	constexpr range(IndexType f, IndexTypeLast l) : first_{f}, last_{l}{}
	constexpr explicit range(IndexType f) : range{f, f + 1}{}
	class const_iterator 
		: public boost::multi::iterator_facade<const_iterator, 
			value_type, std::random_access_iterator_tag, 
			const_reference, difference_type
		>
	{
		typename const_iterator::value_type curr_;
		typename const_iterator::reference dereference() const{return curr_;}
		constexpr void increment(){++curr_;} void decrement(){--curr_;}
		constexpr void advance(typename const_iterator::difference_type n){curr_+=n;}
		constexpr bool equal(const_iterator const& y) const{return curr_ == y.curr_;}
		constexpr auto distance_to(const_iterator const& z) const{return z.curr_-curr_;}
		constexpr explicit const_iterator(value_type current) : curr_(current){}
		friend class range;
	 //   friend class boost::iterator_core_access;
	public:
		using difference_type = std::ptrdiff_t;
		const_iterator() = default;
		constexpr auto operator==(const_iterator const& y) const{return curr_ == y.curr_;}
		constexpr auto operator<(const_iterator const& y) const{return curr_ < y.curr_;}
		constexpr const_iterator& operator++(){++curr_; return *this;}
		constexpr const_iterator& operator--(){--curr_; return *this;}
		constexpr const_iterator& operator-=(typename const_iterator::difference_type n){curr_-=n; return *this;}
		constexpr const_iterator& operator+=(typename const_iterator::difference_type n){curr_+=n; return *this;}
		constexpr auto operator-(const_iterator const& y) const{return curr_ - y.curr_;}
		constexpr const_iterator operator-(typename const_iterator::difference_type n) const{return curr_ - n;}
		constexpr typename const_iterator::reference operator*() const{return curr_;}
		constexpr auto operator[](typename const_iterator::difference_type n) const{return *((*this)+n);}
	};
	using iterator = const_iterator;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	constexpr const_reference first() const{return first_;}
	constexpr const_reference last()  const{return last_;}
	constexpr const_reference operator[](difference_type p) const{return first() + p;}
	constexpr const_reference front() const{return first();}
	constexpr const_reference back()  const{return last() - 1;}
	constexpr const_iterator cbegin() const{return const_iterator{first_};}
	constexpr const_iterator cend()   const{return const_iterator{last_};}
	constexpr reverse_iterator rbegin() const{return reverse_iterator{end()};}
	constexpr reverse_iterator rend() const{return reverse_iterator{begin()};}
	constexpr const_iterator begin() const{return cbegin();}
	constexpr const_iterator end() const{return cend();}
	constexpr bool empty() const{return first_ == last_;}
	friend constexpr bool empty(range const& s){return s.empty();}
	constexpr size_type size() const noexcept{return last_ - first_;}
	friend constexpr size_type size(range const& s){return s.size();}
	friend std::ostream& operator<<(std::ostream& os, range const& s){
		return s.empty()?os<<"[)":os <<"["<< s.first() <<", "<< s.last() <<")";
	}
	friend constexpr const_iterator begin(range const& self){return self.begin();}
	friend constexpr const_iterator end(range const& self){return self.end();}
//	constexpr range& operator=(range const&) = default;
	friend constexpr auto operator==(range const& a, range const& b){
		return (a.empty() and b.empty()) or (a.first_==b.first_ and a.last_==b.last_);
	}
	friend constexpr bool operator!=(range const& r1, range const& r2){return not(r1 == r2);}
	constexpr range::const_iterator find(value_type const& value) const{
		if(value >= last_ or value < first_) return end();
		return begin() + (value - front());
	}
	template<class K>
	constexpr bool contains(K const& k) const{return (k>=first_) and (k<last_);}
	template<class K>
	constexpr auto count(K const& k) const{return contains(k);}
	friend constexpr auto intersection(range const& r1, range const& r2){
		using std::max; using std::min;
		auto f = max(r1.first(), r2.first()); 
		auto l = min(r1.last() , r2.last() );
		f = min(f, l);
		return range<decltype(f), decltype(l)>{f, l};
	}
	constexpr auto contains(value_type const& v) const{return v>=first_ and v<last_;}//?true:false;}
};

template<class IndexType = std::true_type, typename IndexTypeLast = IndexType>
constexpr range<IndexType, IndexTypeLast> make_range(IndexType first, IndexTypeLast last){
	return {first, last};
}

template<class IndexType = std::ptrdiff_t>
class intersecting_range{
	range<IndexType> impl_{std::numeric_limits<IndexType>::min(), std::numeric_limits<IndexType>::max()};
	intersecting_range() = default;
//	constexpr intersecting_range(IndexType first, IndexType last) = delete;//: impl_{first, last}{}
	static constexpr intersecting_range make(IndexType first, IndexType last){
		intersecting_range ret; ret.impl_ = range<IndexType>{first, last}; return ret;
	}
	friend constexpr auto intersection(intersecting_range const& self, range<IndexType> const& other){
		return intersection(self.impl_, other);
	}
	friend constexpr auto intersection(range<IndexType> const& other, intersecting_range const& self){
		return intersection(other, self.impl_);
	}
	friend constexpr auto operator<(intersecting_range const& self, IndexType end){
		return intersecting_range::make(self.impl_.first(), end);
	}
	friend constexpr auto operator<=(IndexType first, intersecting_range const& self){
		return intersecting_range::make(first, self.impl_.last());
	}
public:
	constexpr intersecting_range const& operator*() const&{return *this;}
	static constexpr intersecting_range all(){return {};}
};

MAYBE_UNUSED constexpr intersecting_range<> const all = intersecting_range<>::all();
MAYBE_UNUSED constexpr intersecting_range<> const _   = all;
MAYBE_UNUSED constexpr intersecting_range<> const __  = all;
MAYBE_UNUSED constexpr intersecting_range<> const U   = all;

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
struct extension_t : public range<IndexType, IndexTypeLast>{
	using range<IndexType, IndexTypeLast>::range;
	constexpr extension_t(IndexType f, IndexTypeLast l) noexcept : range<IndexType, IndexTypeLast>{f, l}{}
	// cppcheck-suppress noExplicitConstructor
	constexpr extension_t(IndexType last) noexcept : range<IndexType, IndexTypeLast>(0, last){}
	constexpr extension_t() noexcept : range<IndexType, IndexTypeLast>(){}
	friend constexpr typename extension_t::size_type size(extension_t const& s){return s.size();}
	friend std::ostream& operator<<(std::ostream& os, extension_t const& self){
		if(self.empty()) return os << static_cast<range<IndexType> const&>(self);
		if(self.first() == 0) return os <<"["<< self.last() <<"]";
		return os << static_cast<range<IndexType> const&>(self);
	}
	constexpr IndexType start () const{return this->first();}
	constexpr IndexType finish() const{return this->last ();}
	friend constexpr auto operator==(extension_t const& a, extension_t const& b){return static_cast<range<IndexType> const&>(a)==static_cast<range<IndexType> const&>(b);}
	friend constexpr auto operator!=(extension_t const& a, extension_t const& b){return not(a==b);}
	friend constexpr extension_t intersection(extension_t const& r1, extension_t const& r2){
		using std::max; using std::min;
		auto f = max(r1.first(), r2.first()); 
		auto l = min(r1.last() , r2.last() );
		f = min(f, l);
		return extension_t{f, l};
	}
};

template<class IndexType = std::ptrdiff_t, class IndexTypeLast = decltype(std::declval<IndexType>() + 1)>
constexpr extension_t<IndexType, IndexTypeLast> make_extension_t(IndexType f, IndexTypeLast l){return {f, l};}

template<class IndexTypeLast = std::ptrdiff_t>
constexpr auto make_extension_t(IndexTypeLast l){return make_extension_t(IndexTypeLast{0}, l);}

}}

////////////////////////////////////////////////////////////////////////////////

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

int main(int, char*[]){

}

#endif
#endif

