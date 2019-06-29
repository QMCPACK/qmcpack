#ifdef compile_instructions
(echo "#include\""$0"\"" > $0x.cpp) && time clang++ -O3 -std=c++14 -Wfatal-errors -Wall -D_TEST_MULTI_INDEX_RANGE $0x.cpp -o $0x.x && $0x.x $@ && rm -rf $0x.cpp $0x.x; exit
#endif

#ifndef MULTI_INDEX_RANGE_HPP
#define MULTI_INDEX_RANGE_HPP

//#include<boost/iterator/iterator_facade.hpp>

#include<iterator> // std::random_iterator_tag // std::reverse_iterator

namespace boost{

namespace multi{

template<class Self, typename ValueType, class AccessCategory, typename Reference = ValueType&,  typename DifferenceType = typename std::pointer_traits<ValueType*>::difference_type, typename Pointer = ValueType*>
class iterator_facade{
	using self_type = Self;
	self_type& self(){return *this;}
	self_type const& self() const{return static_cast<Self const&>(*this);}
public:
	using value_type = ValueType;
	using reference = Reference;
	using pointer = Pointer;
	using difference_type = DifferenceType;
	using iterator_category = AccessCategory;
	auto operator!=(self_type const& o) const
//	->decltype(not(o == self()))
	{	return not(o == self());}
//	Self& operator++(){return ++self(); return *this;}
	self_type operator+(difference_type n) const{self_type r = self(); r += n; return r;}
	self_type operator-(difference_type n) const{self_type r = self(); r -= n; return r;}
};

//class iterator_core_access{};
}

namespace multi{

template<class IndexType>
class range{
	IndexType first_;
	IndexType last_;
public:
	using value_type = IndexType;
	using difference_type = std::make_signed_t<value_type>;
	using size_type = difference_type;
	using const_reference = value_type const /*&*/;
	using reference = const_reference;
	using const_pointer = value_type;
	using pointer = value_type;
	range() : first_{}, last_{first_}{}
	template<class Range, typename = std::enable_if_t<std::is_same<std::decay_t<Range>, value_type>{}> >
	constexpr range(Range&& o) : first_(o.first()), last_(o.last()){}
	constexpr range(value_type fl) : first_{fl}, last_{fl + 1}{}
	constexpr range(value_type f, value_type l) : first_{f}, last_{l}{}
	class const_iterator 
		: public boost::multi::iterator_facade<const_iterator, 
			value_type, std::random_access_iterator_tag, 
			const_reference, difference_type
		>
	{
		typename const_iterator::value_type curr_;
		typename const_iterator::reference dereference() const{return curr_;}
		void increment(){++curr_;} void decrement(){--curr_;}
		void advance(typename const_iterator::difference_type n){curr_+=n;}
		bool equal(const_iterator const& y) const{return curr_ == y.curr_;}
		auto distance_to(const_iterator const& z) const{return z.curr_-curr_;}
		constexpr const_iterator(value_type current) : curr_(current){}
		friend class range;
	 //   friend class boost::iterator_core_access;
	public:
		auto operator==(const_iterator const& y) const{return curr_ == y.curr_;}
		const_iterator& operator++(){++curr_; return *this;}
		const_iterator& operator--(){--curr_; return *this;}
		typename const_iterator::reference operator*() const{return curr_;}
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
	friend const_iterator begin(range const& self){return self.begin();}
	friend const_iterator end(range const& self){return self.end();}
	range& operator=(range const&) = default;
	friend constexpr bool operator==(range const& a, range const& b){
		return(a.empty()&& b.empty())||(a.first_==b.first_ && a.last_==b.last_);
	}
	friend constexpr 
	bool operator!=(range const& r1, range const& r2){return not(r1 == r2);}
	size_type count(value_type const& value) const{
		if(value >= last_ or value < first_) return 0;
		return 1;
	}
	range::const_iterator find(value_type const& value) const{
		auto first = begin();
		if(value >= last_ or value < first_) return end();
		return first += value - *first;
	}
	friend range intersection(range const& r1, range const& r2){
		using std::max; using std::min;
		auto first = max(r1.first(), r2.first()); 
		auto last = min(r1.last(), r2.last());
		if(first < last) return {first, last};
		return {};
	}
	bool contains(value_type const& v) const{return (v >= first() and v < last())?true:false;}
};

//using index_range = range<index>;
/*
class strided_index_range : index_range{	
	index stride_;
public:
	strided_index_range(index first, index last, index stride = 1) : index_range{first, last}, stride_{stride}{}
//	explicit operator index_range() const{return *this;}
	index stride() const{return stride_;}
	using index_range::front;
	index back() const{return front() + size()*stride();}
	size_type size() const{return (this->last_ - first_) / stride_;}
	friend std::ostream& operator<<(std::ostream& os, strided_index_range const& self){
		if(empty() 
		if(self.first_ == self.last_) return os << "[)" << '\n';
		return os << '[' << self.first_ << ", " << self.last_ << ')';
	}
};*/

template<class IndexType>
class extension_t : public range<IndexType>{
	using range<IndexType>::range;
	public:
	constexpr extension_t() noexcept : range<IndexType>(0, 0){}
	constexpr extension_t(typename extension_t::value_type last) noexcept : range<IndexType>(0, last){}
	friend constexpr typename extension_t::size_type size(extension_t const& s){return s.size();}
	friend std::ostream& operator<<(std::ostream& os, extension_t const& self){
		if(self.empty()) return os << static_cast<range<IndexType> const&>(self);
		if(self.first() == 0) return os <<"["<< self.last() <<"]";
		return os << static_cast<range<IndexType> const&>(self);
	}
	friend bool operator==(extension_t const& a, extension_t const& b){return static_cast<range<IndexType> const&>(a)==static_cast<range<IndexType> const&>(b);}
	friend bool operator!=(extension_t const& a, extension_t const& b){return not(a==b);}
};

}}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef _TEST_MULTI_INDEX_RANGE

//#include <boost/spirit/include/karma.hpp>

#include<cassert>
#include<iostream>

namespace multi = boost::multi;

using std::cout;
using std::cerr;

int main(){

	cout << multi::range<int>{5, 5} <<'\n';
	cout << multi::extension_t<int>{5} <<'\n';
	multi::extension_t<int> ee{5};
	assert( ee == 5 );
	assert( multi::extension_t<int>{5} == 5 );
	cout << multi::extension_t<int>{5, 7} <<'\n';
	
	assert(( multi::extension_t<int>{5, 12}.count(10) == 1 ));

	assert( size(multi::range<int>{5, 5}) == 0 );
	assert( empty(multi::range<int>{5, 5}) );
	
	assert( size(multi::range<int>{}) == 0 );
	assert( empty(multi::range<int>{}) );
	
	for(auto const& i : multi::range<int>{5, 12}) cout<< i <<' ';
	cout <<'\n';
	
	cout << intersection(multi::range<int>{5, 12}, multi::range<int>{14, 16}) << '\n';
//	for(auto const& i : intersection(multi::range<int>{5, 12}, multi::range<int>{8, 16})) cout<< i <<' ';
//	cout <<'\n';
	
	multi::range<int> rr{5, 12};
	assert( rr.contains(6) );
	assert( not rr.contains(12) );
	for(auto it = rr.begin(); it != rr.end(); ++it) cout<< *it <<' ';
	cout<<'\n';
	for(auto it = rr.rbegin(); it != rr.rend(); ++it) cout<< *it <<' ';
	cout<<'\n';

//	cout<< *rr.rbegin() <<'\n';
//	for(auto it = rr.rbegin(); it != rr.rend(); ++it) cout<< *it <<' ';
//	cout <<'\n';
	
//	multi::extension<int> ei{5, 10}; 
#if 0
	std::iterator_traits<multi::range<multi::index>::const_iterator>::value_type p = multi::index{4};

	{
		multi::index_range ir{5, 10};
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		std::vector<multi::index_range::value_type> v(5);
		copy(begin(ir), end(ir), begin(v));
		assert(v[0] == 5);
		for(auto& i : ir) cout << i << ' ';
		cout << '\n';
		auto f = ir.find(6);
		cerr << "*f " << *f << '\n';
		assert(*f == 6);
		using std::find;
		auto f2 = find(ir.begin(), ir.end(), 12);
		assert(f2 == ir.end());
		auto f3 = find(ir.begin(), ir.end(), 2);
		assert(f3 == ir.end());
	}
/*	{
		multi::strided_index_range ir{6, 12, 2};
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		std::vector<multi::index_range::value_type> v(5);
		copy(begin(ir), end(ir), begin(v));
		assert( v[0] == 6 );
		assert( v[1] == 8 );
		for(auto& i : ir) cout << i <<' ';
		cout <<'\n';
	}*/
	{
		multi::index_range ir(5);
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		assert(*begin(ir) == 5);
		assert(ir.front() == 5);
		assert(ir.back() == 5);
	}
	{
		multi::index_range ir; // partially formed
		ir = multi::index_range{8, 8};
		assert(ir.empty());
	}
	{
		multi::index_range ir = {};
		assert(ir.empty());
	}
	{
		multi::index_extension ie(5);
		cout << ie << " = {" << format(index_ % ", ", ie) << "}";
	}
#endif
}

#endif
#endif

