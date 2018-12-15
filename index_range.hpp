#ifdef compile_instructions
(echo "#include\""$0"\"" > $0x.cpp) && time clang++ -O3 -std=c++14 -Wfatal-errors -Wall -D_TEST_MULTI_INDEX_RANGE $0x.cpp -o $0x.x && $0x.x $@ && rm -rf $0x.cpp $0x.x; exit
#endif

#ifndef MULTI_INDEX_RANGE_HPP
#define MULTI_INDEX_RANGE_HPP

//#include "types.hpp"

#include<cassert>
#include<iterator> // std::random_iterator_tag
#include<limits>

namespace boost{
namespace multi{

using index = std::ptrdiff_t;
using size_type = index;

class index_range{
protected:
	index first_;
	index last_;
public:
	using value_type = index;
	index_range() = default; // partially formed
	constexpr index_range(index fl) noexcept : first_(fl), last_(fl + 1){}
	constexpr index_range(index f, index l) : first_(f), last_(l){assert(first_ <= last_);}
	class const_iterator{// TODO use iterator_facade?
		index_range::value_type current_;
		friend class index_range;
		constexpr const_iterator(index current) : current_(current){}
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type = index_range::value_type;
		using difference_type = index;
		using pointer = void;
		using reference = index const&;
		reference operator*() const{return current_;}
		const_iterator& operator+=(std::ptrdiff_t s){current_ += s; return *this;}
		const_iterator& operator++(){++current_; return *this;}
		const_iterator operator+(std::ptrdiff_t d) const{return const_iterator{*this}+=d;}
		friend constexpr 
		bool operator==(const_iterator const& it1, const_iterator const& it2) noexcept{
			return it1.current_ == it2.current_;
		}
		friend constexpr 
		bool operator!=(const_iterator const& it1, const_iterator const& it2) noexcept{
			return not (it1 == it2);
		}
		friend constexpr 
		difference_type operator-(const_iterator const& it1, const_iterator const& it2) noexcept{
			return it1.current_ - it2.current_;
		}
	};
//	index operator[](std::ptrdiff_t p) const{return first_ + p;}
	constexpr const_iterator begin() const noexcept{return const_iterator{first_};}
	constexpr const_iterator end() const noexcept{return const_iterator{last_};}
	constexpr index front() const noexcept{return first_;}
	constexpr index back() const noexcept{return last_ - 1;}
	constexpr index first() const noexcept{return first_;}
	constexpr index last() const noexcept{return last_;}
	constexpr bool empty() const noexcept{return first_ == last_;}
	friend constexpr bool empty(index_range const& s){return s.empty();}
	constexpr size_type size() const noexcept{return last_ - first_;}
	friend constexpr size_type size(index_range const& s){return s.size();}
	friend std::ostream& operator<<(std::ostream& os, index_range const& self){
		if(self.empty()) return os <<"[)\n";
		return os <<'['<< self.first_ <<", "<< self.last_ <<')';
	}
	friend const_iterator begin(index_range const& self){return self.begin();}
	friend const_iterator end(index_range const& self){return self.end();}
	friend constexpr 
	bool operator==(index_range const& ir1, index_range const& ir2){
		return 
			(ir1.empty() and ir2.empty()) or 
			(ir1.first_ == ir2.first_ and ir1.last_ == ir2.last_)
		;
	}
	friend constexpr 
	bool operator!=(index_range const& ir1, index_range const& ir2){
		return not (ir1 == ir2);
	}
	index_range::const_iterator find(index value) const{
		auto first = begin();
		if(value >= last_ or value < first_) return end();
		first += value - *first;
		return first;		
	}
};

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

class index_extension : public index_range{
	using index_range::index_range;
	public:
	constexpr index_extension() noexcept : index_range(0, 0){}
	constexpr index_extension(index zerobased) noexcept : index_range(0, zerobased){}
	friend std::ostream& operator<<(std::ostream& os, index_extension const& self){
		if(self.empty()) return os << "[)";
		if(self.first() == 0) return os << "[" << self.last() << "]";
		return os << '[' << self.first() << ", " << self.last() << ')';
	}
};

}}

#ifdef _TEST_MULTI_INDEX_RANGE

#include <boost/spirit/include/karma.hpp>

#include<algorithm>
#include<cassert>
#include<iostream>
#include<vector>

namespace multi = boost::multi;
auto index_ = boost::spirit::karma::int_generator<multi::index>{};

using std::cout;
using std::cerr;

int main(){

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
}

#endif
#endif
