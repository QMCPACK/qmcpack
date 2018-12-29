#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ -O3 `#-fconcepts` -std=c++14 -Wall -Wextra `#-fmax-errors=2` -Wfatal-errors -lboost_timer -I${HOME}/prj -D_TEST_BOOST_MULTI_ARRAY_REF $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP

//#include "../multi/index_range.hpp"
//#include "../multi/ordering.hpp"

#include "./index_range.hpp"
#include "./ordering.hpp"

#include<iostream>
#include<algorithm> // transform
#include<array>
#include<cassert>
#include<memory>
#include<vector>

namespace boost{
namespace multi{

using difference_type = std::ptrdiff_t;
using size_type = index;
using index = difference_type;

}}

namespace std{
template<class T>
boost::multi::index_range extension(std::vector<T> const& v){return boost::multi::index_range{0, static_cast<boost::multi::index>(v.size())};}
}

namespace boost{
namespace multi{

using iextension = index_extension;
using irange = index_range;

namespace detail{
	template<typename T, typename... As>
	inline void construct_from_initializer_list(T* p, As&&... as){
		::new(static_cast<void*>(p)) T(std::forward<As>(as)...);
	}
}

template<class T, std::size_t N> 
index_range extension(T(&)[N]){return {0, N};}

template<typename T, std::size_t N>
std::array<T, N-1> tail(std::array<T, N> const& a){
	std::array<T, N-1> ret;
	std::copy(a.begin() + 1, a.end(), ret.begin());
	return ret;
}

namespace detail{
	template<class Tuple>
	auto head(Tuple t)->decltype(std::get<0>(t)){return std::get<0>(t);}
	template<typename Tuple, std::size_t... Ns>
	auto tail_impl(std::index_sequence<Ns...> , Tuple&& t){
		return std::make_tuple(std::get<Ns+1u>(std::forward<Tuple>(t))...);
	}
	template<typename... Ts>
	auto tail(std::tuple<Ts...> const& t){
		return tail_impl( std::make_index_sequence<sizeof...(Ts) - 1u>() , t );
	}
}

template<dimensionality_type D>
struct extents_t{
	extents_t<D-1> sub;
	index_range extent_;
	extents_t<D+1> operator[](index_extension ie) const{
		return extents_t<D+1>{*this, ie};
	}
	friend index_range const& head(extents_t const& self){return head(self.sub);}
	friend extents_t<D-1> const& tail(extents_t const& self){return self.sub;}
};

template<>
struct extents_t<0u>{
	extents_t<1u> operator[](index_extension ie) const;//{return {*this, ie};}
};

template<>
struct extents_t<1u>{
	extents_t<0u> sub;
	index_range extent_;
	friend index_range const& head(extents_t const& self){return self.extent_;}
	friend extents_t<0u> const& tail(extents_t const& self){return self.sub;}
	extents_t<2u> operator[](index_extension ie) const{
		return {*this, ie};
	}
};

inline extents_t<1u> extents_t<0u>::operator[](index_extension ie) const{
	return {*this, ie};
}

static constexpr extents_t<0u> extents;

template<dimensionality_type D>
struct layout_t{
	static constexpr dimensionality_type dimensionality = D;
	layout_t<D-1> sub;
	index stride_;
	index offset_;
	index nelems_;
	using extensions_type = std::array<index_extension, D>;
	auto operator()(index i) const{return i*stride_ - offset_;}
	auto origin() const{return -offset_ + sub.origin();}
	constexpr layout_t(index_extension ie, layout_t<D-1> const& s) : 
		sub{s},
		stride_{ie.size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{ie.size()*sub.num_elements()}
	{}
	constexpr layout_t(extensions_type e = {}) : 
		sub{tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{}
	template<class Extensions, typename = std::enable_if_t<not std::is_base_of<layout_t, std::decay_t<Extensions>>{}> >
	constexpr layout_t(Extensions&& e) :
		sub{detail::tail(std::forward<Extensions>(e))},
		stride_{index_extension(std::get<0>(std::forward<Extensions>(e))).size()*sub.num_elements()!=0?sub.size()*sub.stride():1},
		offset_{0},
		nelems_{index_extension(std::get<0>(std::forward<Extensions>(e))).size()*sub.num_elements()}
	{}
	constexpr layout_t(extents_t<D> e) :
		stride_{head(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{head(e).size()*sub.num_elements()}
	{}	
	constexpr auto nelems() const{return nelems_;}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr bool operator==(layout_t const& o) const{
		return sub==o.sub and stride_==o.stride_ and offset_==o.offset_ and nelems_==o.nelems_;
	}
	constexpr bool operator!=(layout_t const& o) const{return not(*this==o);}
	constexpr size_type num_elements() const{return size()*sub.num_elements();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return not nelems_;}
	friend bool empty(layout_t const& s){return s.empty();}
	constexpr size_type size() const{
		if(nelems_ == 0) return 0;
		assert(stride_ != 0 and nelems_%stride_ == 0 );
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& l){return l.size();}
	layout_t& rotate(){
		using std::swap;
		swap(stride_, sub.stride_);
		swap(offset_, sub.offset_);
		swap(nelems_, sub.nelems_);
		sub.rotate();
		return *this;
	}
	layout_t& unrotate(){
		sub.unrotate();
		using std::swap;
		swap(stride_, sub.stride_);
		swap(offset_, sub.offset_);
		swap(nelems_, sub.nelems_);
		return *this;
	}
	layout_t& rotate(dimensionality_type r){
		if(r >= 0) while(r){rotate(); --r;}
		else return rotate(D - r);
		return *this;
	}
	size_type size(dimensionality_type d) const{return d?sub.size(d-1):size();}
	index stride(dimensionality_type d = 0) const{
		return d?sub.stride(d-1):stride_;
	}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	void offsets_aux(index* it) const{
		*it = offset();
		sub.offsets_aux(++it);
	}
	auto offsets() const{
		std::array<index, D> ret;
		offsets_aux(ret.begin());
		return ret;
	}
	constexpr index offset() const{return offset_;}
	constexpr index offset(dimensionality_type d) const{
		return d?sub.offset(d-1):offset_;
	}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	decltype(auto) shape() const{return sizes();}
	friend decltype(auto) shape(layout_t const& self){return self.shape();}
	auto sizes() const{
		std::array<size_type, D> ret;
		sizes_aux(ret.begin());
		return ret;
	}
	friend constexpr auto sizes(layout_t const& self){return self.sizes();}
	void sizes_aux(size_type* it) const{
		*it = size(); 
		sub.sizes_aux(++it);
	}
	auto strides() const{
		std::array<index, D> ret;
		strides_aux(ret.begin());
		return ret;
	}
	friend constexpr auto strides(layout_t const& self){return self.strides();}
	void strides_aux(size_type* it) const{
		*it = stride();
		sub.strides_aux(++it);
	}
	constexpr index_extension extension() const{
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	constexpr index_extension extension_aux() const{
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	template<dimensionality_type DD = 0>
	constexpr index_extension extension(dimensionality_type d) const{
		return d?sub.extension(d-1):extension();
	}
	constexpr auto extensions() const{
		extensions_type ret;
		extensions_aux(ret.begin());
		return ret;
	}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
	void extensions_aux(index_extension* it) const{
		*it = extension();
		++it;
		sub.extensions_aux(it);
	}
};

template<>
struct layout_t<0u>{
	static constexpr dimensionality_type dimensionality = 0u;
};

template<>
struct layout_t<1u>{
	static constexpr dimensionality_type dimensionality = 1u;
	using extensions_type = std::array<index_extension, 1>;
	index stride_;
	index offset_;
	index nelems_;
	layout_t() = default;
	constexpr layout_t(index_extension ie, layout_t<0> const&) : 
		stride_{1}, 
		offset_{0}, 
		nelems_{ie.size()}
	{}
	template<class Extensions, typename = decltype(std::get<0>(Extensions{}))>
	constexpr layout_t(Extensions e) : 
		stride_{1*1}, offset_{0}, nelems_{std::get<0>(e).size()*1}
	{}
	constexpr auto offset() const{return offset_;}
	constexpr auto offset(dimensionality_type d) const{
		(void)d;
		assert(d==0);
		return offset_;
	}
	constexpr auto nelems() const{return nelems_;}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr size_type size() const{
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return nelems_/stride_;
	}
	constexpr size_type size(dimensionality_type d) const{
		(void)d;
		assert(d == 0 and stride_ != 0 and nelems_%stride_ == 0);
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& self){return self.size();}
	auto sizes() const{
		std::array<size_type, 1> ret;
		sizes_aux(ret.begin());
		return ret;
	}
	void sizes_aux(size_type* it) const{*it = size();}
	constexpr index stride(dimensionality_type d = 0) const{(void)d; assert(d == 0); return stride_;}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	auto strides() const{
		std::array<index, 1> ret;
		strides_aux(ret.begin());
		return ret;
	}
	constexpr auto offsets() const{return std::array<index, 1>{{offset()}};}
	void offsets_aux(index* it) const{*it = offset();}
	void strides_aux(size_type* it) const{*it = stride();}
	constexpr size_type num_elements() const{return size();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return not nelems_;}
	friend bool empty(layout_t const& s){return s.empty();}
	constexpr index_extension extension(dimensionality_type d = 0) const{
		(void)d;
		assert(d == 0);
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	constexpr extensions_type extensions() const{return {extension()};}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
	void extensions_aux(index_extension* it) const{*it = extension();}
	constexpr auto operator()(index i) const{return i*stride_ - offset_;}
	constexpr auto origin() const{return -offset_;}
	constexpr bool operator==(layout_t const& other) const{
		return stride_==other.stride_ and offset_==other.offset_ and nelems_==other.nelems_;
	} 
	constexpr bool operator!=(layout_t const& other) const{return not(*this==other);}
	layout_t& rotate(){return *this;}
	layout_t& unrotate(){return *this;}
	decltype(auto) shape() const{return sizes();}
};

template<typename T, dimensionality_type D, class Alloc = std::allocator<T>>
struct array;

/*
template<class Ptr>
auto default_allocator_of(Ptr const&)
->decltype(typename pointer_traits<Ptr>::default_allocator_type{}){return {};}

template<class Pointer, 
	typename = decltype(default_allocator_of(std::declval<Pointer const&>()))
>
struct pointer_traits : std::pointer_traits<Pointer>{
	using default_allocator_type = 
		decltype(default_allocator_of(std::declval<Pointer const&>()))
//		std::allocator<typename std::pointer_traits<Pointer>::element_type>
	;
};*/


template<class Pointer>
struct pointer_traits : std::pointer_traits<Pointer>{
//	using default_allocator_type = 
//		decltype(default_allocator_of(std::declval<Pointer const&>()));
	using allocator_type = std::allocator<std::decay_t<typename pointer_traits::element_type>>;
	static allocator_type allocator_of(Pointer){return {};}
};

template<class T>
auto default_allocator_of(T*){
	return std::allocator<std::decay_t<typename std::pointer_traits<T*>::element_type>>{};
}

template<class Ptr>
auto default_allocator_of(Ptr){
	return std::allocator<std::decay_t<typename std::pointer_traits<Ptr>::element_type>>{};
}

/*
template<class Pointer>
struct pointer_traits : std::pointer_traits<Pointer>{
	using default_allocator_type = 
		decltype(default_allocator_of(std::declval<Pointer>()))
//		std::allocator<typename std::pointer_traits<Pointer>::element_type>
	;
};*/


template<
	typename T, 
	dimensionality_type D, 
	typename ElementPtr,// = T const*, 
	class Layout = layout_t<D>//, class Allocator = std::allocator<T>
> 
struct basic_array : Layout{
	using layout_t = Layout;
	using index = boost::multi::index;
	using element = T;
	using element_ptr = typename std::decay<ElementPtr>::type;
	using element_const_ptr = typename std::pointer_traits<element_ptr>::template rebind<element const>;
	using value_type      = multi::array<element, D-1>;
	using decay_type = multi::array<element, D, typename pointer_traits<element_ptr>::allocator_type>;
// multi::array<element, D, typename multi::pointer_traits<element_ptr>::default_allocator_type>;
	using difference_type = index;
//	using allocator_type = Allocator;
	using sub_element_ptr = decltype(std::declval<element_ptr>() + std::declval<Layout>().operator()(0));
	using const_reference = basic_array<element, D-1, element_const_ptr>;
	using reference       = basic_array<element, D-1, element_ptr>;
	friend struct basic_array<element, D+1, element_ptr>;
	friend struct basic_array<element, D+1, element_ptr&>;
	friend struct basic_array<element, D, typename std::pointer_traits<element_ptr>::template rebind<element>>;
protected:
	using initializer_list = std::initializer_list<typename basic_array<T, D-1, ElementPtr>::initializer_list>;
	ElementPtr data_; // TODO call it base_ ?
	basic_array(ElementPtr data, Layout layout) : Layout{layout}, data_{data}{}
	// note here that I am protecting the cctor and making the inplace new a friend, this seems to be a way to allow initializer_list for a move-only  
public:
	basic_array(basic_array const& o) : Layout{o.layout()}, data_{o.data_}{}
	template<class TT, class...As> 
	friend void detail::construct_from_initializer_list(TT*, As&&...);
	element_ptr        base(){return data_;}
	element_ptr        base() const{return data_;}
	element_const_ptr cbase() const{return base();}
	Layout const& layout() const{return *this;}
	friend Layout const& layout(basic_array const& self){return self.layout();}
	basic_array(basic_array&& other) : Layout{other.layout()}, data_{other.data_}{}
	operator basic_array<element, D, element_const_ptr>() const{
		return basic_array<element, D, element_const_ptr>(data_, layout());
	}
	element_ptr origin() const{return data_ + Layout::origin();}
	element_const_ptr corigin() const{return origin();}
	friend decltype(auto) origin(basic_array const& self){return self.origin();}
	friend decltype(auto) corigin(basic_array const& self){return self.corigin();}
	reference operator[](index i) const{ assert(i<this->extension().last() and i>=this->extension().first());
		return reference{data_ + Layout::operator()(i), Layout::sub};
	}
	basic_array sliced(index first, index last) const{
		layout_t new_layout = *this;
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {data_ + Layout::operator()(first), new_layout};
	}
	basic_array sliced(index first, index last, index stride) const{
		return sliced(first, last).strided(stride);
	}
	auto range(index_range const& ir) const{
		return sliced(ir.front(), ir.front() + ir.size());
	}
	auto range(index_range const& ir, dimensionality_type n) const{
		return rotated(n).range(ir).rotated(-n);
	}
	auto strided(index s) const{
		layout_t new_layout = *this;
		new_layout.stride_*=s;
		return basic_array{data_ + Layout::operator()(this->extension().front()), new_layout};		
	}
	auto operator()(index_range const& a) const{return range(a);}
	auto operator()(index i) const{return operator[](i);}
	decltype(auto) paren_aux() const{return *this;}
	template<class... As>
	auto paren_aux(index_range a, As&&... as) const{return operator()(a).rotated().paren_aux(as...);}
	auto rotated() const{
		layout_t new_layout = *this; 
		new_layout.rotate();
		return basic_array<T, D, ElementPtr>{data_, new_layout};
	}
	auto unrotated() const{
		layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{data_, new_layout};
	}
	auto rotated(dimensionality_type i) const{
		layout_t new_layout = *this; 
		new_layout.rotate(i);
		return basic_array<T, D, ElementPtr>{data_, new_layout};
	}
	template<class... As>
	auto paren_aux(index i, As&&... as) const{return operator[](i).rotated().paren_aux(as...);}
	template<class... As>
	auto operator()(index_range a, As&&... as) const{
		return paren_aux(a, as...).rotated(-(1 + sizeof...(As)));
	}
	template<class... As>
	auto operator()(index i, As&&... as) const{
		return operator[](i)(as...);//.rotated(-(1 + sizeof...(As)));
	}
	auto operator()(index_range const& a, index_range const& ir1) const{
		return paren_aux(a, ir1).rotated(-2);
	}
	auto operator()(index_range const& a, index       const& i1) const{
		return paren_aux(a, i1).rotated(-2);
	}
	auto operator()(index const& a, index_range const& ir1) const{
		return operator[](a).sliced(ir1.first(), ir1.last());
	}
	decltype(auto) front(){return *begin();}
	decltype(auto) back(){return *(begin() + (Layout::size() - 1));}
	class const_iterator : private const_reference{
		index stride_;
		const_iterator(const_reference const& cr, index stride) : 
			const_reference{cr}, stride_{stride}
		{}
		friend struct basic_array;
		explicit operator bool() const{return this->data_;}
	public:
		using difference_type = basic_array::difference_type; //	using value_type = typename basic_array<T, D, ElementPtr>::value_type;
		using pointer = void;
		using reference = const_reference;
		using iterator_category = std::random_access_iterator_tag;
		using value_type = multi::array<element, D-1>;
		const_iterator(std::nullptr_t = nullptr) : const_reference{}{}
		const_iterator(const_iterator const&) = default;
		const_iterator& operator=(const_iterator const& other) = default;
		const_reference const& operator*() const{assert(operator bool()); return *this;}
		const_reference const* operator->() const{return static_cast<const_reference const*>(this);}
		const_reference operator[](index i) const{return {this->data_ + i*stride_, Layout::sub};}
		const_iterator& operator++(){this->data_ += stride_; return *this;}
		const_iterator& operator--(){this->data_ -= stride_; return *this;}
		const_iterator& operator+=(index d){this->data_ += stride_*d; return *this;}
		const_iterator& operator-=(index d){this->data_ -= stride_*d; return *this;}
		const_iterator operator+(index d){const_iterator ret = *this; return ret += d;}
		const_iterator operator-(index d){const_iterator ret = *this; return ret -= d;}
		ptrdiff_t operator-(const_iterator const& other) const{
			assert(stride_ != 0 and (this->data_ - other.data_)%stride_ == 0);
			assert( this->layout() == other.layout() );
			return (this->data_ - other.data_)/stride_;
		}
		bool operator==(const_iterator const& other) const{
			return this->data_==other.data_ and this->stride_==other.stride_;
		}
		bool operator!=(const_iterator const& other) const{return not((*this)==other);}
		bool operator<(const_iterator const& other) const{
			if(stride_ < 0) return other.data_ - this->data_ < 0;
			return this->data_ - other.data_ > 0;
		}
	};
	struct iterator : private reference{
		index stride_;
		iterator(basic_array<T, D-1, element_ptr> const& cr, index stride) : 
			basic_array<T, D-1, element_ptr>{cr}, stride_{stride}
		{}
		friend struct basic_array;
		explicit operator bool() const{return this->data_;}
	public:
		operator const_iterator() const{
			return const_iterator{static_cast<reference const&>(*this), stride_};
		}
		using difference_type = typename basic_array<T, D, ElementPtr>::difference_type;
		using value_type = typename basic_array<T, D, ElementPtr>::value_type;
		using pointer = void*;
		using reference = typename basic_array<T, D, element_ptr /*ElementPtr*/>::reference; // careful with shadowing reference (there is another one level up)
		iterator(std::nullptr_t = nullptr) : basic_array<T, D-1, element_ptr>{}{}
		using iterator_category = std::random_access_iterator_tag;
		iterator& operator=(iterator const& other) = default;
		template<class It>
		bool operator==(It const& other) const{
			return this->data_ == other.data_ and this->stride_ == other.stride_;
		}
		template<class It>
		bool operator!=(It const& other) const{return not((*this)==other);}
		reference& operator*() const{return *const_cast<iterator*>(this);}
		reference* operator->(){return static_cast<reference*>(this);}
		iterator& operator++(){this->data_ += stride_; return *this;}
		auto operator++(int){
           iterator ret{*this};   // make a copy for result
           ++(*this);              // Now use the prefix version to do the work
           return ret;      
		}
		iterator& operator--(){this->data_ -= stride_; return *this;}
		iterator  operator+(index d){iterator ret = *this; return ret += d;}
		iterator  operator-(index d){iterator ret = *this; return ret -= d;}
		iterator& operator+=(index d){this->data_ += stride_*d; return *this;}
		iterator& operator-=(index d){this->data_ -= stride_*d; return *this;}
		ptrdiff_t operator-(iterator const& other) const{
			assert(stride_ != 0 and (this->data_ - other.data_)%stride_ == 0);
			return (this->data_ - other.data_)/stride_;
		}
		bool operator<(iterator const& other) const{
			if(stride_ < 0) return other.data_ - this->data_ < 0;
			return this->data_ - other.data_ < 0;
		}
		friend void iter_swap(iterator i1, iterator i2){i1->swap(*i2);}
	};
	friend size_type size(basic_array const& self){return self.size();}
	friend auto sizes(basic_array const& self){return self.sizes();}
	friend auto strides(basic_array const& self){return self.strides();}
	friend auto extension(basic_array const& self){return self.extension();}
	friend auto extensions(basic_array const& self){return self.extensions();}
	const_iterator begin(index i) const{
		Layout new_layout = *this;
		new_layout.rotate(i);
		return const_iterator{
			const_reference{data_ + new_layout(0), new_layout.sub},
			new_layout.stride_
		};
	}
	const_iterator end(index i) const{
		Layout new_layout = *this;
		new_layout.rotate(i);
		return const_iterator{
			const_reference{data_ + new_layout(new_layout.size()), new_layout.sub},
			new_layout.stride_
		};
	}
	iterator begin(index i){
		Layout new_layout = *this;
		new_layout.rotate(i);
		return iterator{
			reference{data_ + new_layout(0), new_layout.sub},
			new_layout.stride_
		};
	}
	iterator end(index i){
		Layout new_layout = *this;
		new_layout.rotate(i);
		return iterator{
			reference{data_ + new_layout(new_layout.size()), new_layout.sub},
			new_layout.stride_
		};
	}
	const_iterator cbegin()  const{
		iterator r = begin();
		const_iterator cr = r;
		return r;
	}
	const_iterator cend()    const{return end();}
	iterator begin() const{return {reference{data_ + Layout::operator()(this->extension().first()), Layout::sub}, Layout::stride_};}
	iterator end()   const{return {reference{data_ + Layout::operator()(this->extension().last()), Layout::sub}, Layout::stride_};}

	friend const_iterator begin(basic_array const& self){return self.begin();}
	friend const_iterator end(basic_array const& self){return self.end();}
	friend iterator begin(basic_array& self){return self.begin();}
	friend iterator end(basic_array& self){return self.end();}
	friend const_iterator cbegin(basic_array const& self){return self.begin();}
	friend const_iterator cend(basic_array const& self){return self.end();}
protected:
	template<class It>
	void recursive_assign_(It first, It last){
		auto self_first = this->begin();
		while(first != last){
			self_first->recursive_assign_(first->begin(), first->end());
			++first;
			++self_first;
		}
	}
	template<class M>
	void intersection_assign_(M&& other){
		for(
			auto i = std::max(basic_array::extension(0).first(), other.extension(0).first()); 
			i != std::min(basic_array::extension(0).last(), other.extension(0).last()); ++i
		) operator[](i).intersection_assign_(other[i]);
	}
public:
	basic_array const& operator=(basic_array const& other) const{
		return operator=<basic_array>(other);
	}
	template<class Array>
	basic_array const& operator=(Array const& other) const{
		assert(this->extension() == extension(other));
		using std::copy; using std::begin; using std::end;
		copy(begin(other), end(other), this->begin());
		return *this;
	}
	template<class Array>
	void swap(Array&& other) const{
		assert(this->extension() == extension(other));
		using std::swap_ranges; using std::begin;
		swap_ranges(basic_array::begin(), basic_array::end(), begin(std::forward<Array>(other)));
	}
	template<class Array> 
	bool operator==(Array const& other) const{
		if(this->extension() != other.extension()) return false;
		using std::equal;
		return equal(this->begin(), this->end(), other.begin());
	}
	template<class Array> bool operator!=(Array const& other) const{
		return not((*this) == other);
	}
	template<class Array>
	bool operator<(Array const& other) const{
		using std::lexicographical_compare;
		using std::begin; using std::end;
		return lexicographical_compare(this->begin(), this->end(), begin(other), end(other)); // needs assignable iterator
	}
	template<class Array>
	bool operator>=(Array const& other) const{return not ((*this)<other);}
	template<class Array>
	bool operator<=(Array const& other) const{return (*this)<other or (*this)==other;}
	template<class Array>
	bool operator>(Array const& other) const{return not ((*this)<=other);}
};

namespace adl{
	template<class... As>
	decltype(auto) extension(As&&... as){
		using boost::multi::extension;
		return extension(std::forward<As>(as)...);
	}
}

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{1}, ElementPtr, Layout> : Layout{
	constexpr static dimensionality_type dimensionality = 1;
	using value_type = T;
	using element = value_type;
	using element_ptr = std::decay_t<ElementPtr>;
	using element_const_ptr = typename std::pointer_traits<element_ptr>::template rebind<element const>;
	using layout_t = Layout;
	using const_reference = decltype(*std::declval<element_const_ptr>()); //T const&;
	using reference = decltype(*std::declval<element_ptr>());
protected:
	using initializer_list = std::initializer_list<T>;
	template<class It>
	void recursive_assign_(It first, It last){
		auto self_first = this->begin();
		while(first != last){
			*self_first = *first;
			++first;
			++self_first;
		}
	}
	template<class M>
	void intersection_assign_(M&& other){
		for(auto i = std::max(basic_array::extension(0).first(), other.extension(0).first()); 
			i != std::min(basic_array::extension(0).last(), other.extension(0).last()); ++i
		) operator[](i) = other[i];
	}
	ElementPtr data_;
	basic_array() = delete;
	basic_array(ElementPtr d, layout_t l) : Layout{l}, data_{d}{}
	friend struct basic_array<T, dimensionality_type{2}, element_ptr>;
	friend struct basic_array<T, dimensionality_type{2}, element_ptr&>;
	friend struct basic_array<T, dimensionality_type{2}, typename std::pointer_traits<element_ptr>::template rebind<element>>;
	friend struct basic_array<T, dimensionality_type{1}, element_const_ptr>;
public:
	Layout const& layout() const{return *this;}
	friend Layout const& layout(basic_array const& self){return self.layout();}
	template<class Array>
	basic_array const& operator=(Array const& other) const{
		assert(this->extension() == adl::extension(other));
		using std::copy; using std::begin; using std::end;
		copy(begin(other), end(other), this->begin());
		return *this;
	}
	basic_array const& operator=(basic_array const& other) const{
		assert(this->extension() == extension(other));
		if(Layout::stride() == 1 and other.stride() == 1){
			using std::copy_n;
			copy_n(other.data_, other.num_elements(), data_);
			return *this;
		}else{
			//TODO further optimize for strided case
			return operator=<basic_array>(other);
		}
	}
	basic_array(basic_array<element, 1, typename std::pointer_traits<element_ptr>::template rebind<element>> const& other) : Layout{other}, data_{other.data_}{}
	reference operator[](index i) const{
		assert( i < this->extension().last() and i >= this->extension().first() );
		return *(data_ + Layout::operator()(i)); //data_[Layout::operator()(i)];
	}
	reference operator[](index i){
		assert( i < this->extension().last() and i >= this->extension().first() );
		return *(data_ + Layout::operator()(i));//data_[Layout::operator()(i)];
	}
	basic_array sliced(index first, index last) const{
		layout_t new_layout = *this; 
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {data_ + Layout::operator()(first), new_layout};		
	}
	basic_array sliced(index first, index last, index stride) const{
		return sliced(first, last).strided(stride);
	}
	auto strided(index s) const{
		layout_t new_layout = *this;
		new_layout.stride_*=s;
		return basic_array{data_ + Layout::operator()(this->extension().front()), new_layout};		
	}
	auto range(index_range const& ir) const{return sliced(ir.front(), ir.last());}
	decltype(auto) paren_aux() const{return *this;}
	template<class... As>
	auto paren_aux(index_range a, As&&... as) const{return operator()(a).rotated().paren_aux(as...);}
	template<class... As>
	auto paren_aux(index i, As&&... as) const{return operator[](i).rotated().paren_aux(as...);}
	auto operator()(index_range const& ir) const{return range(ir);}
	auto operator()(index i) const{return operator[](i);}
	decltype(auto) rotated() const{return *this;}
	decltype(auto) unrotated() const{return *this;}	
	decltype(auto) rotated(dimensionality_type) const{return *this;}
	element_ptr origin() const{return data_ + Layout::origin();}
	friend element_ptr origin(basic_array const& self){return self.origin();}
	element_const_ptr corigin() const{return origin();}
	friend element_const_ptr corigin(basic_array const& s){return s.corigin();}
	class const_iterator{
		friend struct basic_array;
		const_iterator(element_ptr d, index s) : data_(d), stride_(s){}
	protected:
		element_ptr data_;
		index stride_;
	public:
		using difference_type = index;
		using value_type = T;
		using pointer = element_const_ptr;//void*;
		using reference = const_reference;
		using iterator_category = std::random_access_iterator_tag;
		const_reference operator*() const{return *data_;}
		element_const_ptr const& operator->() const{return data_;}
		const_reference operator[](index i) const{return *(data_ + i*stride_);}
		const_iterator& operator++(){data_ += stride_; return *this;}
		const_iterator& operator--(){data_ += stride_; return *this;}
		const_iterator& operator+=(ptrdiff_t d){
			this->data_ += stride_*d; 
			return *this;
		}
		const_iterator operator+(ptrdiff_t d) const{
			return const_iterator(*this)+=d;
		}
		ptrdiff_t operator-(const_iterator const& other) const{
			assert(stride_ != 0 and (data_ - other.data_)%stride_ == 0);
			return (data_ - other.data_)/stride_;
		}
		bool operator==(const_iterator const& other) const{
			return data_ == other.data_ and stride_ == other.stride_;
		}
		bool operator!=(const_iterator const& o) const{return not(*this==o);}
	public:
		element_const_ptr get() const{return data_;}
	};
	struct iterator : const_iterator{
		friend struct basic_array;
		using const_iterator::const_iterator;
		reference operator*() const{return *this->data_;}
		element_ptr const& operator->() const{return this->data_;}
		iterator& operator++(){const_iterator::operator++(); return *this;}
		iterator& operator+=(ptrdiff_t d){const_iterator::operator+=(d); return *this;}
		iterator operator+(ptrdiff_t d) const{return iterator(*this)+=d;}
	public:
		element_ptr get() const{return const_iterator::data_;}
	};
	friend size_type num_elements(basic_array const& self){return self.num_elements();}
	friend auto size(basic_array const& self){return self.size();}
	friend index_range extension(basic_array const& self){return self.extension();}
public:
	iterator begin() const{return {data_ + Layout::operator()(this->extension().first()), Layout::stride_};}
	iterator end()   const{return {data_ + Layout::operator()(this->extension().last()), Layout::stride_};}
	friend auto begin(basic_array const& self){return self.begin();}
	friend auto end(basic_array const& self){return self.end();}
	template<class Array>
	bool operator==(Array const& other) const{
		using multi::extension;
		if(this->extension() != extension(other)) return false;
		using std::equal; using std::begin; // using std::end;
		return equal(this->begin(), this->end(), begin(other));
	}
	template<class Array>
	bool operator!=(Array const& other) const{return not((*this)==other);}
	template<class Array>
	bool operator<(Array const& other) const{
		using std::lexicographical_compare;
		using std::begin; using std::end;
		return lexicographical_compare(
			begin(*this), end(*this), 
			begin(other), end(other)
		);
	}
	template<class Array>
	bool operator<=(Array const& other) const{return (*this)==other or (*this)<other;}
	template<class A> bool operator>(A const& other) const{return not((*this)<=other);}
	template<class A> bool operator>=(A const& other) const{return not((*this)<other);}
	template<class Array> void swap(Array&& other){
		assert(this->extension() == extension(other));
		using std::swap_ranges; using std::begin; //using std::end;
		swap_ranges(this->begin(), this->end(), begin(other));
	}
	friend void swap(basic_array& a1, basic_array& a2){return a1.swap(a2);}
	friend void swap(basic_array& a1, basic_array&& a2){return a1.swap(a2);}
	friend void swap(basic_array&& a1, basic_array& a2){return a1.swap(a2);}
	friend void swap(basic_array&& a1, basic_array&& a2){return a1.swap(a2);}
};

template<typename T, dimensionality_type D, typename ElementPtr = T const*>
struct array_cref : basic_array<T, D, ElementPtr>{
	array_cref() = delete;
	array_cref(array_cref const&) = default;
	array_cref(array_cref&&) = default;
	constexpr array_cref(ElementPtr p, typename array_cref::extensions_type e) noexcept
		: basic_array<T, D, ElementPtr>{p, typename array_cref::layout_t{e}}{}
	template<class Extensions> constexpr
	array_cref(ElementPtr p, Extensions e) noexcept 
		: basic_array<T, D, ElementPtr>{p, typename array_cref::layout_t{e}}{}
	protected:
	using basic_array<T, D, ElementPtr>::operator=;
	public:
	void operator=(array_cref const&) = delete;
	typename array_cref::element_ptr        data() const{return array_cref::data_;}
	typename array_cref::element_const_ptr cdata() const{return data();}
	friend decltype(auto)  data(array_cref const& self){return self. data();}
	friend decltype(auto) cdata(array_cref const& self){return self.cdata();}
};

template<typename T, dimensionality_type D, typename ElementPtr = T const*>
struct const_array_ref : basic_array<T, D, ElementPtr>{
	using element_const_ptr = typename const_array_ref::element_const_ptr;
	using element_ptr = typename const_array_ref::element_ptr;
	const_array_ref() = delete; // references must be initialized (bound)
	const_array_ref(const_array_ref const&) = default;
	const_array_ref(const_array_ref&&) = default;
	constexpr const_array_ref(
		ElementPtr p, // typename const_array_ref::element_ptr p, 
		typename const_array_ref::extensions_type e
	) noexcept 
	: basic_array<T, D, ElementPtr>{p, typename const_array_ref::layout_t{e}}
	{}
	template<class Extensions>
	constexpr const_array_ref(ElementPtr p, Extensions e) noexcept 
	: basic_array<T, D, ElementPtr>{p, typename const_array_ref::layout_t{e}}
	{}
	template<class Extensions>
	constexpr const_array_ref(ElementPtr p, layout_t<D> const& lo)
	: basic_array<T, D, ElementPtr>{p, lo}
	{}
	protected:
	using basic_array<T, D, ElementPtr>::operator=;
	public:
	void operator=(const_array_ref const&) = delete;
	element_ptr data() const{return const_array_ref::data_;}
	element_const_ptr cdata() const{return data();}
	friend decltype(auto) data(const_array_ref const& self){return self. data();}
	friend decltype(auto) cdata(const_array_ref const& self){return self.cdata();}
};

}}

namespace boost{
namespace multi{

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref : const_array_ref<T, D, ElementPtr>{
	using const_array_ref<T, D, ElementPtr>::const_array_ref;
	using const_array_ref<T, D, ElementPtr>::operator=;
	array_ref(array_ref const&) = default;
	array_ref(array_ref&&) = default;
	template<class Array, typename = std::enable_if_t<not std::is_base_of<const_array_ref<T, D, ElementPtr>, std::decay_t<Array>>{}> > array_ref& operator=(Array const& other){
		const_array_ref<T, D, ElementPtr>::operator=(other);
		return *this;
	}
	array_ref const& operator=(const_array_ref<T, D, ElementPtr> const& other){
		assert(this->extensions() == other.extensions());
		using std::copy_n;
		copy_n(other.data(), other.num_elements(), array_ref::data());
		return *this;
	}
	array_ref const& operator=(array_ref const& other){
		return operator=(static_cast<const_array_ref<T, D, ElementPtr> const&>(other));
	}
};

template<typename... A> using carray_ref = array_ref<A...>;

}}

namespace std{
template<class It, typename = decltype(It{}.stride_)>
void iter_swap(It it1, It it2){it1->swap(*it2);}
}

#if _TEST_BOOST_MULTI_ARRAY_REF

#include<cassert>
#include<numeric> // iota
#include<iostream>

using std::cout;
namespace multi = boost::multi;

int main(){

#if 0
	{
		double const d2D[4][5] = {{1.,2.},{2.,3.}};
		bm::array_ref<double, 2> d2Rce{&d2D[0][0], {4, 5}};
		assert( &d2Rce[2][3] == &d2D[2][3] );
		assert( d2Rce.size() == 4 );
	//	assert( d2Rce.size<0>() == 4);
	//	assert( d2Rce.size<1>() == 5);
		cout << d2Rce.num_elements() << std::endl;
		assert( d2Rce.num_elements() == 20 );
	}
#endif
	{
		double const dc2D[4][5] = {{1.,2.},{2.,3.}};
		multi::array_cref<double, 2> acrd2D{&dc2D[0][0], {4, 5}};
		assert( &acrd2D[2][3] == &dc2D[2][3] );
		assert( acrd2D.size() == 4);
		assert( acrd2D.size(0) == 4);
		assert( acrd2D.size(1) == 5);
		assert( acrd2D.sizes().size() == 2 );
	//	assert( acrd2D.size<0>() == acrd2D.size() );
	//	assert( acrd2D.size<1>() == 5);
		assert( acrd2D.num_elements() == 20 );

		assert( &acrd2D[2][3] == &dc2D[2][3] );
	
		assert( acrd2D.begin() == begin(acrd2D) );
		assert( acrd2D.begin() != acrd2D.end() );
	}
#if 0
	{
		double* d2p = new double[4*5]; std::iota(d2p, d2p + 4*5, 0);

		bm::array_ref<double, 2> d2R{d2p, 4, 5};
		assert(d2R.size()==4);
	}
	cout << "ddd " << d2R[1][1] << '\n';
	for(int i = 0; i != 4 ||!(cout << '\n'); ++i)
		for(int j = 0; j != 5 ||!(cout << '\n'); ++j)
			cout << d2R[i][j] << ' ';

	for(auto it1 = d2R.begin(); it1 != d2R.end() ||!(cout << '\n'); ++it1)
		for(auto it2 = it1->begin(); it2 != it1->end() ||!(cout << '\n'); ++it2)
			cout << *it2 << ' ';

	for(auto&& row : d2R){
		for(auto&& e : row) cout << e << ' '; cout << '\n';
	} cout << '\n';

	for(auto i : d2R.extension()){
		for(auto j : d2R[i].extension())
			cout << d2R[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';

	for(auto it1 = d2R.begin1(); it1 != d2R.end1() ||!(cout << '\n'); ++it1)
		for(auto it2 = it1->begin(); it2 != it1->end() ||!(cout << '\n'); ++it2)
			cout << *it2 << ' ';

	assert( d2R.begin()[1][1] == 6 );

	assert(d2R.size() == 4);
	auto it = d2R.begin();
	assert((*it)[1] == 1);
	assert( it->operator[](0) == 0 );
	assert( it->operator[](1) == 1 );
	++it;
	assert( it->operator[](0) == 5 );
	assert( it->operator[](1) == 6 );

	assert( *(it->begin()) == 5 );


#if 0
	if(double* d3p = new double[3*4*5]){
		std::iota(d3p, d3p + 3*4*5, 0);
		bm::array_ref<double, 3> d3R{d3p, 3, 4, 5};
		assert(d3R.size() == 3);
		for(int i = 0; i != 3; ++i, cout << '\n')
			for(int j = 0; j != 4; ++j, cout << '\n')
				for(int k = 0; k != 5; ++k)
					cout << d3R[i][j][k] << ' ';
		auto b = d3R.begin();
		
	}
#endif
#endif
}
#endif
#endif

