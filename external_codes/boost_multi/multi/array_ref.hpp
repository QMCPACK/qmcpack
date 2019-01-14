#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ -O3 -std=c++17 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MULTI_ARRAY_REF $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP


#include "utility.hpp"

#include "detail/layout.hpp"
#include "detail/types.hpp" // dimensionality_type
#include "detail/operators.hpp" // random_iterable
#include "detail/memory.hpp" // pointer_traits

#include<algorithm> // copy_n

namespace boost{
namespace multi{

template<typename T, dimensionality_type D, typename ElementPtr, class Layout = layout_t<D>> 
struct basic_array;

template<typename T, dimensionality_type D, class Alloc = std::allocator<T>> 
struct array;

template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
struct array_types : Layout{ //	template<class... As> array_types(As&&... as) : Layout{std::forward<As>(as)...}{}
	using element = T;
	constexpr static dimensionality_type dimensionality = D;
	using element_ptr = ElementPtr;
	using layout_t = Layout;
	using value_type = std::conditional_t<
		dimensionality!=1, 
		array<element, dimensionality-1>, 
		element
	>;
	using decay_type = array<element, dimensionality, typename pointer_traits<element_ptr>::allocator_type>;
	using difference_type = index;
	using reference       = std::conditional_t<
		dimensionality!=1, 
		basic_array<element, dimensionality-1, element_ptr>, 
		//typename std::iterator_traits<element_ptr>::reference 	//	typename pointer_traits<element_ptr>::element_type&
			typename pointer_traits<element_ptr>::element_type&
	>;
	element_ptr        base() const{return base_;} //	element_const_ptr cbase() const{return base();}
	layout_t const& layout() const{return *this;}
	friend layout_t const& layout(array_types const& s){return s.layout();}
	element_ptr            origin() const{return base_+Layout::origin();} //	element_const_ptr     corigin() const{return origin();}
	friend decltype(auto)  origin(array_types const& s){return s.origin();} //	friend decltype(auto) corigin(array_types const& s){return s.corigin();}
protected:
	using derived = basic_array<T, D, ElementPtr, Layout>;
	element_ptr base_;
	array_types() = delete;
	array_types(std::nullptr_t np) : Layout{}, base_{np}{}
	array_types(array_types const&) = default; //	array_types(array_types const& o) : Layout{o}, base_{o.base_}{}
	array_types(layout_t l, element_ptr data) : Layout{l}, base_{data}{}
//	template<class ArrayTypes, typename = decltype(Layout{std::declval<ArrayTypes&&>().layout()}, element_ptr{std::declval<ArrayTypes&&>().base_})>
//	array_types(ArrayTypes&& other) : Layout{other.layout()}, base_{other.base_}{}
public://TODO find why this needs to be public and not protected or friend
	template<typename ElementPtr2>
	array_types(array_types<T, D, ElementPtr2, Layout> const& other) : Layout{other.layout()}, base_{other.base_}{}
	template<class T2, dimensionality_type D2, class E2, class L2> friend struct array_types;
};

template<class Ref, class Layout>	
struct basic_array_ptr : 
	private Ref,
	boost::iterator_facade<
		basic_array_ptr<Ref, Layout>, void, std::random_access_iterator_tag, 
		Ref const&, typename Layout::difference_type
	>,
	boost::totally_ordered<basic_array_ptr<Ref, Layout>>
{
	using difference_type = typename Layout::difference_type;
	using value_type = typename Ref::decay_type;
	using pointer = Ref const*;
	using reference = Ref;//Ref const&;
	using element_type = typename Ref::value_type;
	using iterator_category = std::random_access_iterator_tag;
	basic_array_ptr(std::nullptr_t p = nullptr) : Ref{p}{}
	template<class, class> friend struct basic_array_ptr;
	template<class Other, typename = decltype(typename Ref::types::element_ptr{typename Other::element_ptr{}})> 
	basic_array_ptr(Other const& o) : Ref{o.layout(), o.base()}, stride_{o.stride_}{}
//	template<class Other>//, typename = decltype(typename Ref::types::element_ptr{typename Other::types::element_ptr{}})> 
//	basic_array_ptr(Other const& o) : Ref{o.base(), o.layout()}, stride_{o.stride_}{}
	basic_array_ptr& operator=(basic_array_ptr const& other){
		this->base_ = other.base_;
	//	layout = other.layout;
		static_cast<Layout&>(*this) = other; //impl_) = other.impl_;
		this->stride_ = other.stride_;
		return *this;
	}
//	explicit operator bool() const{return this->base_;}
	Ref const& operator*() const{/*assert(*this);*/ return *this;}
	Ref const* operator->() const{/*assert(*this);*/ return this;}
	Ref        operator[](difference_type n) const{return *(*this + n);}
	template<class O> bool operator==(O const& o) const{return equal(o);}
	bool operator<(basic_array_ptr const& o) const{return distance_to(o) > 0;}
	basic_array_ptr(typename Ref::element_ptr p, Layout l, index stride) : Ref{l, p}, stride_{stride}{}
	template<typename T, dimensionality_type D, typename ElementPtr, class LLayout>
	friend struct basic_array;
	auto stride() const{return stride_;}
	friend index stride(basic_array_ptr const& s){return s.stride();}
	auto base() const{return this->base_;}
	friend auto base(basic_array_ptr const& self){return self.base();}
private:
	index stride_ = {1}; // nice, non-zero default
	using Ref::layout;
	using Ref::base_;
	bool equal(basic_array_ptr const& o) const{return base_==o.base_ && stride_==o.stride_ && layout()==o.layout();}
	void increment(){base_ += stride_;}
	void decrement(){base_ -= stride_;}
	void advance(difference_type n){base_ += stride_*n;}
	difference_type distance_to(basic_array_ptr const& other) const{
		assert( stride_ == other.stride_ and stride_ != 0 and (other.base_ - base_)%stride_ == 0 and layout() == other.layout() );
		return (other.base_ - base_)/stride_;
	}
	friend class boost::iterator_core_access;
};

// TODO remove the T parameter
template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
struct basic_array : 
	boost::multi::partially_ordered2<basic_array<T, D, ElementPtr, Layout>, void>,
	boost::multi::random_iterable<basic_array<T, D, ElementPtr, Layout>>,
	array_types<T, D, ElementPtr, Layout>
{
	using types = array_types<T, D, ElementPtr, Layout>;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr >;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr&>;
//	friend struct basic_array<typename types::element, Layout::rank    , typename std::pointer_traits<typename types::element_ptr>::template rebind<typename types::element>>;
	using types::layout;
	decltype(auto) layout() const{return array_types<T, D, ElementPtr, Layout>::layout();}
protected:
	using types::types;
	template<typename, dimensionality_type, class Alloc> friend struct array;
public:
	using decay_type = array<typename types::element, D, typename pointer_traits<typename types::element_ptr>::allocator_type>;
	decay_type decay() const{return *this;}
	friend auto decay(basic_array const& self){return self.decay();}
	typename types::reference operator[](index i) const{
		assert( this->extension().count(i) );
		return {sub, types::base_ + Layout::operator()(i)};
	}
	basic_array sliced(typename types::index first, typename types::index last) const{
		typename types::layout_t new_layout = *this;
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};
	}
	basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = *this; new_layout.stride_*=s;
		return {new_layout, types::base_ + types::nelems_};
	}
	basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}
	auto range(typename types::index_range const& ir) const{
		return sliced(ir.front(), ir.front() + ir.size());
	}
	auto range(typename types::index_range const& ir, dimensionality_type n) const{
		return rotated(n).range(ir).rotated(-n);
	}
	auto rotated() const{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	friend auto rotated(basic_array const& self){return self.rotated();}
	auto unrotated() const{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	friend auto unrotated(basic_array const& self){return self.unrotated();}
	auto rotated(dimensionality_type i) const{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate(i);
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	basic_array const& operator()() const{return *this;}
	template<class... As>
	auto operator()(index_range a, As... as) const{return range(a).rotated()(as...).unrotated();}
	template<class... As>
	auto operator()(index i, As... as) const{return operator[](i)(as...);}
#define SARRAY2(A1, A2)	auto operator()(A1 a1, A2 a2) const{return operator()<A2>(a1, a2);}
	SARRAY2(index, index ); SARRAY2(irange, index );
	SARRAY2(index, irange); SARRAY2(irange, irange);
#undef SARRAY2
#define SARRAY3(A1, A2, A3) auto operator()(A1 a1, A2 a2, A3 a3) const{return operator()<A2, A3>(a1, a2, a3);}
	SARRAY3(index, index , index ); SARRAY3(irange, index , index );
	SARRAY3(index, index , irange); SARRAY3(irange, index , irange);
	SARRAY3(index, irange, index ); SARRAY3(irange, irange, index );
	SARRAY3(index, irange, irange); SARRAY3(irange, irange, irange);
#undef SARRAY3
#define SARRAY4(A1, A2, A3, A4) auto operator()(A1 a1, A2 a2, A3 a3, A4) const{return operator()<A2, A3, A4>(a1, a2, a3);}
	SARRAY4(index, index, index , index ); SARRAY4(index, irange, index , index );
	SARRAY4(index, index, index , irange); SARRAY4(index, irange, index , irange);
	SARRAY4(index, index, irange, index ); SARRAY4(index, irange, irange, index );
	SARRAY4(index, index, irange, irange); SARRAY4(index, irange, irange, irange);
	SARRAY4(irange, index, index , index ); SARRAY4(irange, irange, index , index );
	SARRAY4(irange, index, index , irange); SARRAY4(irange, irange, index , irange);
	SARRAY4(irange, index, irange, index ); SARRAY4(irange, irange, irange, index );
	SARRAY4(irange, index, irange, irange); SARRAY4(irange, irange, irange, irange);
#undef SARRAY4
private:
	using Layout::nelems_;
	using Layout::stride_;
	using Layout::sub;
public:
	using iterator = basic_array_ptr<typename types::reference, typename types::sub_t>;
private:
	template<class Iterator>
	struct basic_reverse_iterator : 
		std::reverse_iterator<Iterator>,
		boost::totally_ordered<basic_reverse_iterator<Iterator>>
	{
		template<class O, typename = decltype(std::reverse_iterator<Iterator>{base(std::declval<O const&>())})>
		basic_reverse_iterator(O const& o) : std::reverse_iterator<Iterator>{base(o)}{}
		basic_reverse_iterator() : std::reverse_iterator<Iterator>(){}
		explicit basic_reverse_iterator(Iterator it) : std::reverse_iterator<Iterator>(--it){}
		explicit operator Iterator() const{auto ret = this->base(); return ++ret;}
		explicit operator bool() const{auto ret = this->base(); ++ret; return static_cast<bool>(ret);}
		typename Iterator::reference operator*() const{return this->current;}
		typename Iterator::pointer operator->() const{return &this->current;}
		typename Iterator::reference operator[](typename Iterator::difference_type n) const{return *(this->current - n);}
		bool operator<(basic_reverse_iterator const& o) const{return o.base() < this->base();}
	};
public:
	using reverse_iterator = basic_reverse_iterator<iterator>;
	using ptr = basic_array_ptr<basic_array, Layout>;
	ptr operator&() const{return {this->base_, this->layout(), this->nelems_};}
//	friend auto strides(basic_array const& self){return self.strides();}
	iterator begin(index i) const{
		Layout l = static_cast<Layout const&>(*this); l.rotate(i);
		return {types::base_ + l(0       ), l.sub, l.stride_};
	}
	iterator end(index i) const{
		Layout l = static_cast<Layout const&>(*this); l.rotate(i);
		return {types::base_ + l(l.size()), l.sub, l.stride_};
	}
	iterator  begin() const{return {types::base_          , sub, stride_};}
	iterator  end  () const{return {types::base_ + nelems_, sub, stride_};}
//	auto     rbegin() const{return reverse_iterator{end()};}
//	auto     rend  () const{return reverse_iterator{begin()};}
	typename types::reference front() const{return * begin();}
	typename types::reference back()  const{return *this->rbegin();}
protected:
	template<class A>
	void intersection_assign_(A&& other) const{
		for(auto i : intersection(types::extension(), extension(other)))
			operator[](i).intersection_assign_(std::forward<A>(other)[i]);
	}
public:
	template<class It> void assign(It first, It last) const{
		using std::copy; using std::distance;
		assert( distance(first, last) == this->size() );
		copy(first, last, this->begin());
	}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}>>
	basic_array const& operator=(A&& o) const{
		using multi::extension;
		assert(extension(*this) == extension(o));
		using std::begin; using std::end;
		assign(begin(std::forward<A>(o)), end(std::forward<A>(o)));
		return *this;
	}
	basic_array const& operator=(basic_array const& o) const{
		using multi::extension;
		assert(extension(*this) == extension(o));
		using std::begin; using std::end;
		assign(begin(o), end(o));
		return *this;
	}
	template<class Array>
	void swap(Array&& o) const{
		using multi::extension; 
		using std::swap_ranges; 
		using std::begin;
		assert(this->extension() == extension(o));
		swap_ranges(this->begin(), this->end(), begin(std::forward<Array>(o)));
	}
	friend void swap(basic_array const& a1, basic_array const& a2){
		return a1.swap(a2);
	}
	template<class Array> void swap(basic_array const& a1, Array&& a2){return a1.swap(a2);}
	template<class Array> void swap(Array&& a1, basic_array const& a2){return a2.swap(a1);}
	template<class Array> 
	bool operator==(Array const& other) const{
		using multi::extension; using std::equal; using std::begin;
		if(this->extension() != extension(other)) return false;
		return equal(this->begin(), this->end(), begin(other));
	}
private:
	template<class A1, class A2>
	static auto lexicographical_compare(A1 const& a1, A2 const& a2){
		using multi::extension;
		if(extension(a1).first() > extension(a2).first()) return true;
		if(extension(a1).first() < extension(a2).first()) return false;
		using std::lexicographical_compare;
		using std::begin; using std::end;
		return lexicographical_compare(begin(a1), end(a1), begin(a2), end(a2));
	}
public:
	template<class O>
	bool operator<(O const& o) const{return lexicographical_compare(*this, o);}
	template<class O>
	bool operator>(O const& o) const{return lexicographical_compare(o, *this);}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{1}, ElementPtr, Layout> : 
	multi::partially_ordered2<basic_array<T, dimensionality_type(1), ElementPtr, Layout>, void>,
	multi::random_iterable<basic_array<T, dimensionality_type(1), ElementPtr, Layout> >,
	array_types<T, dimensionality_type(1), ElementPtr, Layout>
{
	using types = array_types<T, dimensionality_type{1}, ElementPtr, Layout>;
	using types::types;
	using decay_type = array<typename types::element, dimensionality_type{1}, typename pointer_traits<typename types::element_ptr>::allocator_type>;
	decay_type decay() const{return *this;}
	friend decay_type decay(basic_array const& self){return self.decay();}
protected:
	template<class A>
	void intersection_assign_(A&& other) const{
		for(auto i : intersection(types::extension(), extension(other)))
			operator[](i) = std::forward<A>(other)[i];
	}
protected:
	friend struct basic_array<T, dimensionality_type{types::dimensionality + 1}, typename types::element_ptr>;
	friend struct basic_array<T, dimensionality_type{types::dimensionality + 1}, typename types::element_ptr&>;
//	friend struct basic_array<T, dimensionality_type{types::dimensionality + 1}, typename std::pointer_traits<typename types::element_ptr>::template rebind<typename types::element>>;

public:
	template<class TT, dimensionality_type DD, typename EP, class LLayout>
	friend struct basic_array;
	template<class TT, dimensionality_type DD, class Alloc> friend struct array;
	basic_array(basic_array const&) = default;
	basic_array_ptr<basic_array, Layout> operator&() const{
		return {this->base_, this->layout(), this->nelems_};
	}
	template<class It> void assign(It first, It last) const{
		using std::distance; assert(distance(first, last) == this->size());
		using std::copy; copy(first, last, this->begin());
	}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> >
	basic_array const& operator=(A&& o) const{
		using multi::extension;
		assert(this->extension() == extension(o));
		using std::begin; using std::end;
		this->assign(begin(std::forward<A>(o)), end(std::forward<A>(o)));
		return *this;
	}
	basic_array const& operator=(basic_array const& o) const{
		using multi::extension;
		assert(this->extension() == extension(o));
		using std::begin; using std::end;
		this->assign(begin(o), end(o));
		return *this;
	}
//	basic_array& operator=(basic_array const& o) = delete;
#if 0
	basic_array const& operator=(basic_array const& o) const{
		return operator=<basic_array const&>(o);
	/*	using multi::extension;
		assert(this->extension() == extension(o));
		using std::begin; using std::end;
		this->assign(begin(o), end(o));
		return *this;*/
	}
	basic_array const& operator=(basic_array&& o) const{
	//	return operator=<basic_array&&>(std::move(other));
		using multi::extension;
		assert(this->extension() == extension(o));
		using std::begin; using std::end;
		this->assign(begin(o), end(o));
		return *this;
	}
#endif
	typename types::reference //	decltype(auto) 
	operator[](typename types::index i) const{
//		return types::base_[Layout::operator()(i)];
		return *(types::base_ + Layout::operator()(i)); // less requirements
	}
	basic_array sliced(typename types::index first, typename types::index last) const{
		typename types::layout_t new_layout = *this; 
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};		
	}
	basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = *this;
		new_layout.stride_*=s;
		return {new_layout, types::base_ + Layout::operator()(this->extension().front())};		
	}
	basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}
	auto range(index_range const& ir) const{return sliced(ir.front(), ir.last());}
	decltype(auto) operator()() const{return *this;}
	auto operator()(index_range const& ir) const{return range(ir);}
	auto operator()(typename types::index i) const{return operator[](i);}
	decltype(auto) rotated() const{return *this;}
	decltype(auto) unrotated() const{return *this;}	
	decltype(auto) rotated(dimensionality_type) const{return *this;}
public:
	template<typename Ptr, typename Ref>
	struct basic_iterator : 
		boost::iterator_facade<
			basic_iterator<Ptr, Ref>, 
			typename types::element, std::random_access_iterator_tag, 
			Ref, typename types::difference_type
		>
	{
		template<class Other, typename = decltype(Ptr{typename Other::pointer{}})> 
		basic_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_}{}
		basic_iterator(std::nullptr_t np = nullptr) : data_{np}, stride_{}{}
		explicit operator bool() const{return static_cast<bool>(this->data_);}
		Ref operator[](typename types::difference_type n) const{assert(*this); return *((*this) + n);}
	private:
		basic_iterator(Ptr d, typename basic_array::index s) : data_{d}, stride_{s}{}
		friend struct basic_array;
		Ptr data_;
		typename types::index stride_;
		Ref dereference() const{/*assert(data_);*/ return *data_;}
		bool equal(basic_iterator const& o) const{return data_ == o.data_ and stride_ == o.stride_;}
		void increment(){data_ += stride_;}
		void decrement(){data_ -= stride_;}
		void advance(typename types::difference_type n){data_ += stride_*n;}
		difference_type distance_to(basic_iterator const& other) const{
			assert( stride_ == other.stride_ and stride_ != 0 and (other.data_ - data_)%stride_ == 0);
			return (other.data_ - data_)/stride_;
		}
	    friend class boost::iterator_core_access;
	    auto base() const{return data_;}
	    friend auto base(basic_iterator const& self){return self.base();}
	    auto stride() const{return stride_;}
	    friend auto stride(basic_iterator const& self){return self.stride();}
	};
public:
	using iterator = basic_iterator<typename types::element_ptr, typename types::reference>;//decltype(*std::declval<typename types::element_ptr>())>;
	using reverse_iterator = std::reverse_iterator<iterator>;

	friend size_type num_elements(basic_array const& self){return self.num_elements();}
	friend index_range extension(basic_array const& self){return self.extension();}

	iterator  begin() const{return {types::base_                 , Layout::stride_};}
	iterator  end  () const{return {types::base_ + types::nelems_, Layout::stride_};}
	auto     rbegin() const{return reverse_iterator{end()  };}
	auto     rend  () const{return reverse_iterator{begin()};}

	typename types::reference front() const{return *begin();}
	typename types::reference back () const{return *rbegin();}

	template<class Array> bool operator==(Array const& other) const{
		using multi::extension; using std::equal; using std::begin;
		if(this->extension() != extension(other)) return false;
		return equal(this->begin(), this->end(), begin(other));
	}
	bool operator<(basic_array const& o) const{return operator< <basic_array const&>(o);}
	template<class Array> void swap(Array&& o) const{
		using multi::extension; using std::swap_ranges; using std::begin;
		assert(this->extension() == extension(o));
		swap_ranges(this->begin(), this->end(), begin(std::forward<Array>(o)));
	}
//	friend void swap(basic_array&& a1, basic_array&& a2){return a1.swap(a2);}
	friend void swap(basic_array const& a1, basic_array const& a2){return a1.swap(a2);}
	template<class A> friend void swap(basic_array&& a1, A&& a2){return a1.swap(std::forward<A>(a2));}
	template<class A> friend void swap(A&& a1, basic_array&& a2){return a2.swap(std::forward<A>(a1));}
private:
	template<class A1, class A2>
	static auto lexicographical_compare(A1 const& a1, A2 const& a2){
		using multi::extension;
		if(extension(a1).first() > extension(a2).first()) return true;
		if(extension(a1).first() < extension(a2).first()) return false;
		using std::lexicographical_compare;
		using std::begin; using std::end;
		return lexicographical_compare(begin(a1), end(a1), begin(a2), end(a2));
	}
public:
	template<class O>
	bool operator<(O const& o) const{return lexicographical_compare(*this, o);}
	template<class O>
	bool operator>(O const& o) const{return lexicographical_compare(o, *this);}

};

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref : 
	boost::multi::partially_ordered2<array_ref<T, D, ElementPtr>, void>,
	basic_array<T, D, ElementPtr>
{	
protected:
	constexpr array_ref() noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{}, nullptr}{}
public:
	array_ref(array_ref const&) = default;
//	constexpr array_ref(typename array_ref::extensions_type const& e, ElementPtr p) noexcept
//		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p}{}
	constexpr array_ref(typename array_ref::element_ptr p, typename array_ref::extensions_type const& e) noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p}{}
#if defined(__INTEL_COMPILER) or (__GNUC < 5)
	constexpr array_ref(typename array_ref::element_ptr p, std::initializer_list<typename array_ref::index_extension> il) noexcept : array_ref{p, multi::detail::to_tuple<D, typename array_ref::index_extension>(il)}{}
//	constexpr array_ref(typename array_ref::element_ptr p, std::initializer_list<typename array_ref::index> il) noexcept : array_ref{p, multi::detail::to_tuple<D, typename array_ref::index_extension>(il)}{}
#endif

	using basic_array<T, D, ElementPtr>::operator=;
	using basic_array<T, D, ElementPtr>::operator<;
	using basic_array<T, D, ElementPtr>::operator>;
//	array_ref& operator=(array_ref const&) = delete;
//	template<class Array>//, std::enable_if_t<not std::is_base_of<array_cref<T, D, ElementPtr>, Array>{}>* =0> 
/*	array_ref& operator=(Array&& other) const{
		assert(this->extensions() == extensions(other));
		using std::copy_n;
		copy_n(other.data(), other.num_elements(), array_ref::data());
		return *this;

		array_cref<T, D, ElementPtr>::operator=(other);
		return *this;
	}*/
//	array_ref const& operator=(array_ref<T, D, ElementPtr> const& other)  = default;
/*	array_ref const& operator=(array_ref<T, D, ElementPstr> const& other) const{
		assert(this->extensions() == other.extensions());
		using std::copy_n;
		copy_n(other.data(), other.num_elements(), array_ref::data());
		return *this;
	}*/
//	array_ref const& operator=(array_ref const& o) const{
//		return operator=(static_cast<array_cref<T, D, ElementPtr> const&>(o));
//	}
	typename array_ref::element_ptr data() const{return array_ref::base_;}
	friend typename array_ref::element_ptr data(array_ref const& self){return self.data();}

};

template<class T, dimensionality_type D> 
using array_cref = array_ref<T, D, T const*>;

//template<dimensionality_type D, class P>
//array_ref<typename std::iterator_traits<P>::value_type, D, P> 
//make_array_ref(P p, typename detail::repeat<index_extension, D>::type x){return {p, x};}

template<dimensionality_type D, class P>
array_ref<typename std::iterator_traits<P>::value_type, D, P> 
make_array_ref(P p, index_extensions<D> x){return {p, x};}

template<class P> auto make_array_ref(P p, index_extensions<1> x){return make_array_ref<1>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<2> x){return make_array_ref<2>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<3> x){return make_array_ref<3>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<4> x){return make_array_ref<4>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<5> x){return make_array_ref<5>(p, x);}

//In ICC you need to specify the dimensionality in make_array_ref<D>
#if defined(__INTEL_COMPILER)
template<dimensionality_type D, class P> 
auto make_array_ref(P p, std::initializer_list<index_extension> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
template<dimensionality_type D, class P> 
auto make_array_ref(P p, std::initializer_list<index> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
#endif

#if __cpp_deduction_guides
template<class It> array_ref(It, typename detail::repeat<index_extension, 1>::type)->array_ref<typename std::iterator_traits<It>::value_type, 1, It>;
template<class It> array_ref(It, typename detail::repeat<index_extension, 2>::type)->array_ref<typename std::iterator_traits<It>::value_type, 2, It>;
template<class It> array_ref(It, typename detail::repeat<index_extension, 3>::type)->array_ref<typename std::iterator_traits<It>::value_type, 3, It>;
template<class It> array_ref(It, typename detail::repeat<index_extension, 4>::type)->array_ref<typename std::iterator_traits<It>::value_type, 4, It>;
template<class It> array_ref(It, typename detail::repeat<index_extension, 5>::type)->array_ref<typename std::iterator_traits<It>::value_type, 5, It>;
#endif

}}


#if 0
namespace {
	template<class T, boost::multi::dimensionality_type N, class... Ts> 
	struct rank<boost::multi::array_ref<T, N, Ts...>> 
	: public std::integral_constant<boost::multi::dimensionality_type,  typename boost::multi::array_ref<T, N, Ts...>::rank{}
	>{};

	template<class T, boost::multi::dimensionality_type N, class... Ts> 
	struct rank<boost::multi::basic_array<T, N, Ts...>> 
	: public std::integral_constant<
		boost::multi::dimensionality_type, 
		typename boost::multi::basic_array<T, N, Ts...>::rank{}
	>{};
}
#endif

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


#if _TEST_BOOST_MULTI_ARRAY_REF

#include<cassert>
#include<numeric> // iota
#include<iostream>
#include<vector>

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
		std::string const dc3D[4][2][3] = {
			{{"A0a", "A0b", "A0c"}, {"A1a", "A1b", "A1c"}},
			{{"B0a", "B0b", "B0c"}, {"B1a", "B1b", "B1c"}},
			{{"C0a", "C0b", "C0c"}, {"C1a", "C1b", "C1c"}}, 
			{{"D0a", "D0b", "D0c"}, {"D1a", "D1b", "D1c"}}, 
		};
		multi::array_cref<std::string, 3> A(&dc3D[0][0][0], {4, 2, 3});
		assert( num_elements(A) == 24 and A[2][1][1] == "C1b" );
		
//		auto const& A2 = A.sliced(0, 3).rotated()[1].rotated(-1).rotated(1).sliced(0, 2).rotated(-1);
		auto const& A2 = A.sliced(0, 3).rotated()[1].sliced(0, 2).unrotated();
		assert( multi::rank<std::decay_t<decltype(A2)>>{} == 2 and num_elements(A2) == 6 );
		assert( std::get<0>(sizes(A2)) == 3 and std::get<1>(sizes(A2)) == 2 );
		for(auto i : std::get<0>(extensions(A2)) ){
			for(auto j : std::get<1>(extensions(A2)) ) cout<< A2[i][j] <<' ';
			cout<<'\n';
		}

		auto const& A3 = A({0, 3}, 1, {0, 2});
		assert( multi::rank<std::decay_t<decltype(A3)>>{} == 2 and num_elements(A3) == 6 );
		for(auto i : std::get<0>(extensions(A3)) ){
			for(auto j : std::get<1>(extensions(A3)) ) cout<< A3[i][j] <<' ';
			cout<<'\n';
		}

#if 0
		assert( dimensionality(A(multi::index_range{0, 2}, 1, multi::index_range{0, 2})) == 2 );
		assert( dimensionality(A(multi::irange{0, 2}, 1, multi::irange{0, 2})) == 2 );
		assert( dimensionality(A({0, 2}, 1, multi::irange{0, 2})) == 2 );
//		assert( dimensionality(A({0, 2}, 1, {0, 2})) == 2 ); // not supported, bracket in third arg

		assert( &A(multi::irange{0, 2}, 1, multi::irange{0, 2})[0][0] == &A[0][1][0] );
		assert(  A(multi::irange{0, 2}, 1, multi::irange{0, 2})[0][0] ==  10 );
		assert( &A(multi::irange{0, 2}, 1, multi::irange{0, 2})[1][1] == &A[1][1][1] );
		assert(  A(multi::irange{0, 2}, 1, multi::irange{0, 2})[1][1] == 111 );

		assert( A.sliced(0, 2).rotated(1)[1].rotated(-1).rotated(2).sliced(0, 2).rotated(-2) == A(multi::irange{0, 2}, 1, multi::irange{0, 2}) );
		// .rotated(-1).rotated(2) == .rotated(1) [cancellation]
		assert( A.sliced(0, 2).rotated(1)[1].rotated(1).sliced(0, 2).rotated(-2) == A(multi::irange{0, 2}, 1, multi::irange{0, 2}) );
		// .rotated(1) == .rotated(), .rotated(-2) == .rotated(0) == identity for something with D = 2
		assert( A.sliced(0, 2).rotated()[1].rotated().sliced(0, 2) == A(multi::irange{0, 2}, 1, multi::irange{0, 2}) );
		assert( A({0, 2}, 1)({0, 2}) == A(multi::irange{0, 2}, 1, multi::irange{0, 2}) );

		assert( A.sliced(0, 2).num_elements() == 12 );
		assert( A.sliced(0, 2).rotated()[1].num_elements() == 6 );
		assert( A.sliced(0, 2).rotated(1)[1].rotated(-1).rotated(2).num_elements() == 6 );
		assert( A.sliced(0, 2).rotated()[1].sliced(0, 2).num_elements() == 4 );
		auto const& A2 = A.sliced(0, 2).rotated()[1].sliced(0, 2);
		assert( A2.num_elements() == 4 );
		for(auto i : std::get<0>(extensions(A2)) ){
			for(auto j : std::get<1>(extensions(A2)) ) cout << A2[i][j] << ' ';
			cout << '\n';
		}
	//	assert( A({0,1}, 1)({0,1})[1][1] == A[1][1][1] );
#endif
	}
	return 0;
	{
		double const dc2D[4][5] = {{1.,2.},{2.,3.}};
//		multi::array_cref<double, 2> acrd2D{&dc2D[0][0], {4, 5}};
		multi::array_ref acrd2D(&dc2D[0][0], {4, 5});
		static_assert( decltype(acrd2D)::dimensionality == 2, "!");
		static_assert( acrd2D.dimensionality == 2, "!");
		static_assert( decltype(acrd2D)::rank{} == 2, "!" );
		static_assert( multi::rank<decltype(acrd2D)>{} == 2, "!" );
		assert( &acrd2D[2][3] == &dc2D[2][3] );
		assert( acrd2D.size() == 4);
		assert( size(acrd2D) == 4 );
		assert( acrd2D.size(0) == 4 );
		assert( acrd2D.size(1) == 5 );
		assert( std::tuple_size<decltype(acrd2D.sizes())>{} == acrd2D.dimensionality );
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

	//	boost::totally_ordered<const_iterator>,
	//	boost::additive<const_iterator, basic_array::difference_type>,
	//	boost::unit_steppable<const_iterator>,
	//	public boost::dereferenceable<const_iterator, const_reference const*>,
	//	public boost::indexable<const_iterator, basic_array::difference_type, const_reference>,
	//	const_iterator(const_reference const& cr, index stride) : 
	//		const_reference{cr}, stride_{stride}{}
	//	const_reference const* operator->() const{return std::addressof(operator*());}
	//	using boost::dereferenceable<const_iterator, const_reference const* const&>::operator->;
	//	const_iterator operator+(index d) const{const_iterator ret = *this; return ret += d;}
	//	const_iterator operator-(index d) const{const_iterator ret = *this; return ret -= d;}
	//	const_reference operator[](index i) const{return *(*this+i);}
		//	return {this->data_+i*stride_, this->layout()};
	//	}
#if 0
	struct const_iterator : 
		private const_reference,
		boost::iterator_facade<
			const_iterator, void, std::random_access_iterator_tag, 
			const_reference const&, difference_type
		>,
		boost::dereferenceable<const_iterator, const_reference const*>,
		boost::indexable<const_iterator, difference_type, const_reference>
	{
		using boost::dereferenceable<const_iterator, const_reference const*>::operator->;
		using boost::indexable<const_iterator, difference_type, const_reference>::operator[];
		bool operator==(const_iterator const& o) const{return equal(o);}
		bool operator!=(const_iterator const& o) const{return not equal(o);}
		using typename boost::iterator_facade<
			const_iterator, void, std::random_access_iterator_tag, 
			const_reference const&, difference_type
		>::difference_type;
		using typename boost::iterator_facade<
			const_iterator, void, std::random_access_iterator_tag, 
			const_reference const&, difference_type
		>::reference;
		using value_type = typename const_reference::value_type;
		using pointer = const_reference const*;
		const_iterator(std::nullptr_t p = nullptr) : const_reference{p}{}
		template<class Other> const_iterator(Other const& o) 
			: const_reference{o.data_, o.layout()}, stride_{o.stride_}{}
	//	const_iterator(const_iterator const& o)
	//		: const_reference{o.data_, o.layout()}, stride_{o.stride_}{}
		const_iterator& operator=(const_iterator const& other){
			this->data_ = other.data_;
			static_cast<typename Layout::sub_t&>(*this) = other; //impl_) = other.impl_;
			stride_ = other.stride_;
			return *this;
		}
		explicit operator bool() const{return this->data_;}
		const_reference const& operator*() const{assert(*this); return *this;}
	private:
		index stride_;
		const_iterator(typename const_iterator::element_ptr p, typename Layout::sub_t sub, index stride)
			: const_reference{p, sub}, stride_{stride}{}
		friend const_iterator basic_array::cbegin() const;
		friend const_iterator basic_array::cend() const;
		using const_reference::layout;
		using const_reference::data_;
		bool equal(const_iterator const& o) const{
			return data_==o.data_ && stride_==o.stride_ && layout()==o.layout();
		}
		void increment(){data_ += stride_;}
		void decrement(){data_ -= stride_;}
		void advance(difference_type n){data_ += stride_*n;}
		difference_type distance_to(const_iterator const& other) const{
			assert( stride_ == other.stride_ );
			assert( stride_ != 0 and (other.data_ - data_)%stride_ == 0);
			assert( layout() == other.layout() );
			return (other.data_ - data_)/stride_;
		}
	    friend class boost::iterator_core_access;
	};
#endif
#if 0
	struct iterator : 
		private reference
	{
		index stride_;
		iterator(basic_array<T, D-1, element_ptr> const& cr, index stride) : 
			basic_array<T, D-1, element_ptr>{cr}, stride_{stride}
		{}
		iterator(element_ptr p, typename Layout::sub_t sub, index stride)
			: const_reference{p, sub}, stride_{stride}{}
		friend iterator basic_array::begin() const;
		friend iterator basic_array::end() const;
		friend struct basic_array;
	public:
		iterator(std::nullptr_t np = nullptr) : reference{np}{}
		explicit operator bool() const{return this->data_;}
		using difference_type = typename basic_array<T, D, ElementPtr>::difference_type;
		using value_type = typename basic_array<T, D, ElementPtr>::value_type;
		using pointer = void*;
		using reference = typename basic_array<T, D, element_ptr /*ElementPtr*/>::reference; // careful with shadowing reference (there is another one level up)
		using iterator_category = std::random_access_iterator_tag;
	//	iterator& operator=(iterator const& other) = default;
		iterator& operator=(iterator const& other){
			this->data_ = other.data_;
			static_cast<typename Layout::sub_t&>(*this) = other;
			stride_ = other.stride_;
			return *this;
		}
		template<class It>
		bool operator==(It const& other) const{
			return this->data_ == other.data_ and this->stride_ == other.stride_;
		}
		template<class It>
		bool operator!=(It const& other) const{return not((*this)==other);}
		reference& operator*() const{return *const_cast<iterator*>(this);}
		reference* operator->() const{return static_cast<reference*>(this);}
		reference operator[](index i) const{return {this->data_ + i*stride_, this->layout()};}
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
	private:
	};
#endif
//	template<class Iterator> 
//	struct basic_reverse_iterator : std::reverse_iterator<Iterator>{
//	};
#if 0	
	struct const_reverse_iterator : std::reverse_iterator<const_iterator>{
		template<class Other, typename = decltype(const_iterator{std::declval<Other const&>().base()})> 
		const_reverse_iterator(Other const& o) : std::reverse_iterator<const_iterator>(o.base()){}
		const_reverse_iterator() : std::reverse_iterator<const_iterator>(){}
		explicit const_reverse_iterator(const_iterator it) : std::reverse_iterator<const_iterator>(--it){}
		operator const_iterator() const{auto ret = this->base(); ++ret; return ret;}
		explicit operator bool() const{return static_cast<bool>(this->current);}
		typename const_iterator::reference operator*() const{return this->current;}
		typename const_iterator::pointer operator->() const{return &this->current;}
		typename const_iterator::reference operator[](typename iterator::difference_type n) const{return *(this->current - n);}
	};
	struct reverse_iterator : std::reverse_iterator<iterator>{
		template<class Other, typename = decltype(iterator{std::declval<Other const&>().base()})> 
		reverse_iterator(Other const& o) : std::reverse_iterator<iterator>(o.base()){}
		explicit reverse_iterator(iterator it) : std::reverse_iterator<iterator>(--it){}
		operator iterator() const{auto ret = this->base(); ++ret; return ret;}
		explicit operator bool() const{return static_cast<bool>(this->current);}
		typename iterator::reference operator*() const{return this->current;}
		typename iterator::pointer operator->() const{return &this->current;}
		typename iterator::reference operator[](typename iterator::difference_type n) const{return *(this->current - n);}
	};
#endif
//	using types::dimensionality;
#if 0
	using element = T;
	constexpr static dimensionality_type dimensionality = D;
	using element_ptr = std::decay_t<ElementPtr>;
	using element_const_ptr = typename std::pointer_traits<element_ptr>::template rebind<element const>;
	using layout_t = Layout;
	using index = multi::index;
	using value_type = array<element, dimensionality-1>;
	using decay_type = array<element, dimensionality, typename pointer_traits<element_ptr>::allocator_type>;
	using difference_type = index;
	using const_reference = basic_array<element, dimensionality-1, element_const_ptr>;
	using reference       = basic_array<element, dimensionality-1, element_ptr>;
#endif
//	using element = T;
//	constexpr static dimensionality_type dimensionality = 1;
//	using element_ptr = std::decay_t<ElementPtr>;
//	using element_const_ptr = typename std::pointer_traits<typename types::element_ptr>::template rebind<typename types::element const>;
//	using layout_t = Layout;
//	using index = multi::index;
//	using value_type = T; // array<element, D-1>;
//	using decay_type = array<typename types::element, types::dimensionality, typename pointer_traits<typename types::element_ptr>::allocator_type>;
//	using difference_type = index;
//	using const_reference = typename pointer_traits<typename types::element_ptr>::element_type const&;// basic_array<element, D-1, element_const_ptr>;
//	using reference       = typename pointer_traits<typename types::element_ptr>::element_type&;//basic_array<element, D-1, element_ptr>;

//	constexpr static dimensionality_type dimensionality = 1;
//	using layout_t = Layout;
//	using value_type = T;
//	using element = value_type;
//	using element_ptr = std::decay_t<ElementPtr>;
//	using element_const_ptr = typename std::pointer_traits<element_ptr>::template rebind<element const>;
//	using const_reference = decltype(*std::declval<element_const_ptr>());
//	using reference = decltype(*std::declval<element_ptr>());
//	basic_array(basic_array const& o) : types{o}{}//types{o}, data_{o.data_}{}
//	basic_array(basic_array&& o) noexcept : types{std::move(o)}{}
//	basic_array(basic_array&& o) noexcept:types{std::move(o)}, data_{o.data_}{};
//	typename types::element_ptr data_; //	ElementPtr data_; // TODO call it base_ ?
//	basic_array(typename types::element_ptr data, Layout layout) : types{data, layout}{}//types{layout}, data_{data}{}
#if 0
	struct iterator : const_iterator{
		friend struct basic_array;
		using const_iterator::const_iterator;
		reference operator*() const{return *this->data_;}
		element_ptr const& operator->() const{return this->data_;}
		iterator& operator++(){const_iterator::operator++(); return *this;}
		iterator& operator+=(ptrdiff_t d){const_iterator::operator+=(d); return *this;}
		iterator operator+(ptrdiff_t d) const{return iterator(*this)+=d;}
	public:
	//	element_ptr get() const{return const_iterator::data_;}
		element_ptr base() const{return this->data_;}
	};
#endif
#if 0
	struct iterator : array_cref<T, D, ElementPtr>::iterator{
		using array_cref<T, D, ElementPtr>::iterator::iterator;
//		decltype(auto) cdata() const{
//			return array_cref<T, D, ElementPtr>::iterator::cbase();}
		decltype(auto) data() const{return this->data_;}
		friend decltype(auto) data(iterator const& self){return self.data();}
	//	typename array_ref::reference operator*() const{return *this->data_;}
	//	typename array_ref::element_ptr const& operator->() const{return const_cast<typename array_ref::element_ptr const&>(this->data_);}
		operator typename array_ref::element_ptr() const{return this->data_;}
//		friend decltype(auto) cdata(iterator const& self){return self.cdata();}
	};
#endif
//	typename array_ref::iterator begin() const&{return iterator{basic_array<T, D, ElementPtr>::begin()};}
//	typename array_ref::iterator end() const&{return iterator{basic_array<T, D, ElementPtr>::end()};}
//	typename array_ref::iterator begin() &&{return iterator{basic_array<T, D, ElementPtr>::begin()};}
//	typename array_ref::iterator end() &&{return iterator{basic_array<T, D, ElementPtr>::end()};}
//	struct const_iterator : basic_array<T, D, ElementPtr>::const_iterator{
//		using basic_array<T, D, ElementPtr>::const_iterator::const_iterator;
//		const_iterator(typename basic_array<T, D, ElementPtr>::const_iterator const& o) : basic_array<T, D, ElementPtr>::const_iterator{o}{}
//		decltype(auto) cdata() const{
//			return basic_array<T, D, ElementPtr>::const_iterator::cbase();}
//		decltype(auto) data() const{return cdata();}
//		friend decltype(auto) data(const_iterator const& self){return self.data();}
//		friend decltype(auto) cdata(const_iterator const& self){return self.cdata();}
//	};
//	const_iterator cbegin() const{return {basic_array<T, D, ElementPtr>::cbegin()};}
//	const_iterator cend() const{return {basic_array<T, D, ElementPtr>::cend()};}
//	using iterator = const_iterator;
//	iterator begin() const{return cbegin();}
//	iterator end() const{return cend();}
#if 0
	template<class Ref>	
	struct basic_iterator : 
		private Ref,
		boost::iterator_facade<
			basic_iterator<Ref>, void, std::random_access_iterator_tag, 
			Ref const&, difference_type
		>,
		boost::totally_ordered<basic_iterator<Ref>>
	{
		using difference_type = typename basic_array::difference_type;
		using reference = Ref const&;
		using value_type = typename Ref::decay_type;
		using pointer = Ref const*;
		basic_iterator(std::nullptr_t p = nullptr) : Ref{p}, stride_{0}{}
		template<class Other, typename = decltype(typename Ref::element_ptr{typename Other::element_ptr{}})> 
		basic_iterator(Other const& o) 
			: Ref{o.base_, o.layout()}, stride_{o.stride_}{}
		basic_iterator& operator=(basic_iterator const& other){
			this->base_ = other.base_;
			static_cast<typename Layout::sub_t&>(*this) = other; //impl_) = other.impl_;
			stride_ = other.stride_;
			return *this;
		}
		explicit operator bool() const{return this->base_;}
		Ref const& operator*() const{assert(*this); return *this;}
		Ref const* operator->() const{assert(*this); return this;}
		Ref        operator[](difference_type n) const{return *(*this + n);}
		template<class O> bool operator==(O const& o) const{return equal(o);}
		bool operator<(basic_iterator const& o) const{return distance_to(o) > 0;}
	private:
		index stride_;
		basic_iterator(typename Ref::element_ptr p, typename Layout::sub_t sub, index stride)
			: Ref{p, sub}, stride_{stride}{}
		friend struct basic_array;
		using Ref::layout;
		using Ref::base_;
		bool equal(basic_iterator const& o) const{return base_==o.base_ && stride_==o.stride_ && layout()==o.layout();}
		void increment(){base_ += stride_;}
		void decrement(){base_ -= stride_;}
		void advance(difference_type n){base_ += stride_*n;}
		difference_type distance_to(basic_iterator const& other) const{
			assert( stride_ == other.stride_ and stride_ != 0 and (other.base_ - base_)%stride_ == 0 and layout() == other.layout() );
			return (other.base_ - base_)/stride_;
		}
	    friend class boost::iterator_core_access;
	};
#endif

//template<class I, class T> 
//index_range extension(T const&, I d = {}){assert(d == -1); assert(0); return {};}

