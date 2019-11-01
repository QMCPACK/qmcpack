#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && cuda-g++ -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_MULTI_LAYOUT $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_LAYOUT_HPP
#define MULTI_LAYOUT_HPP

//#include "detail.hpp"
#include "index_range.hpp"
#include "types.hpp"

#include<array>
#include<cassert>
#include<cstddef>
#include<type_traits> // make_signed_t
#include<tuple> //make_tuple
//#include<boost/operators.hpp>
#include "../detail/operators.hpp"

namespace boost{
namespace multi{

namespace detail{

template<typename T, typename... As>
inline void construct_from_initializer_list(T* p, As&&... as){
	::new(static_cast<void*>(p)) T(std::forward<As>(as)...);
}

template<class To, class From, size_t... I>
constexpr auto to_tuple_impl(std::initializer_list<From> il, std::index_sequence<I...>){
	return std::make_tuple(To{il.begin()[I]}...);
}

template<size_t N, class To, class From>
constexpr auto to_tuple(std::initializer_list<From> il){
	return il.size()==N?to_tuple_impl<To>(il, std::make_index_sequence<N>()):throw 0;
}

template<class To, class From, size_t... I>
constexpr auto to_tuple_impl(std::array<From, sizeof...(I)> arr, std::index_sequence<I...>){
	return std::make_tuple(To{std::get<I>(arr)}...);
}

template<class To, size_t N, class From>
constexpr auto to_tuple(std::array<From, N> arr){
	return to_tuple_impl<To, From>(arr, std::make_index_sequence<N>());
}

template <class TT, class Tuple, std::size_t... I>
constexpr std::array<TT, std::tuple_size<std::decay_t<Tuple>>{}> to_array_impl(
	Tuple&& t, std::index_sequence<I...>
){
	return {static_cast<TT>(std::get<I>(std::forward<Tuple>(t)))...};
}
 
template<class T = void, class Tuple, class TT = std::conditional_t<std::is_same<T, void>{}, std::decay_t<decltype(std::get<0>(std::decay_t<Tuple>{}))>, T> >
constexpr std::array<TT, std::tuple_size<std::decay_t<Tuple>>{}> 
to_array(Tuple&& t){
	return 
		to_array_impl<TT>(
			std::forward<Tuple>(t),
			std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple>>{}>{}	
		)
	;
}

template <class Tuple, std::size_t... Ns>
auto tuple_tail_impl(Tuple&& t, std::index_sequence<Ns...>){
   return std::forward_as_tuple(std::forward<decltype(std::get<Ns + 1>(t))>(std::get<Ns + 1>(t))...);
}

template<class Tuple>
auto tuple_tail(Tuple&& t)
->decltype(tuple_tail_impl(t, std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{} - 1>())){//std::tuple<Ts...> t){
	return tuple_tail_impl(t, std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{} - 1>());}

#if 0
template<typename Tuple, size_t... I, size_t S = sizeof...(I)>
auto reverse_impl(Tuple&& t, std::index_sequence<I...>)
    return std::make_tuple(std::get<S - 1 - I>(std::forward<Tuple>(t))...);
}

template<typename Tuple, size_t S = std::tuple_size<std::decay_t<Tuple>>{} >
auto reverse(Tuple&& t){
	return reverse_impl(std::forward<Tuple>(t), std::make_index_sequence<S>());
}
#endif

}

}}

namespace boost{
namespace multi{

#if 0
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
#endif

struct f_tag{};

template<dimensionality_type D>
struct layout_t
	: boost::multi::equality_comparable2<layout_t<D>, void>
{
//	static constexpr dimensionality_type dimensionality = D;
	using dimensionality_type = multi::dimensionality_type;
	using rank = std::integral_constant<dimensionality_type, D>;
	static constexpr dimensionality_type dimensionality(){return rank{};}
	friend constexpr dimensionality_type dimensionality(layout_t const& l){return l.dimensionality();}
	using sub_t = layout_t<D-1>;
	using size_type = multi::size_type;
	using index = multi::index;
	using difference_type = multi::difference_type;
	using index_extension = multi::index_extension;
	using index_range = multi::range<index>;
 	sub_t sub;
	index stride_ = 1;
	index offset_;
	index nelems_;
	struct extensions_type_ : std::decay_t<decltype(tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename sub_t::extensions_type::base_>()))>{
		using base_ = std::decay_t<decltype(tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename sub_t::extensions_type::base_>()))>;
		using base_::base_;
		extensions_type_() = default;
		template<class Array, typename = decltype(std::get<D-1>(std::declval<Array>()))> constexpr extensions_type_(Array const& t) : extensions_type_(t, std::make_index_sequence<D>{}){}
	//	template<class T, std::size_t N> extensions_type_(std::array<T, N> const& t) : extensions_type_(t, std::make_index_sequence<N>{}){}
		base_ const& base() const{return *this;}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
	private:
		template<class Array, std::size_t... I, typename = decltype(base_{std::get<I>(std::declval<Array const&>())...})> constexpr extensions_type_(Array const& t, std::index_sequence<I...>) : base_{std::get<I>(t)...}{}
	//	template<class T, std::size_t N, std::size_t... I> extensions_type_(std::array<T, N> const& t, std::index_sequence<I...>) : extensions_type_{std::get<I>(t)...}{}
	};
	using extensions_type = extensions_type_;
//	using extensions_type = decltype(tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename sub_t::extensions_type>()));
	using strides_type    = decltype(tuple_cat(std::make_tuple(std::declval<index>()), std::declval<typename sub_t::strides_type>()));
//	using extensions_type = typename detail::repeat<index_extension, D>::type;
//	using extensions_io_type = std::array<index_extension, D>;
	auto operator()(index i) const{return i*stride_ - offset_;}
	auto origin() const{return sub.origin() - offset_;}
	constexpr layout_t(index_extension const& ie, layout_t<D-1> const& s) : 
		sub{s},
		stride_{ie.size()*num_elements(sub)!=0?size(sub)*stride(sub):1}, // use .size for nvcc
		offset_{0}, 
		nelems_{ie.size()*num_elements(sub)}                             // use .size fort
	{}
	constexpr layout_t(extensions_type const& e = {}) : 
		sub{detail::tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{}
#if defined(__INTEL_COMPILER) or (defined(__GNUC) && (__GNUC<6))
//	constexpr layout_t(std::initializer_list<index_extension> il) noexcept : layout_t{multi::detail::to_tuple<D, index_extension>(il)}{}
//	constexpr layout_t(std::initializer_list<index> il) noexcept : layout_t{multi::detail::to_tuple<D, index_extension>(il)}{}
	constexpr layout_t(std::array<index_extension, D> x) noexcept : layout_t{multi::detail::to_tuple<index_extension>(x)}{}
#endif
#if 0
	constexpr layout_t(f_tag, extensions_type const& e = {}) : 
		sub{detail::tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{
		reverse();
	}
	void reverse(){
		unrotate();
		sub.reverse();
	}
#endif
/*	template<class Ext>
	constexpr layout_t(Ext const& e) :
		sub{detail::tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{}*/
	template<class StdArray, typename = std::enable_if_t<std::is_same<StdArray, std::array<index_extension, D>>{}> >
	constexpr layout_t(StdArray const& e) : 
		sub{detail::tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{}
//	template<class Extensions//, typename = std::enable_if_t<
	//	not std::is_base_of<layout_t, std::decay_t<Extensions>
	//	and not std::is_same<std::decay_t<Extensions>	
	//>{}> 
//	>
//	constexpr layout_t(Extensions const& e) : layout_t{extensions_type{e}}{}
	//	sub{detail::tail(std::forward<Extensions>(e))},
	//	stride_{index_extension(std::get<0>(std::forward<Extensions>(e))).size()*sub.num_elements()!=0?sub.size()*sub.stride():1},
	//	offset_{0},
	//	nelems_{index_extension(std::get<0>(std::forward<Extensions>(e))).size()*sub.num_elements()}
//	{}
//	constexpr layout_t(extents_t<D> const& e) :
//		stride_{head(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
//		offset_{0}, 
//		nelems_{head(e).size()*sub.num_elements()}
//	{}
	constexpr index nelems() const{return nelems_;}
	friend constexpr index nelems(layout_t const& self){return self.nelems();}
	auto nelems(dimensionality_type d) const{return d?sub.nelems(d-1):nelems_;}
	constexpr bool operator==(layout_t const& o) const{
		return sub==o.sub and stride_==o.stride_ and offset_==o.offset_ and nelems_==o.nelems_;
	}
//	constexpr bool operator!=(layout_t const& o) const{return not(*this==o);}
	constexpr size_type num_elements() const{return size()*sub.num_elements();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return not nelems_;} friend
	constexpr bool empty(layout_t const& s){return s.empty();}
	constexpr size_type size() const {return nelems_/stride_;} 
//	friend constexpr size_type size(layout_t const& l){return l.size();} // removed for gcc6, add "using boost::multi::size"
	friend constexpr size_type size(layout_t const& l){return l.size();}
	size_type size(dimensionality_type d) const{return d?sub.size(d-1):size();}

	constexpr index stride() const{return stride_;}
	index stride(dimensionality_type d) const{return d?sub.stride(d-1):stride();}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	constexpr strides_type strides() const{return tuple_cat(std::make_tuple(stride()), sub.strides());}
	friend constexpr strides_type strides(layout_t const& self){return self.strides();}
	constexpr index offset() const{return offset_;}
	constexpr index offset(dimensionality_type d) const{return d?sub.offset(d-1):offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offsets() const{return tuple_cat(std::make_tuple(offset()), sub.offsets());}
	constexpr auto base_size() const{using std::max; return max(nelems_, sub.base_size());}
	auto is_compact() const{return base_size() == num_elements();}
	decltype(auto) shape() const{return sizes();}
	friend decltype(auto) shape(layout_t const& self){return self.shape();}
	constexpr auto sizes() const{return tuple_cat(std::make_tuple(size()), sub.sizes());}
	template<class T = void>
	constexpr auto sizes_as() const{return detail::to_array<T>(sizes());}
private:
	friend struct layout_t<D+1>; // void layout_t<D+1>::strides_aux(size_type*) const;
public:
	constexpr index_extension extension() const{
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	friend auto extension(layout_t const& self){return self.extension();}
	constexpr index_extension extension_aux() const{
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	template<dimensionality_type DD = 0>
	constexpr index_extension extension(dimensionality_type d) const{
		return d?sub.extension(d-1):extension();
	}
	constexpr extensions_type extensions() const{return tuple_cat(std::make_tuple(extension()), sub.extensions().base());}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
	void extensions_aux(index_extension* it) const{
		*it = extension();
		++it;
		sub.extensions_aux(it);
	}
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

};

template<class T, class Layout>
constexpr auto sizes_as(Layout const& self)
->decltype(self.template sizes_as<T>()){
	return self.template sizes_as<T>();}

template<>
struct layout_t<dimensionality_type{0}>{
	using rank = std::integral_constant<dimensionality_type, 0>;
	static constexpr dimensionality_type dimensionality = 0;
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
	using difference_type = multi::difference_type;
	struct extensions_type_ : std::tuple<>{
		using std::tuple<>::tuple;
		using base_ = std::tuple<>;
		extensions_type_(base_ const& b) : base_(b){}
		extensions_type_() = default;
		base_ const& base() const{return *this;}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
	};
	using extensions_type = extensions_type_;
	extensions_type extensions() const{return extensions_type{};}
	friend auto extensions(layout_t const& self){return self.extensions();}
	auto sizes() const{return std::tuple<>{};}
	friend auto sizes(layout_t const& s){return s.sizes();}
};

template<>
struct layout_t<dimensionality_type{1}>{
	using sub_t = layout_t<dimensionality_type{0}>;
//	typedef boost::multi::dimensionality_type dimensionality_type;
	using rank = std::integral_constant<dimensionality_type, 1>;
	static constexpr dimensionality_type dimensionality = rank{};
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
	using index_extension = multi::index_extension;

//	using extensions_type = std::array<index_extension, 1>;
	using index = multi::index;
	using size_type = multi::size_type;
	using index_range = multi::range<index>;
	using difference_type = multi::difference_type;
	index stride_ = 1;
	index offset_;
	index nelems_;

	struct extensions_type_ : std::tuple<index_extension>{
		using std::tuple<index_extension>::tuple;
		using base_ = std::tuple<index_extension>;
		extensions_type_(index_extension const& ie) : base_{ie}{}
		extensions_type_() = default;
		base_ const& base() const{return *this;}
		template<class Array, typename = decltype(std::get<0>(std::declval<Array const&>()))> extensions_type_(Array const& arr) : base_{std::get<0>(arr)}{}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
	};
//		typedef std::tuple<index_extension> base_;
//		using base_::base_;
//		template<class T, std::size_t N> extensions_type_(std::array<T, N> const& t) : extensions_type_(t, std::make_index_sequence<N>{}){}
//	private:
//		template<class T, std::size_t N, std::size_t... I> extensions_type_(std::array<T, N> const& t, std::index_sequence<I...>) : extensions_type_{std::get<I>(t)...}{}
//	};
	using extensions_type = extensions_type_;
//	using extensions_type = std::tuple<index_extension>;//typename detail::repeat<index_extension, 1u>::type;
	using strides_type = std::tuple<index>;
//	layout_t() = default;
	layout_t() : stride_{1}, offset_{0}, nelems_{0}{}
	constexpr layout_t(index_extension ie, layout_t<0> const&) : 
		stride_{1}, offset_{0}, nelems_{ie.size()}
	{}
#if defined(__INTEL_COMPILER)
	constexpr layout_t(std::initializer_list<index_extension> il) noexcept : layout_t{multi::detail::to_tuple<1, typename layout_t::index_extension>(il)}{}
	constexpr layout_t(std::initializer_list<index> il) noexcept : layout_t{multi::detail::to_tuple<1, typename layout_t::index_extension>(il)}{}
#endif
	template<class Extensions, typename = decltype(std::get<0>(Extensions{}))>
	constexpr layout_t(Extensions e) : 
		stride_{1*1}, offset_{0}, nelems_{std::get<0>(e).size()*1}
	{}
	constexpr auto offset() const{return offset_;}	
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offset(dimensionality_type d) const{return d==0?offset_:throw 0;}
	constexpr auto nelems() const{return nelems_;}
	constexpr auto nelems(dimensionality_type d) const{
		return d==0?nelems_:throw 0;
	}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr size_type size() const{
		return nelems_/stride_; // assert(stride_!=0 and nelems_%stride_ == 0)
	}
	constexpr size_type size(dimensionality_type d) const{
		return d==0?nelems_/stride_:throw 0; // assert(d == 0 and stride_ != 0 and nelems_%stride_ == 0);
	}
	friend constexpr size_type size(layout_t const& self){return self.size();}
	constexpr auto base_size() const{return nelems_;}
	auto is_compact(){return base_size() == num_elements();}
public:
	constexpr auto stride(dimensionality_type d = 0) const{return d?throw 0:stride_;}
//		return d==0?stride_:throw 0;
//	}
	friend constexpr index stride(layout_t const& self){return self.stride();}
public:
	constexpr auto strides() const{return std::make_tuple(stride());}
//	template<class T = void>
//	constexpr auto sizes() const{
//		return detail::to_array<T>(tuple_cat(std::make_tuple(size()), sub.sizes()));
//	}
	constexpr auto sizes() const{return std::make_tuple(size());}
	template<class T = void>
	constexpr auto sizes_as() const{
		return detail::to_array<T>(sizes());
	}
	constexpr auto offsets() const{return std::make_tuple(offset());}

	constexpr size_type num_elements() const{return this->size();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return nelems_ == 0;}
	friend bool empty(layout_t const& s){return s.empty();}
	constexpr index_extension extension() const{
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	friend auto extension(layout_t const& self){return self.extension();}
	constexpr auto extension(dimensionality_type d) const{
		return d==0?index_extension{offset_/stride_, (offset_ + nelems_)/stride_}:throw 0;
	//	assert(stride_ != 0 and nelems_%stride_ == 0);
	//	return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	constexpr extensions_type extensions() const{return extensions_type{extension()};}//std::make_tuple(extension());}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
private:
	friend struct layout_t<2u>;
	void strides_aux(size_type* it) const{*it = stride();}
	void sizes_aux(size_type* it) const{*it = size();}
	void offsets_aux(index* it) const{*it = offset();}
	void extensions_aux(index_extension* it) const{*it = extension();}
public:
	constexpr auto operator()(index i) const{return i*stride_ - offset_;}
	constexpr auto origin() const{return -offset_;}
	constexpr bool operator==(layout_t const& other) const{
		return stride_==other.stride_ and offset_==other.offset_ and nelems_==other.nelems_;
	} 
	constexpr bool operator!=(layout_t const& other) const{return not(*this==other);}
	layout_t& rotate(){return *this;}
	layout_t& unrotate(){return *this;}
//	decltype(auto) shape() const{return sizes();}
};



}}

#if _TEST_MULTI_LAYOUT

#include<cassert>
#include<iostream>
#include<vector>

#include "../../multi/utility.hpp"

using std::cout;
namespace multi = boost::multi;

int main(){
//	auto t = std::make_tuple(1.,2.,3.);
//	auto u = multi::detail::reverse(t);
//	assert( std::get<0>(u) == 3. );
	auto t = multi::detail::to_tuple<3, multi::index_extension>({1,2,3});
	assert( std::get<1>(t) == 2 );
	std::array<multi::index, 3> arr{1,2,3};
	auto u = multi::detail::to_tuple<multi::index_extension>(arr);
	assert( std::get<1>(u) == 2 );
	
 {  multi::layout_t<1> L{}; assert( dimensionality(L)==1 and num_elements(L) == 0 and size(L) == 0 and size(extension(L))==0 and stride(L)!=0 and empty(L) );
}{  multi::layout_t<2> L{}; assert( dimensionality(L)==2 and num_elements(L) == 0 and size(L) == 0 and size(extension(L))==0 and stride(L)!=0 and empty(L) );
}{  multi::layout_t<3> L{}; assert( num_elements(L) == 0 );
}{	multi::layout_t<3> L({{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({{10}, {10}, {10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({10, 10, 10}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({multi::index_extension{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L(std::array<multi::index_extension, 3>{{ {0,10}, {0,10}, {0,10} }}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000);
}{	multi::layout_t<3> L(std::make_tuple(multi::iextension{0, 10}, multi::iextension{0, 10}, multi::iextension{0, 10})); assert(L.num_elements() == 1000);
}{	multi::layout_t<3> L(std::make_tuple(multi::iextension{10}, multi::iextension{10}, multi::iextension{10})); assert( num_elements(L) == 1000);
}{	multi::layout_t<3> L(std::make_tuple(10, 10, multi::iextension{10})); assert( num_elements(L) == 1000 );
}{
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}}); 
	multi::layout_t<3> L2{extensions(L)};

	assert( L.stride() == 20*30 );
	assert( L.offset() == 0 );
	assert( L.nelems() == 10*20*30 );
	
	assert( L.stride(0) == L.stride() );
	assert( L.offset(0) == L.offset() );
	assert( L.nelems(0) == L.nelems() );

	assert( L.stride(1) == 30 );
	assert( L.offset(1) == 0 );
	assert( L.nelems(1) == 20*30 );
	
	assert( L.stride(2) == 1 );
	assert( L.offset(2) == 0 );
	assert( L.nelems(2) == 30 );
	
	assert( L.num_elements() == 10*20*30 );
	assert( L.size() == 10 );
	assert( L.extension().first() == 0 );
	assert( L.extension().last() == 10 );

	assert( L.size(1) == 20 );
	assert( L.extension(1).first() == 0 );
	assert( L.extension(1).last() == 20 );

	assert( L.size(2) == 30 );
	assert( L.extension(2).first() == 0 );
	assert( L.extension(2).last() == 30 );

	assert( std::get<0>(L.strides()) == L.stride(0) );
	assert( std::get<1>(L.strides()) == L.stride(1) );
	assert( std::get<2>(L.strides()) == L.stride(2) );

	auto const& strides = L.strides();
	assert( std::get<0>(strides) == L.stride(0) );
	
}{
	constexpr multi::layout_t<3> L( {{0, 10}, {0, 20}, {0, 30}} );
	static_assert( L.stride() == 20*30, "!");
	using std::get;
	static_assert( get<0>(L.strides()) == 20*30, "!");
	static_assert( get<1>(L.strides()) == 	30, "!");
	static_assert( get<2>(L.strides()) == 	 1, "!");
	static_assert( L.size() == 10, "!");
}{
	std::tuple<int, int, int> ttt = {1,2,3};
	auto arrr = boost::multi::detail::to_array(ttt);
	assert(arrr[1] == 2);
}

}
#endif
#endif

