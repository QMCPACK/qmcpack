#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o$0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
#ifndef MULTI_LAYOUT_HPP
#define MULTI_LAYOUT_HPP
// © Alfredo A. Correa 2018-2020

#include "types.hpp"

#include "../detail/operators.hpp"
#include "../config/NODISCARD.hpp"

#include<type_traits> // make_signed_t

#include<limits>

#ifndef HD
#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD 
#endif
#endif

namespace boost{
namespace multi{

namespace detail{

template<typename T, typename... As>
inline void construct_from_initializer_list(T* p, As&&... as){
	::new(static_cast<void*>(p)) T(std::forward<As>(as)...);
}

template<class To, class From, size_t... I>
constexpr auto to_tuple_impl(std::initializer_list<From> il, std::index_sequence<I...>){
	(void)il;
	return std::make_tuple(To{il.begin()[I]}...);
}

//template<size_t N, class To, class From>
//constexpr auto to_tuple(std::initializer_list<From> il){
//	return il.size()==N?to_tuple_impl<To>(il, std::make_index_sequence<N>()):throw 0;
//}

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

}

}}

namespace boost{
namespace multi{

struct f_tag{};

template<dimensionality_type D, typename SSize=multi::size_type> struct layout_t;

template<typename SSize>
struct layout_t<dimensionality_type{0}, SSize>{
	using size_type = SSize;
	using difference_type = std::make_signed_t<size_type>;
	using index_extension = multi::index_extension;
	using index = difference_type;
	using stride_type=index;
	using offset_type=index;
	using nelems_type=index;
	using index_range = multi::range<index>;
	using rank = std::integral_constant<dimensionality_type, 0>;
	static constexpr dimensionality_type dimensionality = 0;
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
	using strides_type    = std::tuple<>;
	nelems_type nelems_ = 1;//std::numeric_limits<nelems_type>::max(); // 1
	void* stride_ = nullptr;
	void* sub = nullptr;
	struct extensions_type_ : std::tuple<>{
		using std::tuple<>::tuple;
		using base_ = std::tuple<>;
		extensions_type_(base_ const& b) : base_(b){}
		extensions_type_() = default;
		base_ const& base() const{return *this;}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
		operator nelems_type() const{return 1;}
		template<class Archive> void serialize(Archive&, unsigned){}
		static constexpr dimensionality_type dimensionality = 0;
	};
	constexpr layout_t(extensions_type_ const& = {}){}// : nelems_{1}{}
	using extensions_type = extensions_type_;
	constexpr extensions_type extensions() const{return extensions_type{};}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
	constexpr auto sizes() const{return std::tuple<>{};}
	friend auto sizes(layout_t const& s){return s.sizes();}
	nelems_type num_elements() const{return 1;}//nelems_;}
};

template<typename SSize>
struct layout_t<dimensionality_type{1}, SSize>{
	using size_type=SSize; 
	using difference_type=std::make_signed_t<size_type>;
	using index_extension = multi::index_extension;
	using index = difference_type;
	using stride_type=index; 
	using offset_type=index; 
	using nelems_type=index;
	using index_range = multi::range<index>;
	using rank = std::integral_constant<dimensionality_type, 1>;
	static constexpr dimensionality_type dimensionality = rank{};
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
	using sub_t = layout_t<dimensionality_type{0}, SSize>;
protected:
public:
	stride_type stride_ = 1;//std::numeric_limits<stride_type>::max(); 
	offset_type offset_ = 0; 
	nelems_type nelems_ = 0;
	struct extensions_type_ : std::tuple<index_extension>{
		using std::tuple<index_extension>::tuple;
		using base_ = std::tuple<index_extension>;
		extensions_type_(index_extension const& ie) : base_{ie}{}
		extensions_type_() = default;
		base_ const& base() const{return *this;}
		extensions_type_(std::tuple<index_extension> const& t) : std::tuple<index_extension>(t){}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
		template<class Archive>
		void serialize(Archive& ar, unsigned){ar & make_nvp("extension", std::get<0>(*this));}
		static constexpr dimensionality_type dimensionality = 1;
	};
	using extensions_type = extensions_type_;
	using strides_type = std::tuple<index>;
	constexpr layout_t() = default;
	constexpr layout_t(layout_t const&) = default;
	constexpr layout_t(index_extension ie, layout_t<0> const&) :
		stride_{1},
		offset_{ie.first()},
		nelems_{
//			ie.size()<=1?ie.size()*std::numeric_limits<stride_type>::max():ie.size()
			ie.size()
		}{}
	constexpr layout_t(extensions_type e) : layout_t(std::get<0>(e), {}){}
	constexpr layout_t(stride_type stride, offset_type offset, nelems_type nelems) : 
		stride_{stride}, offset_{offset}, nelems_{nelems}
	{}
#if defined(__INTEL_COMPILER)
//	constexpr layout_t(std::initializer_list<index_extension> il) noexcept : layout_t{multi::detail::to_tuple<1, typename layout_t::index_extension>(il)}{}
//	constexpr layout_t(std::initializer_list<index> il) noexcept : layout_t{multi::detail::to_tuple<1, typename layout_t::index_extension>(il)}{}
#endif
	template<class Extensions, typename = decltype(std::get<0>(Extensions{}).size())>
	constexpr layout_t(Extensions const& e)
	: 
	//	stride_{std::get<0>(e).size()<=1?std::numeric_limits<stride_type>::max():1},
	//	nelems_{
		//	std::get<0>(e).size()<=1?std::get<0>(e).size()*std::numeric_limits<stride_type>::max():std::get<0>(e).size()
		//	ie.size()
	//	}{
		nelems_{std::get<0>(e).size()}{
	}
	constexpr auto offset() const{return offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offset(dimensionality_type d) const{assert(d==0); (void)d; return offset_;}
	constexpr auto nelems() const{return nelems_;}
	constexpr auto nelems(dimensionality_type d) const{return d==0?nelems_:throw 0;}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr size_type size() const{//assert(stride_!=0 and nelems_%stride_==0);
		return nelems_/stride_;
	}
	friend constexpr auto size(layout_t const& self){return self.size();}
	constexpr size_type size(dimensionality_type d) const{
		return d==0?nelems_/stride_:throw 0; // assert(d == 0 and stride_ != 0 and nelems_%stride_ == 0);
	}
	constexpr auto base_size() const{return nelems_;}
	auto is_compact() const{return base_size() == num_elements();}
	friend auto is_compact(layout_t const& self){return self.is_compact();}
public:
	constexpr auto stride(dimensionality_type d = 0) const{assert(!d); (void)d; return stride_;}
	friend constexpr index stride(layout_t const& self){return self.stride();}
public:
	constexpr auto strides() const{return std::make_tuple(stride());}
	friend auto strides(layout_t const& self){return self.strides();}
	constexpr auto sizes() const{return std::make_tuple(size());}
	template<class T=void> constexpr auto sizes_as() const{return detail::to_array<T>(sizes());}
	constexpr auto offsets() const{return std::make_tuple(offset());}
	constexpr auto nelemss() const{return std::make_tuple(nelems_);}
	constexpr size_type num_elements() const{return this->size();}
	friend constexpr size_type num_elements(layout_t const& s){return s.num_elements();}
//	constexpr bool empty() const{return nelems_ == 0;}
//	friend constexpr bool empty(layout_t const& s){return s.empty();}
	       constexpr bool is_empty()        const    {return not nelems_;}
	friend constexpr bool is_empty(layout_t const& s){return s.is_empty();}
	       constexpr bool    empty()        const    {return is_empty();}

	constexpr index_extension extension() const{return {offset_/stride_, (offset_+nelems_)/stride_};}
	friend constexpr auto extension(layout_t const& self){return self.extension();}
	constexpr auto extension(dimensionality_type d) const{
		return d==0?index_extension{offset_/stride_, (offset_ + nelems_)/stride_}:throw 0;
	}
	constexpr extensions_type extensions() const{return extensions_type{extension()};}//std::make_tuple(extension());}
	friend constexpr extensions_type extensions(layout_t const& self){return self.extensions();}
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
	template<typename Size>
	layout_t& partition(Size const&){return *this;}
	layout_t& rotate(dimensionality_type = 1){return *this;}
	layout_t& unrotate(dimensionality_type = 1){return *this;}
	constexpr layout_t scale(size_type s) const{return {stride_*s, offset_*s, nelems_*s};}
};

inline typename layout_t<1>::extensions_type_ operator*(layout_t<0>::index_extension const& ie, layout_t<0>::extensions_type_ const&){
	return std::make_tuple(ie);
}

template<dimensionality_type D, typename SSize>
struct layout_t : multi::equality_comparable2<layout_t<D>, void>{
	using dimensionality_type = multi::dimensionality_type;
	using rank = std::integral_constant<dimensionality_type, D>;
	static constexpr dimensionality_type dimensionality(){return rank{};}
	friend constexpr dimensionality_type dimensionality(layout_t const& l){
		return l.dimensionality();
	}
	using sub_type = layout_t<D-1>;
	using size_type = multi::size_type;
	using index = multi::index;
	using difference_type = multi::difference_type;
	using index_extension = multi::index_extension;
	using index_range = multi::range<index>;
	using stride_type = index;
	using offset_type = index;
	using nelems_type = index;
 	sub_type    sub_ = {};
	stride_type stride_ = 1;//std::numeric_limits<stride_type>::max();
	offset_type offset_ = 0;
	nelems_type nelems_ = 0;
	struct extensions_type_ 
		: std::decay_t<decltype(tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename sub_type::extensions_type::base_>()))>
	{
		using base_ = std::decay_t<decltype(tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename sub_type::extensions_type::base_>()))>;
		using base_::base_;
		extensions_type_() = default;
		template<class Array, typename = decltype(std::get<D-1>(std::declval<Array>()))> 
		constexpr extensions_type_(Array const& t) : extensions_type_(t, std::make_index_sequence<static_cast<std::size_t>(D)>{}){}
		extensions_type_(index_extension const& ie, typename layout_t<D-1>::extensions_type_ const& other) : extensions_type_(std::tuple_cat(std::make_tuple(ie), other.base())){}
		base_ const& base() const{return *this;}
		friend decltype(auto) base(extensions_type_ const& s){return s.base();}
		friend typename layout_t<D + 1>::extensions_type_ operator*(index_extension const& ie, extensions_type_ const& self){
			return {std::tuple_cat(std::make_tuple(ie), self.base())};
		}
		template<class Archive, std::size_t... I>
		void serialize_impl(Archive& ar, std::index_sequence<I...>){
			using boost::serialization::make_nvp;
			(void)std::initializer_list<int>{(ar & make_nvp("extension", std::get<I>(*this)),0)...};
		}
		template<class Archive>
		void serialize(Archive& ar, unsigned){
			serialize_impl(ar, std::make_index_sequence<D>{});
		}
		static constexpr dimensionality_type dimensionality = D;
	private:
		template<class Array, std::size_t... I, typename = decltype(base_{std::get<I>(std::declval<Array const&>())...})> constexpr extensions_type_(Array const& t, std::index_sequence<I...>) : base_{std::get<I>(t)...}{}
	//	template<class T, std::size_t N, std::size_t... I> extensions_type_(std::array<T, N> const& t, std::index_sequence<I...>) : extensions_type_{std::get<I>(t)...}{}
	};
	using extensions_type = extensions_type_;
	using strides_type    = decltype(tuple_cat(std::make_tuple(std::declval<index>()), std::declval<typename sub_type::strides_type>()));
//	using extensions_type = typename detail::repeat<index_extension, D>::type;
//	using extensions_io_type = std::array<index_extension, D>;
	HD auto operator()(index i) const{return i*stride_ - offset_;}
	auto origin() const{return sub_.origin() - offset_;}
	constexpr
	layout_t(
		sub_type sub, stride_type stride, offset_type offset, nelems_type nelems
	) : sub_{sub}, stride_{stride}, offset_{offset}, nelems_{nelems}{}
	constexpr 
	layout_t(index_extension const& ie, layout_t<D-1> const& s) : 
		sub_{s},
		stride_{1},//ie.size()*sub_.num_elements()!=0?sub_.size()*sub_.stride():1}, // use .size for nvcc
		offset_{sub_.offset_ + ie.first()*sub_.stride()},
		nelems_{ie.size()*sub_.num_elements()}                             // use .size fort
	{}
	constexpr layout_t() = default;//: sub{}, stride_{1}, offset_{0}, nelems_{0}{} // needs stride != 0 for things to work well in partially formed state
//	constexpr 
	layout_t(extensions_type const& e) :// = {}) : 
		sub_{detail::tail(e)}, 
		stride_{sub_.size()*sub_.stride()},//std::get<0>(e).size()*sub_.num_elements()!=0?sub_.size()*sub_.stride():1}, 
		offset_{sub_.offset_ + std::get<0>(e).first()*sub_.stride()}, 
		nelems_{std::get<0>(e).size()*sub_.num_elements()} 
	{}//assert(0);}
#if defined(__INTEL_COMPILER) or (defined(__GNUC) && (__GNUC<6))
	constexpr 
	layout_t(std::array<index_extension, D> x) noexcept :
		layout_t{multi::detail::to_tuple<index_extension>(x)}
	{}
#endif
	template<class StdArray, typename = std::enable_if_t<std::is_same<StdArray, std::array<index_extension, static_cast<std::size_t>(D)>>{}> >
	constexpr 
	layout_t(StdArray const& e) : 
		sub_{detail::tail(e)}, 
		stride_{1},//std::get<0>(e).size()*sub_.num_elements()!=0?sub_.size()*sub_.stride():1}, 
		offset_{sub_.offset_ + std::get<0>(e).first()*sub_.stride()}, 
		nelems_{std::get<0>(e).size()*sub_.num_elements()}
	{}
	constexpr index nelems() const{return nelems_;}
	friend constexpr index nelems(layout_t const& self){return self.nelems();}
	auto nelems(dimensionality_type d) const{return d?sub_.nelems(d-1):nelems_;}
public:
	constexpr bool operator==(layout_t const& o) const{
		return sub_==o.sub_ and stride_==o.stride_ and offset_==o.offset_ and nelems_==o.nelems_;
	}
	constexpr bool operator!=(layout_t const& o) const{return not((*this)==o);}
public:
	constexpr size_type num_elements() const{return size()*sub_.num_elements();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}

	       constexpr bool is_empty()        const    {return not nelems_;}
	friend constexpr bool is_empty(layout_t const& s){return s.is_empty();}
	NODISCARD(".empty means .is_empty")
	       constexpr bool    empty()        const    {return is_empty();}
	constexpr size_type size() const{
		if(not nelems_) return 0;
		assert(stride_);
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& l){return l.size();}
	size_type size(dimensionality_type d) const{return d?sub_.size(d-1):size();}

	constexpr index stride() const{return stride_;}
	index stride(dimensionality_type d) const{return d?sub_.stride(d-1):stride();}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	constexpr strides_type strides() const{return tuple_cat(std::make_tuple(stride()), sub_.strides());}
	friend constexpr strides_type strides(layout_t const& self){return self.strides();}
	constexpr index offset() const{return offset_;}
	constexpr index offset(dimensionality_type d) const{return d?sub_.offset(d-1):offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offsets() const{return tuple_cat(std::make_tuple(offset()), sub_.offsets());}
	constexpr auto nelemss() const{return tuple_cat(std::make_tuple(nelems()), sub_.nelemss());}

	constexpr auto base_size() const{using std::max; return max(nelems_, sub_.base_size());}
	auto is_compact() const{return base_size() == num_elements();}
	friend auto is_compact(layout_t const& self){return self.is_compact();}
	decltype(auto) shape() const{return sizes();}
	friend decltype(auto) shape(layout_t const& self){return self.shape();}
	constexpr auto sizes() const{return tuple_cat(std::make_tuple(size()), sub_.sizes());}
	template<class T = void>
	constexpr auto sizes_as() const{return detail::to_array<T>(sizes());}
public:
	constexpr index_extension extension() const{
		if(not nelems_) return {};
		assert(stride_);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	friend auto extension(layout_t const& self){return self.extension();}
	constexpr index_extension extension_aux() const{
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	template<dimensionality_type DD = 0>
	constexpr index_extension extension(dimensionality_type d) const{
		return d?sub_.extension(d-1):extension();
	}
	constexpr extensions_type extensions() const{return tuple_cat(std::make_tuple(extension()), sub_.extensions().base());}
	friend constexpr extensions_type extensions(layout_t const& self){return self.extensions();}
	void extensions_aux(index_extension* it) const{
		*it = extension();
		++it;
		sub_.extensions_aux(it);
	}
	template<typename Size>
	layout_t& partition(Size const& s){
		using std::swap;
		stride_ *= s;
	//	offset
		nelems_ *= s;
		sub_.partition(s);
		return *this;
	}
	layout_t& transpose(){
		using std::swap;
		swap(stride_, sub_.stride_);
		swap(offset_, sub_.offset_);
		swap(nelems_, sub_.nelems_);
		return *this;
	}
	layout_t& rotate(){
		using std::swap;
		swap(stride_, sub_.stride_);
		swap(offset_, sub_.offset_);
		swap(nelems_, sub_.nelems_);
		sub_.rotate();
		return *this;
	}
	layout_t& unrotate(){
		sub_.unrotate();
		using std::swap;
		swap(stride_, sub_.stride_);
		swap(offset_, sub_.offset_);
		swap(nelems_, sub_.nelems_);
		return *this;
	}
	layout_t& rotate(dimensionality_type r){
		if(r >= 0) while(r){rotate(); --r;}
		else return rotate(D - r);
		return *this;
	}
	layout_t& unrotate(dimensionality_type r){
		if(r >= 0) while(r){unrotate(); --r;}
		else return unrotate(D-r);
		return *this;
	}
	layout_t scale(size_type s) const{
		return layout_t{sub_.scale(s), stride_*s, offset_*s, nelems_*s};
	}
};

inline typename layout_t<2>::extensions_type_ operator*(layout_t<1>::index_extension const& ie, layout_t<1>::extensions_type_ const& self){
	return layout_t<2>::extensions_type_(ie, self);
}

template<dimensionality_type D>
using extensions_type_ = typename layout_t<D>::extensions_type_;

template<class T, class Layout>
constexpr auto sizes_as(Layout const& self)
->decltype(self.template sizes_as<T>()){
	return self.template sizes_as<T>();}

}}

namespace std{
	template<> struct tuple_size<boost::multi::extensions_type_<0>> : std::integral_constant<boost::multi::dimensionality_type, 0>{};
	template<> struct tuple_size<boost::multi::extensions_type_<1>> : std::integral_constant<boost::multi::dimensionality_type, 1>{};
	template<> struct tuple_size<boost::multi::extensions_type_<2>> : std::integral_constant<boost::multi::dimensionality_type, 2>{};
	template<> struct tuple_size<boost::multi::extensions_type_<3>> : std::integral_constant<boost::multi::dimensionality_type, 3>{};
	template<> struct tuple_size<boost::multi::extensions_type_<4>> : std::integral_constant<boost::multi::dimensionality_type, 4>{};
}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_LAYOUT

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<cassert>
#include<iostream>
#include<vector>

#include "../../multi/utility.hpp"

using std::cout;
namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_layout){
//	auto t = std::make_tuple(1.,2.,3.);
//	auto u = multi::detail::reverse(t);
//	assert( std::get<0>(u) == 3. );
//	auto t = multi::detail::to_tuple<3, multi::index_extension>({1,2,3});
//	assert( std::get<1>(t) == 2 );
//	std::array<multi::index, 3> arr{1,2,3};
//	auto u = multi::detail::to_tuple<multi::index_extension>(arr);
//	assert( std::get<1>(u) == 2 );
 {
	multi::layout_t<0> L; assert( dimensionality(L)==0 and num_elements(L) == 1 ); sizes(L);
}{
	multi::iextensions<0> x{};
	multi::layout_t<0> L(x);
}{  multi::layout_t<1> L{}; assert( dimensionality(L)==1 and num_elements(L) == 0 and size(L) == 0 and size(extension(L))==0 and stride(L)!=0 and is_empty(L) );
}{
	multi::layout_t<2> L({2, 10}); 
	assert( dimensionality(L)==2 );
	assert( num_elements(L) == 20 );
	assert( size(L) == 2 ); 
	assert( size(extension(L))==2 );
	assert( stride(L)==10 );/*std::numeric_limits<std::ptrdiff_t>::max()*/ 
	assert( not is_empty(L) );
}{
	multi::layout_t<1> L(multi::iextensions<1>{20});
	assert( dimensionality(L)==1 and num_elements(L) == 20 and size(L) == 20 );
	assert( stride(L) == 1 );
}{
	multi::layout_t<1> L(multi::iextensions<1>{1});
	assert( dimensionality(L)==1 );
	assert( num_elements(L) == 1 );
	assert( size(L) == 1 );
	assert( stride(L) == 1 );
}{
	multi::layout_t<2> L({1, 10}); 
	assert( dimensionality(L)==2 and num_elements(L) == 10 and size(L) == 1); 
	assert( not is_empty(L) );
	assert( size(extension(L))==1 );
	assert( stride(L)== 10 );//std::numeric_limits<std::ptrdiff_t>::max() );
	using std::get;
	assert( get<0>(strides(L)) == 10);
	assert( get<1>(strides(L)) == 1 );
}{
	multi::layout_t<2> L({10, 1});
	assert( dimensionality(L)==2 );
	assert( num_elements(L) == 10 );
	assert( size(L) == 10 );
	using std::get;
	assert( get<0>(strides(L)) == 1 );
	assert( get<1>(strides(L)) == 1 );
}{  multi::layout_t<2> L{}; assert( dimensionality(L)==2 and num_elements(L) == 0 and size(L) == 0 and size(extension(L))==0 and stride(L)!=0 and is_empty(L) );
}{  multi::layout_t<3> L{}; assert( num_elements(L) == 0 );
}{	multi::layout_t<3> L({{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({{10}, {10}, {10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({10, 10, 10}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({multi::index_extension{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000 );
}{	
//	multi::layout_t<3> L(std::array<multi::index_extension, 3>{{ {0,10}, {0,10}, {0,10} }}); 
//	std::array<multi::index_extension, 3>{{ {0,10}, {0,10}, {0,10} }}; 
//	assert( size(L) == 1 );
//	std::cout << "ne " << num_elements(L) << std::endl;
//	assert( num_elements(L) == 1000 );
}
}

BOOST_AUTO_TEST_CASE(continued){
 {	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); assert( num_elements(L) == 1000);
}{	multi::layout_t<3> L(std::make_tuple(multi::iextension{0, 10}, multi::iextension{0, 10}, multi::iextension{0, 10})); assert(L.num_elements() == 1000);
}{	multi::layout_t<3> L(std::make_tuple(multi::iextension{10}, multi::iextension{10}, multi::iextension{10})); assert( num_elements(L) == 1000);
}{	multi::layout_t<3> L(std::make_tuple(10, 10, multi::iextension{10})); assert( num_elements(L) == 1000 );
}{
	char buffer[sizeof(multi::layout_t<2>)];
	new(buffer) multi::layout_t<2>;
	assert( size( reinterpret_cast<multi::layout_t<2>&>(buffer) ) == 0 );
}{
	multi::layout_t<1> L;
	assert( size(L) == 0 );
}{
	multi::layout_t<2> L;
	assert( size(L) == 0 );
}{
	multi::layout_t<3> L;
	assert( size(L) == 0 );
}{
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}}); 
	multi::layout_t<3> L2{extensions(L)};

	assert( not L.empty() );

	assert( stride(L) == L.stride() );
	assert( offset(L) == L.offset() );
	assert( nelems(L) == L.nelems() );

	assert( stride(L) == 20*30 );
	assert( offset(L) == 0 );
	assert( nelems(L) == 10*20*30 );
	
	assert( L.stride(0) == stride(L) );
	assert( L.offset(0) == offset(L) );
	assert( L.nelems(0) == nelems(L) );

	assert( L.stride(1) == 30 );
	assert( L.offset(1) == 0 );
	assert( L.nelems(1) == 20*30 );
	
	assert( L.stride(2) == 1 );
	assert( L.offset(2) == 0 );
	assert( L.nelems(2) == 30 );

	assert( L.num_elements() == num_elements(L) );
	assert( L.size() == size(L) );
	assert( L.extension() == extension(L) );

	assert( num_elements(L) == 10*20*30 );
	assert( size(L) == 10 );
	assert( extension(L).first() == 0 );
	assert( extension(L).last() == 10 );

	assert( L.size(1) == 20 );
	assert( L.extension(1).first() == 0 );
	assert( L.extension(1).last() == 20 );

	assert( L.size(2) == 30 );
	assert( L.extension(2).first() == 0 );
	assert( L.extension(2).last() == 30 );

	using std::get;
	assert( get<0>(strides(L)) == L.stride(0) );
	assert( get<1>(strides(L)) == L.stride(1) );
	assert( get<2>(strides(L)) == L.stride(2) );

	auto const& strides = L.strides();
	assert( get<0>(strides) == L.stride(0) );
}
{
	multi::layout_t<3> L( {{0, 10}, {0, 20}, {0, 30}} );
	assert( stride(L) == 20*30 );
//	static_assert( L.stride() == 20*30, "!");
//	static_assert( get<0>(L.strides()) == 20*30, "!");
//	static_assert( get<1>(L.strides()) == 	30, "!");
//	static_assert( get<2>(L.strides()) == 	 1, "!");
//	static_assert( size(L) == 10, "!");
//	cerr << size(L) <<'\n';
}
{
	multi::layout_t<1> L({{0, 10}});
	assert( extension(L).first() == 0 );
	assert( extension(L).last() == 10 );
}
{
	multi::layout_t<1> L({{8, 18}});
	assert( extension(L).first() == 8 );
	assert( extension(L).last() == 18 );
}
{
	multi::layout_t<2> L({{0, 10}, {0, 20}});
	assert( extension(L).first() == 0 and extension(L).last() == 10 );
}
{
	multi::layout_t<2> L( {{0, 10}, {11, 31}} );
	BOOST_TEST( size(L) == 10   );
	BOOST_TEST( stride(L) == 20 );
	BOOST_TEST( offset(L) == 0 );
//	auto Lr = L.rotate();
//	assert( extension(L).first() == 0 );
//	assert( extension(L).last() == 10 );
}
{
	multi::layout_t<2> L( {{8, 18}, {0, 20}} );
	assert( size(L) == 10 );
	assert( stride(L) == 20 );
//	std::cout << extension(L).first() << std::endl;
//	assert( extension(L).first() == 8 );
//	assert( extension(L).last() == 18 );
}
{
	multi::layout_t<3> L({{0, 3}, {0, 5}, {10, 17}}); 
	assert( stride(L) == 5*7 );
	BOOST_TEST( stride(L.sub_.sub_) == 1 );
//	BOOST_TEST( offset(L.sub_.sub_) == 10 );
//	assert( nelems(L) == 10*20*30 );
}
{
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}}); 
	assert( stride(L) == 20*30 );
	assert( offset(L) == 0 );
	assert( nelems(L) == 10*20*30 );
}
{
	multi::layout_t<3> L({{10, 20}, {10, 30}, {10, 40}}); 
	assert( stride(L) == 20*30 );
//	BOOST_TEST( offset(L) == 0 );
//	assert( nelems(L) == 10*20*30 );
}
{
	std::tuple<int, int, int> ttt = {1,2,3};
	auto arrr = boost::multi::detail::to_array(ttt);
	assert(arrr[1] == 2);
}

}
#endif
#endif

