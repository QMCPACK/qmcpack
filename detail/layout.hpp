#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
#ifndef MULTI_LAYOUT_HPP
#define MULTI_LAYOUT_HPP
// © Alfredo A. Correa 2018-2020

#include "types.hpp"

#include "../detail/operators.hpp"
#include "../config/NODISCARD.hpp"
#include "../config/ASSERT.hpp"

#include<type_traits> // make_signed_t

#include<limits>

namespace boost{
namespace multi{

namespace detail{

template<typename T, typename... As>
inline constexpr void construct_from_initializer_list(T* p, As&&... as){
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
constexpr auto tuple_tail_impl(Tuple&& t, std::index_sequence<Ns...>){
   return std::forward_as_tuple(std::forward<decltype(std::get<Ns + 1>(t))>(std::get<Ns + 1>(t))...);
}

template<class Tuple>
constexpr auto tuple_tail(Tuple&& t)
->decltype(tuple_tail_impl(t, std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{} - 1>())){//std::tuple<Ts...> t){
	return tuple_tail_impl(t, std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{} - 1>());}

}

}}

namespace boost{
namespace multi{

struct f_tag{};

template<dimensionality_type D, typename SSize=multi::size_type> struct layout_t;

template<dimensionality_type D> struct extensions_t;

template<> struct extensions_t<0> :
		std::tuple<>{
typedef std::tuple<> base_;
	static constexpr dimensionality_type dimensionality = 0;
	using nelems_type = index;
	using std::tuple<>::tuple;
	extensions_t() = default;
	constexpr base_ const& base() const{return *this;}
	friend constexpr decltype(auto) base(extensions_t const& s){return s.base();}
	constexpr operator nelems_type() const{return 1;}
	template<class Archive> void serialize(Archive&, unsigned){}
	constexpr size_type num_elements() const{return 1;}
	constexpr std::tuple<> from_linear(nelems_type n) const{assert(n < num_elements()); (void)n;
		return {};
	}
	friend constexpr std::tuple<> operator%(nelems_type n, extensions_t const& s){return s.from_linear(n);}
	friend constexpr extensions_t intersection(extensions_t const&, extensions_t const&){return {};}
};

template<> struct extensions_t<1> : 
		std::tuple<multi::index_extension>{
typedef std::tuple<multi::index_extension> base_;
	static constexpr dimensionality_type dimensionality = 1;
	using nelems_type = index;
	using index_extension = multi::index_extension;
	using std::tuple<index_extension>::tuple;
	// cppcheck-suppress noExplicitConstructor ; because one wants to interpret a single number as a range from 0
	constexpr extensions_t(index_extension const& ie) : base_{ie}{}
	extensions_t() = default;
	constexpr base_ const& base() const{return *this;}
	// cppcheck-suppress noExplicitConstructor ; I don't know why TODO
	constexpr extensions_t(base_ const& t) : std::tuple<index_extension>(t){}
	friend constexpr decltype(auto) base(extensions_t const& s){return s.base();}
	template<class Archive> void serialize(Archive& ar, unsigned){ar & multi::archive_traits<Archive>::make_nvp("extension", std::get<0>(*this));}
	constexpr size_type num_elements() const{return std::get<0>(*this).size();}
	constexpr std::tuple<multi::index> from_linear(nelems_type n) const{assert(n < num_elements());
		return std::tuple<multi::index>{n};
	}
	friend constexpr std::tuple<multi::index> operator%(nelems_type n, extensions_t const& s){return s.from_linear(n);}
	friend extensions_t intersection(extensions_t const& x1, extensions_t const& x2){
		return extensions_t(std::tuple<index_extension>(intersection(std::get<0>(x1), std::get<0>(x2))));
	}
};

template<dimensionality_type D>
struct extensions_t : 
		std::decay_t<decltype(std::tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename extensions_t<D-1>::base_>()))>{
typedef std::decay_t<decltype(std::tuple_cat(std::make_tuple(std::declval<index_extension>()), std::declval<typename extensions_t<D-1>::base_>()))> base_;
	using base_::base_;
	static constexpr dimensionality_type dimensionality = D;
	extensions_t() = default;
	using nelems_type = multi::index;
	template<class Array, typename = decltype(std::get<D-1>(std::declval<Array>()))> 
	constexpr extensions_t(Array const& t) : extensions_t(t, std::make_index_sequence<static_cast<std::size_t>(D)>{}){}
	constexpr extensions_t(index_extension const& ie, typename layout_t<D-1>::extensions_type const& other) : extensions_t(std::tuple_cat(std::make_tuple(ie), other.base())){}
	constexpr base_ const& base() const{return *this;}
	friend constexpr decltype(auto) base(extensions_t const& s){return s.base();}
	friend constexpr typename layout_t<D + 1>::extensions_type operator*(index_extension const& ie, extensions_t const& self){
		return {std::tuple_cat(std::make_tuple(ie), self.base())};
	}
	constexpr explicit operator bool() const{return not layout_t<D>{*this}.empty();}
	template<class Archive, std::size_t... I>
	void serialize_impl(Archive& ar, std::index_sequence<I...>){
	//	using boost::serialization::make_nvp;
	//	(void)std::initializer_list<int>{(ar & make_nvp("extension", std::get<I>(*this)),0)...};
		(void)std::initializer_list<int>{(ar & multi::archive_traits<Archive>::make_nvp("extension", std::get<I>(*this)),0)...};
	//	(void)std::initializer_list<int>{(ar & boost::serialization::nvp<std::remove_reference_t<decltype(std::get<I>(*this))> >{"extension", std::get<I>(*this)},0)...};
	}
	constexpr auto from_linear(nelems_type n) const{
		auto const sub_extensions = extensions_t<D-1>(detail::tuple_tail(this->base()));
		auto const sub_num_elements = sub_extensions.num_elements();
		return std::tuple_cat(std::make_tuple(n/sub_num_elements), sub_extensions.from_linear(n%sub_num_elements));
	}
	friend constexpr auto operator%(nelems_type n, extensions_t const& s){return s.from_linear(n);}
	template<class Archive>
	void serialize(Archive& ar, unsigned){
		serialize_impl(ar, std::make_index_sequence<D>{});
	}
private:
	template<class Array, std::size_t... I, typename = decltype(base_{std::get<I>(std::declval<Array const&>())...})> 
	constexpr extensions_t(Array const& t, std::index_sequence<I...>) : base_{std::get<I>(t)...}{}
//	template<class T, std::size_t N, std::size_t... I> extensions_type_(std::array<T, N> const& t, std::index_sequence<I...>) : extensions_type_{std::get<I>(t)...}{}
	static constexpr size_type multiply_fold(){return 1;}
	static constexpr size_type multiply_fold(size_type const& a0){return a0;}
	template<class...As> static constexpr size_type multiply_fold(size_type const& a0, As const&...as){return a0*multiply_fold(as...);}
	template<std::size_t... I> constexpr size_type num_elements_impl(std::index_sequence<I...>) const{return multiply_fold(std::get<I>(*this).size()...);}
public:
	constexpr size_type num_elements() const{return num_elements_impl(std::make_index_sequence<D>{});}
	friend constexpr extensions_t intersection(extensions_t const& x1, extensions_t const& x2){
		return extensions_t(
			std::tuple_cat(
				std::tuple<index_extension>(intersection(std::get<0>(x1), std::get<0>(x2))),
				intersection( extensions_t<D-1>(detail::tuple_tail(x1.base())), extensions_t<D-1>(detail::tuple_tail(x2.base())) ).base()
			)
		);
	}
};

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
	using extensions_type = extensions_t<0>;
	constexpr layout_t(extensions_type const& = {}){}// : nelems_{1}{}
	constexpr extensions_type extensions() const{return extensions_type{};}
	friend constexpr auto extensions(layout_t const& self){return self.extensions();}
	constexpr auto sizes() const{return std::tuple<>{};}
	friend constexpr auto sizes(layout_t const& s){return s.sizes();}
	constexpr nelems_type num_elements() const{return 1;}//nelems_;}
	constexpr bool operator==(layout_t const&) const{return true ;}
	constexpr bool operator!=(layout_t const&) const{return false;}
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
	using extensions_type = extensions_t<1>;
	using strides_type = std::tuple<index>;
	layout_t() = default;
	layout_t(layout_t const&) = default;
	constexpr layout_t(index_extension ie, layout_t<0> const&) :
		stride_{1},
		offset_{ie.first()},
		nelems_{
//			ie.size()<=1?ie.size()*std::numeric_limits<stride_type>::max():ie.size()
			ie.size()
		}{}
	constexpr explicit layout_t(extensions_type e) : layout_t(std::get<0>(e), {}){}
	constexpr layout_t(stride_type stride, offset_type offset, nelems_type nelems) : 
		stride_{stride}, offset_{offset}, nelems_{nelems}
	{}
	constexpr auto offset() const{return offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offset(dimensionality_type d) const{assert(d==0); (void)d; return offset_;}
	constexpr auto nelems() const{return nelems_;}
	constexpr auto nelems(dimensionality_type d) const{return d==0?nelems_:throw 0;}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr size_type size() const{//assert(stride_!=0 and nelems_%stride_==0);
		MULTI_ACCESS_ASSERT(stride_);
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& self){return self.size();}
	constexpr size_type size(dimensionality_type d) const{
		return d==0?nelems_/stride_:throw 0; // assert(d == 0 and stride_ != 0 and nelems_%stride_ == 0);
	}
	constexpr layout_t& reindex(index i){offset_ = i*stride_; return *this;}
	constexpr auto base_size() const{return nelems_;}
	       constexpr auto is_compact() const{return base_size() == num_elements();}
	friend constexpr auto is_compact(layout_t const& self){return self.is_compact();}
public:
	constexpr auto stride(dimensionality_type d = 0) const{assert(!d); (void)d; return stride_;}
	friend constexpr index stride(layout_t const& self){return self.stride();}
public:
	constexpr auto strides() const{return std::make_tuple(stride());}
	friend constexpr auto strides(layout_t const& self){return self.strides();}
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

	constexpr index_extension extension()        const&{
		assert(stride_);
		return {offset_/stride_, (offset_+nelems_)/stride_};
	} friend
	constexpr index_extension extension(layout_t const& s){return s.extension();}
	constexpr auto extension(dimensionality_type d) const{
		assert(stride_);
		return d==0?index_extension{offset_/stride_, (offset_ + nelems_)/stride_}:throw 0;
	}
	constexpr extensions_type extensions() const{return extensions_type{extension()};}//std::make_tuple(extension());}
	friend constexpr extensions_type extensions(layout_t const& self){return self.extensions();}
private:
	friend struct layout_t<2u>;
	void constexpr strides_aux(size_type* it) const{*it = stride();}
	void constexpr sizes_aux(size_type* it) const{*it = size();}
	void constexpr offsets_aux(index* it) const{*it = offset();}
	void constexpr extensions_aux(index_extension* it) const{*it = extension();}
public:
	constexpr auto operator()(index i) const{return i*stride_ - offset_;}
	constexpr std::ptrdiff_t at(index i) const{return offset_ + i*stride_;}
	constexpr std::ptrdiff_t operator[](index i) const{return at(i);}
	constexpr auto origin() const{return -offset_;}
	constexpr bool operator==(layout_t const& other) const{
		return stride_==other.stride_ and offset_==other.offset_ and nelems_==other.nelems_;
	}
	constexpr bool operator!=(layout_t const& other) const{return not(*this==other);}
	template<typename Size>
	constexpr layout_t& partition(Size const&){return *this;}
	constexpr layout_t&   rotate(dimensionality_type = 1){return *this;}
	constexpr layout_t& unrotate(dimensionality_type = 1){return *this;}
	constexpr layout_t scale(size_type s) const{return {stride_*s, offset_*s, nelems_*s};}
	constexpr layout_t& reverse(){return *this;}
};

inline constexpr typename layout_t<1>::extensions_type operator*(layout_t<0>::index_extension const& ie, layout_t<0>::extensions_type const&){
	return typename layout_t<1>::extensions_type{std::make_tuple(ie)};
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
	using extensions_type = extensions_t<D>;
	using strides_type    = decltype(tuple_cat(std::make_tuple(std::declval<index>()), std::declval<typename sub_type::strides_type>()));
//	using extensions_type = typename detail::repeat<index_extension, D>::type;
//	using extensions_io_type = std::array<index_extension, D>;
	constexpr auto operator()(index i) const{return i*stride_ - offset_;}
	constexpr auto origin() const{return sub_.origin() - offset_;}
	constexpr sub_type at(index i) const{//assert( this->extension().contains(i) ); see why it gives false positives
		auto ret = sub_;
		ret.offset_ += offset_ + i*stride_;
		return ret;
	}
	constexpr sub_type operator[](index i) const{return at(i);}
	constexpr layout_t(
		sub_type sub, stride_type stride, offset_type offset, nelems_type nelems
	) : sub_{sub}, stride_{stride}, offset_{offset}, nelems_{nelems}{}
	layout_t() = default;
	layout_t(layout_t const&) = default;
	constexpr explicit layout_t(extensions_type const& e) :
		sub_{detail::tail(e)}, 
		stride_{sub_.size()*sub_.stride()},//std::get<0>(e).size()*sub_.num_elements()!=0?sub_.size()*sub_.stride():1}, 
		offset_{std::get<0>(e).first()*stride_}, //sub_.offset_ + std::get<0>(e).first()*sub_.stride()}, //sub_.stride()  offset_ = i*stride_}, 
		nelems_{std::get<0>(e).size()*sub_.num_elements()} 
	{}
	template<class StdArray, typename = std::enable_if_t<std::is_same<StdArray, std::array<index_extension, static_cast<std::size_t>(D)>>{}> >
	constexpr 
	explicit layout_t(StdArray const& e) : 
		sub_{detail::tail(e)}, 
		stride_{1},//std::get<0>(e).size()*sub_.num_elements()!=0?sub_.size()*sub_.stride():1}, 
		offset_{sub_.offset_ + std::get<0>(e).first()*sub_.stride()}, 
		nelems_{std::get<0>(e).size()*sub_.num_elements()}
	{}
	constexpr index nelems() const{return nelems_;}
	friend constexpr index nelems(layout_t const& self){return self.nelems();}
	constexpr auto nelems(dimensionality_type d) const{return d?sub_.nelems(d-1):nelems_;}
public:
	constexpr bool operator==(layout_t const& o) const{
		return sub_==o.sub_ and stride_==o.stride_ and offset_==o.offset_ and nelems_==o.nelems_;
	}
	constexpr bool operator!=(layout_t const& o) const{return not((*this)==o);}
public:
	constexpr layout_t& reindex(index i){offset_ = i*stride_; return *this;}
	template<class... Idx>
	constexpr layout_t& reindex(index i, Idx... is){reindex(i).rotate().reindex(is...).unrotate(); return *this;}
	constexpr size_type num_elements() const{return size()*sub_.num_elements();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}

	       constexpr bool is_empty()        const    {return not nelems_;}
	friend constexpr bool is_empty(layout_t const& s){return s.is_empty();}
	NODISCARD(".empty means .is_empty")
	       constexpr bool    empty()        const    {return is_empty();}
	constexpr size_type size() const{
		if(not nelems_) return 0;
		MULTI_ACCESS_ASSERT(stride_);
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& l){return l.size();}
	constexpr size_type size(dimensionality_type d) const{return d?sub_.size(d-1):size();}

	constexpr index stride() const{return stride_;}
	constexpr index stride(dimensionality_type d) const{return d?sub_.stride(d-1):stride();}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	constexpr strides_type strides() const{return tuple_cat(std::make_tuple(stride()), sub_.strides());}
	friend constexpr strides_type strides(layout_t const& self){return self.strides();}
	constexpr index offset() const{return offset_;}
	constexpr index offset(dimensionality_type d) const{return d?sub_.offset(d-1):offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offsets() const{return tuple_cat(std::make_tuple(offset()), sub_.offsets());}
	constexpr auto nelemss() const{return tuple_cat(std::make_tuple(nelems()), sub_.nelemss());}

	constexpr auto base_size() const{using std::max; return max(nelems_, sub_.base_size());}
	constexpr auto is_compact() const{return base_size() == num_elements();}
	friend constexpr auto is_compact(layout_t const& self){return self.is_compact();}
	constexpr decltype(auto) shape() const{return sizes();}
	friend constexpr decltype(auto) shape(layout_t const& self){return self.shape();}
	constexpr auto sizes() const{return tuple_cat(std::make_tuple(size()), sub_.sizes());}
	template<class T = void>
	constexpr auto sizes_as() const{return detail::to_array<T>(sizes());}
public:
	constexpr index_extension extension()        const&{
		if(not nelems_) return {};
		assert(stride_); 
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	friend constexpr index_extension extension(layout_t const& s){return s.extension();}
	constexpr index_extension extension_aux() const{
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	template<dimensionality_type DD = 0>
	constexpr index_extension extension(dimensionality_type d) const{return d?sub_.extension(d-1):extension();}
	constexpr extensions_type extensions() const{return tuple_cat(std::make_tuple(extension()), sub_.extensions().base());}
	friend constexpr extensions_type extensions(layout_t const& self){return self.extensions();}
	constexpr void extensions_aux(index_extension* it) const{
		*it = extension();
		++it;
		sub_.extensions_aux(it);
	}
	template<typename Size>
	constexpr layout_t& partition(Size const& s){
		using std::swap;
		stride_ *= s;
	//	offset
		nelems_ *= s;
		sub_.partition(s);
		return *this;
	}
	constexpr layout_t& transpose(){
		using std::swap;
		swap(stride_, sub_.stride_);
		swap(offset_, sub_.offset_);
		swap(nelems_, sub_.nelems_);
		return *this;
	}
	constexpr layout_t& reverse(){
		unrotate();
		sub_.reverse();
		return *this;
	}
	constexpr layout_t&   rotate(){transpose(); sub_.  rotate(); return *this;}
	constexpr layout_t& unrotate(){sub_.unrotate(); transpose(); return *this;}
	constexpr layout_t& rotate(dimensionality_type r){
		if(r >= 0) while(r){rotate(); --r;}
		else return rotate(D - r);
		return *this;
	}
	constexpr layout_t& unrotate(dimensionality_type r){
		if(r >= 0) while(r){unrotate(); --r;}
		else return unrotate(D-r);
		return *this;
	}
	constexpr layout_t scale(size_type s) const{
		return layout_t{sub_.scale(s), stride_*s, offset_*s, nelems_*s};
	}
};

inline constexpr typename layout_t<2>::extensions_type operator*(layout_t<1>::index_extension const& ie, layout_t<1>::extensions_type const& self){
	return layout_t<2>::extensions_type(ie, self);
}

//template<dimensionality_type D>
//using extensions_type_ = typename layout_t<D>::extensions_type;

template<class T, class Layout>
constexpr auto sizes_as(Layout const& self)
->decltype(self.template sizes_as<T>()){
	return self.template sizes_as<T>();}

}}

namespace std{
	template<> struct tuple_size<boost::multi::extensions_t<0>> : std::integral_constant<boost::multi::dimensionality_type, 0>{};
	template<> struct tuple_size<boost::multi::extensions_t<1>> : std::integral_constant<boost::multi::dimensionality_type, 1>{};
	template<> struct tuple_size<boost::multi::extensions_t<2>> : std::integral_constant<boost::multi::dimensionality_type, 2>{};
	template<> struct tuple_size<boost::multi::extensions_t<3>> : std::integral_constant<boost::multi::dimensionality_type, 3>{};
	template<> struct tuple_size<boost::multi::extensions_t<4>> : std::integral_constant<boost::multi::dimensionality_type, 4>{};
}

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<cassert>

#include "../utility.hpp"
#include "../array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_layout_with_offset){
	{
		multi::layout_t<1> l1(multi::iextension(2, 5));
		BOOST_TEST( l1.extension().start()  == 2 );
		BOOST_TEST( l1.extension().finish() == 5 );
	}
	{
		boost::multi::layout_t<2, long>::extensions_type x = {multi::iextension(2, 5), multi::iextension(0, 5)};
		multi::layout_t<2> l2(x);
		BOOST_TEST( l2.extension().start()  == std::get<0>(x).start()  );
		BOOST_TEST( l2.extension().finish() == std::get<0>(x).finish() );
	}
	{
		multi::layout_t<2> l2({multi::iextension(0, 3), multi::iextension(2, 7)});
	//	BOOST_TEST( l2.sub_.extension().start() == 2 );
		BOOST_TEST( std::get<1>(l2.extensions()).start()  == 2 );
		BOOST_TEST( std::get<1>(l2.extensions()).finish() == 7 );
	}
}

template<class T> void whats(T&&) = delete;

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
//	multi::iextensions<2> x(2, 3);
//	whats( multi::iextensions<2>(2, 5).operator bool() ); // TODO add operator bool to iextensions<2>
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
	multi::layout_t<1> L({{0, 10}});
	assert( size(L) == 10 );
	assert( extension(L).start () ==  0 );
	assert( extension(L).finish() == 10 );

	L.reindex(1);
	assert( size(L) == 10 );
	assert( extension(L).start () ==  1 );
	assert( extension(L).finish() == 11 );
}{
	multi::layout_t<2> L;
	assert( size(L) == 0 );
}{
	multi::layout_t<2> L({{0, 10}, {0, 20}});
	assert( size(L) == 10 );
	assert( extension(L).start () ==  0 );
	assert( extension(L).finish() == 10 );
	
	L.reindex(1);
	assert( extension(L).start () ==  1 );
	assert( extension(L).finish() == 11 );

	L.rotate().reindex(3).unrotate();
	BOOST_TEST_REQUIRE( extension(L).start () ==  1 );
	BOOST_TEST_REQUIRE( extension(L).finish() == 11 );

	BOOST_TEST_REQUIRE( std::get<0>(extensions(L)).start () == 1 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(L)).start () == 3 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(L)).finish() == 23 );
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

BOOST_AUTO_TEST_CASE(layout_to_offset){
	multi::layout_t<3> L({10, 20, 30});
	multi::array<double, 3> A({10, 20, 30});
	assert( L[0][0][0] == &A[0][0][0] - data_elements(A) );
	assert( L[0][0][1] == &A[0][0][1] - data_elements(A) );
	assert( L[0][0][2] == &A[0][0][2] - data_elements(A) );
	assert( L[0][1][2] == &A[0][1][2] - data_elements(A) );
	assert( L[3][1][2] == &A[3][1][2] - data_elements(A) );
}

BOOST_AUTO_TEST_CASE(layout_to_offset_sub){
	multi::array<double, 3> A({10, 20, 30});
	auto&& s = A({2, 6}, {4, 8}, {10, 20});
	auto l = layout(s);
	assert( l[0][0][0] == &s[0][0][0] - base(s) );
	assert( l[0][0][1] == &s[0][0][1] - base(s) );
	assert( l[0][0][2] == &s[0][0][2] - base(s) );
	assert( l[0][1][2] == &s[0][1][2] - base(s) );
	assert( l[3][1][2] == &s[3][1][2] - base(s) );
}


#endif
#endif

