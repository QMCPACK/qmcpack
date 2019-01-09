#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ `#-DNDEBUG` -std=c++14 -Wall -Wextra -I$HOME/prj -D_TEST_MULTI_LAYOUT $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
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
#include<boost/operators.hpp>

namespace boost{
namespace multi{

namespace detail{

template<typename T, typename... As>
inline void construct_from_initializer_list(T* p, As&&... as){
	::new(static_cast<void*>(p)) T(std::forward<As>(as)...);
}

template<class Tuple>
constexpr auto head(Tuple&& t)
->decltype(std::get<0>(std::forward<Tuple>(t))){
	return std::get<0>(std::forward<Tuple>(t));}
template<typename Tuple, std::size_t... Ns>
constexpr auto tail_impl(std::index_sequence<Ns...> , Tuple&& t){
	return std::make_tuple(std::get<Ns+1u>(std::forward<Tuple>(t))...);
}
template<class Tuple>
constexpr auto tail(Tuple const& t)
->decltype(tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t)){
	return tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t);}

template<typename T, std::size_t N>
std::array<T, N-1> tail(std::array<T, N> const& a){
	std::array<T, N-1> ret;
	std::copy(a.begin() + 1, a.end(), ret.begin());
	return ret;
}

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

template<dimensionality_type D>
struct layout_t
	: boost::equality_comparable<layout_t<D>>
{
//	static constexpr dimensionality_type dimensionality = D;
	using dimensionality_type = multi::dimensionality_type;
	static constexpr dimensionality_type rank = D;
	static constexpr dimensionality_type dimensionality(){return rank;}
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality();}
	using sub_t = layout_t<D-1>;
	using index = multi::index;
	using difference_type = multi::difference_type;
	using index_extension = multi::index_extension;
	using index_range = multi::range<index>;
 	sub_t sub;
	index stride_;
	index offset_;
	index nelems_;
	using extensions_type = typename detail::repeat<index_extension, D>::type;
	using extensions_io_type = std::array<index_extension, D>;
	auto operator()(index i) const{return i*stride_ - offset_;}
	auto origin() const{return -offset_ + sub.origin();}
	constexpr layout_t(index_extension const& ie, layout_t<D-1> const& s) : 
		sub{s},
		stride_{size(ie)*num_elements(sub)!=0?size(sub)*stride(sub):1}, 
		offset_{0}, 
		nelems_{size(ie)*num_elements(sub)}
	{}
	constexpr layout_t(extensions_type const& e = {}) : 
		sub{detail::tail(e)}, 
		stride_{std::get<0>(e).size()*sub.num_elements()!=0?sub.size()*sub.stride():1}, 
		offset_{0}, 
		nelems_{std::get<0>(e).size()*sub.num_elements()} 
	{}
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
	constexpr auto nelems() const{return nelems_;}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	auto nelems(dimensionality_type d) const{return d?sub.nelems(d-1):nelems_;}
	constexpr bool operator==(layout_t const& o) const{
		return sub==o.sub and stride_==o.stride_ and offset_==o.offset_ and nelems_==o.nelems_;
	}
//	constexpr bool operator!=(layout_t const& o) const{return not(*this==o);}
	constexpr size_type num_elements() const{return size()*sub.num_elements();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return not nelems_;}
	friend bool empty(layout_t const& s){return s.empty();}
	constexpr size_type size() const{
	//	assert(nelems_ != 0);
		if(nelems_ == 0) return 0;
		assert(stride_ != 0 and nelems_%stride_ == 0 );
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& l){return l.size();}
	size_type size(dimensionality_type d) const{return d?sub.size(d-1):size();}

	constexpr index stride() const{return stride_;}
	index stride(dimensionality_type d) const{return d?sub.stride(d-1):stride();}
	friend constexpr index stride(layout_t const& self){return self.stride();}
	constexpr auto strides() const{return tuple_cat(std::make_tuple(stride()), sub.strides());}
	friend constexpr decltype(auto) strides(layout_t const& self){return self.strides();}


	void offsets_aux(index* it) const{
		*it = offset();
		sub.offsets_aux(++it);
	}
	constexpr index offset() const{return offset_;}
	constexpr index offset(dimensionality_type d) const{return d?sub.offset(d-1):offset_;}
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offsets() const{return tuple_cat(std::make_tuple(offset()), sub.offsets());}

	decltype(auto) shape() const{return sizes();}
	friend decltype(auto) shape(layout_t const& self){return self.shape();}

	constexpr auto sizes() const{return tuple_cat(std::make_tuple(size()), sub.sizes());}
	friend constexpr auto sizes(layout_t const& self){return self.sizes();}

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
	constexpr auto extensions() const{
		auto e = extension();
		return std::tuple_cat(std::make_tuple(e), sub.extensions());
	}
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

template<>
struct layout_t<dimensionality_type{0}>{
	static constexpr dimensionality_type dimensionality = 0;
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
};

template<>
struct layout_t<dimensionality_type{1}>{
	using sub_t = layout_t<dimensionality_type{0}>;
//	typedef boost::multi::dimensionality_type dimensionality_type;
	static constexpr dimensionality_type rank = 1;
	static constexpr dimensionality_type dimensionality = rank;
	friend constexpr auto dimensionality(layout_t const& l){return l.dimensionality;}
	using index_extension = multi::index_extension;
	using extensions_type = typename detail::repeat<index_extension, 1u>::type;
//	using extensions_type = std::array<index_extension, 1>;
	using index = multi::index;
	using index_range = multi::range<index>;
	using difference_type = multi::difference_type;
	index stride_;
	index offset_;
	index nelems_;
//	layout_t() = default;
	layout_t() : stride_{1}, offset_{0}, nelems_{0}{}
	constexpr layout_t(index_extension ie, layout_t<0> const&) : 
		stride_{1}, offset_{0}, nelems_{ie.size()}
	{}
	template<class Extensions, typename = decltype(std::get<0>(Extensions{}))>
	constexpr layout_t(Extensions e) : 
		stride_{1*1}, offset_{0}, nelems_{std::get<0>(e).size()*1}
	{}
	constexpr auto offset() const{return offset_;}	
	friend constexpr index offset(layout_t const& self){return self.offset();}
	constexpr auto offset(dimensionality_type d) const{
		assert(d==0); (void)d;
		return offset_;
	}
	constexpr auto nelems() const{return nelems_;}
	constexpr auto nelems(dimensionality_type d) const{
		assert(d==0); (void)d;
		return nelems_;
	}
	friend constexpr auto nelems(layout_t const& self){return self.nelems();}
	constexpr size_type size() const{
		if(nelems_ == 0) return 0;
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return nelems_/stride_;
	}
	constexpr size_type size(dimensionality_type d) const{
		(void)d;
		assert(d == 0 and stride_ != 0 and nelems_%stride_ == 0);
		return nelems_/stride_;
	}
	friend constexpr size_type size(layout_t const& self){return self.size();}
public:
	constexpr index stride(dimensionality_type d = 0) const{(void)d; assert(d == 0); return stride_;}
	friend constexpr index stride(layout_t const& self){return self.stride();}
private:
	struct strides_t{
		layout_t const& l_;
		strides_t(strides_t const&) = delete;
		constexpr decltype(auto) operator[](dimensionality_type d) const{
			assert(d == 0); (void)d;
			return l_.stride();
		}
	};
public:
	constexpr auto strides() const{return std::make_tuple(stride());}
	constexpr auto sizes() const{return std::make_tuple(size());}
	constexpr auto offsets() const{return std::make_tuple(offset());}

	constexpr size_type num_elements() const{return this->size();}
	friend size_type num_elements(layout_t const& s){return s.num_elements();}
	constexpr bool empty() const{return nelems_ == 0;}
	friend bool empty(layout_t const& s){return s.empty();}
	constexpr index_extension extension() const{
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	friend auto extension(layout_t const& self){return self.extension();}
	constexpr index_extension extension(dimensionality_type d) const{
		(void)d;
		assert(d == 0);
		assert(stride_ != 0 and nelems_%stride_ == 0);
		return {offset_/stride_, (offset_ + nelems_)/stride_};
	}
	constexpr extensions_type extensions() const{return std::make_tuple(extension());}
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
	decltype(auto) shape() const{return sizes();}
};

}}

#if _TEST_MULTI_LAYOUT

#include<cassert>
#include<iostream>
#include<vector>

using std::cout;
namespace multi = boost::multi;

int main(){

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
}
}
#endif
#endif

