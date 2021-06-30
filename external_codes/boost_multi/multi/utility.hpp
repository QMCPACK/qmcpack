#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2021

#ifndef MULTI_UTILITY_HPP
#define MULTI_UTILITY_HPP

#include "detail/layout.hpp"

#if(__cplusplus >= 201703L)
#include<iterator> // std::size in c++17
#endif

namespace boost{
namespace multi{

template<class Array, typename Reference = void, typename Element = void>
struct array_traits;

template<class Array, typename Reference, typename Element>
struct array_traits{
	using reference = typename Array::reference;
	using element   = typename Array::element;
	using element_ptr = typename Array::element_ptr;
	using decay_type = typename Array::decay_type;
	using default_allocator_type = typename Array::default_allocator_type;
};

/*template<class Element>
struct array_traits<Element, void, void>{
	using reference = Element&;
	using element   = Element;
};*/

template<class T, typename = typename T::rank>
       std::true_type  has_rank_aux(T const&){return {};}
inline std::false_type has_rank_aux(...     ){return {};}

template<class T> struct has_rank : decltype(has_rank_aux(std::declval<T>())){};

template<typename T> struct rank;

template<typename T, typename = std::enable_if_t<has_rank<T>{}> > 
constexpr typename T::rank rank_aux(T const&);

template<typename T, typename = std::enable_if_t<not has_rank<T>{}> > 
constexpr std::integral_constant<size_t, std::rank<T>{}> rank_aux(T const&);

template<typename T> struct rank : decltype(rank_aux(std::declval<T>())){};

#if not defined(__cpp_lib_nonmember_container_access) or __cpp_lib_nonmember_container_access < 201411
template<class Container>
constexpr auto size(Container const& con)
->std::make_signed_t<decltype(con.size())>{
	return con.size();}
#else
#endif

//template <class T>
//constexpr auto size(T const& t)
//->decltype(t.size()){
//	return t.size();}

template<class Pointer, std::enable_if_t<std::is_pointer<Pointer>{}, int> =0> constexpr std::ptrdiff_t stride(Pointer){return 1;}
template<class Pointer, std::enable_if_t<std::is_pointer<Pointer>{}, int> =0> constexpr Pointer base(Pointer d){return d;}

template<class T, class U>
auto reinterpret_pointer_cast(U* other)
->decltype(reinterpret_cast<T*>(other)){
	return reinterpret_cast<T*>(other);}

template <class T, std::size_t N>
constexpr auto size(const T(&)[N]) noexcept{return multi::size_type{N};}

template<class T, typename = typename T::get_allocator>
std::true_type has_get_allocator_aux(T const&);
inline std::false_type has_get_allocator_aux(...);

template <class T, std::size_t N>
constexpr std::allocator<std::decay_t<typename std::remove_all_extents<T[N]>::type>> 
get_allocator(T(&)[N]) noexcept{return {};}

template<class T>
constexpr 
auto get_allocator(T* const&)
->decltype(std::allocator<typename std::iterator_traits<T*>::value_type>{}){
	return std::allocator<typename std::iterator_traits<T*>::value_type>{};}

template<class T>
constexpr 
auto default_allocator_of(T*)
->decltype(std::allocator<typename std::iterator_traits<T*>::value_type>{}){
	return std::allocator<typename std::iterator_traits<T*>::value_type>{};}

template<class T>
constexpr 
auto to_address(T* const& t)
->decltype(t){
	return t;}

template<class Archive> struct archive_traits{ // TODO implemente a poors man nvp that works with boost serialization, is it possible?
	template<class T>
	static constexpr T& make_nvp(char const* /*name*/, T& v){return v;}
};

//template<class It>
//constexpr auto get_allocator(It const& it)
//->decltype(get_allocator(to_address(it))){
//	return get_allocator(to_address(it));}

template<class T>
auto has_get_allocator_aux(T const& t)->decltype(t.get_allocator(), std::true_type {});
//inline auto has_get_allocator_aux(...)->decltype(                   std::false_type{});
template<class T> struct has_get_allocator : decltype(has_get_allocator_aux(std::declval<T>())){};

//template<class T, typename = std::enable_if_t<has_get_allocator<T>>
//decltype(auto) get_allocator(T const& v){return v.get_allocator();}

//template<class It, typename = std::enable_if_t<not has_get_allocator<It>{}>>
//auto get_allocator(It)->std::allocator<typename std::iterator_traits<It>::value_type>{return {};}

template<class T1, class T2, typename Ret = T1>// std::common_type_t<T1, T2>> 
Ret common(T1 const& t1, T2 const& t2){
	return t1==t2?t1:Ret{};}

//template<class It>//, class ss = decltype(get_allocator(typename std::iterator_traits<Iterator>::pointer{}))> 
//auto get_allocator(It const&, void* = 0)
//->decltype(std::allocator<typename std::iterator_traits<It>::value_type>{}){
//	return std::allocator<typename std::iterator_traits<It>::value_type>{};}

//->decltype(get_allocator(typename std::iterator_traits<Iterator>::pointer{})){
//	return get_allocator(typename std::iterator_traits<Iterator>::pointer{});}

//->decltype(std::allocator<typename std::iterator_traits<Iterator>::value_type>{}){
//	return std::allocator<typename std::iterator_traits<Iterator>::value_type>{};}

template<class T>
auto has_num_elements_aux(T const& t)->decltype(t.num_elements(), std::true_type {});
inline auto has_num_elements_aux(...       )->decltype(                  std::false_type{});
template<class T> struct has_num_elements : decltype(has_num_elements_aux(std::declval<T>())){};

template<class A, typename = std::enable_if_t<has_num_elements<A>{}> > 
constexpr auto num_elements(A const& arr)
->std::make_signed_t<decltype(arr.num_elements())>{
	return arr.num_elements();}

template<class T>
auto has_size_aux(T const& t)->decltype(t.size(), std::true_type {});
inline auto has_size_aux(...       )->decltype(                  std::false_type{});
template<class T> struct has_size : decltype(has_size_aux(std::declval<T>())){};

template<class T>
auto has_data_elements_aux(T&& t)->decltype(t.data_elements(), std::true_type {});
auto has_data_elements_aux(...  )->decltype(                   std::false_type{});
template<class T> struct has_data_elements : decltype(has_data_elements_aux(std::declval<T>())){};

template<class T>
auto has_data_aux(T&& t)->decltype(t.data_elements(), std::true_type {});
auto has_data_aux(...  )->decltype(                   std::false_type{});
template<class T> struct has_data : decltype(has_data_aux(std::declval<T>())){};

//template<class T, typename = std::enable_if_t<not has_num_elements<T>{} and not (has_size<T>{} and has_data<T>{})>> 
//constexpr size_type num_elements(T const&){return 1;}

template<class A, typename = std::enable_if_t<not has_num_elements<A>{} and has_size<A>{} and has_data<A>{}> >
constexpr auto num_elements(A const& arr)
->std::make_signed_t<decltype(arr.size())>{
	return arr.size();}

template<class A, typename = std::enable_if_t<has_data_elements<A>{}> > 
constexpr auto data_elements(A const& arr)
->decltype(arr.data_elements()){
	return arr.data_elements();}

//template<class A, typename = std::enable_if_t<(not has_data_elements<A>{}) and (has_data<A>{} and has_size<A>{})>>
//constexpr auto data_elements(A&& arr)
//->decltype(arr.data()){
//{	return arr.data();}

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
[[deprecated("use constexpr data_elements() or base() to extract pointer")]] 
auto data(T& t){return &t;}

template<class T, typename = std::enable_if_t<not std::is_array<T>{} and not has_data_elements<T>{} and not has_data<T>{}>>
constexpr auto data_elements(T& t){return &t;}

//template<class T>
//constexpr size_type num_elements(T&&) noexcept{return 1;} // this should be before the rest of `num_elements` functions?

template<class A> struct num_elements_t: std::integral_constant<size_type, 1>{};

template<class T, std::size_t N> struct num_elements_t<T[N]>: std::integral_constant<size_type, (N*num_elements_t<T>{})>{};

template<class T, std::size_t N> struct num_elements_t<T(&)[N]>: num_elements_t<T[N]>{};

template<class T, std::size_t N>
constexpr auto num_elements(const T(&/*t*/)[N]) noexcept{return num_elements_t<T[N]>{};}


template<class Vector>
constexpr auto num_elements(Vector const& v, std::enable_if_t<std::is_same<typename Vector::pointer, decltype(std::declval<Vector>().data())>{}, int> =0)
->std::make_signed_t<decltype(v.size())>{
	return v.size();}

template<class Vector, typename = std::enable_if_t<std::is_same<typename Vector::pointer, decltype(std::declval<Vector>().data())>{}> >
auto data_elements(Vector const& v)
->decltype(v.data()){
	return v.data();}

template <class T, std::size_t N>
constexpr std::ptrdiff_t stride(const T(&/*t*/)[N]) noexcept{return num_elements_t<T>{};}

template <class T, std::size_t N>
constexpr bool is_compact(const T(&)[N]) noexcept{return true;}

//template <class T, std::size_t N>
//constexpr ptrdiff_t stride(T(&t)[N]) noexcept{return num_elements(t[0]);}

template<class T, std::size_t N>
constexpr std::ptrdiff_t offset(const T(&)[N]) noexcept{return 0;}

template<class T, std::size_t N>
[[deprecated("use data_elements instead")]] // this name is bad because when the element belongs to std:: then std::data is picked up by ADL and the 
constexpr auto data(T(&t)[N]) noexcept{return data(t[0]);}

template<class T, std::size_t N>
constexpr auto data_elements(T(&t)[N]) noexcept{return data_elements(t[0]);}

//template<class T, std::size_t N>
//constexpr auto data(const T(&t)[N]) noexcept{return data(t[0]);}

template<class T>
auto has_dimensionality_aux(T const& t)->decltype(t.dimensionality(), std::true_type {});
inline auto has_dimensionality_aux(...)->decltype(                    std::false_type{});
template<class T> struct has_dimensionality : decltype(has_dimensionality_aux(std::declval<T>())){};

template<class Container, typename = std::enable_if_t<has_dimensionality<Container>{}> >
constexpr auto dimensionality(Container const& con)
->decltype(con.dimensionality()){
	return con.dimensionality();}

template<class T>
auto has_dimensionaliy_member_aux(T const& t)->decltype((size_t(t.dimensionality), std::true_type{}));
inline auto has_dimensionaliy_member_aux(...       )->decltype(                          std::false_type{});
template<class T> struct has_dimensionality_member : decltype(has_dimensionaliy_member_aux(std::declval<T>())){};

template<class C> constexpr auto dimensionality(C const& c)->decltype(c.dimensionality){return c.dimensionality;}

template<class T, typename = std::enable_if_t<not has_dimensionality_member<T>{}>>
constexpr auto dimensionality(T const&, void* = nullptr){return 0;}

template<class T, std::size_t N>
constexpr auto dimensionality(T const(&t)[N]){return 1 + dimensionality(t[0]);}

template<class T, typename = decltype(std::declval<T>().sizes())>
std::true_type has_sizes_aux(T const&);
inline std::false_type has_sizes_aux(...);

template<class T> struct has_sizes : decltype(has_sizes_aux(std::declval<T>())){};

template<class Array, typename = std::enable_if_t<has_sizes<Array>{}> >
constexpr auto sizes(Array const& arr)
->decltype(arr.sizes()){
	return arr.sizes();}

template<class T, typename = std::enable_if_t<not has_sizes<T>{}> >
inline constexpr std::tuple<> sizes(T const&){return {};}

inline decltype(auto) base(std::tuple<> const& a){return a;}

template<class T, std::size_t N>
constexpr auto sizes(const T(&t)[N]) noexcept{
//	using std::size; // this line needs c++17
	using boost::multi::size;
	return std::tuple_cat(std::make_tuple(boost::multi::size(t)), sizes(t[0]));
}

//template<class T>
//constexpr auto origin(T& t){return &t;}
template<class T, std::size_t N>
constexpr auto origin(T(&t)[N]) noexcept{return reinterpret_cast<std::remove_all_extents_t<T[N]>*>(&t);}

//template<class T, std::size_t N>
//constexpr auto base(const T(&t)[N]) noexcept{
//	return reinterpret_cast<std::remove_all_extents_t<T[N]> const*>(&t);
//}

template<class T, std::size_t N>
constexpr auto base(T(&t)[N]) noexcept{
	return data_elements(t);
//	return reinterpret_cast<std::remove_all_extents_t<T[N]>*>(&t);
}


template<class T, std::size_t N>
constexpr auto base(T(*&t)[N]) noexcept{return base(*t);}
//reinterpret_cast<std::remove_all_extents_t<T[N]>*>(&t);}

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
constexpr auto base(T const* t) noexcept{return t;}

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
constexpr auto base(T* t) noexcept{return t;}

//template<class T, typename = std::enable_if_t<not std::is_array<T>{}> > 
//constexpr std::ptrdiff_t stride(T*&/*t*/) noexcept{
//	return 1;
//}
//template<class T, typename = std::enable_if_t<not std::is_array<T>{}> > 
//constexpr std::ptrdiff_t stride(T const*&/*t*/) noexcept{
//	return 1;
//}

//inline auto base(double& d){return &d;}
//inline auto base(float& f){return &f;}
//inline auto base(std::complex<double>& c){return &c;}
//inline auto base(std::complex<float>& z){return &z;}

template<class T, std::enable_if_t<std::is_standard_layout<T>{} and std::is_trivial<T>{}, int> = 0>
auto base(T& t){return &t;}

//template<class T, std::enable_if_t<std::is_pod<std::decay_t<T>>{}, int> = 0>
//auto stride(T& t) = delete;

//template<class T> constexpr std::ptrdiff_t stride(T const*/*t*/) noexcept{return 1;}

//template<typename T, std::size_t N>
//constexpr std::ptrdiff_t stride(T(*&/*t*/)[N]) noexcept{return N;}

template<class T>
constexpr auto corigin(const T& t){return &t;}
template<class T, std::size_t N>
constexpr auto corigin(const T(&t)[N]) noexcept{return corigin(t[0]);}

template<class T, typename = decltype(std::declval<T>().extension())>
       std::true_type  has_extension_aux(T const&);
inline std::false_type has_extension_aux(...     );
template<class T> struct has_extension : decltype(has_extension_aux(std::declval<T>())){};

template<class Container, class=std::enable_if_t<not has_extension<Container>{}>>
auto extension(Container const& c) // TODO consider "extent"
->decltype(multi::extension_t<std::make_signed_t<decltype(size(c))>>(0, size(c))){
	return multi::extension_t<std::make_signed_t<decltype(size(c))>>(0, size(c));}

template<class T, typename = decltype(std::declval<T>().extensions())>
       std::true_type  has_extensions_aux(T const&);
inline std::false_type has_extensions_aux(...     );

template<class T> struct has_extensions : decltype(has_extensions_aux(std::declval<T>())){};

template<class T, typename = std::enable_if_t<has_extensions<T>{}> >
auto extensions(T const& t)
->decltype(t.extensions()){
	return t.extensions();}

template<class T, typename = std::enable_if_t<not has_extensions<T>{}> >
constexpr /*std::tuple<>*/ multi::layout_t<0>::extensions_type extensions(T const&){return {};}

template<class T, size_t N>
constexpr auto extensions(T(&t)[N]){
//	return tuple_cat(std::make_tuple(index_extension(N)), extensions(t[0]));
	return index_extension(N)*extensions(t[0]);
}

template<dimensionality_type D>
struct extensions_aux{
	template<class T>
	static auto call(T const& t){
		return std::tuple_cat(std::make_tuple(t.extension()), extensions<D-1>(t));
	}
};

template<> struct extensions_aux<0>{
	template<class T> static auto call(T const& ){return std::make_tuple();}
};

template<dimensionality_type D, class T>
auto extensions(T const& t){
	return extensions_aux<D>::call(t);
//	if constexpr(D != 0)
//		return std::tuple_cat(std::make_tuple(t.extension()), extensions<D-1>(t));
//	else
//		return std::make_tuple();
}

template<class T1> struct extensions_t_aux;

template<class T1, class T2> auto extensions_(T2 const& t2){
	return extensions_t_aux<T1>::call(t2);
}
/*template<class T1, class T2> auto extensions(T2 const& t2){
	if constexpr(std::is_same<T1, T2>{})
		return std::make_tuple();
	else
		return std::tuple_cat(std::make_tuple(t2.extension()), extensions<T1>(*begin(t2)));
}*/

template<class T1> struct extension_t_aux{
	static auto call(T1 const&){return std::make_tuple();}
	template<class T2>
	static auto call(T2 const& t2){return std::tuple_cat(std::make_tuple(t2.extension()), extensions_<T1>(*begin(t2)));}
};

template<class T, typename = decltype(std::declval<T const&>().layout())>
std::true_type has_layout_member_aux(T const&);
inline std::false_type has_layout_member_aux(...);

template<class T>
struct has_layout_member : decltype(has_layout_member_aux(std::declval<T const&>())){};

template<class T, typename = std::enable_if_t<has_layout_member<T const&>{}> >
auto layout(T const& t)
->decltype(t.layout()){
	return t.layout();}

template<class T, typename = std::enable_if_t<not has_layout_member<T const&>{}> >
layout_t<0> layout(T const&){return {};}

template<class T, std::size_t N>
constexpr auto layout(T(&t)[N]){
//	return multi::layout_t<multi::dimensionality(t)>(multi::extensions(t));
	return multi::layout_t<std::rank<T[N]>{}>(multi::extensions(t));
}

template<class T, std::size_t N>
constexpr auto strides(T(&t)[N]){return layout(t).strides();}

}}


// std::array specializations
namespace boost{
namespace multi{

template<class T, std::size_t N>
struct array_traits<std::array<T, N>>{
	static constexpr dimensionality_type dimensionality = 1;
	using reference = T&;
	using value_type = std::decay_t<T>;
	using pointer = T*;
	using element = value_type;
	using element_ptr = pointer;
	using decay_type = std::array<value_type, N>;
};

template<class T, std::size_t N, std::size_t M>
struct array_traits<std::array<std::array<T, M>, N>>{
	static constexpr dimensionality_type dimensionality = 1 + array_traits<std::array<T, M>>::dimensionality;
	using reference = std::array<T, M>&;
	using value_type = std::array<std::decay_t<T>, M>;
	using pointer = std::array<T, M>*;
	using element = typename array_traits<std::array<T, M>>::element;
	using element_ptr = typename array_traits<std::array<T, M>>::element;
	using decay_type = std::array<value_type, M>;
};

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N>&       arr) noexcept{return arr.data();}
template<class T, std::size_t M, std::size_t N> constexpr auto data_elements(std::array<std::array<T, M>, N>& arr) noexcept{return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N> const&       arr) noexcept{return arr.data();}
template<class T, std::size_t M, std::size_t N> constexpr auto data_elements(std::array<std::array<T, M>, N> const& arr) noexcept{return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N>&&       arr) noexcept{return arr.data();}
template<class T, std::size_t M, std::size_t N> constexpr auto data_elements(std::array<std::array<T, M>, N>&& arr) noexcept{return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr std::ptrdiff_t num_elements(std::array<T, N> const& /*unused*/) noexcept{return N;}
template<class T, std::size_t M, std::size_t N> 
constexpr std::ptrdiff_t num_elements(std::array<std::array<T, M>, N> const& arr){return N*num_elements(arr[0]);}

template<class T, std::size_t N> 
constexpr dimensionality_type dimensionality(std::array<T, N> const& /*unused*/){return 1;}

template<class T, std::size_t M, std::size_t N> 
constexpr dimensionality_type dimensionality(std::array<std::array<T, M>, N> const& arr){return 1 + dimensionality(arr[0]);}

#if (__cplusplus < 201703L)
// this conflicts with std::size in nvcc 11 and c++17
template<class T, std::size_t N>
constexpr auto size(std::array<T, N> const& /*arr*/){
	return multi::size_type{N};
}
#endif

template<class T, std::size_t N>
constexpr auto extensions(std::array<T, N> const& /*arr*/){
	return multi::extensions_t<1>{{0, N}};
}

template<class T, std::size_t N, std::size_t M>
auto extensions(std::array<std::array<T, N>, M> const& arr){
	return multi::iextension(M)*extensions(arr[0]);
}

template<class T, std::size_t N>
constexpr auto stride(std::array<T, N> const& /*arr*/){
	return multi::size_type{1}; // multi::stride_type?
}

template<class T, std::size_t N, std::size_t M>
constexpr auto stride(std::array<std::array<T, N>, M> const& arr){
	return num_elements(arr[0]);
}

template<class T, std::size_t N>
constexpr auto layout(std::array<T, N> const& arr){
	return multi::layout_t<multi::array_traits<std::array<T, N>>::dimensionality>{multi::extensions(arr)};
}

}}


namespace boost{
namespace serialization{
//	template<class Archive, template<class,  std::size_t> class ArrayRef, class E,  std::size_t D>
//	inline void serialize(Archive&, ArrayRef<E, D>&, const unsigned int){
//		assert(0);
//	}
}}

#if not __INCLUDE_LEVEL__ // TEST BELOW

#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi zero dimensionality"
#include<boost/test/unit_test.hpp>

#include<cassert>
#include<iostream>
#include<vector>
//#include<cmath>

using std::cout;
namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_utility_test){

	static_assert( std::is_same<std::iterator_traits<double const*>::value_type, double>{}, "!");

	using multi::corigin;
	using multi::dimensionality;
	using multi::extension;
	using multi::extensions;
	using multi::origin;
	using multi::size;
	using multi::sizes;
	using multi::num_elements;
{
	double A[4] = {1.,2.,3.,4.};
	assert( dimensionality(A) == 1 );
	assert( extension(A).first() == 0 );
	assert( extension(A).last() == 4 );
	assert( origin(A) == &A[0] );
	assert( size(A) == 4 );
	assert( std::get<0>(sizes(A)) == size(A) );
	using multi::get_allocator;

	static_assert(std::is_same<decltype(get_allocator(A)), std::allocator<double> >{}, "!");

	using std::addressof;
//	using multi::data;
	using multi::data_elements;
	static_assert( std::is_same<decltype(data_elements(A)), double*>{} , "!");
//	assert( data(A) == addressof(A[0]) );
	assert( data_elements(A) == addressof(A[0]) );
}{
	double const A[4] = {1.,2.,3.,4.};
	f(A);
}{
	double A[2][3] = {{1.,2.,3.},{4.,5.,6.}};
	assert( dimensionality(A) == 2 );
	assert( extension(A).first() == 0 );
	assert( extension(A).last() == 2 );
//	int a = extensions(A);
	assert( origin(A) == &A[0][0] );
	*origin(A) = 99.;
	assert( A[0][0] == 99. );	
	assert( corigin(A) == &A[0][0] );
	assert( size(A) == 2 );
	assert( std::get<0>(sizes(A)) == size(A) );
	BOOST_REQUIRE( num_elements(A) == 6 );
	static_assert( num_elements(A) == 6 , "!" );

}{
	double const A[2][3] = {{1.,2.,3.},{4.,5.,6.}};
	assert( origin(A) == &A[0][0] );
	assert( origin(A) == &A[0][0] );
//	*origin(A) = 99.;
}{
//	static_assert( multi::rank<std::array<double, 10>>{} == 1, "!" );
//	static_assert( multi::rank<std::array<std::array<double, 2>, 10>>{} == 2, "!" );
//	std::array<std::array<double, 2>, 10> a;
//	auto x = multi::extensions(a);
//	assert( std::get<0>(x) == 10 );
//	assert( std::get<1>(x) == 2 );
//	std::cout << multi::num_elements(a) << std::endl;
//	assert( multi::num_elements(a) == 20 );
}{
	
}

}

#endif
#endif

