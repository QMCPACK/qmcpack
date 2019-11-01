#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -std=c++17 -Wall `#-Wextra` -I$HOME/prj -D_TEST_MULTI_UTILITY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_UTILITY_HPP
#define MULTI_UTILITY_HPP

//#include "detail/index_range.hpp"
#include "detail/layout.hpp"

#include<tuple> 

namespace boost{
namespace multi{

//template<class Array, typename Reference = void> struct array_traits;

template<class Array>///, typename std::enable_if_t< Array::reference>
struct array_traits{
	using reference = typename Array::reference;
};

template<class T, size_t N>
struct array_traits<T[N]>{
	using reference = T&;
};


template<class T, typename = typename T::rank>
std::true_type has_rank_aux(T const&){return {};}
inline std::false_type has_rank_aux(...){return {};}

template<class T> struct has_rank : decltype(has_rank_aux(std::declval<T>())){};

template<typename T> struct rank;

template<typename T, size_t N> 
constexpr std::integral_constant<size_t, 1 + multi::rank<T>{}> rank_aux(std::array<T, N> const&);

template<typename T, typename = std::enable_if_t<has_rank<T>{}> > 
constexpr typename T::rank rank_aux(T const&);

template<typename T, typename = std::enable_if_t<not has_rank<T>{}> > 
constexpr std::integral_constant<size_t, std::rank<T>{}> rank_aux(T const&);

template<typename T> struct rank : decltype(rank_aux(std::declval<T>())){};

#if __cpp_lib_nonmember_container_access < 201411
template<class Container>
constexpr auto size(Container const& con)
->decltype(con.size()){
	return con.size();}
#else
#endif

template <class T, std::size_t N>
constexpr auto size(const T(&)[N]) noexcept{return multi::size_type{N};}

template<class T>
auto has_num_elements_aux(T const& t)->decltype(t.num_elements(), std::true_type {});
inline auto has_num_elements_aux(...       )->decltype(                  std::false_type{});
template<class T> struct has_num_elements : decltype(has_num_elements_aux(std::declval<T>())){};

template<class A, typename = std::enable_if_t<has_num_elements<A>{}> > 
constexpr auto num_elements(A const& arr)
->decltype(arr.num_elements()){
	return arr.num_elements();}
	
template<class T, typename = std::enable_if_t<!has_num_elements<T>{}>> 
constexpr size_type num_elements(T const&){return 1;}

template<class T>
auto has_data_elements_aux(T&& t)->decltype(t.data_elements(), std::true_type {});
auto has_data_elements_aux(...  )->decltype(                   std::false_type{});
template<class T> struct has_data_elements : decltype(has_data_elements_aux(std::declval<T>())){};

template<class A, typename = std::enable_if_t<has_data_elements<A>{}> > 
constexpr auto data_elements(A const& arr)
->decltype(arr.data_elements()){
	return arr.data_elements();}

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
[[deprecated("use constexpr data_elements")]] auto data(T& t){return &t;}

template<class T, typename = std::enable_if_t<not std::is_array<T>{} and not has_data_elements<T>{}>>
constexpr auto data_elements(T& t){return &t;}

template<class T, std::size_t N>
constexpr auto num_elements(const T(&t)[N]) noexcept{return N*num_elements(t[0]);}

template<class T, size_t N>
constexpr auto num_elements(std::array<T, N> arr){return N*num_elements(arr[0]);}

template <class T, std::size_t N>
constexpr auto stride(const T(&t)[N]) noexcept{return num_elements(t[0]);}

template<class T, std::size_t N>
constexpr std::ptrdiff_t offset(const T(&)[N]) noexcept{return 0;}

template<class T, std::size_t N>
[[deprecated("use data_elements instead")]] // this name is bad because when the element belongs to std:: then std::data is picked up by ADL and the 
constexpr auto data(T(&t)[N]) noexcept{return data(t[0]);}

template<class T, std::size_t N>
constexpr auto data_elements(T(&t)[N]) noexcept{return data_elements(t[0]);}

//template<class T, std::size_t N>
//constexpr auto data(const T(&t)[N]) noexcept{return data(t[0]);}

template<class Container>
auto extension(Container const& c) // TODO consider "extent"
->decltype(range<decltype(size(c))>{0, size(c)}){
	return range<decltype(size(c))>{0, size(c)};}

template<class T>
auto has_dimensionaliy_aux(T const& t)->decltype(t.dimensionality(), std::true_type {});
inline auto has_dimensionaliy_aux(...       )->decltype(                    std::false_type{});
template<class T> struct has_dimensionality : decltype(has_dimensionaliy_aux(std::declval<T>())){};

template<class Container, typename = std::enable_if_t<has_dimensionality<Container>{}> >
constexpr auto dimensionality(Container const& con)
->decltype(con.dimensionality()){
	return con.dimensionality();}

template<class T>
auto has_dimensionaliy_member_aux(T const& t)->decltype(size_t(t.dimensionality), std::true_type {});
inline auto has_dimensionaliy_member_aux(...       )->decltype(                          std::false_type{});
template<class T> struct has_dimensionality_member : decltype(has_dimensionaliy_member_aux(std::declval<T>())){};

template<class C> constexpr auto dimensionality(C const& c)->decltype(c.dimensionality){return c.dimensionality;}

template<class T, typename = std::enable_if_t<not has_dimensionality_member<T>{}>>
constexpr auto dimensionality(T const&, void* = 0){return 0;}

template<class T, std::size_t N>
constexpr auto dimensionality(T const(&t)[N]){return 1 + dimensionality(t[0]);}

template<class T, std::size_t N>
constexpr auto dimensionality(std::array<T, N> const&){return 1 + dimensionality<T>();}

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

template<class T>
constexpr auto corigin(const T& t){return &t;}
template<class T, std::size_t N>
constexpr auto corigin(const T(&t)[N]) noexcept{return corigin(t[0]);}

template<class T, typename = decltype(std::declval<T>().extensions())>
std::true_type has_extensions_aux(T const&);
inline std::false_type has_extensions_aux(...);

template<class T> struct has_extensions : decltype(has_extensions_aux(std::declval<T>())){};

template<class T, typename = std::enable_if_t<has_extensions<T>{}> >
auto extensions(T const& t)
->decltype(t.extensions()){
	return t.extensions();}

template<class T, typename = std::enable_if_t<not has_extensions<T>{}> >
constexpr std::tuple<> extensions(T const&){return {};}

template<class T, size_t N>
constexpr auto extensions(T(&t)[N]){return tuple_cat(std::make_tuple(index_extension(N)), extensions(t[0]));}

template<class T, size_t N>
constexpr auto extensions(std::array<T, N> const& t){return tuple_cat(std::make_tuple(index_extension(N)), extensions(t[0]));}

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

template<class T1, class T2> auto extensions_t(T2 const& t2){
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
	static auto call(T2 const& t2){return std::tuple_cat(std::make_tuple(t2.extension()), extensions_t<T1>(*begin(t2)));}
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


}}

#if _TEST_MULTI_UTILITY

#include<cassert>
#include<iostream>
#include<vector>
#include<cmath>

using std::cout;
namespace multi = boost::multi;

template<class T>
void f(T&& t){
	using multi::dimensionality;
	std::cout<< dimensionality(t) <<'\n';
}

template<class T> void f();
int main(){

	using T4 = typename multi::array_traits<typename std::remove_reference_t<double(&)[4][4]>>::reference;
	f<T4>();

	using multi::corigin;
	using multi::dimensionality;
	using multi::extension;
	using multi::origin;
	using multi::size;
	using multi::sizes;
	using multi::num_elements;
;{
	double A[4] = {1.,2.,3.,4.};
	assert( dimensionality(A) == 1 );
	assert( extension(A).first() == 0 );
	assert( extension(A).last() == 4 );
	assert( origin(A) == &A[0] );
	assert( size(A) == 4 );
	assert( std::get<0>(sizes(A)) == size(A) );
	using std::addressof;
//	using multi::data;
	using multi::data_elements;
	static_assert( std::is_same<decltype(data_elements(A)), double*>{} );
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
	assert( origin(A) == &A[0][0] );
	*origin(A) = 99.;
	assert( A[0][0] == 99. );	
	assert( corigin(A) == &A[0][0] );
	assert( size(A) == 2 );
	assert( std::get<0>(sizes(A)) == size(A) );
	assert( num_elements(A) == 6. );
}{
	double const A[2][3] = {{1.,2.,3.},{4.,5.,6.}};
	assert( origin(A) == &A[0][0] );
	assert( origin(A) == &A[0][0] );
//	*origin(A) = 99.;
}{

	static_assert( multi::rank<std::array<double, 10>>{} == 1, "!" );
	static_assert( multi::rank<std::array<std::array<double, 2>, 10>>{} == 2, "!" );
	std::array<std::array<double, 2>, 10> a;
	auto x = multi::extensions(a);
	assert( std::get<0>(x) == 10 );
	assert( std::get<1>(x) == 2 );
	std::cout << multi::num_elements(a) << std::endl;
	assert( multi::num_elements(a) == 20 );
}

}

#endif
#endif

