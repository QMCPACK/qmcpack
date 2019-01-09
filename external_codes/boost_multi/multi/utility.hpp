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

#if __cplusplus < 201703L
template<class Container>
constexpr auto size(Container const& con)
->decltype(con.size()){
	return con.size();}
#else
#endif

template <class T, std::size_t N>
constexpr auto size(const T(&)[N]) noexcept{return multi::size_type{N};}

template<class A> 
constexpr auto num_elements(A const& arr)
->decltype(arr.num_elements()){
	return arr.num_elements();}

constexpr size_type num_elements(...){return 1;}
template<class T, std::size_t N>
constexpr auto num_elements(const T(&t)[N]) noexcept{return size(t)*num_elements(t[0]);}

template <class T, std::size_t N>
constexpr auto stride(const T(&t)[N]) noexcept{return num_elements(t[0]);}

template <class T, std::size_t N>
constexpr std::ptrdiff_t offset(const T(&)[N]) noexcept{return 0;}

template<class Container>
auto extension(Container const& c) // TODO consider "extent"
->decltype(range<decltype(size(c))>{0, size(c)}){
	return range<decltype(size(c))>{0, size(c)};}

template<class Container>
constexpr auto dimensionality(Container const& con)
->decltype(con.dimensionality()){
	return con.dimensionality();}
	
template<class Container>
constexpr auto dimensionality(Container const&)
->decltype(Container::dimensionality){
	return Container::dimensionality;}

constexpr auto dimensionality(...){return 0;}
	
template<class T, std::size_t N>
constexpr auto dimensionality(const T(&t)[N]){return 1 + dimensionality(t[0]);}

template<class Array>
constexpr auto sizes(Array const& arr)
->decltype(arr.sizes()){
	return arr.sizes();}

template<class T>
constexpr std::tuple<> sizes(...){return {};}

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

template<class T>
auto extensions(T const& t)
->decltype(t.extensions()){
	return t.extensions();}

inline std::array<index_extension, 0> 
extensions(...){return {};}

template<dimensionality_type> struct extension_aux;

template<dimensionality_type D, class T>
auto extensions(T const& t);

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
std::false_type has_layout_member_aux(...);

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

using std::cout;
namespace multi = boost::multi;

template<class T>
void f(T&& t){
	using multi::dimensionality;
	std::cout<< dimensionality(t) <<'\n';
}

int main(){

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
//	extensions(A[
	assert( origin(A) == &A[0] );
	assert( size(A) == 4 );
	assert( std::get<0>(sizes(A)) == size(A) );
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

};

}

#endif
#endif

