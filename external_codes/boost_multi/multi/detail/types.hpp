#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra -D_TEST_MULTI_DETAIL_TYPES $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_DETAIL_TYPES_HPP
#define MULTI_DETAIL_TYPES_HPP

//#include "detail.hpp"
#include "index_range.hpp"

#include<tuple> // make_tuple
#include<array>
#include<cassert>
#include<cstddef>
#include<type_traits> // make_signed_t

namespace boost{
namespace multi{

namespace detail{

template<typename, typename>
struct append_to_type_seq{};

template<typename T, typename... Ts, template<typename...> class TT>
struct append_to_type_seq<T, TT<Ts...>>{
    using type = TT<Ts..., T>;
};

template<typename T, unsigned int N, template<typename...> class TT = std::tuple> 
struct repeat{
    using type = typename
        append_to_type_seq<
            T,
            typename repeat<T, N-1, TT>::type
        >::type;
};

template<typename T, template<typename...> class TT>
struct repeat<T, 0, TT>{
	using type = TT<>;
};

template<class T, std::size_t N>
auto array_size_impl(const std::array<T, N>&) 
    -> std::integral_constant<std::size_t, N>;

template<class... T>
auto array_size_impl(const std::tuple<T...>&) 
    -> std::integral_constant<std::size_t, std::tuple_size<std::tuple<T...>>{}>;

template<class Array>
using array_size = decltype(array_size_impl(std::declval<const Array&>()));

template<class Array>
constexpr auto static_size() -> decltype(array_size<Array>::value){
    return array_size<Array>::value;
}
template<class Array>
constexpr auto static_size(Array const&) -> decltype(static_size<Array>()){
    return static_size<Array>();
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
->decltype(tail_impl(std::make_index_sequence<(static_size<Tuple>())-1>(), t)){
	return tail_impl(std::make_index_sequence<(static_size<Tuple>())-1>(), t);}
//->decltype(tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t)){
//	return tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t);}

template<typename T, std::size_t N>
std::array<T, N-1> tail(std::array<T, N> const& a){
	std::array<T, N-1> ret;
	std::copy(a.begin() + 1, a.end(), ret.begin());
	return ret;
}

}

using size_type = std::make_signed_t<std::size_t>;

using index               = std::make_signed_t<size_type>;
using difference_type     = std::make_signed_t<index>;
using index_range         = range<index>;
using index_extension     = extension_t<index>;
using dimensionality_type = index;

using iextension = index_extension;
using irange     = index_range;

template<dimensionality_type D> using index_extensions = typename detail::repeat<index_extension, D>::type;
//template<dimensionality_type D> using iextensions = index_extensions<D>;

template<dimensionality_type D> 
struct iextensions : detail::repeat<index_extension, D>::type{
	using base_ = typename detail::repeat<index_extension, D>::type;
	using base_::base_;
//	template<class... Args, typename = std::enable_if_t<sizeof...(Args)==D>>
//	iextensions(Args... args) : detail::repeat<index_extension, D>::type{args...}{}
	template<class T>
	iextensions(std::array<T, D> const& arr) : iextensions(arr, std::make_index_sequence<D>{}){}//detail::repeat<index_extension, D>::type{as_tuple(arr)}{}
	iextensions(std::array<iextension, D> const& arr) : iextensions(arr, std::make_index_sequence<D>{}){}//detail::repeat<index_extension, D>::type{as_tuple(arr)}{}
	base_ const& base() const{return *this;}
	friend decltype(auto) base(iextensions const& s){return s.base();}
private:
	template <class T, size_t... Is> 
	iextensions(std::array<T, D> const& arr, std::index_sequence<Is...>) : iextensions{arr[Is]...}{}
};

//template<dimensionality_type D>
//using extensions_t = iextensions<D>;

#if __cpp_deduction_guides
template<class... Args> iextensions(Args...) -> iextensions<sizeof...(Args)>;
#endif

template<dimensionality_type D, class Tuple>
auto contains(index_extensions<D> const& ie, Tuple const& tp){
//	using detail::head;
//	using detail::tail;
	return contains(head(ie), head(tp)) and contains(tail(ie), tail(tp));
}

}}

#if _TEST_MULTI_DETAIL_TYPES

#include<cassert>
#include<iostream>
#include<vector>

using std::cout;
namespace multi = boost::multi;

int main(){}
#endif
#endif

