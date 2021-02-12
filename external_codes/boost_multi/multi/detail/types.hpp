#ifdef COMPILATION_INSTRUCTIONS
$CXXX $CXXFLAGS $0 -o $0$X&&$0$X&&rm $0$X;exit
#endif
//  © Alfredo A. Correa 2018-2020

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

using size_type = std::make_signed_t<std::size_t>;

using index               = std::make_signed_t<size_type>;
using difference_type     = std::make_signed_t<index>;
using index_range         = range<index>;
using index_extension     = extension_t<index>;
using dimensionality_type = index;

using iextension = index_extension;
using irange     = index_range;

namespace detail{

template<typename, typename>
struct append_to_type_seq{};

template<typename T, typename... Ts, template<typename...> class TT>
struct append_to_type_seq<T, TT<Ts...>>{
    using type = TT<Ts..., T>;
};

template<typename T, dimensionality_type N, template<typename...> class TT = std::tuple> 
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
constexpr auto array_size_impl(const std::array<T, N>&) 
    -> std::integral_constant<std::size_t, N>;

template<class... T>
constexpr auto array_size_impl(const std::tuple<T...>&) 
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

//TODO consolidate with tuple_tail defined somewhere else
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
//->decltype(tail_impl(std::make_index_sequence<(static_size<Tuple>())-1>(), t)){
{	return tail_impl(std::make_index_sequence<(static_size<Tuple>())-1>(), t);}
//->decltype(tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t)){
//	return tail_impl(std::make_index_sequence<(std::tuple_size<Tuple>{})-1>(), t);}

template<typename T, std::size_t N>
constexpr std::array<T, N-1> tail(std::array<T, N> const& a){
	std::array<T, N-1> ret;
	std::copy(a.begin() + 1, a.end(), ret.begin());
	return ret;
}

}

template<typename T, dimensionality_type D>
struct initializer_list{
	using type = std::initializer_list<typename initializer_list<T, D-1>::type>;
};
template<typename T>
struct initializer_list<T, 1>{using type = std::initializer_list<T>;};

template<typename T, dimensionality_type D>
using initializer_list_t = typename initializer_list<T, D>::type;

template<dimensionality_type D> using index_extensions = typename detail::repeat<index_extension, D>::type;
//template<dimensionality_type D> using iextensions = index_extensions<D>;

template<dimensionality_type D> 
struct iextensions : detail::repeat<index_extension, D>::type{
	static constexpr dimensionality_type dimensionality = D;
	using base_ = typename detail::repeat<index_extension, D>::type;
	using base_::base_;
//	template<class... Args, typename = std::enable_if_t<sizeof...(Args)==D>>
//	iextensions(Args... args) : detail::repeat<index_extension, D>::type{args...}{}
	iextensions() = default;
	template<class T>
	constexpr iextensions(std::array<T, static_cast<std::size_t>(D)> const& arr) : iextensions(arr, std::make_index_sequence<static_cast<std::size_t>(D)>{}){}//detail::repeat<index_extension, D>::type{as_tuple(arr)}{}
	constexpr iextensions(std::array<iextension, static_cast<std::size_t>(D)> const& arr) : iextensions(arr, std::make_index_sequence<static_cast<std::size_t>(D)>{}){}//detail::repeat<index_extension, D>::type{as_tuple(arr)}{}
	constexpr base_ const& base() const{return *this;}
	friend constexpr decltype(auto) base(iextensions const& s){return s.base();}
private:
//	template<std::size_t... Ns> bool bool_aux(std::index_sequence<Ns...>) const{return (not std::get<0>(*this).empty() and tail(*this));}
public:
//	explicit operator bool() const{return bool_aux(std::make_index_sequence<D>());}
//	bool is_empty() const{return bool_aux(std::make_index_sequence<D>());}
//	bool empty() const{return is_empty();}
private:
	template <class T, size_t... Is> 
	constexpr iextensions(std::array<T, static_cast<std::size_t>(D)> const& arr, std::index_sequence<Is...>) : iextensions{arr[Is]...}{}
};

#if defined(__cpp_deduction_guides) and __cpp_deduction_guides >= 201703
template<class... Args> iextensions(Args...) -> iextensions<sizeof...(Args)>;
#endif

template<dimensionality_type D, class Tuple>
constexpr auto contains(index_extensions<D> const& ie, Tuple const& tp){
//	using detail::head;
//	using detail::tail;
	return contains(head(ie), head(tp)) and contains(tail(ie), tail(tp));
}

}}

#if defined(__cpp_structured_bindings) and __cpp_structured_bindings>=201606
namespace std{ // this is for structured binding
	template<boost::multi::dimensionality_type D> struct tuple_size<boost::multi::iextensions<D>> : std::integral_constant<size_t, D> { };
	template<size_t N, boost::multi::dimensionality_type D> struct tuple_element<N, boost::multi::iextensions<D>> : tuple_element<N, typename boost::multi::iextensions<D>::base_>{};
}
#endif

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#include<range/v3/begin_end.hpp>
//#include<range/v3/utility/concepts.hpp>

#include<cassert>
#include<iostream>
#include<numeric> // accumulate
#include<vector>


using std::cout;
namespace multi = boost::multi;


template<class T> T what(T&&) = delete;

int main(){

	multi::index_extension x(10);

	assert( *begin(x) == 0 );
	assert( size(x) == 10 );
	assert( x[0] == 0 );
	assert( x[1] == 1 );
	assert( x[9] == 9 );
	
	auto b = begin(x);
	assert( b[0] == x[0] );
	assert( b[1] == x[1] );

//	static_assert( ranges::forward_iterator< std::decay_t<decltype(b)> > , "!");

	assert( std::accumulate( begin(x), end(x), 0) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

	std::iterator_traits<std::decay_t<decltype(begin(x))>>::difference_type d; (void)d;
//	for(auto i : x) std::cout << i << std::endl;

	{
		multi::iextensions<3> ies({{0, 3}, {0, 4}, {0, 5}});
		assert( std::get<1>(ies).size() == 4 );
		auto [is, js, ks] = ies;
		assert( is.size() == 3 );
	}

}


#endif
#endif

