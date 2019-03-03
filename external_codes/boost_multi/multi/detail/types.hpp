#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra -D_TEST_MULTI_DETAIL_TYPES $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_DETAIL_TYPES_HPP
#define MULTI_DETAIL_TYPES_HPP

//#include "detail.hpp"
#include "index_range.hpp"

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
template<dimensionality_type D> using iextensions = index_extensions<D>;

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

