#ifndef BOOST_MPI3_DETAIL_VALUE_TRAITS_HPP
#define BOOST_MPI3_DETAIL_VALUE_TRAITS_HPP

//#include "../../mpi3/detail/package_archive.hpp"
//#include "../../mpi3/communicator.hpp"

#include "./datatype.hpp"
#include "./iterator.hpp"
//#include "./just.hpp"

#include<iterator>
#include<type_traits>

namespace boost {
namespace mpi3 {
namespace detail {

/*
template<class T> 
struct is_memcopyable : std::integral_constant<bool, std::is_trivially_copyable<T>{}>{};

template<class T, class... Ts> 
struct is_memcopyable<std::tuple<T, Ts...>> : 
    std::integral_constant<bool, 
        is_memcopyable<T>{} and is_memcopyable<std::tuple<Ts...>>{}
    >
{};

template<class T1, class T2> 
struct is_memcopyable<std::pair<T1, T2>> : 
    std::integral_constant<bool, 
        is_memcopyable<T1>{} and is_memcopyable<T2>{}
    >
{};
*/
/*
template<
	class T, 
	typename = decltype(std::declval<mpi3::detail::package_oarchive&>() << std::declval<T>()),
	typename = decltype(std::declval<mpi3::detail::package_iarchive&>() >> std::declval<T>())
>
std::true_type  is_serializable_aux(T);
std::false_type is_serializable_aux(...);
*/

template<class T>
struct is_serializable : decltype(is_serializable_aux(std::declval<T>())){};

struct value_unspecified_tag{};
struct nonmemcopyable_serializable_tag : value_unspecified_tag{
	using base = value_unspecified_tag;
};
struct memcopyable_tag : value_unspecified_tag{
	using base = value_unspecified_tag;
};
struct basic_tag : memcopyable_tag{using base = memcopyable_tag;};

//template<class V, typename = std::enable_if_t<is_memcopyable<V>{} and not is_basic<V>{}>>
//memcopyable_tag value_category_aux(V&&);
template<class V, typename = std::enable_if_t<is_basic<V>{}>>
basic_tag value_category_aux(V const&);

template<class V, std::size_t N>
value_unspecified_tag value_category_aux(V[N]);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : legacy compatibility
value_unspecified_tag value_category_aux(... );

template<class V>
struct value_category{using type = decltype(value_category_aux(std::declval<V>()));};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class V>
using value_category_t = typename value_category<V>::type;

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_DETAIL_VALUE_TRAITS

//#include "../../mpi3/environment.hpp"
//#include "../../mpi3/main.hpp"

//#include<deque>
//#include<list>
//#include<vector>

//#include<iostream>

//namespace mpi3 = boost::mpi3;
//using std::cout;

//template<class It>
//std::string f(It it, mpi3::detail::value_unspecified_tag){
//	return "memcopyable_tag";
//};

//template<class It>
//std::string f(It it, mpi3::detail::memcopyable_tag){
//	return "memcopyable_tag";
//};

//template<class It>
//std::string f(It it, mpi3::detail::basic_tag const&){
//	return "basic_tag";
//};

//template<class It> std::string f(It&& it){
//	return f(
//		std::forward<It>(it),
//		typename boost::mpi3::detail::value_category<typename std::iterator_traits<It>::value_type>::type{}
//	);
//}

//int mpi3::main(int, char*[], mpi3::communicator){

//	{
//		std::list<std::tuple<double, double>> l1;
//	//	assert( f(begin(l1)) == "memcopyable_tag" );
//	//	std::list<double> l2;
//	//	assert( f(begin(l2)) == "basic_tag" );
//	//	static_assert(std::is_trivially_copyable<std::complex<double>>{}, "complex is not trivially copyable?");
//	//	assert( f(nullptr, mpi3::detail::value_category_t<std::tuple<double, double>>{}) == "memcopyable_tag" );
//	//	assert( f(nullptr, mpi3::detail::value_category_t<std::string>{}) == "memcopyable_tag" );
//	}

//	return 0;
//}
//#endif
#endif

