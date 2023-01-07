//#if COMPILATION_INSTRUCTIONS
//(echo "#include\""$0"\"">$0x.cpp) && mpic++ -O3 -std=c++17 -Wfatal-errors -D_BOOST_MPI3_MAIN_ENVIRONMENT -D_TEST_MPI3_DETAIL_ITERATOR $0x.cpp -o $0x.x && time mpirun -n 1 $0x.x $@ && rm -f $0x.cpp; exit
//#endif
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MPI3_DETAIL_ITERATOR_HPP
#define MPI3_DETAIL_ITERATOR_HPP

#include<iterator> // iterator_traits 

namespace boost{
namespace mpi3{
namespace detail{

template<class T> struct is_declared_contiguous : std::false_type{};

#if 0
template<
	class It, 
	class V = typename std::iterator_traits<It>::value_type,
	typename = std::enable_if_t<
		/**/std::is_convertible<It, typename std::vector<V>::const_iterator>{}
		or  std::is_convertible<It, typename std::basic_string<V>::const_iterator>{} 
		or  std::is_convertible<It, typename boost::container::vector<V>::const_iterator>{}
		or  std::is_pointer<It>{}
		or  std::is_constructible<typename std::iterator_traits<It>::pointer, It>{}
		or  std::is_convertible<It, typename std::iterator_traits<It>::pointer>{}
		or  is_declared_contiguous<It>{}
	>
>
typename std::iterator_traits<It>::pointer
data(It const& it){return std::addressof(*it);}
#endif

template<int N = 10> struct priority : priority<N-1>{};
template<> struct priority<0>{};

template<class It, typename = decltype(std::declval<It>().base())> 
typename std::iterator_traits<It>::pointer
data(It it, priority<5> /*priority*/){return it.base();}

template<class It, typename = decltype(std::declval<It>().get_ptr())> 
typename std::iterator_traits<It>::pointer
data(It it, priority<0> /*priority*/){return it.get_ptr();}

template<class It> auto data(It it) -> decltype(data(it, priority<>{})){
	return data(it, priority<>{});
}

template<class T> auto data(T* t){return t;}

template<class It>
typename std::iterator_traits<It>::value_type const*
cdata(It it){return data(it);}

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_DETAIL_ITERATOR

//#include "../../mpi3/allocator.hpp"
//#include "../../mpi3/main.hpp"
//#include "../../mpi3/vector.hpp"

//#include<boost/container/static_vector.hpp>
//#include<boost/container/flat_set.hpp>
//#include<boost/container/flat_map.hpp>
//#include<boost/container/small_vector.hpp>
//#include<boost/container/stable_vector.hpp>

//#include<boost/array.hpp>

//#include<array>
//#include<list>
//#include<queue>
//#include<set>
//#include<valarray>
//#include<vector>

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int, char*[], mpi3::environment&){

//	using mpi3::detail::data;
//	using mpi3::detail::cdata;

//	{
//		std::vector<double> v(30);
//		assert( data(begin(v)) == &v[0] );
//		assert( data(begin(v) + 10) == &v[10] );
//	}
//	{
//		std::array<int, 4> arr;
//		assert( data(begin(arr)) == &arr[0] );
//		assert( data(begin(arr) + 2) == &arr[2] );
//	}
//	{
//		std::valarray<double> varr(30);
//		varr[0] = 11.;
//		assert( data(begin(varr)) == &varr[0] );
//		assert( data(begin(varr) + 4) == &varr[4] );
//	}
//	{
//		std::set<double> s; 
//		s.insert(13.);
//	//	assert( *data(begin(s)) == 13. ); // data fails here
//	}
//	{
//		std::vector<double> v(30);
//		v[0] = 10.;
//		assert( *data(v.cbegin()) == 10. );
//	}
//	{
//		std::array<int, 4> arr;
//		assert( data(arr.cbegin()) == &arr[0] );
//		*data(arr.begin()) = 12.;
//		assert( arr[0] == 12. );
//	}
//	{
//		std::array<int, 4> arr;
//		arr[0] = 10.;
//		assert( *cdata(arr.begin()) == 10. );
//	}
//	{
//		std::string s = "hola";
//		assert( *data(begin(s)) == 'h' );
//		assert( data(begin(s) + 2) == &s[2] );
//	}
//	// all boost.containers iterator that are contiguous are particular cases of vector
//	{
//		boost::container::static_vector<double, 30> const sv(20);
//		assert( *data(sv.begin()) == 0. );
//	}
//	{
//		boost::container::flat_set<double> fs = {1.,2.,3.}; // flat_set::iterator is the same type as static_vector::iterator
//		assert( data(begin(fs))[1] == 2. );
//	}
//	{
//		boost::container::flat_map<double, int> fm = {{1., 10},{2.,20},{3.,30}};		assert( data(begin(fm))[1] == (std::pair<double, int>(2., 20)) );
//	}
//	{
//		boost::container::small_vector<double, 5> sv(3); sv[1] = 2.;// = {1.,2.,3.};
//		assert( data(sv.begin())[1] == 2. );
//	}
//	{
//		boost::container::stable_vector<double> v(10); v[1] = 5.;
//	//	assert ( data(v.begin())[1] == 5. ); // stable vector is not contiguous
//	}
//	{
//		std::vector<bool> v(30, false); v[2] = true;
//		assert(begin(v)[2] == true);
//	//	assert( data(v.begin())[1] == true ); // vector<bool>::iterator is not contiguous
//	}
//	{
//		assert(mpi3::detail::is_contiguous<std::vector<double>::iterator>{});
//	}
//	{
//		std::vector<double, mpi3::allocator<double>> v(100);
//		assert( data(begin(v)) == &v[0] );
//		assert( mpi3::detail::is_contiguous<std::vector<double>::iterator>{} );
//		assert((mpi3::detail::is_contiguous<std::vector<double, mpi3::uallocator<double>>::iterator>{}) );
//		assert( mpi3::detail::is_contiguous<std::string::iterator>{} );
//		assert( mpi3::detail::is_basic<std::string::iterator::value_type>{} );
//	}
//	return 0;
//}
//#endif
#endif
