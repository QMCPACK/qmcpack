#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++11 -Wfatal-errors -D_TEST_BOOST_MPI3_DETAIL_DATATYPE -lboost_serialization $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp $0x.x; exit
#endif
#ifndef BOOST_MPI3_DETAIL_DATATYPE_HPP
#define BOOST_MPI3_DETAIL_DATATYPE_HPP

#include<mpi.h>

#include<boost/serialization/strong_typedef.hpp>

#include<type_traits>
#include<utility> // pair
#include<tuple>
#include<complex>

namespace boost{
namespace mpi3{
namespace detail{

using float_int       = std::pair<float, int>;
using long_int        = std::pair<long, int>;
using double_int      = std::pair<double, int>;
using short_int       = std::pair<short, int>;
using int_int         = std::pair<int, int>;
using long_double_int = std::pair<long double, int>;

using float_float             = std::pair<float, float>;
using double_double           = std::pair<double, double>;
using long_double_long_double = std::pair<long double, long double>;

using cxx_float_complex       = std::complex<float>;
using cxx_double_complex      = std::complex<double>;
using cxx_long_double_complex = std::complex<long double>;

using cxx_bool = bool;

using wchar = wchar_t;

#if(__cplusplus >= 201703L)
using byte = std::byte;
#endif

BOOST_STRONG_TYPEDEF(unsigned char, packed);

template<class T> struct basic_datatype;

#define MPI3_DECLARE_DATATYPE(TypE, MpiiD) \
template<> struct basic_datatype<TypE>{ \
/*	constexpr*/ operator MPI_Datatype() const{return MpiiD;} \
/*	static constexpr MPI_Datatype value = MpiiD;*/ \
}

// basic data types http://beige.ucs.indiana.edu/I590/node100.html
MPI3_DECLARE_DATATYPE(char, MPI_CHAR);
MPI3_DECLARE_DATATYPE(unsigned char, MPI_UNSIGNED_CHAR);

#if(__cplusplus >= 201703L)
MPI3_DECLARE_DATATYPE(byte, MPI_BYTE);
#endif
MPI3_DECLARE_DATATYPE(wchar, MPI_WCHAR);

MPI3_DECLARE_DATATYPE(short, MPI_SHORT);
MPI3_DECLARE_DATATYPE(unsigned short, MPI_UNSIGNED_SHORT);
MPI3_DECLARE_DATATYPE(int, MPI_INT);
MPI3_DECLARE_DATATYPE(unsigned int, MPI_UNSIGNED);
MPI3_DECLARE_DATATYPE(long, MPI_LONG);
MPI3_DECLARE_DATATYPE(unsigned long, MPI_UNSIGNED_LONG);
MPI3_DECLARE_DATATYPE(float, MPI_FLOAT);
MPI3_DECLARE_DATATYPE(double, MPI_DOUBLE);
MPI3_DECLARE_DATATYPE(long double, MPI_LONG_DOUBLE);
MPI3_DECLARE_DATATYPE(long long int, MPI_LONG_LONG_INT);

MPI3_DECLARE_DATATYPE(cxx_float_complex, MPI_CXX_FLOAT_COMPLEX);
MPI3_DECLARE_DATATYPE(cxx_double_complex, MPI_CXX_DOUBLE_COMPLEX);
MPI3_DECLARE_DATATYPE(cxx_long_double_complex, MPI_CXX_DOUBLE_COMPLEX);

MPI3_DECLARE_DATATYPE(float_float, MPI_CXX_FLOAT_COMPLEX);
	static_assert(sizeof(std::pair<float, float>) == sizeof(std::complex<float>), "checking that complex mem layout maps to pair");
MPI3_DECLARE_DATATYPE(double_double, MPI_CXX_DOUBLE_COMPLEX);
	static_assert(sizeof(std::pair<double, double>) == sizeof(std::complex<double>), "checking that complex mem layout maps to pair");
MPI3_DECLARE_DATATYPE(decltype(std::tuple<double,double>{}), MPI_CXX_DOUBLE_COMPLEX);
MPI3_DECLARE_DATATYPE(long_double_long_double, MPI_CXX_DOUBLE_COMPLEX);
	static_assert(sizeof(std::pair<long double, long double>) == sizeof(std::complex<long double>), "checking that complex mem layout maps to pair");

MPI3_DECLARE_DATATYPE(float_int, MPI_FLOAT_INT);
MPI3_DECLARE_DATATYPE(long_int, MPI_LONG_INT);
MPI3_DECLARE_DATATYPE(double_int, MPI_DOUBLE_INT);
MPI3_DECLARE_DATATYPE(short_int, MPI_SHORT_INT);
MPI3_DECLARE_DATATYPE(int_int, MPI_2INT);
MPI3_DECLARE_DATATYPE(long_double_int, MPI_LONG_DOUBLE_INT);

//BOOST_MPI3_DECLARE_DATATYPE(std::intptr_t, MPI_AINT);
//BOOST_MPI3_DECLARE_DATATYPE(std::size_t, MPI_AINT);
MPI3_DECLARE_DATATYPE(void*, MPI_AINT);

//BOOST_MPI3_DECLARE_DATATYPE(bool, MPI_INT);
MPI3_DECLARE_DATATYPE(bool, MPI_CXX_BOOL);

MPI3_DECLARE_DATATYPE(packed, MPI_PACKED);

// LB 
// UB

#undef BOOST_MPI3_DECLARE_DATATYPE

template<class T, class = decltype(basic_datatype<T>{})>
std::true_type  is_basic_aux(T  );
std::false_type is_basic_aux(...);

template<class T> 
struct is_basic : decltype(is_basic_aux(std::declval<T>())){};

}}}

#if _TEST_BOOST_MPI3_DETAIL_DATATYPE

#include<iostream>
#include<string>
#include<boost/type_index.hpp>

namespace mpi3 = boost::mpi3;
using std::cout;

int main(){

	using mpi3::detail::is_basic;

	static_assert( is_basic<int>{}, "");
	static_assert( is_basic<double>{}, "");
	static_assert( is_basic<mpi3::detail::float_int>{}, "");

	static_assert( not is_basic<std::string>{}, "");

}

#endif
#endif

