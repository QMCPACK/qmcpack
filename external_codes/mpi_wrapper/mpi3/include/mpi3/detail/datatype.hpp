// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2017-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_DETAIL_DATATYPE_HPP
#define BOOST_MPI3_DETAIL_DATATYPE_HPP

// #define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#if defined(__NVCC__)
#include <thrust/complex.h>
#endif

#include<cassert>
#include<complex>
#include<cstddef> // std::byte
#include<tuple>
#include<type_traits>
#include<utility> // pair

namespace boost {
namespace mpi3 {

template<class T>
struct vlp {  // value location pair
	T value;  // NOLINT(misc-non-private-member-variables-in-classes)
	int location;  // NOLINT(misc-non-private-member-variables-in-classes)
	bool operator<(vlp<T> const& other) const {
	//  if(value == other.value) {return location < other.location;}  // partial order on purpose?
		return value < other.value;
	}
};

namespace detail {

using float_int       = std::pair<float      , int>;
using long_int        = std::pair<long       , int>;  // NOLINT(google-runtime-int) : long <-> int64
using double_int      = std::pair<double     , int>;
using short_int       = std::pair<short      , int>;  // NOLINT(google-runtime-int) : short <-> int16
using int_int         = std::pair<int        , int>;

using long_double_int = std::pair<long double, int>;

using float_float             = std::pair<float, float>;
using double_double           = std::pair<double, double>;
using long_double_long_double = std::pair<long double, long double>;

using cxx_float_complex       = std::complex<float>;
using cxx_double_complex      = std::complex<double>;
using cxx_long_double_complex = std::complex<long double>;

// using cxx_2double_complex     = std::pair<std::complex<double>, std::complex<double>>;

using cxx_bool = bool;

using wchar = wchar_t;

#if(__cplusplus >= 201703L)
using byte = std::byte;
#endif

class packed {
	unsigned char t{};

 public:
    explicit packed(unsigned char t_) : t{t_} {};
    packed() = default;

	packed& operator=(unsigned char const& rhs) {t = rhs; return *this;}

    explicit operator const unsigned char&() const {return t;}
    explicit operator       unsigned char&()       {return t;}

    bool operator==(packed const& rhs) const {return t == rhs.t;}
    bool operator< (packed const& rhs) const {return t < rhs.t;}
};

template<class T> struct basic_datatype;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MPI3_DECLARE_DATATYPE(TypE, MpiiD) \
template<> struct basic_datatype<TypE> { \
/*	constexpr*/ operator MPI_Datatype() const { \
	assert(MPI_DOUBLE_COMPLEX != MPI_DATATYPE_NULL );  /* NOLINT(cert-dcl03-c,hicpp-static-assert,misc-static-assert)*/ \
	assert( (MpiiD) != MPI_DATATYPE_NULL );            /* NOLINT(cert-dcl03-c,hicpp-static-assert,misc-static-assert) in some MPI distros this is not constexpr */ /*this system doesn't support this type*/ \
		return MpiiD; \
	} \
	auto get() const -> MPI_Datatype {return MpiiD;} \
/*	static constexpr MPI_Datatype value = MpiiD;*/ \
}

// basic data types http://beige.ucs.indiana.edu/I590/node100.html
MPI3_DECLARE_DATATYPE(char                   , MPI_CHAR);
MPI3_DECLARE_DATATYPE(unsigned char          , MPI_UNSIGNED_CHAR);

#if(__cplusplus >= 201703L)
MPI3_DECLARE_DATATYPE(byte                   , MPI_BYTE);
#endif
MPI3_DECLARE_DATATYPE(wchar                  , MPI_WCHAR);

MPI3_DECLARE_DATATYPE(short                  , MPI_SHORT);
MPI3_DECLARE_DATATYPE(unsigned short         , MPI_UNSIGNED_SHORT);
MPI3_DECLARE_DATATYPE(int                    , MPI_INT);
MPI3_DECLARE_DATATYPE(unsigned int           , MPI_UNSIGNED);
MPI3_DECLARE_DATATYPE(long                   , MPI_LONG);
MPI3_DECLARE_DATATYPE(unsigned long          , MPI_UNSIGNED_LONG);
MPI3_DECLARE_DATATYPE(float                  , MPI_FLOAT);
MPI3_DECLARE_DATATYPE(double                 , MPI_DOUBLE);
MPI3_DECLARE_DATATYPE(long double            , MPI_LONG_DOUBLE);
MPI3_DECLARE_DATATYPE(long long int          , MPI_LONG_LONG_INT);

MPI3_DECLARE_DATATYPE(bool                   , MPI_C_BOOL);  // C++ binding not used MPI_CXX_BOOL);

// MPI_INT8_T	int8_t
// MPI_INT16_T	int16_t
// MPI_INT32_T	int32_t
// MPI_INT64_T	int64_t
// MPI_UINT8_T	uint8_t
// MPI_UINT16_T	uint16_t
// MPI_UINT32_T	uint32_t
// MPI_UINT64_T	uint64_t

MPI3_DECLARE_DATATYPE(cxx_float_complex      , MPI_C_FLOAT_COMPLEX);
MPI3_DECLARE_DATATYPE(cxx_double_complex     , MPI_C_DOUBLE_COMPLEX);
MPI3_DECLARE_DATATYPE(cxx_long_double_complex, MPI_C_LONG_DOUBLE_COMPLEX);

// MPI3_DECLARE_DATATYPE(cxx_2double_complex    , MPI_2DOUBLE_COMPLEX);  // not available in mpich

// TODO(correaa) these types below probably don't behave correctly for reductions with multiplication

MPI3_DECLARE_DATATYPE(float_float            , MPI_COMPLEX);  static_assert(sizeof(std::pair<float, float>) == sizeof(std::complex<float>), "checking that complex mem layout maps to pair");
MPI3_DECLARE_DATATYPE(double_double          , MPI_DOUBLE_COMPLEX); static_assert(sizeof(std::pair<double, double>) == sizeof(std::complex<double>), "checking that complex mem layout maps to pair");
MPI3_DECLARE_DATATYPE(decltype(std::tuple<double,double>{}), MPI_DOUBLE_COMPLEX);
MPI3_DECLARE_DATATYPE(long_double_long_double, MPI_DOUBLE_COMPLEX); static_assert(sizeof(std::pair<long double, long double>) == sizeof(std::complex<long double>), "checking that complex mem layout maps to pair");

#if defined(__NVCC__)
MPI3_DECLARE_DATATYPE(thrust::complex<double>, MPI_DOUBLE_COMPLEX);
#endif

MPI3_DECLARE_DATATYPE(float_int              , MPI_FLOAT_INT);
MPI3_DECLARE_DATATYPE(long_int               , MPI_LONG_INT);
MPI3_DECLARE_DATATYPE(double_int             , MPI_DOUBLE_INT);
MPI3_DECLARE_DATATYPE(short_int              , MPI_SHORT_INT);
MPI3_DECLARE_DATATYPE(int_int                , MPI_2INT);
MPI3_DECLARE_DATATYPE(long_double_int        , MPI_LONG_DOUBLE_INT);

MPI3_DECLARE_DATATYPE(vlp<float>             , MPI_FLOAT_INT);
MPI3_DECLARE_DATATYPE(vlp<long>              , MPI_LONG_INT);
MPI3_DECLARE_DATATYPE(vlp<double>            , MPI_DOUBLE_INT);
MPI3_DECLARE_DATATYPE(vlp<short>             , MPI_SHORT_INT);
MPI3_DECLARE_DATATYPE(vlp<int>               , MPI_2INT);
MPI3_DECLARE_DATATYPE(vlp<long double>       , MPI_LONG_DOUBLE_INT);

//BOOST_MPI3_DECLARE_DATATYPE(std::intptr_t  , MPI_AINT);
//BOOST_MPI3_DECLARE_DATATYPE(std::size_t    , MPI_AINT);
MPI3_DECLARE_DATATYPE(void*                  , MPI_AINT);

MPI3_DECLARE_DATATYPE(packed                 , MPI_PACKED);

// LB
// UB

#undef BOOST_MPI3_DECLARE_DATATYPE

template<class T, class = decltype(basic_datatype<T>{})>
std::true_type  is_basic_aux(T  );
std::false_type is_basic_aux(...);

template<class T>
struct is_basic : decltype(is_basic_aux(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class T> constexpr bool is_basic_v = is_basic<T>::value;

}  // end namespace detail

template<class T> class datatype;

template<class T> class default_datatype {
 public:
	template<class TT = T, class R = decltype(boost::mpi3::detail::basic_datatype<TT>{}.get())>
	R operator()() const { return boost::mpi3::detail::basic_datatype<T>{}.get(); }
	template<class TT = T, class R = decltype(boost::mpi3::detail::basic_datatype<TT>{}.get())>
	R get() const { return boost::mpi3::detail::basic_datatype<T>{}.get(); }
};

template<class T>
auto datatype_detect(...) -> default_datatype<T>;

template<class T, class U, class RealType = decltype(U{}.real()), class = decltype(U{}.imag())>
auto datatype_detect(U const&) -> default_datatype<std::complex<RealType>>;

template<class T> class datatype :  public decltype(datatype_detect<T>(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg) detection idiom

template<class T, class = decltype(datatype<T>{}())>
std::true_type  has_datatype_aux(T  );
std::false_type has_datatype_aux(...);

template<class T>
struct has_datatype : decltype(has_datatype_aux(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class T> constexpr bool has_datatype_v = has_datatype<T>::value;

}  // end namespace mpi3
}  // end namespace boost
#endif
