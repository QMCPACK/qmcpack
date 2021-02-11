#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0.$X `pkg-config --cflags --libs cudart-11.0 cublas-11.0 blas` -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_ADAPTORS_CUDA_CUBLAS_ERROR_HPP
#define MULTI_ADAPTORS_CUDA_CUBLAS_ERROR_HPP

#include<cublas_v2.h> // cublasStatus_t

#include<string>
#include<system_error> // std::error_category
#include<type_traits> // std::underlying_type

namespace boost{
namespace multi::cuda::cublas{

enum class error : typename std::underlying_type<cublasStatus_t>::type{
	success               = CUBLAS_STATUS_SUCCESS,
	not_initialized       = CUBLAS_STATUS_NOT_INITIALIZED,
	allocation_failed     = CUBLAS_STATUS_ALLOC_FAILED,
	invalid_value         = CUBLAS_STATUS_INVALID_VALUE,
	architecture_mismatch = CUBLAS_STATUS_ARCH_MISMATCH,
	mapping_error         = CUBLAS_STATUS_MAPPING_ERROR,
	execution_failed      = CUBLAS_STATUS_EXECUTION_FAILED,
	internal_error        = CUBLAS_STATUS_INTERNAL_ERROR,
	not_supported         = CUBLAS_STATUS_NOT_SUPPORTED,
	license_error         = CUBLAS_STATUS_LICENSE_ERROR
};

std::string inline error_string(enum cublas::error err){ //https://stackoverflow.com/questions/13041399/equivalent-of-cudageterrorstring-for-cublas
	switch(err){
	case cublas::error::success              : return "CUBLAS_STATUS_SUCCESS"         ;
	case cublas::error::not_initialized      : return "CUBLAS_STATUS_NOT_INITIALIZED" ;
	case cublas::error::allocation_failed    : return "CUBLAS_STATUS_ALLOC_FAILED"    ;
	case cublas::error::invalid_value        : return "CUBLAS_STATUS_INVALID_VALUE"   ;
	case cublas::error::architecture_mismatch: return "CUBLAS_STATUS_ARCH_MISMATCH"   ;
	case cublas::error::mapping_error        : return "CUBLAS_STATUS_MAPPING_ERROR"   ;
	case cublas::error::execution_failed     : return "CUBLAS_STATUS_EXECUTION_FAILED";
	case cublas::error::internal_error       : return "CUBLAS_STATUS_INTERNAL_ERROR"  ;
	case cublas::error::not_supported        : return "CUBLAS_STATUS_NOT_SUPPORTED"   ;
	case cublas::error::license_error        : return "CUBLAS_STATUS_LICENSE_ERROR"   ;
	}
	return "cublas status <unknown>";
}

struct error_category : std::error_category{
	char const* name() const noexcept override{return "cublas wrapper";}
	std::string message(int err) const override{return error_string(static_cast<enum cublas::error>(err));}
	static error_category& instance(){static cublas::error_category instance; return instance;}
};

inline std::error_code make_error_code(cublas::error err) noexcept{
	return std::error_code(int(err), cublas::error_category::instance());
}

}
}

namespace std{
	template<> struct is_error_code_enum<::boost::multi::cuda::cublas::error> : true_type{};
}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_CUDA

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../array.hpp"
//#include "../../utility.hpp"

//#include "../../adaptors/cuda.hpp"
//#include "../../adaptors/blas.hpp"
//#include "../../adaptors/blas/cuda.hpp"

#include<cassert>

namespace multi = boost::multi;


BOOST_AUTO_TEST_CASE(multi_cublas_error){

	BOOST_CHECK_THROW(
		throw (std::system_error{multi::cuda::cublas::make_error_code(multi::cuda::cublas::error::not_initialized), "error test"}), 
		std::system_error
	);

}

#endif
#endif

