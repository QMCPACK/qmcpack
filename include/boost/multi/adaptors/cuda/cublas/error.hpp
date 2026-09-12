// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_CUDA_CUBLAS_ERROR_HPP
#define BOOST_MULTI_ADAPTORS_CUDA_CUBLAS_ERROR_HPP

#if !defined(MULTI_USE_HIP)
#include<cublas_v2.h> // cublasStatus_t
#else
#include<hipblas/hipblas.h> // cublasStatus_t
#endif

#include<string>
#include<system_error> // std::error_category
#include<type_traits> // std::underlying_type

#if !defined(MULTI_USE_HIP)
#define hicup(name) cuda##name
#define hicu(name) cu##name
#define HICU(name) CU##name
#define HICUP(name) CU##name
#else
#define hicup(name) hip##name
#define hicu(name) hip##name
#define HICU(name) HIP##name
#define HICUP(name) HIP##name
#endif

namespace boost::multi::cuda::cublas{

enum class error : typename std::underlying_type<hicu(blasStatus_t)>::type{
	success               = HICUP(BLAS_STATUS_SUCCESS),
	not_initialized       = HICUP(BLAS_STATUS_NOT_INITIALIZED),
	allocation_failed     = HICUP(BLAS_STATUS_ALLOC_FAILED),
	invalid_value         = HICUP(BLAS_STATUS_INVALID_VALUE),
	architecture_mismatch = HICUP(BLAS_STATUS_ARCH_MISMATCH),
	mapping_error         = HICUP(BLAS_STATUS_MAPPING_ERROR),
	execution_failed      = HICUP(BLAS_STATUS_EXECUTION_FAILED),
	internal_error        = HICUP(BLAS_STATUS_INTERNAL_ERROR),
	not_supported         = HICUP(BLAS_STATUS_NOT_SUPPORTED),
//  license_error         = HICUP(BLAS_STATUS_LICENSE_ERROR),  // not supported by hip
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
//  case cublas::error::license_error        : return "CUBLAS_STATUS_LICENSE_ERROR"   ;
	}
	return "cublas status <unknown>";
}

struct error_category : std::error_category{
	char const* name() const noexcept override{return "cublas wrapper";}
	std::string message(int err) const override{return error_string(static_cast<enum cublas::error>(err));}
	static error_category& instance(){static cublas::error_category instance; return instance;}
};

inline std::error_code make_error_code(cublas::error err) noexcept{
	return {int(err), cublas::error_category::instance()};
}

}

namespace std {
	template<> struct is_error_code_enum<::boost::multi::cuda::cublas::error> : true_type{};
}

#undef hicu
#undef hicup
#undef HICU
#undef HICUP

#endif
