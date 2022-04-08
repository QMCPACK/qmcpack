#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_LAPACK_CUDA $0.cpp -o $0x `pkg-config --libs blas` -lcudart -lcublas -lcusolver&&$0x&&rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_ADAPTORS_LAPACK_CUDA_HPP
#define MULTI_ADAPTORS_LAPACK_CUDA_HPP

#include "../../memory/adaptors/cuda/ptr.hpp"
#include "../../memory/adaptors/cuda/managed/ptr.hpp"
#include "../../memory/adaptors/cuda/managed/allocator.hpp"

#include "../../adaptors/cuda.hpp"

//#include<cublas_v2.h>
#include <cusolverDn.h>

//#include<iostream> // debug

#include <boost/log/trivial.hpp>

#include<complex>
#include<memory>

#include "../blas/filling.hpp"

namespace boost{
namespace multi{

namespace cusolver{

enum class status : typename std::underlying_type<cusolverStatus_t>::type{
	success                   = CUSOLVER_STATUS_SUCCESS, // "The operation completed successfully."
	not_initialized           = CUSOLVER_STATUS_NOT_INITIALIZED, // "The cuSolver library was not initialized. This is usually caused by the lack of a prior call, an error in the CUDA Runtime API called by the cuSolver routine, or an error in the hardware setup. To correct: call cusolverCreate() prior to the function call; and check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed."
	allocation_failed         = CUSOLVER_STATUS_ALLOC_FAILED, // "Resource allocation failed inside the cuSolver library. This is usually caused by a cudaMalloc() failure. To correct: prior to the function call, deallocate previously allocated memory as much as possible."
	invalid_value             = CUSOLVER_STATUS_INVALID_VALUE, // "An unsupported value or parameter was passed to the function (a negative vector size, for example). To correct: ensure that all the parameters being passed have valid values."
	architecture_mismatch     = CUSOLVER_STATUS_ARCH_MISMATCH, // "The function requires a feature absent from the device architecture; usually caused by the lack of support for atomic operations or double precision. To correct: compile and run the application on a device with compute capability 2.0 or above."
	execution_failed          = CUSOLVER_STATUS_EXECUTION_FAILED, // "The GPU program failed to execute. This is often caused by a launch failure of the kernel on the GPU, which can be caused by multiple reasons. To correct: check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed."
	internal_error            = CUSOLVER_STATUS_INTERNAL_ERROR, // "An internal cuSolver operation failed. This error is usually caused by a cudaMemcpyAsync() failure. To correct: check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed. Also, check that the memory passed as a parameter to the routine is not being deallocated prior to the routine’s completion."
	matrix_type_not_supported = CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED // "The matrix type is not supported by this function. This is usually caused by passing an invalid matrix descriptor to the function. To correct: check that the fields in descrA were set correctly."
};

std::string inline status_string(enum status s){ //https://stackoverflow.com/questions/13041399/equivalent-of-cudageterrorstring-for-cublas
	switch(s){
		case status::success                   : return "The operation completed successfully.";
		case status::not_initialized           : return "The cuSolver library was not initialized. This is usually caused by the lack of a prior call, an error in the CUDA Runtime API called by the cuSolver routine, or an error in the hardware setup. To correct: call cusolverCreate() prior to the function call; and check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed.";
		case status::allocation_failed         : return "Resource allocation failed inside the cuSolver library. This is usually caused by a cudaMalloc() failure. To correct: prior to the function call, deallocate previously allocated memory as much as possible.";
		case status::invalid_value             : return "An unsupported value or parameter was passed to the function (a negative vector size, for example). To correct: ensure that all the parameters being passed have valid values.";
		case status::architecture_mismatch     : return "The function requires a feature absent from the device architecture; usually caused by the lack of support for atomic operations or double precision. To correct: compile and run the application on a device with compute capability 2.0 or above.";
		case status::execution_failed          : return "The GPU program failed to execute. This is often caused by a launch failure of the kernel on the GPU, which can be caused by multiple reasons. To correct: check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed.";
		case status::internal_error            : return "An internal cuSolver operation failed. This error is usually caused by a cudaMemcpyAsync() failure. To correct: check that the hardware, an appropriate version of the driver, and the cuSolver library are correctly installed. Also, check that the memory passed as a parameter to the routine is not being deallocated prior to the routine’s completion.";
		case status::matrix_type_not_supported : return "The matrix type is not supported by this function. This is usually caused by passing an invalid matrix descriptor to the function. To correct: check that the fields in descrA were set correctly.";
	}
	return "cublas status <unknown>";
}
struct error_category : std::error_category{
	char const* name() const noexcept override{return "cusolver wrapper";}
	std::string message(int err) const override{return cusolver::status_string(static_cast<enum cusolver::status>(err));}
	static error_category& instance(){static cusolver::error_category instance; return instance;}
};
inline std::error_code make_error_code(cusolver::status s) noexcept{
	return std::error_code(int(s), cusolver::error_category::instance());
}

struct version_t{
	int major = -1, minor =-1, patch=-1;
	friend std::ostream& operator<<(std::ostream& os, version_t const& self){
		return os<< self.major <<'.'<< self.minor <<'.'<<self.patch <<'\n';
	}
};

auto version(){
	version_t ret;
	cusolverGetProperty(MAJOR_VERSION, &ret.major);
	cusolverGetProperty(MINOR_VERSION, &ret.minor);
	cusolverGetProperty(PATCH_LEVEL, &ret.patch);
	return ret;
}

namespace dense{

using blas::filling;

struct context : std::unique_ptr<std::decay_t<decltype(*cusolverDnHandle_t{})>, decltype(&cusolverDnDestroy)>{
	context()  : std::unique_ptr<std::decay_t<decltype(*cusolverDnHandle_t{})>, decltype(&cusolverDnDestroy)>{
		[]{
			cusolverDnHandle_t h; 
			auto s=cusolverDnCreate(&h); assert(CUSOLVER_STATUS_SUCCESS==s and h);
			return h;
		}(), &cusolverDnDestroy
	}{}
	template<class A> int potrf_buffer_size(filling uplo, A const& a);
};

template<typename T>
struct cusolverDn;//{
//	template<class... ArgsA3, class... ArgsB2>
//	static auto potrf_bufferSize(ArgsA3... argsa3, T* ptr, ArgsB2... argsb2);
//};

template<>
struct cusolverDn<float>{
	template<class... A3, class... B2>
	static auto potrf_bufferSize(A3... a3, double* ptr, B2... b2)
	->decltype(cusolverDnSpotrf_bufferSize(a3..., ptr, b2...)){
		return cusolverDnSpotrf_bufferSize(a3..., ptr, b2...);}
	template<class... A> static auto syev(A... a)
	->decltype(cusolverDnSsyevd_bufferSize(a...)){
		return cusolverDnSsyevd_bufferSize(a...);}
	template<class... A> static auto syev(A... a)
	->decltype(cusolverDnSsyevd(a...)){
		return cusolverDnSsyevd(a...);}
};

template<>
struct cusolverDn<double>{
	template<class... A3, class... B2>
	static auto potrf_bufferSize(A3... a3, double* ptr, B2... b2)
	->decltype(cusolverDnDpotrf_bufferSize(a3..., ptr, b2...)){
		return cusolverDnDpotrf_bufferSize(a3..., ptr, b2...);}
	template<class... A> static auto syevd_bufferSize(A... a)
//	->decltype(cusolverDnDsyevd_bufferSize(a...)){
	{	return cusolverDnDsyevd_bufferSize(a...);}
	template<class... A> static auto syevd(A... a)
	->decltype(cusolverDnDsyevd(a...))
	{	return cusolverDnDsyevd(a...);}
};

template<>
struct cusolverDn<std::complex<float>>{
	template<class... A3, class... B2>
	static auto potrf_bufferSize(A3... a3, std::complex<double>* ptr, B2... b2)
	->decltype(cusolverDnCpotrf_bufferSize(a3..., reinterpret_cast<cuComplex*>(ptr), b2...)){
		return cusolverDnCpotrf_bufferSize(a3..., reinterpret_cast<cuComplex*>(ptr), b2...);}
};

template<>
struct cusolverDn<std::complex<double>>{
	static auto translate(std::complex<double>* p){return reinterpret_cast<cuDoubleComplex*>(p);}
	template<class T> static auto translate(T t){return t;}
	template<class... A> static auto potrf_bufferSize(A... a)
	->decltype(cusolverDnZpotrf_bufferSize(translate(a)...)){
		return cusolverDnZpotrf_bufferSize(translate(a)...);}
	template<class... A> static auto potrf(A... a)
	->decltype(cusolverDnZpotrf(translate(a)...)){
		return cusolverDnZpotrf(translate(a)...);}
};

}
}

namespace memory{
namespace cuda{

template<class UL, class S, class PtrT, typename T = typename std::pointer_traits<PtrT>::element_type>
void potrf(UL ul, S n, PtrT A, S incx, int& info){
	boost::multi::cusolver::dense::context ctx; //BOOST_LOG_TRIVIAL(trace)<<"cuda::potrf called on size/stride "<< n <<' '<< incx <<'\n';
	int lwork = -1;
	{
		auto s = cusolver::dense::cusolverDn<T>::potrf_bufferSize(ctx.get(), ul=='U'?CUBLAS_FILL_MODE_UPPER:CUBLAS_FILL_MODE_LOWER, n, raw_pointer_cast(A), incx, &lwork);
		assert(s == CUSOLVER_STATUS_SUCCESS); assert(lwork >= 0);
	}
	multi::cuda::array<T, 1> work(lwork);
	multi::cuda::static_array<int, 0> devInfo;
	auto s = cusolver::dense::cusolverDn<T>::potrf(ctx.get(), ul=='U'?CUBLAS_FILL_MODE_UPPER:CUBLAS_FILL_MODE_LOWER, n, raw_pointer_cast(A), incx, raw_pointer_cast(base(work)), lwork, raw_pointer_cast(base(devInfo)) );
	assert(s == CUSOLVER_STATUS_SUCCESS);
	cudaDeviceSynchronize();
	info = devInfo();
}

// https://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-syevd
template<class S, class PtrT, typename T = typename std::pointer_traits<PtrT>::element_type>
void syev(char jobz, char uplo, S n, PtrT a, S lda, PtrT w, PtrT /*work*/, S /*lwork*/, int& info){
	boost::multi::cusolver::dense::context ctx;
	int lwork_needed = -1;
	{
		auto s = cusolver::dense::cusolverDn<T>::syevd_bufferSize(
			ctx.get(), jobz=='V'?CUSOLVER_EIG_MODE_VECTOR:CUSOLVER_EIG_MODE_NOVECTOR, 
			uplo=='U'?CUBLAS_FILL_MODE_UPPER:CUBLAS_FILL_MODE_LOWER, 
			n,
			raw_pointer_cast(a),
			lda,
			raw_pointer_cast(w),
			&lwork_needed
		);
		assert(s == CUSOLVER_STATUS_SUCCESS); assert(lwork_needed >= 0);
	}
	multi::cuda::array<T, 1> tmp_work(lwork_needed); // buffers needs no-managed memory!
	multi::cuda::static_array<int, 0> devInfo;
	auto s = cusolver::dense::cusolverDn<T>::syevd(
		ctx.get(), jobz=='V'?CUSOLVER_EIG_MODE_VECTOR:CUSOLVER_EIG_MODE_NOVECTOR, 
		uplo=='U'?CUBLAS_FILL_MODE_UPPER:CUBLAS_FILL_MODE_LOWER, 
		n,
		raw_pointer_cast(a),
		lda,
		raw_pointer_cast(w),
		raw_pointer_cast(tmp_work.data_elements()),
		tmp_work.size(),
		raw_pointer_cast(base(devInfo))
	);
	if( s != CUSOLVER_STATUS_SUCCESS ) throw std::system_error{cusolver::make_error_code(static_cast<cusolver::status>(s)), "cannot call cusolver function "};
//	cudaDeviceSynchronize();
	info = devInfo();
}

namespace managed{
	template<class UL, class S, class PtrT, typename T = typename std::pointer_traits<PtrT>::element_type>
	auto potrf(UL ul, S n, PtrT A, S incx, int& info)
	->decltype(cuda::potrf(ul, n, cuda::ptr<T>(A), incx, info)){
		return cuda::potrf(ul, n, cuda::ptr<T>(A), incx, info);}

	template<class S, class PtrT, class P2, typename T = typename std::pointer_traits<PtrT>::element_type>
	auto syev(char jobz, char uplo, S n, PtrT a, S lda, PtrT w, P2 work, S lwork, int& info)
	->decltype(cuda::syev(jobz, uplo, n, cuda::ptr<T>(a), lda, cuda::ptr<T>(w), cuda::ptr<T>(work), lwork, info)){
		return cuda::syev(jobz, uplo, n, cuda::ptr<T>(a), lda, cuda::ptr<T>(w), cuda::ptr<T>(work), lwork, info);}

}

}
}

}}

namespace std{template<> struct is_error_code_enum<::boost::multi::cusolver::status> : true_type{};}

///////////////////////////////////////////////////////////////////////////////

#if _TEST_MULTI_ADAPTORS_LAPACK_CUDA

#include "../../array.hpp"
#include "../../utility.hpp"
#include<cassert>

namespace multi = boost::multi;

int main(){
	std::cout << "cusolver version " << multi::cusolver::version() << std::endl;
	multi::cusolver::dense::context c;
//	multi::cublas_context c;
//	assert( c.version() >= 10100 );
}

#endif
#endif

