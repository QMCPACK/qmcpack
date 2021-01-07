#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS -include"boost/log/trivial.hpp" -D'MULTI_MARK_SCOPE(MsG)=BOOST_LOG_TRIVIAL(trace)<<MsG' -DBOOST_LOG_DYN_LINK $0 -o $0x `pkg-config --cflags --libs cudart-11.0 cublas-11.0 blas` -lboost_unit_test_framework -lboost_log -lboost_thread -lboost_system -lboost_log_setup -lpthread&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_CUDA_HPP
#define MULTI_ADAPTORS_BLAS_CUDA_HPP

#include "../blas/../../config/MARK.hpp" // MULTI_MARK_SCOPE

#include "../../memory/adaptors/cuda/ptr.hpp"
#include "../../memory/adaptors/cuda/managed/ptr.hpp"
#include "../../memory/adaptors/cuda/managed/allocator.hpp"

#include<cublas_v2.h>

#include "../cublas/error.hpp"

#include<thrust/complex.h>

#define DECLRETURN(ExpR) ->decltype(ExpR){return ExpR;}
#define JUSTRETURN(ExpR)                 {return ExpR;}

#include<complex>

///////////////////

#include<system_error>

#define CUBLAS_CALL(CodE) \
	MULTI_MARK_SCOPE("multi::cublas "#CodE); \
	auto s = static_cast<enum cublas::error>(CodE); \
	cudaDeviceSynchronize(); /*TODO make this more specific to mananged ptr and specific handle*/ \
	if(s != cublas::error::success) throw std::system_error{make_error_code(s), "cannot call cublas function "#CodE };

cublasStatus_t cublasZdot (cublasHandle_t handle, int n,
                           const double2          *x, int incx,
                           const double2          *y, int incy,
                           double2          *result) = delete;

namespace boost{
namespace multi{

namespace cublas{
	using Complex = cuComplex;
	using DoubleComplex = cuDoubleComplex;
	namespace {
		template<class T> struct complex_t;
		template<> struct complex_t<float>{using type = Complex;};
		template<> struct complex_t<double>{using type = DoubleComplex;};
	}
	template<class T> using complex = typename complex_t<T>::type;

// 2.2.7. cublasPointerMode_t https://docs.nvidia.com/cuda/cublas/index.html#cublaspointermode_t
	enum class pointer_mode : std::underlying_type<cublasPointerMode_t>::type{
		host   = CUBLAS_POINTER_MODE_HOST, 
		device = CUBLAS_POINTER_MODE_DEVICE
	};
	template<class T> enum pointer_mode scalar_kind(memory::cuda::ptr<T>){return pointer_mode::device;}
	template<class T> enum pointer_mode scalar_kind(T*){return pointer_mode::host;}
}

using v = void;
using S = float;
using D = double;
using C = cublas::complex<float>;
using Z = cublas::complex<double>;

template<class T = void> struct cublas1{};
template<class T = void> struct cublas2{};
template<class T = void> struct cublas3{};

#define DEFINE_CUBLAS1(UppeR, LowR) \
	template<> struct cublas1<UppeR>{ \
		template<class...As> static auto iamax(As...as){return cublasI##LowR##amax(as...);} \
		/*amin */ \
		template<class...As> static auto asum (As...as){return cublas##UppeR##asum (as...);} \
		/*axpy */ \
		template<class...As> static auto copy (As...as){return cublas##UppeR##copy (as...);} \
		template<class...As> static auto dot  (As...as){return cublas##UppeR##dot  (as...);} \
		template<class...As> static auto dotu (As...as){return cublas##UppeR##dotu (as...);} \
		template<class...As> static auto dotc (As...as){return cublas##UppeR##dotc (as...);} \
		template<class...As> static auto nrm2 (As...as){return cublas##UppeR##nrm2 (as...);} \
		/*rot  */ \
		/*rotg */ \
		/*rotmg*/ \
		template<class...As> static auto scal (As...as){return cublas##UppeR##scal (as...);} \
		/*swap */ \
	}

DEFINE_CUBLAS1(S, s);
DEFINE_CUBLAS1(D, d);

#define DEFINE_CUBLAS1_COMPLEX(UppeR, LowR, ReaLUppeR, ReaLLowR) \
	template<> struct cublas1<UppeR>{ \
		template<class...As> static auto iamax(As...as){return cublasI##LowR##amax(as...);}    \
		/*amin */ \
		template<class...As> static auto asum (As...as){return cublas##ReaLUppeR##LowR##asum (as...);}   \
		/*axpy */ \
		template<class...As> static auto copy (As...as){return cublas##UppeR##copy (as...);}   \
		template<class...As> static auto dot  (As...as){return cublas##UppeR##dotu  (as...);} \
		template<class...As> static auto dotu (As...as){return cublas##UppeR##dotu (as...);}   \
		template<class...As> static auto dotc (As...as){return cublas##UppeR##dotc (as...);}   \
		template<class...As> static auto nrm2 (As...as){return cublas##UppeR##nrm2 (as...);}   \
		/*rot  */ \
		/*rotg */ \
		/*rotmg*/ \
		template<class...As> static auto scal (As...as){return cublas##UppeR##scal (as...);}   \
		/*swap */ \
	}

DEFINE_CUBLAS1_COMPLEX(C, c, S, s);
DEFINE_CUBLAS1_COMPLEX(Z, z, D, d);

template<class T> struct nrm2_result;//{using type = T;};
template<> struct nrm2_result<S>{using type = S;};
template<> struct nrm2_result<D>{using type = D;};
template<> struct nrm2_result<C>{using type = S;};
template<> struct nrm2_result<Z>{using type = D;};

template<> struct cublas1<void>{
// 2.5.1. cublasI<t>amax() https://docs.nvidia.com/cuda/cublas/index.html#cublasi-lt-t-gt-amax
	template<class T> static cublasStatus_t iamax(cublasHandle_t handle, int n, const T* x, int incx, int *result   ){return cublas1<T>::iamax(handle, n, x, incx, result);}
// 2.5.3. cublas<t>asum() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-asum
	template<class T1, class T2> static cublasStatus_t asum (cublasHandle_t handle, int n, T1 const* x, int incx, T2* result     ){return cublas1<T1>::asum(handle, n, x, incx, result);}
// 2.5.5. cublas<t>copy() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-copy
	template<class T> static cublasStatus_t copy (cublasHandle_t handle, int n, const T* x, int incx, T* y, int incy){return cublas1<T>::copy(handle, n, x, incx, y, incy);}
// 2.5.6. cublas<t>dot() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-dot
	template<class T> static auto dot(cublasHandle_t handle, int n, const T* x, int incx, const T* y, int incy, T* result)
	->decltype(cublas1<T>::dot(handle, n, x, incx, y, incy, result)){MULTI_MARK_SCOPE("function dot");
		return cublas1<T>::dot(handle, n, x, incx, y, incy, result);}
	template<class T> static auto dotu(cublasHandle_t handle, int n, const T* x, int incx, const T* y, int incy, T* result)
	->decltype(cublas1<T>::dotu(handle, n, x, incx, y, incy, result)){MULTI_MARK_SCOPE("function dotu");
		return cublas1<T>::dotu(handle, n, x, incx, y, incy, result);}
	template<class T> static auto dotc(cublasHandle_t handle, int n, const T* x, int incx, const T* y, int incy, T* result)
	->decltype(cublas1<T>::dotc(handle, n, x, incx, y, incy, result)){MULTI_MARK_PRETTY_FUNCTION;
		return cublas1<T>::dotc(handle, n, x, incx, y, incy, result);}
// 2.5.7. cublas<t>nrm2() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-nrm2
	template<class T> static auto nrm2(cublasHandle_t handle, int n,
                            const T           *x, int incx, typename nrm2_result<T>::type  *result){return cublas1<T>::nrm2(handle, n, x, incx, result);}
// 2.5.12. cublas<t>scal()	https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-scale
	template<class T> static cublasStatus_t scal(cublasHandle_t handle, int n,
                            const T           *alpha,
                            T           *x, int incx){return cublas1<T>::scal(handle, n, alpha, x, incx);}
};

template<> struct cublas2<void>{
// 2.6.16. cublas<t>trsv() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-trsv
	template<class T> static cublasStatus_t trsv(cublasHandle_t handle, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           int n, const T           *A, int lda,
                           T           *x, int incx){return cublas2<T>::trsv(handle, uplo, trans, diag, n, A, lda, x, incx);}
};

template<> struct cublas2<S>{template<class...A> static auto trsv(A...a){return cublasStrsv(a...);}};
template<> struct cublas2<D>{template<class...A> static auto trsv(A...a){return cublasDtrsv(a...);}};
template<> struct cublas2<C>{template<class...A> static auto trsv(A...a){return cublasCtrsv(a...);}};
template<> struct cublas2<Z>{template<class...A> static auto trsv(A...a){return cublasZtrsv(a...);}};

template<> struct cublas3<S>{
	template<class...As> static auto gemm (As...as){CUBLAS_CALL(cublasSgemm(as...));}
	template<class...As> static auto syrk (As...as){CUBLAS_CALL(cublasSsyrk(as...));}
//	template<class...As> static auto herk (As...as){return CUBLAS_CALL(cublasSherk)(as...);}
	template<class...As> static auto trsm (As...as){CUBLAS_CALL(cublasStrsm(as...));}
};
template<> struct cublas3<D>{
	template<class...As> static auto gemm (As...as){ CUBLAS_CALL(cublasDgemm(as...));}
	template<class...As> static auto syrk (As...as){ CUBLAS_CALL(cublasDsyrk(as...));}
//	template<class...As> static auto herk (As...as){return cublas_call(cublasDherk)(as...);}
	template<class...As> static auto trsm (As...as){ CUBLAS_CALL(cublasDtrsm(as...));}
};
template<> struct cublas3<C>{
	template<class...As> static auto gemm (As...as){ CUBLAS_CALL(cublasCgemm(as...));}
	template<class...As> static auto syrk (As...as){ CUBLAS_CALL(cublasCsyrk(as...));}
	template<class...As> static auto herk (As...as){ CUBLAS_CALL(cublasCherk(as...));}
	template<class...As> static auto trsm (As...as){ CUBLAS_CALL(cublasCtrsm(as...));}
};
template<> struct cublas3<Z>{
	template<class...As> static auto gemm (As...as){ CUBLAS_CALL(cublasZgemm(as...));}
	template<class...As> static auto syrk (As...as){ CUBLAS_CALL(cublasZsyrk(as...));}
	template<class...As> static auto herk (As...as){ CUBLAS_CALL(cublasZherk(as...));}
	template<class...As> static auto trsm (As...as){ CUBLAS_CALL(cublasZtrsm(as...));}
};

template<class T> struct herk_scalar;
template<> struct herk_scalar<C>{using type = S;};
template<> struct herk_scalar<Z>{using type = D;};

template<class T> struct asum_scalar;
template<> struct asum_scalar<C>{using type = S;};
template<> struct asum_scalar<Z>{using type = D;};

template<class T> using herk_scalar_t = typename herk_scalar<T>::type;

template<> struct cublas3<void>{
// 2.7.1. cublas<t>gemm() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
	template<class T> static auto gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           int m, int n, int k,
                           const T           *alpha,
                           const T           *A, int lda,
                           const T           *B, int ldb,
                           const T           *beta,
                           T           *C, int ldc){MULTI_MARK_PRETTY_FUNCTION; return cublas3<T>::gemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);}
// 2.7.6. cublas<t>syrk() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-syrk
	template<class T> static auto syrk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           int n, int k,
                           const T           *alpha,
                           const T           *A, int lda,
                           const T           *beta,
                           T           *C, int ldc){return cublas3<T>::syrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);}
// 2.7.13. cublas<t>herk() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-herk
	template<class T2, class T3> static auto herk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           int n, int k,
                           const herk_scalar_t<T2> *alpha,
                           const T2       *A, int lda,
                           const herk_scalar_t<T2> *beta,
                           T3       *C, int ldc){return cublas3<T2>::herk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);}
// 2.7.10. cublas<t>trsm() https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-trsm
	template<class T> static auto trsm(cublasHandle_t handle,
                           cublasSideMode_t side, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           int m, int n,
                           std::add_const_t<T>           *alpha,
                           std::add_const_t<T>           *A, int lda,
                           T           *B, int ldb){return cublas3<T>::trsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);}
};

namespace cublas{

template<class T, std::enable_if_t<not std::is_integral<T>{}, int> =0> decltype(auto) translate(T t){return t;}
template<class T, std::enable_if_t<not std::is_copy_constructible<std::decay_t<T>>{}, int> =0> T& translate(T& t){return t;}

auto translate(std::complex<float> const * t){return reinterpret_cast<cublas::complex<float>  const*>(t);}
auto translate(std::complex<float>       * t){return reinterpret_cast<cublas::complex<float>       *>(t);}
auto translate(std::complex<double> const* t){return reinterpret_cast<cublas::complex<double> const*>(t);}
auto translate(std::complex<double>      * t){return reinterpret_cast<cublas::complex<double>      *>(t);}

auto translate(thrust::complex<double> const* t){return reinterpret_cast<cublas::complex<double> const*>(t);}
auto translate(thrust::complex<double>      * t){return reinterpret_cast<cublas::complex<double>      *>(t);}

template<class T> auto translate(memory::cuda::ptr<T>          p) DECLRETURN(translate(raw_pointer_cast(p)))
template<class T> auto translate(memory::cuda::managed::ptr<T> p) DECLRETURN(translate(raw_pointer_cast(p)))

//auto translate(context& c){return c;}

template<class T, std::enable_if_t<std::is_integral<T>{},int> = 0> 
auto translate(T n){
	assert(n <= +static_cast<T>(std::numeric_limits<int>::max()));
	assert(n >  -static_cast<T>(std::numeric_limits<int>::max()));
	return static_cast<T>(n);
}

auto translate(char O)->cublasOperation_t{
	switch(O){case 'N': return CUBLAS_OP_N; case 'T': return CUBLAS_OP_T; case 'C': return CUBLAS_OP_C;} assert(0); 
	return CUBLAS_OP_N;
}

struct context : std::unique_ptr<std::decay_t<decltype(*cublasHandle_t{})>, decltype(&cublasDestroy)>{
	context()  : std::unique_ptr<std::decay_t<decltype(*cublasHandle_t{})>, decltype(&cublasDestroy)>(
		[]{MULTI_MARK_SCOPE("multi::cublas::create context"); cublasHandle_t h; cublasCreate(&h); return h;}(), &cublasDestroy
	){}
	int version() const{
		int ret; cublasGetVersion(get(), &ret); return ret;
	}
	context(context&& other) noexcept = default;
	~context() noexcept = default;
// 2.4.7. cublasGetPointerMode()
	auto get_pointer_mode() const{
		cublasPointerMode_t ret; cublasGetPointerMode(get(), &ret);
		return static_cast<enum pointer_mode>(ret);
	}
// 2.4.8. cublasSetPointerMode() https://docs.nvidia.com/cuda/cublas/index.html#cublassetpointermode
	context& set_pointer_mode(enum pointer_mode m){
		cublasSetPointerMode(get(), static_cast<cublasPointerMode_t>(m)); return *this;
	}
	//set_stream https://docs.nvidia.com/cuda/cublas/index.html#cublassetstream
	//get_stream https://docs.nvidia.com/cuda/cublas/index.html#cublasgetstream
	//get_pointer_mode https://docs.nvidia.com/cuda/cublas/index.html#cublasgetpointermode
	//set_pointer_mode https://docs.nvidia.com/cuda/cublas/index.html#cublasgetpointermode
	template<class...As> auto iamax(As...as) const DECLRETURN(cublas1<>::iamax(get(), translate(as)...))
	template<class...As> auto asum (As...as) const DECLRETURN(cublas1<>::asum (get(), translate(as)...))
	template<class...As> auto scal (As...as) const DECLRETURN(cublas1<>::scal (get(), translate(as)...))
	template<class...As> auto dot  (As...as) const DECLRETURN(cublas1<>::dot  (get(), translate(as)...))
	template<class...As> auto dotu (As...as) const DECLRETURN(cublas1<>::dotu (get(), translate(as)...))
	template<class...As> auto dotc (As...as) const DECLRETURN(cublas1<>::dotc (get(), translate(as)...))
	template<class S, class Ptr, class T>
	auto nrm2(S n, Ptr p, S incx, memory::cuda::ptr<T> result) // no const because the method is not thread safe
	->decltype(cublas1<>::nrm2 (get(), translate(n), translate(p), translate(incx), translate(result))){set_pointer_mode(pointer_mode::device);
		auto r=cublas1<>::nrm2 (get(), translate(n), translate(p), translate(incx), translate(result)); set_pointer_mode(pointer_mode::host);
		return r;
	}
	template<class S, class Ptr, class T>
	auto nrm2(S n, Ptr p, S incx, T* result) const{
		return cublas1<>::nrm2 (get(), translate(n), translate(p), translate(incx), translate(result));
	}
	template<class...As> auto copy (As...as) const DECLRETURN(cublas1<>::copy (get(), translate(as)...))
	template<class...As> auto trsv (As...as) const{return cublas2<>::trsv(get(), translate(as)...);}

	template<typename... As> auto gemm(As... as) DECLRETURN(cublas3<>::gemm(get(), translate(as)...))

	template<class...As> auto syrk (As...as) const{return cublas3<>::syrk(get(), translate(as)...);}
	template<class...As> auto herk (As...as) const{return cublas3<>::herk(get(), translate(as)...);}
	template<class...As> auto trsm (As...as) const{return cublas3<>::trsm(get(), translate(as)...);}
};
}

}}

namespace boost{
namespace multi{

namespace blas{

template<class T> boost::multi::cublas::context default_context_of(memory::cuda::ptr<T>){return {};}
template<class T> boost::multi::cublas::context default_context_of(memory::cuda::managed::ptr<T>){return {};}

}

namespace memory{namespace cuda{
	using boost::multi::blas::default_context_of; // to please nvcc 'default_context_of' should be declared prior to the call site or in namespace 'boost::multi::memory::cuda'
}}

}}

namespace boost{
namespace multi{

namespace memory{
namespace cuda{

template<class... As>
auto iamax(As... as)
->decltype(cublas::context{}.iamax(as..., std::declval<int*>()), int()){
	int r; cublas::context{}.iamax(as..., &r); return r-1;}

template<class ComplexTconst, typename S>//, typename T = typename std::decay_t<ComplexTconst>::value_type>
auto asum(S n, cuda::ptr<ComplexTconst> x, S incx){
	decltype(std::abs(ComplexTconst{})) r;
	cublas::context{}.asum(n, raw_pointer_cast(x), incx, &r);
	return r;
}

template<class...As> auto copy(As... as) DECLRETURN(cublas::context{}.copy(as...))
template<class...As> auto scal(As... as) DECLRETURN(cublas::context{}.scal(as...))
//template<class...As> auto dot (As... as) DECLRETURN(cublas::context{}.dot (as...))
template<class...As> auto dotu(As... as) DECLRETURN(cublas::context{}.dotu(as...))
template<class...As> auto dotc(As... as) DECLRETURN(cublas::context{}.dotc(as...))
template<class...As> auto nrm2(As... as) DECLRETURN(cublas::context{}.nrm2(as...))

template<class S, class Tconst, class T>
auto trsv(char ul, char transA, char a_diag, S n, memory::cuda::ptr<Tconst> A, S lda, memory::cuda::ptr<T> X, S ldc){
	cublasFillMode_t uplo = [ul](){
		switch(ul){
			case 'U': return CUBLAS_FILL_MODE_UPPER;
			case 'L': return CUBLAS_FILL_MODE_LOWER;
		} assert(0); return CUBLAS_FILL_MODE_UPPER;
	}();
	cublasOperation_t cutransA = [transA](){
		switch(transA){
			case 'N': return CUBLAS_OP_N;
			case 'T': return CUBLAS_OP_T;
			case 'C': return CUBLAS_OP_C;
		} assert(0); return CUBLAS_OP_N;
	}();
	auto cudiag = a_diag=='N'?CUBLAS_DIAG_NON_UNIT:CUBLAS_DIAG_UNIT;
	return cublas::context{}.trsv(uplo, cutransA, cudiag, n, A, lda, X, ldc);
}

template<class... As>
auto gemm(As... as)
->decltype(cublas::context{}.gemm(as...)){
	return cublas::context{}.gemm(as...);}

template<class Tconst, class T, class UL, class C, class S, class Real>
void syrk(UL ul, C transA, S n, S k, Real alpha, multi::memory::cuda::ptr<Tconst> A, S lda, Real beta, multi::memory::cuda::ptr<T> CC, S ldc){
	cublasFillMode_t uplo = [ul](){
		switch(ul){
			case 'U': return CUBLAS_FILL_MODE_UPPER;
			case 'L': return CUBLAS_FILL_MODE_LOWER;
		} assert(0); return CUBLAS_FILL_MODE_UPPER;
	}();
	cublasOperation_t cutransA = [transA](){
		switch(transA){
			case 'N': return CUBLAS_OP_N;
			case 'T': return CUBLAS_OP_T;
			case 'C': return CUBLAS_OP_C;
		} assert(0); return CUBLAS_OP_N;
	}();
	return cublas::context{}.syrk(uplo, cutransA, n, k, &alpha, static_cast<T const*>(A), lda, &beta, static_cast<T*>(CC), ldc);
}

template<class Tconst, class T, class UL, class C, class S, class Real>
auto herk(UL ul, C transA, S n, S k, Real alpha, memory::cuda::ptr<Tconst> A, S lda, Real beta, memory::cuda::ptr<T> CC, S ldc){
	cublasFillMode_t uplo = [ul](){
		switch(ul){
			case 'U': return CUBLAS_FILL_MODE_UPPER;
			case 'L': return CUBLAS_FILL_MODE_LOWER;
		} assert(0); return CUBLAS_FILL_MODE_UPPER;
	}();
	cublasOperation_t cutransA = [transA](){
		switch(transA){
			case 'N': return CUBLAS_OP_N;
			case 'T': return CUBLAS_OP_T;
			case 'C': return CUBLAS_OP_C;
		} assert(0); return CUBLAS_OP_N;
	}();
	return cublas::context{}.herk(uplo, cutransA, n, k, &alpha, raw_pointer_cast(A), lda, &beta, raw_pointer_cast(CC), ldc);
}

template<class Side, class Fill, class Trans, class Diag, typename Size, class Tconst, class T/*, class Alpha*/>
auto trsm(Side /*cublasSideMode_t*/ side, /*cublasFillMode_t*/ Fill uplo, /*cublasOperation_t*/ Trans trans, /*cublasDiagType_t*/ Diag diag,
                           Size m, Size n, T alpha, cuda::ptr<Tconst> A, Size lda, cuda::ptr<T> B, Size ldb)
->decltype(cublas::context{}.trsm(
		side=='L'?CUBLAS_SIDE_LEFT:CUBLAS_SIDE_RIGHT, uplo=='L'?CUBLAS_FILL_MODE_LOWER:CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, diag=='N'?CUBLAS_DIAG_NON_UNIT:CUBLAS_DIAG_UNIT, m, n, &alpha, raw_pointer_cast(A), lda, raw_pointer_cast(B), ldb))
{
	cublasOperation_t trans_cu = [&]{
		switch(trans){
			case 'N': return CUBLAS_OP_N;
			case 'T': return CUBLAS_OP_T;
			case 'C': return CUBLAS_OP_C;
		} __builtin_unreachable();
	}();
//	T alpha_{alpha};
	return cublas::context{}.trsm(
		side=='L'?CUBLAS_SIDE_LEFT:CUBLAS_SIDE_RIGHT, uplo=='L'?CUBLAS_FILL_MODE_LOWER:CUBLAS_FILL_MODE_UPPER, trans_cu, diag=='N'?CUBLAS_DIAG_NON_UNIT:CUBLAS_DIAG_UNIT, m, n, &alpha, raw_pointer_cast(A), lda, raw_pointer_cast(B), ldb);
}

}}}}

namespace boost{namespace multi{namespace memory{namespace cuda{namespace managed{

using cuda::iamax;
using cuda::asum;
using cuda::copy;
using cuda::scal;
//using cuda::dot;
using cuda::dotu;
using cuda::dotc;
using cuda::nrm2;

template<class S, class Tconst, class T>
auto trsv(char ul, char transA, char a_diag, S n, multi::memory::cuda::managed::ptr<Tconst> A, S lda, cuda::managed::ptr<T> X, S ldc){
	cuda::trsv(ul, transA, a_diag, n, cuda::ptr<Tconst>(A), lda, cuda::ptr<T>(X), ldc);
}

using cuda::gemm;
using cuda::syrk;
using cuda::herk;

template<class Side, class Fill, class Trans, class Diag, typename Size, class Tconst, class T>
auto trsm(Side /*cublasSideMode_t*/ side, /*cublasFillMode_t*/ Fill uplo, /*cublasOperation_t*/ Trans trans, /*cublasDiagType_t*/ Diag diag,
                           Size m, Size n, T alpha, cuda::managed::ptr<Tconst> A, Size lda, cuda::managed::ptr<T> B, Size ldb){
	return trsm(side, uplo, trans, diag, m, n, alpha, cuda::ptr<Tconst>(A), lda, cuda::ptr<T>(B), ldb);
}

}}}}}

///////////////////////////////////////////////////////////////////////////////

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_CUDA

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include "../../adaptors/cuda.hpp"
#include "../../adaptors/blas.hpp"
#include "../../adaptors/blas/cuda.hpp"

#include<cassert>

namespace multi = boost::multi;

#if 0
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_version){
	multi::cublas::context c;
	BOOST_REQUIRE( c.version() >= 10100 );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_iamax){
	using complex = std::complex<double>;
	complex const I{0,1};
	{
		multi::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::managed::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
}
#endif

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_dot){
	using complex = std::complex<double>;
	complex const I{0,1};
	multi::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
	multi::array<complex, 1> const B = {2. + 3.*I, 4., 5. + 6.*I, 7.};
	namespace blas = multi::blas;
	{
		multi::cuda::array<complex, 1> const A_gpu = A, B_gpu = B;
		using blas::dot;
		BOOST_REQUIRE( dot(blas::C(A_gpu), B_gpu) == dot(blas::C(A), B) );
	}
	{
		multi::cuda::managed::array<complex, 1> const A_mng = A, B_mng = B;
		using blas::dot;
		BOOST_REQUIRE( dot(blas::C(A_mng), A_mng) == dot(blas::C(A), A) );
	}
}


#endif
#endif

