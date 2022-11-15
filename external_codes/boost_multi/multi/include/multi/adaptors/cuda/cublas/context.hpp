// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa
#pragma once

#include "../../../config/MARK.hpp"
#include "../../../adaptors/cuda/cublas/call.hpp"

#include "../../../adaptors/blas/traits.hpp"
#include "../../../adaptors/blas/core.hpp"

#include "../../../memory/adaptors/cuda/ptr.hpp"
#include "../../../memory/adaptors/cuda/managed/ptr.hpp"

#include <thrust/system/cuda/memory.h>

#include<mutex>

namespace boost {
namespace multi::cuda::cublas {

class operation {
	cublasOperation_t impl_;

 public:
	explicit operation(char trans) : impl_{[=]{
		switch(trans) {
			case 'N': return CUBLAS_OP_N;
			case 'T': return CUBLAS_OP_T;
			case 'C': return CUBLAS_OP_C;
			default : assert(0);
		}
		return cublasOperation_t{};
	}()} {}
	operator cublasOperation_t() const{return impl_;}
};

class side {
	cublasSideMode_t impl_;

 public:
	explicit side(char trans) : impl_{[=] {
		switch(trans) {
			case 'L': return CUBLAS_SIDE_LEFT;
			case 'R': return CUBLAS_SIDE_RIGHT;
		}
		assert(0); return cublasSideMode_t{};
	}()} {}
	operator cublasSideMode_t() const {return impl_;}
};

class filling {
	cublasFillMode_t impl_;

 public:
	explicit filling(char trans) : impl_{[=] {
		switch(trans) {
			case 'L': return CUBLAS_FILL_MODE_LOWER;
			case 'U': return CUBLAS_FILL_MODE_UPPER;
		}
		assert(0); return cublasFillMode_t{};
	}()} {}
	operator cublasFillMode_t() const {return impl_;}
};

class diagonal {
	cublasDiagType_t impl_;

 public:
	explicit diagonal(char trans) : impl_{[=] {
		switch(trans) {
			case 'N': return CUBLAS_DIAG_NON_UNIT;
			case 'U': return CUBLAS_DIAG_UNIT;
		}
		assert(0); return cublasDiagType_t{};
	}()} {}
	operator cublasDiagType_t() const {return impl_;}
};

using blas::is_z;
using blas::is_d;
using std::is_assignable;
using std::is_assignable_v;
using std::is_convertible_v;

class context : private std::unique_ptr<std::decay_t<decltype(*cublasHandle_t{})>, decltype(&cublasDestroy)> {
	using pimpl_t = std::unique_ptr<std::decay_t<decltype(*cublasHandle_t{})>, decltype(&cublasDestroy)>;
	cudaStream_t stream() const {cudaStream_t streamId; cublas::call<cublasGetStream>(this->get(), &streamId); return streamId;}
	template<auto Function, class... Args>
	void sync_call(Args... args) {
		call<Function>(this->get(), args...);
		this->synchronize();
	}

 public:
	using pimpl_t::get;
	static context& get_instance() {
		thread_local context ctxt;
		return ctxt;
	};
	context() : pimpl_t{[] {cublasHandle_t h; cublasCreate(&h); return h;}(), &cublasDestroy} {}
	using ssize_t = int;
	static int version() {int ret; cublas::call<cublasGetVersion>(nullptr, &ret); return ret;}
	void synchronize() {
		cudaError_t	e = cudaDeviceSynchronize();
		//cudaError_t e = cudaStreamSynchronize(stream());
		if(e != cudaSuccess) {throw std::runtime_error{"cannot synchronize stream in cublas context"};}
	}
	template<class ALPHA, class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type,
		std::enable_if_t<is_d<X>{} and is_d<Y>{}, int> = 0
	//	std::enable_if_t<is_d<X>{} and is_d<Y>{} and is_assignable<Y&, ALPHA{}*X{} + Y{}>{} and is_convertible_v<XP, thrust::cuda::pointer<X>> and is_convertible_v<YP, thrust::cuda::pointer<Y>>, int> = 0
	>
	void axpy(ssize_t n, ALPHA const* alpha, XP x, ssize_t incx, YP y, ssize_t incy) {
		sync_call<cublasDaxpy>(
			n,
			(double const*)alpha,  // TODO(correaa) use static_cast
			(double const*)raw_pointer_cast(x), incx,
			(double*)raw_pointer_cast(y), incy
		);
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename std::pointer_traits<CCP>::element_type,
		std::enable_if_t<
			is_z<AA>{} and is_z<BB>{} and is_z<CC>{} and is_z<ALPHA>{} and is_z<BETA>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{})>{} and
			std::is_convertible_v<AAP, memory::cuda::ptr<AA>> and std::is_convertible_v<BBP, memory::cuda::ptr<BB>> and std::is_convertible_v<CCP, memory::cuda::ptr<CC>>
		,int> =0
	>
	void gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {
		MULTI_MARK_SCOPE("cublasZgemm");
		sync_call<cublasZgemm>(cublas::operation{transA}, cublas::operation{transB}, m, n, k, (cuDoubleComplex const*)alpha, (cuDoubleComplex const*)raw_pointer_cast(aa), lda, (cuDoubleComplex const*)raw_pointer_cast(bb), ldb, (cuDoubleComplex const*)beta, (cuDoubleComplex*)raw_pointer_cast(cc), ldc);
	}
	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename std::pointer_traits<CCP>::element_type,
		std::enable_if_t<
			is_d<AA>{} and is_d<BB>{} and is_d<CC>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{})>{} and
			std::is_convertible_v<AAP, memory::cuda::ptr<AA>> and std::is_convertible_v<BBP, memory::cuda::ptr<BB>> and std::is_convertible_v<CCP, memory::cuda::ptr<CC>>
		,int> =0
	>
	void gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {
		MULTI_MARK_SCOPE("cublasDgemm");
		sync_call<cublasDgemm>(cublas::operation{transA}, cublas::operation{transB}, m, n, k, (double const*)alpha, (double const*)raw_pointer_cast(aa), lda, (double const*)raw_pointer_cast(bb), ldb, (double const*)beta, (double*)raw_pointer_cast(cc), ldc);
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename std::pointer_traits<CCP>::element_type,
		std::enable_if_t<
			is_z<AA>{} and is_z<BB>{} and is_z<CC>{} and is_z<ALPHA>{} and is_z<BETA>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{})>{} and
			std::is_convertible_v<AAP, ::thrust::cuda::universal_pointer<AA>> and std::is_convertible_v<BBP, ::thrust::cuda::universal_pointer<BB>> and std::is_convertible_v<CCP, ::thrust::cuda::universal_pointer<CC>>
		,int> =0
	>
	void gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {
		MULTI_MARK_SCOPE("cublasZgemm");
		sync_call<cublasZgemm>(cublas::operation{transA}, cublas::operation{transB}, m, n, k, (cuDoubleComplex const*)alpha, (cuDoubleComplex const*)raw_pointer_cast(aa), lda, (cuDoubleComplex const*)raw_pointer_cast(bb), ldb, (cuDoubleComplex const*)beta, (cuDoubleComplex*)raw_pointer_cast(cc), ldc);
	}
	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename std::pointer_traits<CCP>::element_type,
		std::enable_if_t<
			is_d<AA>{} and is_d<BB>{} and is_d<CC>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{})>{} and
			std::is_convertible_v<AAP, ::thrust::cuda::universal_pointer<AA>> and std::is_convertible_v<BBP, ::thrust::cuda::universal_pointer<BB>> and std::is_convertible_v<CCP, ::thrust::cuda::universal_pointer<CC>>
		,int> =0
	>
	void gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {
		MULTI_MARK_SCOPE("cublasDgemm");
		sync_call<cublasDgemm>(cublas::operation{transA}, cublas::operation{transB}, m, n, k, (double const*)alpha, (double const*)raw_pointer_cast(aa), lda, (double const*)raw_pointer_cast(bb), ldb, (double const*)beta, (double*)raw_pointer_cast(cc), ldc);
	}

	template<class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BBP, class BB = typename pointer_traits<BBP>::element_type,
		std::enable_if_t<
			is_z<AA>{} and is_z<BB>{} and is_assignable<BB&, decltype(AA{}*BB{}/ALPHA{})>{} and is_assignable<BB&, decltype(ALPHA{}*BB{}/AA{})>{} and 
			is_convertible_v<AAP, memory::cuda::ptr<AA>> and is_convertible_v<BBP, memory::cuda::ptr<BB>>
		,int> =0
	>
	void trsm(char side, char ul, char transA, char diag, ssize_t m, ssize_t n, ALPHA alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb) {
		sync_call<cublasZtrsm>(cublas::side{side}, cublas::filling{ul}, cublas::operation{transA}, cublas::diagonal{diag}, m, n, (cuDoubleComplex const*)&alpha, (cuDoubleComplex const*)raw_pointer_cast(aa), lda, (cuDoubleComplex*)raw_pointer_cast(bb), ldb);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_d<XX>{} and is_d<YY>{} and is_d<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, memory::cuda::ptr<XX>> and is_convertible_v<YYP, memory::cuda::ptr<YY>> and is_convertible_v<RRP, RR*>
		, int> =0
	>
	void dot(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		cublasPointerMode_t mode;
		auto s = cublasGetPointerMode(get(), &mode); assert( s == CUBLAS_STATUS_SUCCESS );
		assert( mode == CUBLAS_POINTER_MODE_HOST );
		sync_call<cublasDdot>(n, raw_pointer_cast(xx), incx, raw_pointer_cast(yy), incy, rr);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_d<XX>{} and is_d<YY>{} and is_d<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, ::thrust::cuda::universal_pointer<XX>> and is_convertible_v<YYP, ::thrust::cuda::universal_pointer<YY>> and is_convertible_v<RRP, RR*>
		, int> =0
	>
	void dot(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		cublasPointerMode_t mode;
		auto s = cublasGetPointerMode(get(), &mode); assert( s == CUBLAS_STATUS_SUCCESS );
		assert( mode == CUBLAS_POINTER_MODE_HOST );
		sync_call<cublasDdot>(n, raw_pointer_cast(xx), incx, raw_pointer_cast(yy), incy, rr);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_z<XX>{} and is_z<YY>{} and is_z<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, memory::cuda::ptr<XX>> and is_convertible_v<YYP, memory::cuda::ptr<YY>> and is_convertible_v<RRP, RR*>
		, int> =0
	>
	void dotc(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		cublasPointerMode_t mode;
		auto s = cublasGetPointerMode(get(), &mode); assert( s == CUBLAS_STATUS_SUCCESS );
		assert( mode == CUBLAS_POINTER_MODE_HOST );
		sync_call<cublasZdotc>(n, (cuDoubleComplex const*)raw_pointer_cast(xx), incx, (cuDoubleComplex const*)raw_pointer_cast(yy), incy, (cuDoubleComplex*)rr);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_z<XX>{} and is_z<YY>{} and is_z<RR>{} and std::is_assignable_v<RR&, decltype(XX{}*YY{})> and
			is_convertible_v<XXP, ::thrust::cuda::universal_pointer<XX>> and is_convertible_v<YYP, ::thrust::cuda::universal_pointer<YY>> and is_convertible_v<RRP, RR*>
		, int> =0
	>
	void dotc(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		cublasPointerMode_t mode;
		auto s = cublasGetPointerMode(get(), &mode); assert( s == CUBLAS_STATUS_SUCCESS );
		assert( mode == CUBLAS_POINTER_MODE_HOST );
		sync_call<cublasZdotc>(n, (cuDoubleComplex const*)raw_pointer_cast(xx), incx, (cuDoubleComplex const*)raw_pointer_cast(yy), incy, (cuDoubleComplex*)rr);
	}
};

}  // end namespace multi::cuda::cublas
}  // end namespace boost

namespace boost::multi::blas {

	template<> struct is_context<boost::multi::cuda::cublas::context > : std::true_type {};
	template<> struct is_context<boost::multi::cuda::cublas::context&> : std::true_type {};

	template<class Ptr, class T = typename std::pointer_traits<Ptr>::element_type, std::enable_if_t<std::is_convertible<Ptr, multi::memory::cuda::ptr<T>>{}, int> =0>
	boost::multi::cuda::cublas::context* default_context_of(Ptr const&) {
		namespace multi = boost::multi;
		return &multi::cuda::cublas::context::get_instance();
	}

	template<class T>
	boost::multi::cuda::cublas::context* default_context_of(boost::multi::memory::cuda::managed::ptr<T> const&) {
		namespace multi = boost::multi;
		return &multi::cuda::cublas::context::get_instance();
	}

	template<class T, class R>
	boost::multi::cuda::cublas::context* default_context_of(::thrust::pointer<T, ::thrust::cuda_cub::tag, R> const&) {
		namespace multi = boost::multi;
		return &multi::cuda::cublas::context::get_instance();
	}

}
