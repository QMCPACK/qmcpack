// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#pragma once

// #include <multi/config/MARK.hpp>
#include <boost/multi/adaptors/cuda/cublas/call.hpp>

#include <boost/multi/adaptors/blas/traits.hpp>
#include <boost/multi/adaptors/blas/core.hpp>

#if !defined(MULTI_USE_HIP)
#include <thrust/system/cuda/memory.h>  // for thrust::cuda::pointer
#else
#include <thrust/system/hip/memory.h>  // for thrust::cuda::pointer
#include <hipblas/hipblas.h>
#endif

#if not defined(MULTI_USE_HIP)
#define hicup(name) cuda##name
#define hicu(name) cu##name
#define HICU(name) CU##name
#else
#define hicup(name) hip##name
#define hicu(name) hip##name
#define HICU(name) HIP##name
#endif

namespace boost {
namespace multi::cuda::cublas {

class operation {
	hicu(blasOperation_t) impl_;

 public:
	explicit operation(char trans) : impl_{[=]{
		switch(trans) {
			case 'N': return HICU(BLAS_OP_N);
			case 'T': return HICU(BLAS_OP_T);
			case 'C': return HICU(BLAS_OP_C);
			default : assert(0);
		}
		return hicu(blasOperation_t){};
	}()} {}
	operator hicu(blasOperation_t)() const{return impl_;}
};

class side {
	hicu(blasSideMode_t) impl_;

 public:
	explicit side(char trans) : impl_{[=] {
		switch(trans) {
			case 'L': return HICU(BLAS_SIDE_LEFT);
			case 'R': return HICU(BLAS_SIDE_RIGHT);
			default: assert(0);

		}
		assert(0); return hicu(blasSideMode_t){};
	}()} {}
	operator hicu(blasSideMode_t)() const {return impl_;}
};

class filling {
	hicu(blasFillMode_t) impl_;

 public:
	explicit filling(char trans) : impl_{[=] {
		switch(trans) {
			case 'L': return HICU(BLAS_FILL_MODE_LOWER);
			case 'U': return HICU(BLAS_FILL_MODE_UPPER);
		}
		assert(0); return hicu(blasFillMode_t){};
	}()} {}
	operator hicu(blasFillMode_t)() const {return impl_;}
};

class diagonal {
	hicu(blasDiagType_t) impl_;

 public:
	explicit diagonal(char trans) : impl_{[=] {
		switch(trans) {
			case 'N': return HICU(BLAS_DIAG_NON_UNIT);
			case 'U': return HICU(BLAS_DIAG_UNIT);
		}
		assert(0); return hicu(blasDiagType_t){};
	}()} {}
	operator hicu(blasDiagType_t)() const {return impl_;}
};

using blas::is_s;
using blas::is_d;
using blas::is_c;
using blas::is_z;

using std::is_assignable;
using std::is_assignable_v;
using std::is_convertible_v;

// enum class type {S, D, C, Z};

// template<class T>
// constexpr auto type_of(T const& = {}) -> cublas::type {
//  static_assert(is_s<T>{} || is_d<T>{} || is_c<T>{} || is_z<T>{});
//       if(is_s<T>{}) {return type::S;}
//  else if(is_d<T>{}) {return type::D;}
//  else if(is_c<T>{}) {return type::C;}
//  else if(is_z<T>{}) {return type::Z;}
// }

#if defined(__HIP_PLATFORM_NVIDIA__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
	using Complex = hipblasComplex;
	using DoubleComplex = hipblasDoubleComplex;
#else  //  __CUDA__  __NVCC__  or clang cuda
	using Complex = cuComplex;
	using DoubleComplex = cuDoubleComplex;
#endif

template<class T, std::enable_if_t<is_s<T>{}, int> =0> constexpr auto data_cast(T      * p) {return reinterpret_cast<float        *>(p);}
template<class T, std::enable_if_t<is_d<T>{}, int> =0> constexpr auto data_cast(T      * p) {return reinterpret_cast<double       *>(p);}
template<class T, std::enable_if_t<is_c<T>{}, int> =0> constexpr auto data_cast(T      * p) {return reinterpret_cast<Complex      *>(p);}
template<class T, std::enable_if_t<is_z<T>{}, int> =0> constexpr auto data_cast(T      * p) {return reinterpret_cast<DoubleComplex*>(p);}

template<class T, std::enable_if_t<is_s<T>{}, int> =0> constexpr auto data_cast(T const* p) {return reinterpret_cast<float        const*>(p);}
template<class T, std::enable_if_t<is_d<T>{}, int> =0> constexpr auto data_cast(T const* p) {return reinterpret_cast<double       const*>(p);}
template<class T, std::enable_if_t<is_c<T>{}, int> =0> constexpr auto data_cast(T const* p) {return reinterpret_cast<Complex      const*>(p);}
template<class T, std::enable_if_t<is_z<T>{}, int> =0> constexpr auto data_cast(T const* p) {return reinterpret_cast<DoubleComplex const*>(p);}

class context : private std::unique_ptr<typename std::pointer_traits<hicu(blasHandle_t)>::element_type, decltype(&hicu(blasDestroy))> {
	using pimpl_t = std::unique_ptr<typename std::pointer_traits<hicu(blasHandle_t)>::element_type, decltype(&hicu(blasDestroy))>;
	hicup(Stream_t) stream() const {hicup(Stream_t) streamId; cuda::cublas::call<hicu(blasGetStream)>(this->get(), &streamId); return streamId;}
	template<auto Function, class... Args>
	void sync_call(Args... args) const {
		call<Function>(const_cast<context*>(this)->get(), args...);
		this->synchronize();
	}
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
	context() : pimpl_t{[] {hicu(blasHandle_t) h; hicu(blasCreate)(&h); return h;}(), &hicu(blasDestroy)} {}
	using ssize_t = int;
	// static int version() {int ret; cuda::cublas::call<hicu(blasGetVersion)>(nullptr, &ret); return ret;} // no hipblasGetVersion available
	void synchronize() const {
		// cudaError_t e = cudaDeviceSynchronize();
		auto s = stream();
		if(s != 0) {throw std::logic_error("CUBLAS stream expected to be zero");}
		hicup(Error_t) e = hicup(StreamSynchronize)(s);
		if(e != hicup(Success)) {throw std::runtime_error{"cannot synchronize stream in cublas context"};}
	}

	template<
		class XP, class X = typename std::pointer_traits<XP>::element_type,
		class YP, class Y = typename std::pointer_traits<YP>::element_type,
		class = decltype(std::swap(std::declval<X&>(), std::declval<Y&>())),
		std::enable_if_t<std::is_convertible_v<XP, ::thrust::hicup()::pointer<X>>, int> = 0
	>
	void swap(ssize_t n, XP x, ssize_t incx, YP y, ssize_t incy) const {
		if(is_s<X>{}) {sync_call<hicu(blasSswap)>(n, (float        *)raw_pointer_cast(x), incx, (float        *)raw_pointer_cast(y), incy);}
		if(is_d<X>{}) {sync_call<hicu(blasDswap)>(n, (double       *)raw_pointer_cast(x), incx, (double       *)raw_pointer_cast(y), incy);}
		if(is_c<X>{}) {sync_call<hicu(blasCswap)>(n, (Complex      *)raw_pointer_cast(x), incx, (Complex      *)raw_pointer_cast(y), incy);}
		if(is_z<X>{}) {sync_call<hicu(blasZswap)>(n, (DoubleComplex*)raw_pointer_cast(x), incx, (DoubleComplex*)raw_pointer_cast(y), incy);}
	}

	template<
		class XP, class X = typename std::pointer_traits<XP>::element_type,
		class YP, class Y = typename std::pointer_traits<YP>::element_type,
		class = decltype(std::declval<Y&>() = std::declval<X&>()),
		std::enable_if_t<std::is_convertible_v<XP, ::thrust::hicup()::pointer<X>>, int> = 0
	>
	void copy(ssize_t n, XP x, ssize_t incx, YP y, ssize_t incy) const {
		if(is_s<X>{}) {sync_call<hicu(blasScopy)>(n, (float         const*)raw_pointer_cast(x), incx, (float        *)raw_pointer_cast(y), incy);}
		if(is_d<X>{}) {sync_call<hicu(blasDcopy)>(n, (double        const*)raw_pointer_cast(x), incx, (double       *)raw_pointer_cast(y), incy);}
		if(is_c<X>{}) {sync_call<hicu(blasCcopy)>(n, (Complex       const*)raw_pointer_cast(x), incx, (Complex      *)raw_pointer_cast(y), incy);}
		if(is_z<X>{}) {sync_call<hicu(blasZcopy)>(n, (DoubleComplex const*)raw_pointer_cast(x), incx, (DoubleComplex*)raw_pointer_cast(y), incy);}
	}

	template<class ALPHA, class XP, class X = typename std::pointer_traits<XP>::element_type,
		class = decltype(std::declval<X&>() *= ALPHA{}),
		std::enable_if_t<std::is_convertible_v<XP, ::thrust::hicup()::pointer<X>>, int> = 0
	>
	void scal(ssize_t n, ALPHA const& alpha, XP x, ssize_t incx) const {
		if(is_s<X>{}) {sync_call<hicu(blasSscal)>(n, (float         const*)alpha, (float        *)::thrust::raw_pointer_cast(x), incx);}
		if(is_d<X>{}) {sync_call<hicu(blasDscal)>(n, (double        const*)alpha, (double       *)::thrust::raw_pointer_cast(x), incx);}
		if(is_c<X>{}) {sync_call<hicu(blasCscal)>(n, (Complex       const*)alpha, (Complex      *)::thrust::raw_pointer_cast(x), incx);}
		if(is_z<X>{}) {sync_call<hicu(blasZscal)>(n, (DoubleComplex const*)alpha, (DoubleComplex*)::thrust::raw_pointer_cast(x), incx);}
	}

	template<class ALPHA, class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type,
		typename = decltype(std::declval<Y&>() = ALPHA{}*X{} + Y{}),
		std::enable_if_t<std::is_convertible_v<XP, ::thrust::hicup()::pointer<X>> and std::is_convertible_v<YP, ::thrust::hicup()::pointer<Y>>, int> = 0
	>
	void axpy(ssize_t n, ALPHA const* alpha, XP x, ssize_t incx, YP y, ssize_t incy) {
		if(is_d<X>{}) {sync_call<hicu(blasDaxpy)>(n, (double        const*)alpha, (double        const*)raw_pointer_cast(x), incx, (double       *)raw_pointer_cast(y), incy);}
		if(is_z<X>{}) {sync_call<hicu(blasZaxpy)>(n, (DoubleComplex const*)alpha, (DoubleComplex const*)raw_pointer_cast(x), incx, (DoubleComplex*)raw_pointer_cast(y), incy);}
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class XXP, class XX = typename std::pointer_traits<XXP>::element_type, class BETA, class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		typename = decltype(std::declval<YY&>() = ALPHA{}*(AA{}*XX{} + AA{}*XX{})),
		std::enable_if_t<std::is_convertible_v<AAP, ::thrust::hicup()::pointer<AA>> and std::is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> and std::is_convertible_v<YYP, ::thrust::hicup()::pointer<YY>>, int> = 0
	>
	auto gemv(char transA, ssize_t m, ssize_t n, ALPHA const* alpha, AAP aa, ssize_t lda, XXP xx, ssize_t incx, BETA const* beta, YYP yy, ssize_t incy) {
		if(is_d<AA>{}) {sync_call<hicu(blasDgemv)>(operation{transA}, m, n, (double        const*)alpha, (double        const*)::thrust::raw_pointer_cast(aa), lda, (double          const*)::thrust::raw_pointer_cast(xx), incx, (double          const*)beta, (double         *)::thrust::raw_pointer_cast(yy), incy);}
		if(is_z<AA>{}) {sync_call<hicu(blasZgemv)>(operation{transA}, m, n, (DoubleComplex const*)alpha, (DoubleComplex const*)::thrust::raw_pointer_cast(aa), lda, (DoubleComplex const*)::thrust::raw_pointer_cast(xx), incx, (DoubleComplex const*)beta, (DoubleComplex*)::thrust::raw_pointer_cast(yy), incy);}
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename std::pointer_traits<CCP>::element_type,
		typename = decltype(std::declval<CC&>() = ALPHA{}*(AA{}*BB{} + AA{}*BB{})),
		class = std::enable_if_t<std::is_convertible_v<AAP, ::thrust::hicup()::pointer<AA>> and std::is_convertible_v<BBP, ::thrust::hicup()::pointer<BB>> and std::is_convertible_v<CCP, ::thrust::hicup()::pointer<CC>>>
	>
	void gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {
		/*MULTI_MARK_SCOPE("cublasXgemm");*/
		if(is_d<AA>{}) {sync_call<hicu(blasDgemm)>(cuda::cublas::operation{transA}, cuda::cublas::operation{transB}, m, n, k, (double        const*)alpha, (double        const*)::thrust::raw_pointer_cast(aa), lda, (double             const*)::thrust::raw_pointer_cast(bb), ldb, (double          const*)beta, (double         *)::thrust::raw_pointer_cast(cc), ldc);}
		if(is_z<AA>{}) {sync_call<hicu(blasZgemm)>(cuda::cublas::operation{transA}, cuda::cublas::operation{transB}, m, n, k, (DoubleComplex const*)alpha, (DoubleComplex const*)::thrust::raw_pointer_cast(aa), lda, (DoubleComplex const*)::thrust::raw_pointer_cast(bb), ldb, (DoubleComplex const*)beta, (DoubleComplex*)::thrust::raw_pointer_cast(cc), ldc);}
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type,
		std::enable_if_t<
			is_z<AA>{} and is_z<BB>{} and is_assignable<BB&, decltype(AA{}*BB{}/ALPHA{})>{} and is_assignable<BB&, decltype(ALPHA{}*BB{}/AA{})>{} and 
			is_convertible_v<AAP, ::thrust::hicup()::pointer<AA>> and is_convertible_v<BBP, ::thrust::hicup()::pointer<BB>>
		,int> =0
	>
	void trsm(char side, char ul, char transA, char diag, ssize_t m, ssize_t n, ALPHA alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb) {
		sync_call<hicu(blasZtrsm)>(cuda::cublas::side{side}, cuda::cublas::filling{ul}, cuda::cublas::operation{transA}, cuda::cublas::diagonal{diag}, m, n, (DoubleComplex const*)&alpha, (DoubleComplex*)raw_pointer_cast(aa), lda, (DoubleComplex*)raw_pointer_cast(bb), ldb);
	}

	template<class ALPHA, class AAP, class AA = typename std::pointer_traits<AAP>::element_type, class BBP, class BB = typename std::pointer_traits<BBP>::element_type,
		std::enable_if_t<
			is_d<AA>{} and is_d<BB>{} and is_assignable<BB&, decltype(AA{}*BB{}/ALPHA{})>{} and is_assignable<BB&, decltype(ALPHA{}*BB{}/AA{})>{} and 
			is_convertible_v<AAP, ::thrust::hicup()::pointer<AA>> and is_convertible_v<BBP, ::thrust::hicup()::pointer<BB>>
		,int> =0
	>
	void trsm(char side, char ul, char transA, char diag, ssize_t m, ssize_t n, ALPHA alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb) {
		sync_call<hicu(blasDtrsm)>(
			cuda::cublas::side{side},
			cuda::cublas::filling{ul},
			cuda::cublas::operation{transA},
			cuda::cublas::diagonal{diag}, 
			m, n, (double const*)&alpha, (double const*)raw_pointer_cast(aa), lda, (double*)raw_pointer_cast(bb), ldb
		);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_d<XX>{} and is_d<YY>{} and is_d<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> and is_convertible_v<YYP, ::thrust::hicup()::pointer<YY>>
			and is_convertible_v<RRP, RR*>
		, int> =0
	>
	void dot(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		hicu(blasPointerMode_t) mode;
		auto s = hicu(blasGetPointerMode)(get(), &mode); assert( s == HICU(BLAS_STATUS_SUCCESS) );
		assert( mode == HICU(BLAS_POINTER_MODE_HOST) );
		sync_call<hicu(blasDdot)>(n, ::thrust::raw_pointer_cast(xx), incx, ::thrust::raw_pointer_cast(yy), incy, rr);
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_z<XX>{} and is_z<YY>{} and is_z<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> and is_convertible_v<YYP, ::thrust::hicup()::pointer<YY>>
			and (is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>> or is_convertible_v<RRP, RR*>)
		, int> =0
	>
	void dotc(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		hicu(blasPointerMode_t) mode;
		auto s = hicu(blasGetPointerMode)(get(), &mode); assert( s == HICU(BLAS_STATUS_SUCCESS) );
		assert( mode == HICU(BLAS_POINTER_MODE_HOST) );
	//  cublasSetPointerMode(get(), CUBLAS_POINTER_MODE_DEVICE);
		if constexpr(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {
			sync_call<hicu(blasZdotc)>(n, (DoubleComplex const*)::thrust::raw_pointer_cast(xx), incx, (DoubleComplex const*)::thrust::raw_pointer_cast(yy), incy, (DoubleComplex*)::thrust::raw_pointer_cast(rr) );
		} else {
			sync_call<hicu(blasZdotc)>(n, (DoubleComplex const*)::thrust::raw_pointer_cast(xx), incx, (DoubleComplex const*)::thrust::raw_pointer_cast(yy), incy, (DoubleComplex*)rr);
		}
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_z<XX>{}  and is_d<RR>{} and is_assignable<RR&, decltype(XX{}.real())>{} and
			is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> and (is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>> or is_convertible_v<RRP, RR*>)
		, int> =0
	>
	void asum(int n, XXP xx, int incx, RRP rr) {
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_DEVICE));}
		if constexpr(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {
			sync_call<hicu(blasDzasum)>(n, (DoubleComplex const*)::thrust::raw_pointer_cast(xx), incx, (double*)::thrust::raw_pointer_cast(rr) );
		} else {
			sync_call<hicu(blasDzasum)>(n, (DoubleComplex const*)::thrust::raw_pointer_cast(xx), incx, (double*)                           rr  );
		}
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_HOST));}
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type
		,
		std::enable_if_t<
			is_z<XX>::value
		, int> =0
		// ,
		// std::enable_if_t<
		//  is_z<XX>{}  && is_d<RR>{} && is_assignable<RR&, decltype(XX{}.real())>{} &&
		//  is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> && (is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>> or is_convertible_v<RRP, RR*>)
		// , int> =0
	>
	void nrm2(int n, XXP xx, int incx, RRP rr) {
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_DEVICE));}
		if constexpr(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {
			sync_call<hicu(blasDznrm2)>(n, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(xx)), incx, reinterpret_cast<double*>(::thrust::raw_pointer_cast(rr)) );
		} else {
			sync_call<hicu(blasDznrm2)>(n, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(xx)), incx, reinterpret_cast<double*>(                           rr ) );
		}
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_HOST));}
	}

	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type
		,
		std::enable_if_t<
			is_d<XX>::value
		, int> =0
		// ,
		// std::enable_if_t<
		//  is_z<XX>{}  && is_d<RR>{} && is_assignable<RR&, decltype(XX{}.real())>{} &&
		//  is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> && (is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>> or is_convertible_v<RRP, RR*>)
		// , int> =0
	>
	void nrm2(int n, XXP xx, int incx, RRP rr) {
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_DEVICE));}
		if constexpr(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {
			sync_call<hicu(blasDnrm2)>(n, ::thrust::raw_pointer_cast(xx), incx, reinterpret_cast<double*>(::thrust::raw_pointer_cast(rr)) );
		} else {
			sync_call<hicu(blasDnrm2)>(n, ::thrust::raw_pointer_cast(xx), incx, reinterpret_cast<double*>(                           rr ) );
		}
		if(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {hicu(blasSetPointerMode)(get(), HICU(BLAS_POINTER_MODE_HOST));}
	}


	template<
		class XXP, class XX = typename std::pointer_traits<XXP>::element_type,
		class YYP, class YY = typename std::pointer_traits<YYP>::element_type,
		class RRP, class RR = typename std::pointer_traits<RRP>::element_type,
		std::enable_if_t<
			is_z<XX>{} and is_z<YY>{} and is_z<RR>{} and is_assignable<RR&, decltype(XX{}*YY{})>{} and
			is_convertible_v<XXP, ::thrust::hicup()::pointer<XX>> and is_convertible_v<YYP, ::thrust::hicup()::pointer<YY>>
			and (is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>> or is_convertible_v<RRP, RR*>)
		, int> =0
	>
	void dotu(int n, XXP xx, int incx, YYP yy, int incy, RRP rr) {
		hicu(blasPointerMode_t) mode;
		auto s = hicu(blasGetPointerMode)(get(), &mode); assert( s == HICU(BLAS_STATUS_SUCCESS) );
		assert( mode == HICU(BLAS_POINTER_MODE_HOST) );
	//  cublasSetPointerMode(get(), CUBLAS_POINTER_MODE_DEVICE);
		if constexpr(is_convertible_v<RRP, ::thrust::hicup()::pointer<RR>>) {
			sync_call<hicu(blasZdotu)>(n, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(xx)), incx, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(yy)), incy, reinterpret_cast<DoubleComplex*>(::thrust::raw_pointer_cast(rr)) );
		} else {
			sync_call<hicu(blasZdotu)>(n, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(xx)), incx, reinterpret_cast<DoubleComplex const*>(::thrust::raw_pointer_cast(yy)), incy, reinterpret_cast<DoubleComplex*>(rr));
		}
	//  cublasSetPointerMode(get(), CUBLAS_POINTER_MODE_HOST);
	}
};

}  // end namespace multi::cuda::cublas
}  // end namespace boost

namespace boost::multi::blas {

	template<> struct is_context<boost::multi::cuda::cublas::context > : std::true_type {};
	template<> struct is_context<boost::multi::cuda::cublas::context&> : std::true_type {};

	template<class Ptr, class T = typename std::pointer_traits<Ptr>::element_type, std::enable_if_t<std::is_convertible<Ptr, ::thrust::hicup()::pointer<T>>{}, int> =0>
	boost::multi::cuda::cublas::context* default_context_of(Ptr const&) {
		namespace multi = boost::multi;
		return &multi::cuda::cublas::context::get_instance();
	}

	template<class T, class R>
	boost::multi::cuda::cublas::context*
	#if defined(__HIPCC__)
	default_context_of(::thrust::pointer<T, ::thrust::hip::tag, R> const&) {
	#else  //  __NVCC__
	default_context_of(::thrust::pointer<T, ::thrust::cuda_cub::tag, R> const&) {
	#endif
		namespace multi = boost::multi;
		return &multi::cuda::cublas::context::get_instance();
	}
}

#undef hicup
#undef hicu
#undef HICU
