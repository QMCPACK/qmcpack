// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_CUSOLVER_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_CUSOLVER_HPP
#pragma once

// cuSOLVER (CUDA) / hipSOLVER (ROCm) backend for the Multi LAPACK adaptor.
// Compiling with -DMULTI_USE_HIP (or including lapack/hipsolver.hpp instead) selects hipSOLVER,
// in the same way blas is handled by adaptors/cuda/cublas/context.hpp.
// It provides potrf/syev overloads, found by ADL, for thrust device (and universal) pointers,
// so the generic wrappers in lapack/potrf.hpp and lapack/syev.hpp dispatch to the GPU
// when used with thrust-based arrays.
// This header supersedes the (legacy) lapack/cuda.hpp, which hooks on pre-thrust pointer types;
// do not include both.

#include <boost/multi/adaptors/blas/traits.hpp>  // for is_s, is_d, is_c, is_z

#if !defined(MULTI_USE_HIP)
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <thrust/system/cuda/memory.h>  // for thrust::cuda::pointer
#else
#include <hip/hip_runtime.h>
#include <hipsolver/hipsolver.h>
#include <thrust/system/hip/memory.h>  // for thrust::hip::pointer
#endif

#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <memory>       // for unique_ptr, pointer_traits
#include <stdexcept>    // for runtime_error
#include <string>       // for to_string
#include <type_traits>  // for enable_if_t, is_convertible_v

#if !defined(MULTI_USE_HIP)
#define hicusolver(name) cusolver##name      // NOLINT(cppcoreguidelines-macro-usage) e.g. cusolverStatus_t   / hipsolverStatus_t
#define hicusolverDn(name) cusolverDn##name  // NOLINT(cppcoreguidelines-macro-usage) e.g. cusolverDnCreate   / hipsolverCreate (regular hipSOLVER API mirrors the cusolverDn API)
#define HICUSOLVER(name) CUSOLVER##name      // NOLINT(cppcoreguidelines-macro-usage) e.g. CUSOLVER_STATUS_SUCCESS / HIPSOLVER_STATUS_SUCCESS
#define hicup(name) cuda##name               // NOLINT(cppcoreguidelines-macro-usage) runtime API, e.g. cudaMalloc / hipMalloc
#ifndef thrust_hicup
#define thrust_hicup thrust::cuda  // NOLINT(cppcoreguidelines-macro-usage) also defined by adaptors/cuda/cublas/context.hpp
#endif
#else
#define hicusolver(name) hipsolver##name    // NOLINT(cppcoreguidelines-macro-usage)
#define hicusolverDn(name) hipsolver##name  // NOLINT(cppcoreguidelines-macro-usage)
#define HICUSOLVER(name) HIPSOLVER##name    // NOLINT(cppcoreguidelines-macro-usage)
#define hicup(name) hip##name               // NOLINT(cppcoreguidelines-macro-usage)
#ifndef thrust_hicup
#define thrust_hicup thrust::hip  // NOLINT(cppcoreguidelines-macro-usage)
#endif
#endif

namespace boost::multi::cusolver {

using blas::is_s;
using blas::is_d;
using blas::is_c;
using blas::is_z;

#if !defined(MULTI_USE_HIP)
using Complex       = cuComplex;
using DoubleComplex = cuDoubleComplex;
#else
using Complex       = hipFloatComplex;
using DoubleComplex = hipDoubleComplex;
#endif

template<class T, std::enable_if_t<is_s<T>{}, int> = 0> constexpr auto data_cast(T* ptr) { return reinterpret_cast<float*>(ptr); }          // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
template<class T, std::enable_if_t<is_d<T>{}, int> = 0> constexpr auto data_cast(T* ptr) { return reinterpret_cast<double*>(ptr); }         // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
template<class T, std::enable_if_t<is_c<T>{}, int> = 0> constexpr auto data_cast(T* ptr) { return reinterpret_cast<Complex*>(ptr); }        // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
template<class T, std::enable_if_t<is_z<T>{}, int> = 0> constexpr auto data_cast(T* ptr) { return reinterpret_cast<DoubleComplex*>(ptr); }  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

inline void solver_check(hicusolver(Status_t) stat, char const* what) {
	if(stat != HICUSOLVER(_STATUS_SUCCESS)) {
		throw std::runtime_error{std::string{what} + " failed with solver status " + std::to_string(static_cast<int>(stat))};
	}
}

inline void runtime_check(hicup(Error_t) err, char const* what) {
	if(err != hicup(Success)) {
		throw std::runtime_error{std::string{what} + " failed with runtime error " + std::to_string(static_cast<int>(err))};
	}
}

#if !defined(MULTI_USE_HIP)
inline auto fill_mode(char uplo) {
	assert(uplo == 'U' || uplo == 'L');
	return uplo == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
}
#else
inline auto fill_mode(char uplo) {
	assert(uplo == 'U' || uplo == 'L');
	return uplo == 'U' ? HIPSOLVER_FILL_MODE_UPPER : HIPSOLVER_FILL_MODE_LOWER;
}
#endif

inline auto eig_mode(char jobz) {
	assert(jobz == 'V' || jobz == 'N');
	return jobz == 'V' ? HICUSOLVER(_EIG_MODE_VECTOR) : HICUSOLVER(_EIG_MODE_NOVECTOR);
}

// scratch device memory allocated directly through the runtime,
// so this header stays usable from host-only translation units (no kernel launches involved)
class device_buffer {
	void* ptr_ = nullptr;

 public:
	explicit device_buffer(std::size_t bytes) { runtime_check(hicup(Malloc)(&ptr_, bytes), "device malloc"); }
	device_buffer(device_buffer const&)                    = delete;
	device_buffer(device_buffer&&)                         = delete;
	auto operator=(device_buffer const&) -> device_buffer& = delete;
	auto operator=(device_buffer&&) -> device_buffer&      = delete;
	~device_buffer() { [[maybe_unused]] auto const err = hicup(Free)(ptr_); }  // cannot throw from a destructor

	template<class T = void> auto data() const { return static_cast<T*>(ptr_); }
};

inline auto device_info_to_host(int const* dinfo) -> int {
	int info = -1;
	runtime_check(hicup(Memcpy)(&info, dinfo, sizeof(int), hicup(MemcpyDeviceToHost)), "device info memcpy");
	return info;
}

class context : private std::unique_ptr<typename std::pointer_traits<hicusolverDn(Handle_t)>::element_type, decltype(&hicusolverDn(Destroy))> {
	using pimpl_t = std::unique_ptr<typename std::pointer_traits<hicusolverDn(Handle_t)>::element_type, decltype(&hicusolverDn(Destroy))>;

 public:
	using pimpl_t::get;
	context() : pimpl_t{
		[] {
			hicusolverDn(Handle_t) handle{};
			solver_check(hicusolverDn(Create)(&handle), "solver handle create");
			return handle;
		}(),
		&hicusolverDn(Destroy)
	} {}
	static auto get_instance() -> context& {
		thread_local context ctxt;
		return ctxt;
	}
};

// the LAPACK char conventions (and the row-major <-> column-major flip) are handled by the callers
// in lapack/potrf.hpp and lapack/syev.hpp, exactly as for the CPU bindings in lapack/core.hpp

template<class T>
void potrf(char uplo, int n, T* aa, int lda, int& info) {  // NOLINT(readability-identifier-length) conventional lapack name
	static_assert(is_s<T>{} || is_d<T>{} || is_c<T>{} || is_z<T>{}, "potrf supports float/double real and complex elements");

	auto& ctx = context::get_instance();

	int lwork = -1;
	/**/ if constexpr(is_s<T>{}) { solver_check(hicusolverDn(Spotrf_bufferSize)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, &lwork), "spotrf buffer size"); }
	else if constexpr(is_d<T>{}) { solver_check(hicusolverDn(Dpotrf_bufferSize)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, &lwork), "dpotrf buffer size"); }
	else if constexpr(is_c<T>{}) { solver_check(hicusolverDn(Cpotrf_bufferSize)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, &lwork), "cpotrf buffer size"); }
	else if constexpr(is_z<T>{}) { solver_check(hicusolverDn(Zpotrf_bufferSize)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, &lwork), "zpotrf buffer size"); }
	assert(lwork >= 0);

	device_buffer work{sizeof(T) * static_cast<std::size_t>(lwork)};
	device_buffer dinfo{sizeof(int)};

	/**/ if constexpr(is_s<T>{}) { solver_check(hicusolverDn(Spotrf)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, data_cast(work.data<T>()), lwork, dinfo.data<int>()), "spotrf"); }
	else if constexpr(is_d<T>{}) { solver_check(hicusolverDn(Dpotrf)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, data_cast(work.data<T>()), lwork, dinfo.data<int>()), "dpotrf"); }
	else if constexpr(is_c<T>{}) { solver_check(hicusolverDn(Cpotrf)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, data_cast(work.data<T>()), lwork, dinfo.data<int>()), "cpotrf"); }
	else if constexpr(is_z<T>{}) { solver_check(hicusolverDn(Zpotrf)(ctx.get(), fill_mode(uplo), n, data_cast(aa), lda, data_cast(work.data<T>()), lwork, dinfo.data<int>()), "zpotrf"); }

	runtime_check(hicup(DeviceSynchronize)(), "device synchronize");
	info = device_info_to_host(dinfo.data<int>());  // 0: success, >0: order of first non-positive-definite minor, as in LAPACK
}

// implemented in terms of the divide-and-conquer version (syevd), like the cusolverDn examples;
// the workspace is sized and allocated internally, the caller-provided lapack-style work array is ignored
template<class T>
void syev(char jobz, char uplo, int n, T* aa, int lda, T* ww, int& info) {  // NOLINT(readability-identifier-length) conventional lapack name
	static_assert(is_s<T>{} || is_d<T>{}, "syev(d) supports real float/double elements; complex Hermitian (heev) is not implemented yet");

	auto& ctx = context::get_instance();

	int lwork = -1;
	if constexpr(is_s<T>{}) { solver_check(hicusolverDn(Ssyevd_bufferSize)(ctx.get(), eig_mode(jobz), fill_mode(uplo), n, data_cast(aa), lda, data_cast(ww), &lwork), "ssyevd buffer size"); }
	else /*  is_d<T>{}  */ { solver_check(hicusolverDn(Dsyevd_bufferSize)(ctx.get(), eig_mode(jobz), fill_mode(uplo), n, data_cast(aa), lda, data_cast(ww), &lwork), "dsyevd buffer size"); }
	assert(lwork >= 0);

	device_buffer work{sizeof(T) * static_cast<std::size_t>(lwork)};
	device_buffer dinfo{sizeof(int)};

	if constexpr(is_s<T>{}) { solver_check(hicusolverDn(Ssyevd)(ctx.get(), eig_mode(jobz), fill_mode(uplo), n, data_cast(aa), lda, data_cast(ww), data_cast(work.data<T>()), lwork, dinfo.data<int>()), "ssyevd"); }
	else /*  is_d<T>{}  */ { solver_check(hicusolverDn(Dsyevd)(ctx.get(), eig_mode(jobz), fill_mode(uplo), n, data_cast(aa), lda, data_cast(ww), data_cast(work.data<T>()), lwork, dinfo.data<int>()), "dsyevd"); }

	runtime_check(hicup(DeviceSynchronize)(), "device synchronize");
	info = device_info_to_host(dinfo.data<int>());
}

}  // end namespace boost::multi::cusolver

namespace thrust {

// ADL hooks: the generic wrappers in lapack/potrf.hpp and lapack/syev.hpp make unqualified calls with
// the arrays' element pointers; for thrust device (or universal) pointers those calls land here.
// The enable_if keeps other thrust systems (omp, tbb, cpp) away from the GPU backend.

template<
	class S1, class S2, class Ptr, class T = typename std::pointer_traits<Ptr>::element_type,
	std::enable_if_t<std::is_convertible_v<Ptr, ::thrust_hicup::pointer<T>>, int> = 0
>
void potrf(char uplo, S1 n, Ptr aa, S2 lda, int& info) {  // NOLINT(readability-identifier-length) conventional lapack name
	::boost::multi::cusolver::potrf(uplo, static_cast<int>(n), ::thrust::raw_pointer_cast(aa), static_cast<int>(lda), info);
}

template<
	class S, class PtrA, class PtrW, class PtrWork,
	class T = typename std::pointer_traits<PtrA>::element_type,
	std::enable_if_t<
		std::is_convertible_v<PtrA, ::thrust_hicup::pointer<T>> &&
		std::is_convertible_v<PtrW, ::thrust_hicup::pointer<T>>, int
	> = 0
>
void syev(char jobz, char uplo, S n, PtrA aa, S lda, PtrW ww, PtrWork /*work*/, S /*lwork*/, int& info) {  // NOLINT(readability-identifier-length) conventional lapack name
	::boost::multi::cusolver::syev(jobz, uplo, static_cast<int>(n), ::thrust::raw_pointer_cast(aa), static_cast<int>(lda), ::thrust::raw_pointer_cast(ww), info);
}

}  // end namespace thrust

#undef hicusolver
#undef hicusolverDn
#undef HICUSOLVER
#undef hicup

#endif  // BOOST_MULTI_ADAPTORS_LAPACK_CUSOLVER_HPP
