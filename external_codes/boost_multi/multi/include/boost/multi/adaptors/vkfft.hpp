// Copyright 2024-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_VKFFT_HPP
#define BOOST_MULTI_ADAPTORS_VKFFT_HPP

#include <boost/multi/array.hpp>

#include <thrust/memory.h>  // for raw_pointer_cast

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

// VkFFT selects its GPU backend at compile time via the VKFFT_BACKEND macro
//   0 - Vulkan, 1 - CUDA, 2 - HIP, 3 - OpenCL, 4 - Level Zero, 5 - Metal
// If the user did not set it, pick a sensible default from the compiler in use.
#if !defined(VKFFT_BACKEND)
#if defined(__HIPCC__)
#define VKFFT_BACKEND 2
#elif defined(__CUDACC__) || defined(__NVCC__)
#define VKFFT_BACKEND 1
#else
#define VKFFT_BACKEND 0
#endif
#endif

#if !defined(VKFFT_MAX_FFT_DIMENSIONS)
#define VKFFT_MAX_FFT_DIMENSIONS 12
#endif

#if(VKFFT_BACKEND == 1)
#include <cuda.h>  // CUDA driver API (VkFFT CUDA backend uses the driver API + NVRTC)
#elif(VKFFT_BACKEND == 2)
#include <hip/hip_runtime.h>  // VkFFT HIP backend (uses HIPRTC)
#endif

// VkFFT is third-party and header-only, and it is not warning-clean under nvcc's
// `-Werror all-warnings` (which CMake's CMAKE_COMPILE_WARNING_AS_ERROR enables):
// e.g. `uint64_t initPageSize = -1;` in vkFFT_UpdateBuffers.h raises diagnostic
// #68-D ("integer conversion resulted in a change of sign").  Silence VkFFT's own
// diagnostics across its include only; our code keeps the full warning set.
#if defined(__NVCC__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress integer_sign_change  // #68-D, from `uint64_t initPageSize = -1;`
#pragma nv_diag_suppress set_but_not_used     // #550-D, e.g. `temp_int` in vkFFT_CodeGen/.../vkFFT_Zeropad.h
#pragma nv_diag_suppress declared_but_not_referenced  // #177-D, e.g. `raderFFTKernel` in vkFFT_CodeGen/.../vkFFT_Constants.h
#endif

#include <vkFFT.h>  // external VkFFT library (header-only, multi-header since v1.3)

#if defined(__NVCC__)
#pragma nv_diagnostic pop
#endif

namespace boost::multi::vkfft {

// This slice of the adaptor implements the CUDA and HIP backends (they share the
// same `void**` device-pointer path) for interleaved complex<double> (C2C)
// transforms, mirroring the scope of boost/multi/adaptors/cufft.hpp.  Metal and
// Vulkan need a buffer-handle abstraction and are not covered yet.
static_assert(VKFFT_BACKEND == 1 || VKFFT_BACKEND == 2, "boost::multi::vkfft currently only supports the CUDA (VKFFT_BACKEND=1) and HIP (VKFFT_BACKEND=2) backends");

// Direction of the transform.  VkFFTAppend takes -1 for the forward transform
// (exponent sign -1, same convention as FFTW/cuFFT) and +1 for the inverse.
class sign {
	int impl_ = 0;

 public:
	sign()                            = default;
	constexpr explicit sign(int impl) : impl_{impl} {}
	constexpr explicit operator int() const { return impl_; }

	friend constexpr auto operator==(sign a, sign b) { return a.impl_ == b.impl_; }
	friend constexpr auto operator!=(sign a, sign b) { return a.impl_ != b.impl_; }
};

constexpr sign forward{-1};
constexpr sign none{0};
constexpr sign backward{+1};

static_assert(forward != none && none != backward && backward != forward);

namespace detail {

#if(VKFFT_BACKEND == 1)
using device_type = CUdevice;
#elif(VKFFT_BACKEND == 2)
using device_type = hipDevice_t;
#endif

// A single, stable device handle.  VkFFTConfiguration stores the *pointer* we
// hand it and dereferences it again at append time, so it must outlive every
// plan; a function-local static is the simplest object with that lifetime.
inline auto device_ptr() -> device_type* {
	static device_type const dev = [] {
		device_type d{};
#if(VKFFT_BACKEND == 1)
		static std::once_flag once;
		std::call_once(once, [] { cuInit(0); });
		if(cuCtxGetDevice(&d) != CUDA_SUCCESS) {  // no current context (pure runtime-API user): fall back to device 0
			cuDeviceGet(&d, 0);
		}
#elif(VKFFT_BACKEND == 2)
		int ordinal = 0;
		if(hipGetDevice(&ordinal) == hipSuccess) {  // hipDevice_t is the plain ordinal
			static_cast<void>(hipDeviceGet(&d, ordinal));
		}
#endif
		return d;
	}();
	return const_cast<device_type*>(&dev);  // NOLINT(cppcoreguidelines-pro-type-const-cast) VkFFT wants a non-const pointer it does not write through
}

template<class Tuple>
auto tuple_to_u64_array(Tuple const& tup) {
	return std::apply(
		[](auto... elems) {
			return std::array<std::uint64_t, sizeof...(elems)>{static_cast<std::uint64_t>(elems)...};
		},
		tup
	);
}

inline void check(VkFFTResult res, char const* what) {
	if(res != VKFFT_SUCCESS) {
		throw std::runtime_error{std::string{"vkfft: "} + what + " failed (VkFFTResult " + std::to_string(static_cast<int>(res)) + ")"};
	}
}

}  // end namespace detail

// Environment is kept for API parity with fftw::environment.  VkFFT itself has
// no global state to set up or tear down; this only forces device initialization.
struct environment {
	environment() { detail::device_ptr(); }

	environment(environment const&) = delete;
	environment(environment&&)      = delete;
	auto operator=(environment const&) -> environment& = delete;
	auto operator=(environment&&) -> environment&      = delete;
	~environment()                                     = default;
};

template<dimensionality_type D>
class plan {
	static_assert(D >= 1 && D <= VKFFT_MAX_FFT_DIMENSIONS, "VkFFT dimensionality out of range (raise -DVKFFT_MAX_FFT_DIMENSIONS)");

	VkFFTApplication app_{};
	std::uint64_t    in_bytes_{};
	std::uint64_t    out_bytes_{};

	struct axis_ {
		std::uint64_t n;   // extent
		std::uint64_t is;  // input  stride, in elements
		std::uint64_t os;  // output stride, in elements
		bool          tr;  // is this axis transformed?
	};

 public:
	using size_type = std::uint64_t;

	// Non-movable: VkFFT keeps the addresses of `in_bytes_`/`out_bytes_` (via
	// VkFFTConfiguration::bufferSize) and dereferences them again at append time,
	// so a plan must stay put once initialized.  `cached_plan` builds it in place.
	plan()            = delete;
	plan(plan const&) = delete;
	plan(plan&&)      = delete;

	auto operator=(plan const&) -> plan& = delete;
	auto operator=(plan&&) -> plan&      = delete;

	template<
		class ILayout, class OLayout,
		std::enable_if_t<static_cast<dimensionality_type>(std::decay_t<ILayout>::dimensionality) == D && static_cast<dimensionality_type>(std::decay_t<OLayout>::dimensionality) == D, int> = 0>
	plan(std::array<bool, D> which, ILayout const& in_layout, OLayout const& out_layout) {
		assert(in_layout.sizes() == out_layout.sizes());

		auto const ns   = detail::tuple_to_u64_array(in_layout.sizes());
		auto const istr = detail::tuple_to_u64_array(in_layout.strides());
		auto const ostr = detail::tuple_to_u64_array(out_layout.strides());

		if(std::none_of(which.begin(), which.end(), [](bool b) { return b; })) {
			throw std::runtime_error{"vkfft: no transformed axis is not supported"};
		}

		// Order the axes the way VkFFT expects them: fastest (smallest stride)
		// first.  VkFFT axis 0 carries an implicit unit stride (the config has no
		// field for it), so whichever multi axis is contiguous -- transformed or
		// batched, and regardless of its position in the layout -- must land at
		// VkFFT position 0.  multi's `.strides()` are *not* guaranteed sorted (a
		// subarray / rotated view can put the unit stride in the middle), so we
		// sort here rather than assume axis D-1 is the fast one.
		std::array<axis_, static_cast<std::size_t>(D)> ax{};
		for(std::size_t i = 0; i != static_cast<std::size_t>(D); ++i) { ax[i] = axis_{ns[i], istr[i], ostr[i], which[i]}; }
		std::stable_sort(ax.begin(), ax.end(), [](axis_ const& a, axis_ const& b) { return a.is < b.is; });  // this is the natural order for vfFFT

		assert(ax.front().is == 1 && "vkfft: the input needs a unit-stride (contiguous) axis");
		assert(ax.front().os == 1 && "vkfft: the output's contiguous axis must be the same one as the input's");

		VkFFTConfiguration cfg{};  // MUST be zero-initialized (C struct)
		cfg.FFTdim = D;
		for(dimensionality_type j = 0; j != D; ++j) {
			cfg.size[j]          = ax[static_cast<std::size_t>(j)].n;
			cfg.omitDimension[j] = ax[static_cast<std::size_t>(j)].tr ? 0U : 1U;  // 0 = transform this axis, 1 = skip it (batch)
		}
		// bufferStride[j] is the element stride to step once along VkFFT axis j+1.
		for(dimensionality_type j = 0; j + 1 < D; ++j) {
			cfg.inputBufferStride[j] = ax[static_cast<std::size_t>(j) + 1].is;
			cfg.bufferStride[j]      = ax[static_cast<std::size_t>(j) + 1].os;
		}
		cfg.inputBufferStride[D - 1] = ax.back().is * ax.back().n;  // batch stride; inert while numberBatches stays 1
		cfg.bufferStride[D - 1]      = ax.back().os * ax.back().n;

		cfg.doublePrecision = 1;
		cfg.device          = detail::device_ptr();

		// Always configure an out-of-place ("formatted input") transform: the main
		// `buffer` is the destination, `inputBuffer` the source.  A caller doing an
		// in-place transform simply passes the same device pointer as both at launch
		// (VkFFT reads the input into registers before writing the buffer).  We cannot
		// tell in-place from out-of-place here because only layouts, not pointers, are
		// known at plan-construction time.
		cfg.isInputFormatted = 1;

		in_bytes_           = ax.back().is * ax.back().n * 2 * sizeof(double);
		out_bytes_          = ax.back().os * ax.back().n * 2 * sizeof(double);
		cfg.bufferSize      = &out_bytes_;
		cfg.inputBufferSize = &in_bytes_;

		detail::check(initializeVkFFT(&app_, cfg), "initializeVkFFT");
	}

	~plan() { deleteVkFFT(&app_); }

	template<class IPtr, class OPtr>
	auto execute(IPtr idata, OPtr odata, sign dir) const
		-> decltype((void)(::thrust::raw_pointer_cast(idata), ::thrust::raw_pointer_cast(odata))) {
		static_assert(sizeof(*::thrust::raw_pointer_cast(idata)) == 2 * sizeof(double), "vkfft first slice handles interleaved complex<double> only");

		auto* in_ptr  = ::thrust::raw_pointer_cast(idata);
		auto* out_ptr = ::thrust::raw_pointer_cast(odata);

		void* in_raw  = const_cast<void*>(static_cast<void const*>(in_ptr));  // NOLINT(cppcoreguidelines-pro-type-const-cast) VkFFT does not write the input in C2C out-of-place
		void* out_raw = static_cast<void*>(out_ptr);

		VkFFTLaunchParams lp{};
		lp.buffer      = &out_raw;
		lp.inputBuffer = &in_raw;  // may alias out_raw for an in-place transform

		detail::check(VkFFTAppend(const_cast<VkFFTApplication*>(&app_), static_cast<int>(dir), &lp), "VkFFTAppend");  // NOLINT(cppcoreguidelines-pro-type-const-cast) VkFFT append mutates internal scratch only
	}

	template<class IPtr, class OPtr>
	auto execute_forward(IPtr idata, OPtr odata) const -> decltype(execute(idata, odata, forward)) {
		return execute(idata, odata, forward);
	}
	template<class IPtr, class OPtr>
	auto execute_backward(IPtr idata, OPtr odata) const -> decltype(execute(idata, odata, backward)) {
		return execute(idata, odata, backward);
	}

	template<class IPtr, class OPtr>
	void operator()(IPtr idata, OPtr odata, sign dir) const { execute(idata, odata, dir); }
};

template<dimensionality_type D>
class cached_plan {
	typename std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D>>::iterator it_;

 public:
	cached_plan(cached_plan const&) = delete;
	cached_plan(cached_plan&&)      = delete;
	auto operator=(cached_plan const&) -> cached_plan& = delete;
	auto operator=(cached_plan&&) -> cached_plan&      = delete;
	~cached_plan()                                     = default;

	cached_plan(std::array<bool, D> which, multi::layout_t<D, multi::ssize_t> in, multi::layout_t<D, multi::ssize_t> out) {
		thread_local std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D>>& LEAKY_cache =
			*new std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D>>;  // NOLINT(cppcoreguidelines-owning-memory) intentional leak, mirrors cufft::cached_plan

		auto const key = std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>{which, in, out};
		it_            = LEAKY_cache.find(key);
		if(it_ == LEAKY_cache.end()) {
			it_ = LEAKY_cache.try_emplace(key, which, in, out).first;  // constructs plan<D> in place (it is non-movable)
		}
	}

	template<class IPtr, class OPtr>
	auto execute(IPtr idata, OPtr odata, sign dir) -> decltype(std::declval<plan<D> const&>().execute(idata, odata, dir)) {
		return it_->second.execute(idata, odata, dir);
	}
};

template<class In, class Out, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft(std::array<bool, +D> which, In const& in, Out&& out, sign dir)
	-> decltype(cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), dir), std::forward<Out>(out)) {
	cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), dir);
	return std::forward<Out>(out);
}

template<class In, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft(std::array<bool, +D> which, In&& in, sign dir)
	-> decltype(dft(which, in, in, dir), std::forward<In>(in)) {
	dft(which, in, in, dir);
	return std::forward<In>(in);
}

template<class In, class R = typename std::decay_t<In>::decay_type>
auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In const& in, sign dir) -> R {
	static_assert(std::is_trivially_default_constructible_v<typename std::decay_t<In>::element>);
	R ret(in.extensions(), get_allocator(in));
	vkfft::dft(which, in, ret, dir);
	return ret;
}

template<class In, class Out, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft_forward(std::array<bool, +D> which, In const& in, Out&& out)
	-> decltype(vkfft::dft(which, in, std::forward<Out>(out), forward)) {
	return vkfft::dft(which, in, std::forward<Out>(out), forward);
}

template<class In, class Out, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft_backward(std::array<bool, +D> which, In const& in, Out&& out)
	-> decltype(vkfft::dft(which, in, std::forward<Out>(out), backward)) {
	return vkfft::dft(which, in, std::forward<Out>(out), backward);
}

template<class In, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft_forward(std::array<bool, +D> which, In&& in) -> decltype(vkfft::dft(which, std::forward<In>(in), forward)) {
	return vkfft::dft(which, std::forward<In>(in), forward);
}

template<class In, dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft_backward(std::array<bool, +D> which, In&& in) -> decltype(vkfft::dft(which, std::forward<In>(in), backward)) {
	return vkfft::dft(which, std::forward<In>(in), backward);
}

}  // end namespace boost::multi::vkfft

#endif  // BOOST_MULTI_ADAPTORS_VKFFT_HPP
