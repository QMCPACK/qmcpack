// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_CUFFT_HPP
#define BOOST_MULTI_ADAPTORS_CUFFT_HPP

#include "boost/multi/array.hpp"
#include "boost/multi/detail/config/NODISCARD.hpp"
#include "boost/multi/utility.hpp"

#include <thrust/memory.h>  // for raw_pointer_cast

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <tuple>
#include <type_traits>

#if !defined(__HIP_ROCclr__)
#include <cufft.h>
#include <cufftXt.h>
#endif

namespace boost::multi::cufft {

// cuFFT API errors
// static auto cuda_get_error_enum(cufftResult error) -> char const* {
// 	switch(error) {
// 	case CUFFT_SUCCESS: return "CUFFT_SUCCESS";

// 	case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
// 	case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
// #ifdef CUFFT_INCOMPLETE_PARAMETER_LIST
// 	case CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
// #endif
// 	case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
// 	case CUFFT_INVALID_DEVICE: return "CUFFT_INVALID_DEVICE";
// 	case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
// 	case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
// 	case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
// 	case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
// 	case CUFFT_NO_WORKSPACE: return "CUFFT_NO_WORKSPACE";
// 	case CUFFT_NOT_IMPLEMENTED: return "CUFFT_NOT_IMPLEMENTED";
// 	case CUFFT_NOT_SUPPORTED: return "CUFFT_NOT_SUPPORTED";
// 	// #if !defined(__HIP_PLATFORM_NVIDIA__)
// 	// case CUFFT_PARSE_ERROR:    return "CUFFT_PARSE_ERROR";
// 	// #endif
// 	case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
// 	case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
// 	// #if !defined(__HIP_PLATFORM_NVIDIA__)
// 	// case CUFFT_LICENSE_ERROR:  return "CUFFT_LICENSE_ERROR";
// 	// #endif
// 	default: assert(0);
// 	}
// 	return "<unknown>";
// }

#define cufftSafeCall(err) implcufftSafeCall(err, __FILE__, __LINE__)
inline void implcufftSafeCall(cufftResult err, const char* /*file*/, const int /*line*/) {
	if(CUFFT_SUCCESS != err) {
		// std::cerr << "CUFFT error in file " << file << ", line " << line << "\nerror " << err << ": " << cuda_get_error_enum(err) << "\n";
		// fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n", __FILE__, __LINE__, err,
		//                         _cudaGetErrorEnum(err));
		cudaDeviceReset() == cudaSuccess ? void() : assert(0);
	}
}

class sign {
	int impl_ = 0;

 public:
	sign() = default;
	constexpr explicit sign(int impl) : impl_{impl} {}
	constexpr operator int() const { return impl_; }

	constexpr auto operator==(sign const& other) const { return impl_ == other.impl_; }
	constexpr auto operator!=(sign const& other) const { return impl_ != other.impl_; }
};

constexpr sign forward{CUFFT_FORWARD};
constexpr sign none{0};
constexpr sign backward{CUFFT_INVERSE};
// constexpr sign backward{CUFFT_BACKWARD};

static_assert(forward != none && none != backward && backward != forward);

struct cufft_iodim64 {
	std::ptrdiff_t n;
	std::ptrdiff_t is;
	std::ptrdiff_t os;
};  // based on fftw_iodim64

template<dimensionality_type DD, class Alloc = void*>
class plan {
	Alloc                                              alloc_;
	::size_t                                           workSize_ = 0;
	void*                                              workArea_{};
	cufftHandle                                        h_{};  // TODO(correaa) put this in a unique_ptr
	std::array<std::pair<bool, cufft_iodim64>, DD + 1> which_iodims_{};
	int                                                first_howmany_{};

	// mutable bool used_ = false;

	using complex_type = cufftDoubleComplex;

 public:
	using allocator_type = Alloc;

	template<
		class ILayout, class OLayout, dimensionality_type D = std::decay_t<ILayout>::dimensionality,
		class = std::enable_if_t<D == std::decay_t<OLayout>::dimensionality>>
	plan(std::array<bool, +D> which, ILayout const& in, OLayout const& out) : plan(which, in, out, allocator_type{}) {}

	plan()            = delete;
	plan(plan const&) = delete;

	plan(plan&& other) noexcept
	: alloc_{std::move(other.alloc_)},
	  workSize_{std::exchange(other.workSize_, {})},
	  workArea_{std::exchange(other.workArea_, {})},
	  h_{std::exchange(other.h_, {})},
	  which_iodims_{std::exchange(other.which_iodims_, {})},
	  first_howmany_{std::exchange(other.first_howmany_, {})} {
		// other.used_ = true;  // moved-from object cannot be used
		// used_       = false;
	}

	auto operator=(plan const&) = delete;
	auto operator=(plan&&)      = delete;

	template<
		class ILayout, class OLayout, dimensionality_type D = std::decay_t<ILayout>::dimensionality,
		class = std::enable_if_t<D == std::decay_t<OLayout>::dimensionality>>
	plan(std::array<bool, +D> which, ILayout const& in, OLayout const& out, allocator_type const& alloc) : alloc_{alloc} {
		// used_ = false;
		assert(in.sizes() == out.sizes());

		auto const sizes_tuple   = in.sizes();
		auto const istride_tuple = in.strides();
		auto const ostride_tuple = out.strides();

		using boost::multi::detail::get;
		auto which_iodims = std::apply([](auto... elems) {
			return std::array<std::pair<bool, cufft_iodim64>, sizeof...(elems) + 1>{
  // TODO(correaa) added one element to avoid problem with gcc 13 static analysis (out-of-bounds)
				std::pair<bool, cufft_iodim64>{
											   get<0>(elems),
											   cufft_iodim64{get<1>(elems), get<2>(elems), get<3>(elems)}
				}
				 ...,
				std::pair<bool, cufft_iodim64>{}
			};
		},
									   boost::multi::detail::tuple_zip(which, sizes_tuple, istride_tuple, ostride_tuple));

		std::stable_sort(which_iodims.begin(), which_iodims.end() - 1, [](auto const& alpha, auto const& omega) { return get<1>(alpha).is > get<1>(omega).is; });

		auto const part = std::stable_partition(which_iodims.begin(), which_iodims.end() - 1, [](auto elem) { return std::get<0>(elem); });

		std::array<cufft_iodim64, D> dims{};
		auto const                   dims_end = std::transform(which_iodims.begin(), part, dims.begin(), [](auto elem) { return elem.second; });

		// std::array<cufftw_iodim64, D> howmany_dims{};
		// auto const howmany_dims_end = std::transform(part, which_iodims.end() -1, howmany_dims.begin(), [](auto elem) {return elem.second;});

		which_iodims_  = which_iodims;
		first_howmany_ = part - which_iodims.begin();

		////////////////////////////////////////////////////////////////////////

		std::array<int, D> istrides{};
		std::array<int, D> ostrides{};
		std::array<int, D> ion{};

		auto const istrides_end = std::transform(dims.begin(), dims_end, istrides.begin(), [](auto elem) { return static_cast<int>(elem.is); });
		auto const ostrides_end = std::transform(dims.begin(), dims_end, ostrides.begin(), [](auto elem) { return static_cast<int>(elem.os); });
		auto const ion_end      = std::transform(dims.begin(), dims_end, ion.begin(), [](auto elem) { return static_cast<int>(elem.n); });

		int  istride = *(istrides_end - 1);
		auto inembed = istrides;
		inembed.fill(0);
		int  ostride = *(ostrides_end - 1);
		auto onembed = ostrides;
		onembed.fill(0);

		for(std::ptrdiff_t idx = 1; idx != ion_end - ion.begin(); ++idx) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) replace with algorithm
			assert(ostrides[idx - 1] >= ostrides[idx]);
			assert(ostrides[idx - 1] % ostrides[idx] == 0);
			onembed[idx] = ostrides[idx - 1] / ostrides[idx];
			assert(istrides[idx - 1] % istrides[idx] == 0);
			inembed[idx] = istrides[idx - 1] / istrides[idx];
		}

		if(dims_end == dims.begin()) {
			throw std::runtime_error{"no ffts in any dimension is not supported"};
		}

		while(first_howmany_ < D - 1) {  // NOLINT(altera-id-dependent-backward-branch) TODO(correaa) replace with algorithm
			int nelems = 1;

			for(int idx = first_howmany_ + 1; idx != D; ++idx) {
				nelems *= which_iodims_[idx].second.n;
			}  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) replace with algorithm
			if(
				which_iodims_[first_howmany_].second.is == nelems && which_iodims_[first_howmany_].second.os == nelems
			) {
				which_iodims_[first_howmany_ + 1].second.n *= which_iodims_[first_howmany_].second.n;
				++first_howmany_;
			} else {
				break;
			}
		}

		if(first_howmany_ == D) {
			if constexpr(std::is_same_v<Alloc, void*>) {
				assert(dims_end - dims.begin() < 4);  // cufft cannot do 4D FFT
				cufftSafeCall(::cufftPlanMany(
					/*cufftHandle *plan*/ &h_,
					/*int rank*/ dims_end - dims.begin(),
					/*int *n*/ ion.data(),
					/*int *inembed*/ inembed.data(),
					/*int istride*/ istride,
					/*int idist*/ 1,  // stride(first),
					/*int *onembed*/ onembed.data(),
					/*int ostride*/ ostride,
					/*int odist*/ 1,  // stride(d_first),
					/*cufftType type*/ CUFFT_Z2Z,
					/*int batch*/ 1  // BATCH
				));
			} else {
				cufftSafeCall(cufftCreate(&h_));
				cufftSafeCall(cufftSetAutoAllocation(h_, false));
				cufftSafeCall(cufftMakePlanMany(
					/*cufftHandle *plan*/ h_,
					/*int rank*/ dims_end - dims.begin(),
					/*int *n*/ ion.data(),
					/*int *inembed*/ inembed.data(),
					/*int istride*/ istride,
					/*int idist*/ 1,  // stride(first),
					/*int *onembed*/ onembed.data(),
					/*int ostride*/ ostride,
					/*int odist*/ 1,  // stride(d_first),
					/*cufftType type*/ CUFFT_Z2Z,
					/*int batch*/ 1,  // BATCH
					/*size_t **/ &workSize_
				));
				cufftSafeCall(cufftGetSize(h_, &workSize_));
				workArea_ = ::thrust::raw_pointer_cast(alloc_.allocate(workSize_));
				static_assert(sizeof(Alloc) == 1000);
				// auto s = cudaMalloc(&workArea_, workSize_);
				// if(s != cudaSuccess) {throw std::runtime_error{"L212"};}
				cufftSafeCall(cufftSetWorkArea(h_, workArea_));
			}
			if(!h_) {
				throw std::runtime_error{"cufftPlanMany null"};
			}
			return;
		}

		std::sort(which_iodims_.begin() + first_howmany_, which_iodims_.begin() + D, [](auto const& alpha, auto const& omega) { return get<1>(alpha).n > get<1>(omega).n; });

		if(first_howmany_ <= D - 1) {
			if constexpr(std::is_same_v<Alloc, void*>) {  // NOLINT(bugprone-branch-clone) workaround bug in DeepSource
				cufftSafeCall(::cufftPlanMany(
					/*cufftHandle *plan*/ &h_,
					/*int rank*/ dims_end - dims.begin(),
					/*int *n*/ ion.data(),
					/*int *inembed*/ inembed.data(),
					/*int istride*/ istride,
					/*int idist*/ which_iodims_[first_howmany_].second.is,
					/*int *onembed*/ onembed.data(),
					/*int ostride*/ ostride,
					/*int odist*/ which_iodims_[first_howmany_].second.os,
					/*cufftType type*/ CUFFT_Z2Z,
					/*int batch*/ which_iodims_[first_howmany_].second.n
				));
			} else {
				cufftSafeCall(cufftCreate(&h_));
				cufftSafeCall(cufftSetAutoAllocation(h_, false));
				cufftSafeCall(cufftMakePlanMany(
					/*cufftHandle *plan*/ h_,
					/*int rank*/ dims_end - dims.begin(),
					/*int *n*/ ion.data(),
					/*int *inembed*/ inembed.data(),
					/*int istride*/ istride,
					/*int idist*/ which_iodims_[first_howmany_].second.is,
					/*int *onembed*/ onembed.data(),
					/*int ostride*/ ostride,
					/*int odist*/ which_iodims_[first_howmany_].second.os,
					/*cufftType type*/ CUFFT_Z2Z,
					/*int batch*/ which_iodims_[first_howmany_].second.n,
					/*size_t **/ &workSize_
				));
				cufftSafeCall(cufftGetSize(h_, &workSize_));
				workArea_ = ::thrust::raw_pointer_cast(alloc_.allocate(workSize_));
				cufftSafeCall(cufftSetWorkArea(h_, workArea_));
			}
			if(!h_) {
				throw std::runtime_error{"cufftPlanMany null"};
			}
			++first_howmany_;
			return;
		}
		// throw std::runtime_error{"cufft not implemented yet"};
	}

 private:
	template<typename = void>
	void ExecZ2Z_(complex_type const* idata, complex_type* odata, int direction) const {
		// used_ = true;
		cufftSafeCall(cufftExecZ2Z(h_, const_cast<complex_type*>(idata), odata, direction));  // NOLINT(cppcoreguidelines-pro-type-const-cast) wrap legacy interface
																							  // cudaDeviceSynchronize();
	}

 public:
	template<class IPtr, class OPtr>
	auto execute(IPtr idata, OPtr odata, int direction) const
		-> decltype((void)(reinterpret_cast<complex_type const*>(::thrust::raw_pointer_cast(idata)),
						   reinterpret_cast<complex_type*>(::thrust::raw_pointer_cast(odata)))) {  // TODO(correaa) make const
		// used_ = true;
		if(first_howmany_ == DD) {
			ExecZ2Z_(reinterpret_cast<complex_type const*>(::thrust::raw_pointer_cast(idata)), reinterpret_cast<complex_type*>(::thrust::raw_pointer_cast(odata)), direction);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) wrap a legacy interface
			return;
		}
		if(first_howmany_ == DD - 1) {
			if(which_iodims_[first_howmany_].first) {
				throw std::runtime_error{"logic error"};
			}

			// std::cout << "doing an inefficient loop of " << which_iodims_[first_howmany_].second.n << std::endl;

			for(int idx = 0; idx != which_iodims_[first_howmany_].second.n; ++idx) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
				cufftExecZ2Z(
					h_,
					const_cast<complex_type*>(reinterpret_cast<complex_type const*>(::thrust::raw_pointer_cast(idata + idx * which_iodims_[first_howmany_].second.is))),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
					reinterpret_cast<complex_type*>(::thrust::raw_pointer_cast(odata + idx * which_iodims_[first_howmany_].second.os)),                                   // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
					direction
				);
			}
			return;
		}
		if(first_howmany_ == DD - 2) {
			if(which_iodims_[first_howmany_ + 0].first) {
				throw std::runtime_error{"logic error0"};
			}
			if(which_iodims_[first_howmany_ + 1].first) {
				throw std::runtime_error{"logic error1"};
			}

			// std::cout << "doing an inefficient 2D loop of " << which_iodims_[first_howmany_].second.n << " " <<  which_iodims_[first_howmany_ + 1].second.n << std::endl;

			for(int idx = 0; idx != which_iodims_[first_howmany_].second.n; ++idx) {          // NOLINT(altera-unroll-loops,altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use an algorithm
				for(int jdx = 0; jdx != which_iodims_[first_howmany_ + 1].second.n; ++jdx) {  // NOLINT(altera-unroll-loops,altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use an algorithm
					cufftExecZ2Z(
						h_,
						const_cast<complex_type*>(reinterpret_cast<complex_type const*>(::thrust::raw_pointer_cast(idata + idx * which_iodims_[first_howmany_].second.is + jdx * which_iodims_[first_howmany_ + 1].second.is))),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
						reinterpret_cast<complex_type*>(::thrust::raw_pointer_cast(odata + idx * which_iodims_[first_howmany_].second.os + jdx * which_iodims_[first_howmany_ + 1].second.os)),                                   // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
						direction
					);
				}
			}
			return;
		}
		throw std::runtime_error{"error2"};
	}

	template<class IPtr, class OPtr>
	void execute_forward(IPtr idata, OPtr odata) {  // TODO(correaa) make const
		execute(idata, odata, cufft::forward);
	}
	template<class IPtr, class OPtr>
	void execute_backward(IPtr idata, OPtr odata) {  // TODO(correaa) make const
		execute(idata, odata, cufft::backward);
	}

	template<class IPtr, class OPtr>
	void operator()(IPtr idata, OPtr odata, int direction) const {
		// used_ = true;
		ExecZ2Z_(reinterpret_cast<complex_type const*>(::thrust::raw_pointer_cast(idata)), reinterpret_cast<complex_type*>(::thrust::raw_pointer_cast(odata)), direction);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
	}
	template<class I, class O>
	auto execute_dft(I&& in, O&& out, int direction) const -> O&& {
		// used_ = true;
		ExecZ2Z_(
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(in))),   // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) legay interface
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(out))),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) legay interface
			direction
		);
		return std::forward<O>(out);
	}

	~plan() {
		if constexpr(!std::is_same_v<Alloc, void*>) {
			if(workSize_ > 0) {
				alloc_.deallocate(typename std::allocator_traits<Alloc>::pointer(reinterpret_cast<char*>(workArea_)), workSize_);
			}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) legacy interface
		}
		if(h_ != 0) {
			cufftSafeCall(cufftDestroy(h_));
		}
		// if(!used_) {
		//  std::cerr <<"Warning: cufft plan was never used\n";
		//  std::terminate();
		// }
	}

	using size_type  = int;
	using ssize_type = int;
};

template<dimensionality_type D, class Alloc = void*>
class cached_plan {
	typename std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc>>::iterator it_;

 public:
	cached_plan(cached_plan const&) = delete;
	cached_plan(cached_plan&&)      = delete;

	auto operator=(cached_plan const&) -> cached_plan& = delete;
	auto operator=(cached_plan&&) -> cached_plan&      = delete;

	~cached_plan() = default;

	cached_plan(std::array<bool, D> which, boost::multi::layout_t<D, boost::multi::ssize_t> in, boost::multi::layout_t<D, boost::multi::ssize_t> out, Alloc const& alloc = {}) {  // NOLINT(fuchsia-default-arguments-declarations)
		thread_local std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc>>& LEAKY_cache = *new std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc>>;
		it_                                                                                                                         = LEAKY_cache.find(std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>{which, in, out});
		if(it_ == LEAKY_cache.end()) {
			it_ = LEAKY_cache.insert(std::make_pair(std::make_tuple(which, in, out), plan<D, Alloc>(which, in, out, alloc))).first;
		}
	}
	template<class IPtr, class OPtr>
	auto execute(IPtr idata, OPtr odata, int direction)
		-> decltype((void)(std::declval<
							   typename std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc>>::iterator&>()
							   ->second.execute(idata, odata, direction))) {
		// assert(it_ != LEAKY_cache.end());
		it_->second.execute(idata, odata, direction);
	}
};

// template<typename In, class Out, dimensionality_type D = In::rank::value, std::enable_if_t<!multi::has_get_allocator<In>::value, int> =0, typename = decltype(::thrust::raw_pointer_cast(std::declval<In const&>().base()))>
// auto dft(std::array<bool, +D> which, In const& in, Out&& out, int sgn)
// ->decltype(cufft::cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), sgn), std::forward<Out>(out)) {
// 	return cufft::cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), sgn), std::forward<Out>(out);
// }

template<typename In, class Out, dimensionality_type D = In::dimensionality>  // , std::enable_if_t<    multi::has_get_allocator<In>::value, int> =0, typename = decltype(raw_pointer_cast(std::declval<In const&>().base()))>
auto dft(std::array<bool, +D> which, In const& in, Out&& out, int sgn)
	-> decltype(cufft::cached_plan<D /*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/>{which, in.layout(), out.layout() /*, i.get_allocator()*/}.execute(in.base(), out.base(), sgn), std::forward<Out>(out)) {
	if constexpr(D == 4) {
		if(which == std::array<bool, D>{true, true, true, true}) {
			auto const [is, js, ks, ls] = in.extents();
			for(auto i : is)
				for(auto j : js) {
					cufft::dft({true, true}, in[i][j], out[i][j], sgn);
				}
			for(auto k : ks)
				for(auto l : ls) {
					cufft::dft({true, true}, out.rotated().rotated()[k][l], out.rotated().rotated()[k][l], sgn);
				}
			return std::forward<Out>(out);
		}
	}
	return cufft::cached_plan<D /*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/>{which, in.layout(), out.layout() /*, i.get_allocator()*/}.execute(in.base(), out.base(), sgn), std::forward<Out>(out);
}

template<typename In, class Out, dimensionality_type D = In::dimensionality>  //, std::enable_if_t<not multi::has_get_allocator<In>::value, int> =0>
auto dft_forward(std::array<bool, +D> which, In const& in, Out&& out) -> Out&& {
	//->decltype(cufft::plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::forward), std::forward<Out>(o)) {
	return cufft::cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), cufft::forward), std::forward<Out>(out);
}

template<typename In, class Out, dimensionality_type D = In::dimensionality>  //, std::enable_if_t<not multi::has_get_allocator<In>::value, int> =0>
auto dft_backward(std::array<bool, +D> which, In const& in, Out&& out) -> Out&& {
	//->decltype(cufft::plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o)) {
	return cufft::cached_plan<D>{which, in.layout(), out.layout()}.execute(in.base(), out.base(), cufft::backward), std::forward<Out>(out);
}

// template<typename In, class Out, dimensionality_type D = In::rank::value, class = typename In::allocator_type, std::enable_if_t<    multi::has_get_allocator<In>::value, int> =0>
// auto dft_backward(std::array<bool, +D> which, In const& i, Out&& o) -> Out&& {
// //->decltype(cufft::plan<D, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char> >{which, i.layout(), o.layout(), i.get_allocator()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o)) {
//  return cufft::cached_plan<D/*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/>{which, i.layout(), o.layout()/*, i.get_allocator()*/}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o); }

template<typename In, typename R = multi::array<typename In::element, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
BOOST_MULTI_NODISCARD("when first argument is const")
auto dft(In const& in, int sgn) -> R {
	static_assert(std::is_trivially_default_constructible<typename In::element>{});
	R ret(extents(in), get_allocator(in));
	cufft::dft(in, ret, sgn);
	// if(cudaDeviceSynchronize() != cudaSuccess) throw std::runtime_error{"Cuda error: Failed to synchronize"};
	return ret;
}

template<class Array, std::size_t... Ns>
constexpr auto array_tail_impl(Array const& arr, std::index_sequence<Ns...> /*unused*/) {
	return std::array<typename Array::value_type, std::tuple_size<Array>{} - 1>{std::get<Ns + 1>(arr)...};
}

template<class Array>
constexpr auto array_tail(Array const& arr)
	-> decltype(array_tail_impl(arr, std::make_index_sequence<std::tuple_size<Array>{} - 1>())) {
	return array_tail_impl(arr, std::make_index_sequence<std::tuple_size<Array>{} - 1>());
}

template<typename In, std::size_t D = In::dimensionality>
BOOST_MULTI_NODISCARD("when passing a const argument")
auto dft(std::array<bool, D> which, In const& in, int sign) -> std::decay_t<decltype(dft(which, in, typename In::decay_type(extents(in), get_allocator(in)), sign))> { return dft(which, in, typename In::decay_type(extents(in), get_allocator(in)), sign); }

template<typename In, std::size_t D = In::dimensionality>  // TODO(correaa) check that the type of In a decay_type (otherwise there is no ::dimensionality)
auto dft(std::array<bool, D> which, In&& in, int sign)
	-> decltype(dft(which, in, in, sign), std::forward<In>(in)) {
	return dft(which, in, in, sign), std::forward<In>(in);
}

template<typename Array, typename A> BOOST_MULTI_NODISCARD("when passing a const argument")
auto dft_forward(Array arr, A const& in)
	-> decltype(cufft::dft(arr, in, cufft::forward)) {
	return cufft::dft(arr, in, cufft::forward);
}

// template<typename Array, dimensionality_type D> NODISCARD("when passing a const argument")
// auto dft_forward(Array arr, multi::cuda::array<std::complex<double>, D>&& a)
// ->decltype(cufft::dft(arr, a, cufft::forward), multi::cuda::array<std::complex<double>, D>{}){//assert(0);
//  return cufft::dft(arr, a, cufft::forward), std::move(a);}

template<typename A> BOOST_MULTI_NODISCARD("when passing a const argument")
auto dft_forward(A const& arr)
	-> decltype(cufft::dft(arr, cufft::forward)) {
	return cufft::dft(arr, cufft::forward);
}

template<typename... As> auto dft_backward(As&&... as)
	-> decltype(cufft::dft(std::forward<As>(as)..., cufft::backward)) {
	return cufft::dft(std::forward<As>(as)..., cufft::backward);
}

template<typename Array, typename A> BOOST_MULTI_NODISCARD("when passing a const argument")
auto dft_backward(Array arr, A const& in)
	-> decltype(cufft::dft(arr, in, cufft::backward)) {
	return cufft::dft(arr, in, cufft::backward);
}

template<typename A> BOOST_MULTI_NODISCARD("when passing a const argument")
auto dft_backward(A const& arr)
	-> decltype(cufft::dft(arr, cufft::backward)) {
	return cufft::dft(arr, cufft::backward);
}

}  // end namespace boost::multi::cufft

#endif
