// Copyright 2020-2024 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_FFT_HPP
#define BOOST_MULTI_ADAPTORS_FFT_HPP

#include "../adaptors/fftw.hpp"

#if defined(__CUDA__) || defined(__NVCC__)
#include "../adaptors/cufft.hpp"
#elif defined(__HIPCC__)
#include "../adaptors/hipfft.hpp"
#endif

#define BOOST_MULTI_DECLRETURN_(ExpR) -> decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost {
namespace multi {
namespace fft {

	static inline constexpr int forward = static_cast<int>(fftw::forward);
	static inline constexpr int none = static_cast<int>(fftw::none);
	static inline constexpr int backward = static_cast<int>(fftw::backward);

	static_assert( forward != none and none != backward and backward != forward );

	template<std::size_t I> struct priority : std::conditional_t<I==0, std::true_type, struct priority<I-1>>{};

	template<class... Args> auto dft_aux_(priority<0>, Args&&... args) BOOST_MULTI_DECLRETURN_(  fftw::dft_backward(std::forward<Args>(args)...))
	template<class... Args> auto dft_aux_(priority<1>, Args&&... args) BOOST_MULTI_DECLRETURN_(cufft ::dft_backward(std::forward<Args>(args)...))
	template<class... Args> auto dft(Args&&... args) BOOST_MULTI_DECLRETURN_(dft_backward_aux_(priority<1>{}, std::forward<Args>(args)...))
	template<class In, class... Args> auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In const& in, Args&&... args) -> decltype(auto) {return dft_aux_(priority<1>{}, which, in, std::forward<Args>(args)...);}

	template<class... Args> auto dft_forward_aux_(priority<0>, Args&&... args) BOOST_MULTI_DECLRETURN_(  fftw::dft_forward(std::forward<Args>(args)...))
	template<class... Args> auto dft_forward_aux_(priority<1>, Args&&... args) BOOST_MULTI_DECLRETURN_(cufft ::dft_forward(std::forward<Args>(args)...))
	template<class In, class... Args> auto dft_forward(std::array<bool, In::dimensionality> which, In const& in, Args&&... args) -> decltype(auto) {return dft_forward_aux_(priority<1>{}, which, in, std::forward<Args>(args)...);}

	template<class... Args> auto dft_backward_aux_(priority<0>, Args&&... args) BOOST_MULTI_DECLRETURN_(  fftw::dft_backward(std::forward<Args>(args)...))
	template<class... Args> auto dft_backward_aux_(priority<1>, Args&&... args) BOOST_MULTI_DECLRETURN_(cufft ::dft_backward(std::forward<Args>(args)...))
	template<class In, class... Args> auto dft_backward(std::array<bool, In::dimensionality> which, In const& in, Args&&... args) -> decltype(auto) {return dft_backward_aux_(priority<1>{}, which, in, std::forward<Args>(args)...);}

}}}

#undef BOOST_MULTI_DECLRETURN_

#endif  // BOOST_MULTI_ADAPTORS_FFT_HPP
