// Copyright 2020-2025 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_FFT_HPP
#define BOOST_MULTI_ADAPTORS_FFT_HPP

#include "boost/multi/adaptors/fftw.hpp"
#include "boost/multi/array.hpp"  // for extents_t, get, array

#if defined(__CUDA__) || defined(__NVCC__)
#include "boost/multi/adaptors/cufft.hpp"
#elif defined(__HIPCC__)
#include "boost/multi/adaptors/hipfft.hpp"
#endif

#include <array>        // for array
#include <cstddef>      // for size_t
#include <iterator>     // for random_access_iterator_tag
#include <tuple>        // for apply
#include <type_traits>  // for decay_t, conditional_t, true_type
#include <utility>      // for forward
// IWYU pragma: no_include <variant>                        // for get

#define BOOST_MULTI_DECLRETURN_(ExpR) \
	->decltype(ExpR) { return ExpR; }  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define BOOST_MULTI_JUSTRETURN_(ExpR) \
	->decltype(auto) { return ExpR; }  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost::multi::fft {

static inline int const forward  = static_cast<int>(fftw::forward);
static inline int const none     = static_cast<int>(fftw::none);
static inline int const backward = static_cast<int>(fftw::backward);

// static_assert( forward != none && none != backward && backward != forward );

template<std::size_t I> struct priority : std::conditional_t<I == 0, std::true_type, struct priority<I - 1>> {};

template<class... Args> auto dft_aux(priority<0> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(fftw::dft(std::forward<Args>(args)...))
#if defined(__CUDA__) || defined(__NVCC__) || defined(__HIPCC__)
	template<class... Args>
	auto dft_aux(priority<1> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(::boost::multi::cufft ::dft(std::forward<Args>(args)...))
#endif
		template<class... Args>
		auto dft(Args&&... args) BOOST_MULTI_DECLRETURN_(dft_aux_(priority<1>{}, std::forward<Args>(args)...))

			template<class In, class... Args>
			auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In&& in, Args&&... args) BOOST_MULTI_DECLRETURN_(dft_aux(priority<1>{}, which, std::forward<In>(in), std::forward<Args>(args)...))

				template<class... Args>
				auto dft_forward_aux(priority<0> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(fftw::dft_forward(std::forward<Args>(args)...))
#if defined(__CUDA__) || defined(__NVCC__) || defined(__HIPCC__)
					template<class... Args>
					auto dft_forward_aux(priority<1> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(cufft ::dft_forward(std::forward<Args>(args)...))
#endif
						template<class In, class... Args>
						auto dft_forward(std::array<bool, std::decay_t<In>::dimensionality> which, In&& in, Args&&... args) BOOST_MULTI_DECLRETURN_(dft_forward_aux(priority<1>{}, which, std::forward<In>(in), std::forward<Args>(args)...))

							template<class... Args>
							auto dft_backward_aux(priority<0> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(fftw::dft_backward(std::forward<Args>(args)...))
#if defined(__CUDA__) || defined(__NVCC__) || defined(__HIPCC__)
								template<class... Args>
								auto dft_backward_aux(priority<1> /*unused*/, Args&&... args) BOOST_MULTI_DECLRETURN_(cufft ::dft_backward(std::forward<Args>(args)...))
#endif
									template<class In, class... Args>
									auto dft_backward(std::array<bool, In::dimensionality> which, In const& in, Args&&... args) -> decltype(auto) {
	return dft_backward_aux(priority<1>{}, which, in, std::forward<Args>(args)...);
}

template<class In, class Direction>
class dft_range {
 public:
	static constexpr auto dimensionality = std::decay_t<In>::dimensionality;
	auto                  operator+() const {
        multi::array<typename std::decay_t<In>::element, dimensionality> ret = *this;
        return ret;
	}

 private:
	std::array<bool, dimensionality> which_;
	In                               in_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
	Direction                        dir_;

	struct const_iterator : private std::decay_t<In>::const_iterator {  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
		// static constexpr auto dimensionality = In::const_iterator::dimensionality;

	 private:
		bool                                 do_;
		std::array<bool, dimensionality - 1> sub_which_;
		Direction                            dir_;

	 public:
		const_iterator() = default;  // cppcheck-suppress [uninitMemberVar];
		const_iterator(
			typename std::decay_t<In>::const_iterator it,
			bool doo, std::array<bool, std::decay_t<In>::dimensionality - 1> sub_which,
			Direction dir
		) : std::decay_t<In>::const_iterator{it}, do_{doo}, sub_which_{sub_which}, dir_{dir} {}

		using typename std::decay_t<In>::const_iterator::difference_type;
		using typename std::decay_t<In>::const_iterator::value_type;
		using pointer           = void;  // void*;
		using reference         = dft_range<typename std::decay_t<In>::const_iterator::reference, Direction>;
		using iterator_category = std::random_access_iterator_tag;

		auto        operator+(difference_type n) const { return const_iterator{static_cast<typename std::decay_t<In>::const_iterator const&>(*this) + n, do_, sub_which_, dir_}; }
		friend auto operator-(const_iterator const& lhs, const_iterator const& rhs) {
			return static_cast<typename std::decay_t<In>::const_iterator const&>(lhs) - static_cast<typename std::decay_t<In>::const_iterator const&>(rhs);
		}

		auto operator++() -> const_iterator&;
		auto operator--() -> const_iterator&;

		auto operator++(int) -> const_iterator;
		auto operator--(int) -> const_iterator;

		auto operator==(const_iterator const& other) const -> bool;
		auto operator!=(const_iterator const& other) const -> bool;

		// auto operator*() const -> reference;

		class fake_array {
			multi::extents_t<dimensionality - 1> extensions_;

			public:
			explicit fake_array(multi::extents_t<dimensionality - 1> ext) : extensions_{ext} {}
			[[deprecated]] auto extensions() const { return extensions_; }
			auto extents() const { return extensions_; }

			//  multi::size_t size_;
			[[deprecated]] auto extension() const {
				using std::get;
				return get<0>(extents());
			}
			[[nodiscard]] auto extent() const {
				using std::get;
				return get<0>(extents());
			}

			auto size() const { return extent().size(); }
		};

		auto operator*() const -> fake_array {
			fake_array fa{(*static_cast<typename std::decay_t<In>::const_iterator const&>(*this)).extents()};
			return fa;
		}

	 private:
		template<class It>
		auto copy_(const_iterator const& last, It const& first_d) const -> decltype(auto) {
			auto const count = last - *this;
			dft(
				std::apply([doo = do_](auto... es) { return std::array{doo, es...}; }, sub_which_),
				multi::const_subarray<typename std::decay_t<In>::const_iterator::element, dimensionality, typename std::decay_t<In>::const_iterator::element_ptr>(
					static_cast<typename std::decay_t<In>::const_iterator const&>(*this),
					static_cast<typename std::decay_t<In>::const_iterator const&>(last)
				),
				multi::subarray<typename std::decay_t<In>::element, dimensionality, typename It::element_ptr>(
					first_d, first_d + count
				),
				fftw::sign{dir_}
			);
			return first_d + count;
		}

	 public:
		template<class It>
		auto capy(const_iterator const& last, It const& first_d) const -> decltype(auto) {
			return copy_(last, first_d);
		}

		template<class It>
		friend auto copy(const_iterator const& first, const_iterator const& last, It const& first_d) -> decltype(auto) {
			return first.copy_(last, first_d);
		}

		template<class It, class Size>
		friend auto copy_n(const_iterator const& first, Size const& count, It const& first_d) -> decltype(auto) {
			return first.copy_(first + count, first_d);
		}

		template<class It, class Size>
		friend auto uninitialized_copy_n(const_iterator const& first, Size const& count, It const& first_d) -> decltype(auto) {
			return copy_n(first, count, first_d);
		}

		template<class It>
		friend auto uninitialized_copy(const_iterator const& first, const_iterator const& last, It const& first_d) -> decltype(auto) {
			return first.copy_(last, first_d);
		}
	};

 public:
	template<class In2>
	dft_range(std::array<bool, std::decay_t<In>::dimensionality> which, In2&& in, Direction dir) : which_{which}, in_(std::forward<In2>(in)), dir_{dir} {}
	auto begin() const {
		return const_iterator(in_.begin(), which_[0], std::apply([](auto /*e0*/, auto... es) { return std::array<bool, dimensionality - 1>{es...}; }, which_), dir_);
	}
	auto end() const {
		return const_iterator(in_.end(), which_[0], std::apply([](auto /*e0*/, auto... es) { return std::array<bool, dimensionality - 1>{es...}; }, which_), dir_);
	}

	auto extensions() const { return in_.extensions(); }
	auto extents() const { return in_.extents(); }

	auto size() const { return in_.size(); }
};

template<class In, class Direction>
auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In&& in, Direction dir) {
	return dft_range<In&&, Direction>(which, std::forward<In>(in), dir);
}

template<class In>
auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In&& in) {
	return dft(which, std::forward<In>(in), fft::forward);
}

// #if defined(__clang__)
// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wunused-value"
// #elif defined(__GNUC__)
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wunused-value"
// #endif
template<class In>
auto dft_all(In&& in) {
	auto const all_true = std::apply([](auto... es) { return std::array{((void)es, true)...}; }, std::array<bool, std::decay_t<In>::dimensionality>{});
	return dft(all_true, std::forward<In>(in), fft::forward);
}

template<class In>
auto idft_all(In&& in) {
	auto const all_true = std::apply([](auto... es) { return std::array{(es, true)...}; }, std::array<bool, std::decay_t<In>::dimensionality>{});
	return dft(all_true, std::forward<In>(in), fft::backward);
}
// #if defined(__clang__)
// #pragma clang diagnostic pop
// #elif defined(__GNUC__)
// #pragma GCC diagnostic pop
// #endif
template<class In, class Direction>
auto idft(std::array<bool, In::dimensionality> which, In&& in) {
	return dft(which, std::forward<In>(in), fft::forward);
}

template<class In>
auto dft_forward(std::array<bool, std::decay_t<In>::dimensionality> which, In&& in) {
	return dft(which, std::forward<In>(in), fft::forward);
}

template<class In>
auto dft_backward(std::array<bool, In::dimensionality> which, In&& in) {
	return dft(which, std::forward<In>(in), fft::backward);
}

}  // end namespace boost::multi::fft

#undef BOOST_MULTI_DECLRETURN_
#undef BOOST_MULTI_JUSTRETURN_

#endif  // BOOST_MULTI_ADAPTORS_FFT_HPP
