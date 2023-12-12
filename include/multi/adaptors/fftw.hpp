// Copyright 2018-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP
#pragma once

#include <multi/array.hpp>

#include <multi/adaptors/fftw/memory.hpp>

#include <algorithm>  // sort
#include <chrono>
#include <complex>
#include <numeric>  // accumulate

#if HAVE_FFTW3_THREADS
#include <thread>
#endif

#include <fftw3.h>  // external fftw3 library

namespace boost::multi {
namespace fftw {

using std::as_const;

struct flags {
	using underlying_type = decltype(FFTW_PRESERVE_INPUT);  // NOLINT(hicpp-signed-bitwise) : macro definition in external library

 private:
	underlying_type underlying_;

 public:
	constexpr explicit flags(underlying_type underlying) : underlying_{underlying} {}
	constexpr explicit    operator underlying_type() const { return underlying_; }
	friend constexpr auto operator|(flags f1, flags f2) { return flags{f1.underlying_ | f2.underlying_}; }
};

constexpr flags estimate{FFTW_ESTIMATE};  // NOLINT(hicpp-signed-bitwise) : defined in an external lib 1U << 6
constexpr flags measure{FFTW_MEASURE};

constexpr flags preserve_input{FFTW_PRESERVE_INPUT};  // NOLINT(hicpp-signed-bitwise) : defined in an external lib 1U << 4
// // NOLINT(): this is a defect in FFTW https://github.com/FFTW/fftw3/issues/246

}  // end namespace fftw

#if 0
template<typename Size>
auto fftw_plan_dft_1d(
	Size N, 
	std::complex<double> const* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
#ifndef NDEBUG
	auto check = in[N/3]; // check that const data will not been overwritten
#endif
	assert( fftw::alignment_of(in) == fftw::alignment_of(out) );
	auto ret=::fftw_plan_dft_1d(N, (fftw_complex*)in, (fftw_complex*)out, sign, flags | FFTW_PRESERVE_INPUT );
	assert(check == in[N/3]); // check that const data has not been overwritten
	return ret;
}

template<typename Size>
auto fftw_plan_dft_1d(
	Size N, 
	std::complex<double>* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert( fftw::alignment_of(in) == fftw::alignment_of(out) );
	return ::fftw_plan_dft_1d(N, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_2d(
	Size N1, Size N2, 
	std::complex<double> const* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert( fftw::alignment_of(in) == fftw::alignment_of(out) );
#ifndef NDEBUG
	auto check = in[N1*N2/3]; // check that const data will not been overwritten
#endif
	auto ret = ::fftw_plan_dft_2d(N1, N2, (fftw_complex*)in, (fftw_complex*)out, sign, flags | FFTW_PRESERVE_INPUT);
	assert( check == in[N1*N2/3] ); // check that const data has not been overwritten
	return ret;
}

template<typename Size>
auto fftw_plan_dft_2d(
	Size N1, Size N2, 
	std::complex<double>* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_2d(N1, N2, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_3d(
	Size N1, Size N2, Size N3, 
	std::complex<double>* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename Size>
auto fftw_plan_dft_3d(
	Size N1, Size N2, Size N3, 
	std::complex<double> const* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags | FFTW_PRESERVE_INPUT);
}
#endif

#if 0
template<typename Rank>
auto fftw_plan_dft(
	Rank r, int* ns, 
	std::complex<double>* in, std::complex<double>* out, 
	int sign, unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft(r, ns, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename RankType>
auto fftw_plan_dft(
	RankType r, int* ns, 
	std::complex<double> const* in, std::complex<double>* out, 
	int sign, unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
){
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw::alignment_of(in) == fftw::alignment_of(out));
#ifndef NDEBUG
	size_t ne = 1; for(RankType i = 0; i != r; ++i) ne*=ns[i];
	auto check = in[ne/3]; // check that const data will not been overwritten
#endif
	auto ret=::fftw_plan_dft(r, ns, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
	assert(check == in[ne/3]); // check that const data has not been overwritten
	return ret;
}
#endif

#if 0
template<typename In, typename Out>
auto fftw_plan_dft_1d(
	In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE
){
	static_assert(in.dimensionality == 1, "!"); assert(size(in) == size(out));
	assert( in.is_compact() ); assert( out.is_compact() );
	return multi::fftw_plan_dft_1d(size(in), data_elements(in), data_elements(out), sign, flags);
}

template<class In, class Out>
auto fftw_plan_dft_2d(
	In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE
){
	static_assert(in.dimensionality == 2, "!"); assert(in.sizes() == out.sizes());
	assert( in.is_compact() ); assert( out.is_compact() );
	return multi::fftw_plan_dft_2d(
		sizes(in)[0], sizes(in)[1], 
		data_elements(in), data_elements(out), sign, flags
	);
}

template<class In, class Out>
auto fftw_plan_dft_3d(
	In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE
){
	static_assert(in.dimensionality == 3, "!"); assert(in.sizes() == out.sizes());
	assert( in.is_compact() ); assert( out.is_compact() );
	return multi::fftw_plan_dft_3d(
		sizes(in)[0], sizes(in)[1], sizes(in)[2],
		data(in), data(out),
		sign, flags
	);
}
#endif

#if 0
template<class T, class Tpl>
constexpr auto to_array(Tpl const& tpl) {
	return std::apply(
		[](auto const&... elems) { return std::array<T, std::tuple_size<Tpl>::value>{static_cast<T>(elems)...}; },
		tpl
	);
}

template<
	typename It1, class It2,
	std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign, fftw::flags flags)
	-> fftw_plan {

	static_assert(sizeof(*base(first)) == sizeof((*base(first)).real()) + sizeof((*base(first)).imag()), "input  must have complex pod layout");
	static_assert(sizeof(*base(first)) == sizeof(fftw_complex), "input  must have complex pod layout");
	static_assert(sizeof(*base(d_first)) == sizeof((*base(d_first)).real()) + sizeof((*base(d_first)).imag()), "output must have complex pod layout");
	static_assert(sizeof(*base(d_first)) == sizeof(fftw_complex), "output must have complex pod layout");

	assert(strides(*first) == strides(*last));
	assert(sizes(*first) == sizes(*d_first));

	auto const ssn_tuple = multi::detail::tuple_zip(strides(*first), strides(*d_first), sizes(*first));
	auto       ssn       = std::apply([](auto... ssn) {
                                                                                                                                                                                                      using boost::multi::detail::get;
                                                                                                                                                                                                      return std::array<boost::multi::detail::tuple<int, int, int>, sizeof...(ssn)>{
                                                                                                                                                                                                                                                                                                         boost::multi::detail::mk_tuple(static_cast<int>(get<0>(ssn)), static_cast<int>(get<1>(ssn)), static_cast<int>(get<2>(ssn)))...};
                                                                                                   },
                                                                                                                         ssn_tuple);
	std::sort(ssn.begin(), ssn.end(), std::greater<>{});

	auto const istrides = [&]() {
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> istrides{};
		using boost::multi::detail::get;
		std::transform(ssn.begin(), ssn.end(), istrides.begin(), [](auto elem) { return get<0>(elem); });
		return istrides;
	}();

	auto const ostrides = [&]() {
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> ostrides{};
		using boost::multi::detail::get;
		std::transform(ssn.begin(), ssn.end(), ostrides.begin(), [](auto elem) { return get<1>(elem); });
		return ostrides;
	}();
	assert(std::is_sorted(ostrides.begin(), ostrides.end(), std::greater<>{}));  // otherwise ordering is incompatible

	auto const ion = [&]() {
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> ion{};
		using boost::multi::detail::get;
		std::transform(ssn.begin(), ssn.end(), ion.begin(), [](auto elem) { return get<2>(elem); });
		return ion;
	}();

	auto const inembed = [&]() {
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value + 1> inembed{};
		std::adjacent_difference(
			istrides.rbegin(), istrides.rend(), inembed.rbegin(), [](auto alpha, auto omega) {assert(omega != 0 and alpha%omega == 0); return alpha/omega; }
		);
		return inembed;
	}();

	auto const onembed = [&]() {
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value + 1> onembed{};
		std::adjacent_difference(
			ostrides.rbegin(), ostrides.rend(), onembed.rbegin(), [](auto alpha, auto omega) {assert(omega != 0 and alpha%omega == 0); return alpha/omega; }
		);
		return onembed;
	}();

	auto ret = ::fftw_plan_many_dft(
		/*int           rank    */ ion.size(),
		/*const int*    n       */ ion.data(),
		/*int           howmany */ last - first,
		/*fftw_complex* in      */ reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(static_cast<std::complex<double> const*>(base(first)))),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) input data
		/*const int*    inembed */ inembed.data(),
		/*int           istride */ istrides.back(),
		/*int           idist   */ stride(first),
		/*fftw_complex* out     */ reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(d_first))),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) adapt types
		/*const int*    onembed */ onembed.data(),
		/*int           ostride */ ostrides.back(),
		/*int           odist   */ stride(d_first),
		/*int                   */ sign,
		/*unsigned              */ static_cast<unsigned>(flags)
	);
	assert(ret);  // if you get null here it could be because your library doesn't support this fftw call mode
	return ret;
}

template<
	typename It1, class It2,
	std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign)
	-> fftw_plan {
	return fftw_plan_many_dft(first, last, d_first, sign, fftw::estimate);
}
#endif

template<class InPtr, class In, class OutPtr, class Out, dimensionality_type D = In::rank_v>
auto fftw_plan_dft(std::array<bool, +D> which, InPtr in_base, In const& in_layout, OutPtr out_base, Out&& out_layout, int sign, fftw::flags /*flags*/) -> fftw_plan {
	assert(in_layout.extensions() == out_layout.extensions());

	auto const sizes_tuple = in_layout.sizes();

	auto const istride_tuple = in_layout.strides();
	auto const ostride_tuple = out_layout.strides();

	using boost::multi::detail::get;
	auto which_iodims = std::apply(
		[](auto... elems) {
			return std::array<std::pair<bool, fftw_iodim64>, sizeof...(elems) + 1>{  // added one element to avoid problem with gcc 13 static analysis (out-of-bounds)
  // clang-format off
				std::pair{get<0>(elems), fftw_iodim64{get<1>(elems), get<2>(elems), get<3>(elems)}}..., {},  // added one element to avoid problem with gcc 13 static analysis (out-of-bounds)
  // clang-format on
			};
		},
		boost::multi::detail::tuple_zip(which, sizes_tuple, istride_tuple, ostride_tuple)
	);
	auto const part = std::stable_partition(which_iodims.begin(), which_iodims.end() - 1, [](auto elem) { return std::get<0>(elem); });

	std::array<fftw_iodim64, D> dims{};
	std::array<fftw_iodim64, D> howmany_dims{};

	auto const                  dims_end = std::transform(which_iodims.begin(), part, dims.begin(), [](auto elem) { return elem.second; });
	auto const                  howmany_dims_end = std::transform(part, which_iodims.end() - 1, howmany_dims.begin(), [](auto elem) { return elem.second; });

	assert(in_base);
	assert(out_base);

	assert((sign == -1) or (sign == +1));

	fftw_plan ret = fftw_plan_guru64_dft(
		/*int                 rank         */ dims_end - dims.begin(),
		/*const fftw_iodim64 *dims         */ dims.data(),
		/*int                 howmany_rank */ howmany_dims_end - howmany_dims.begin(),
		/*const fftw_iodim   *howmany_dims */ howmany_dims.data(),
		/*fftw_complex       *in           */ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(/*static_cast<std::complex<double> const *>*/ (in_base))),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) FFTW is taken as non-const while it is really not touched
		/*fftw_complex       *out          */ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(/*static_cast<std::complex<double>       *>*/ (out_base))),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast)
		sign, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
	);

	assert(ret && "fftw lib returned a null plan, if you are using MKL check the limitations of their fftw interface");
	// https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/appendix-d-fftw-interface-to-intel-math-kernel-library/fftw3-interface-to-intel-math-kernel-library/using-fftw3-wrappers.html
	return ret;
}

// template<
//  class In, class Out, dimensionality_type D = std::decay_t<In>::rank_v,
//  class = std::enable_if_t<D == std::decay_t<Out>::rank_v>//,
// //  class = decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/ (base(std::declval<Out&>()))))  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
//  >
// auto fftw_plan_dft(std::array<bool, +D> which, In&& in, Out&& out, int sign) -> fftw_plan {
//  return fftw_plan_dft(which, in.base(), in.layout(), out.base(), out.layout(), sign, fftw::estimate);
// }

// template<dimensionality_type D, class PtrIn, class PtrOut, typename = decltype(reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(std::declval<PtrOut&>())))>  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
// auto fftw_plan_dft(multi::layout_t<D> const& in_layout, PtrIn in_base, multi::layout_t<D> const& out_layout, PtrOut out_base, int dir, fftw::flags flags) {
//  using multi::sizes;
//  using multi::strides;

//  assert(in_layout.sizes() == out_layout.sizes());

//  auto const dims = std::apply([](auto... elems) {
//      using boost::multi::detail::get;
//      return std::array<fftw_iodim64, sizeof...(elems)>{
//          fftw_iodim64{get<0>(elems), get<1>(elems), get<2>(elems)}
//                                                                                                                                                                                                                                                                                                          ...
//                                                                                                                                                                                                       };
//  },
//                               boost::multi::detail::tuple_zip(in_layout.sizes(), in_layout.strides(), out_layout.strides()));

//  auto ret = fftw_plan_guru64_dft(
//      /*int                 rank         */ dir ? D : 0,
//      /*const fftw_iodim64 *dims         */ dims.data(),
//      /*int                 howmany_rank */ 0,
//      /*const fftw_iodim   *howmany_dims */ nullptr,  // howmany_dims.data(),
//      /*fftw_complex       *in           */ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(in_base))),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) : interact with legacy code
//      /*fftw_complex       *out          */ reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(out_base)),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
//      dir, static_cast<unsigned>(flags)
//  );
//  assert(ret);
//  return ret;
// }

// template<dimensionality_type D>  //, typename = decltype(reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(std::declval<PtrOut&>())))>  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
// auto fftw_plan_dft(multi::layout_t<D> const& in_layout, multi::layout_t<D> const& out_layout, int dir, fftw::flags flags) {
//  return fftw_plan_dft(in_layout, nullptr, out_layout, nullptr, dir, flags | fftw::estimate);
// }

// template<class In, class Out, dimensionality_type D = In::rank::value, typename = decltype(reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(base(std::declval<Out&>()))))>  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
// auto fftw_plan_dft(In const& in, Out&& out, int dir, fftw::flags flags) {
//  return fftw_plan_dft(in.layout(), in.base(), out.layout(), out.base(), dir, flags);
// }

template<class In, class Out, dimensionality_type D = In::rank_v, typename = decltype(reinterpret_cast<fftw_complex*>(detail::implicit_cast<std::complex<double>*>(base(std::declval<Out&>()))))>  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
auto fftw_plan_dft(In const& in, Out&& out, int dir) {
	return fftw_plan_dft(in, out, dir, fftw::estimate);
}

namespace fftw {

#if HAVE_FFTW3_THREADS
inline void initialize_threads() {
	int good = fftw_init_threads();
	assert(good);
	(void)good;
}
#else
inline void initialize_threads() {}
#endif

enum sign : decltype(FFTW_FORWARD) {
	backward = FFTW_BACKWARD,
	none     = 0,
	forward  = FFTW_FORWARD,
};

static_assert(forward != none and none != backward and backward != forward);

enum class direction : decltype(FFTW_FORWARD) {
	backward = FFTW_BACKWARD,
	none     = 0,
	forward  = FFTW_FORWARD,
};

class plan;

class environment {
	static void cleanup() { ::fftw_cleanup(); }
	static void set_timelimit(std::chrono::duration<double> limit) {
		::fftw_set_timelimit(limit.count());
	}
	static void unset_timelimit() { ::fftw_set_timelimit(FFTW_NO_TIMELIMIT); }

 public:
	environment()                      = default;
	environment(environment const&)    = delete;
	environment(environment&&)         = delete;
	auto operator=(environment const&) = delete;
	auto operator=(environment&&)      = delete;

	template<class In, class Out>
	auto make_plan_forward(std::array<bool, +In::rank_v> which, In const& in, Out&& out);
	template<class In, class Out>
	auto make_plan_backward(std::array<bool, +In::rank_v> which, In const& in, Out&& out);

	~environment() { cleanup(); }
};

class plan {
	plan() : impl_{nullptr, &fftw_destroy_plan} {}
	// std::shared_ptr<std::remove_pointer_t<fftw_plan> const> impl_;
	std::unique_ptr<std::remove_pointer_t<fftw_plan>, decltype(&fftw_destroy_plan)> impl_;

 public:
	plan(plan const&) = delete;
	plan(plan&&)      = delete;
	~plan()           = default;

	template<class InPtr, class In, class OutPtr, class Out>
	explicit plan(
		std::array<bool, In::rank_v> which,
		InPtr in_base, In in_layout,
		OutPtr out_base, Out out_layout, sign ss
	) : impl_{fftw_plan_dft(which, in_base, in_layout, out_base, out_layout, ss, fftw::estimate), &fftw_destroy_plan} {
		assert(impl_);
	}

	template<class InPtr, class In, class OutPtr, class Out>
	static auto forward(std::array<bool, In::rank_v> which, InPtr in_base, In in_layout, OutPtr out_base, Out out_layout) {
		return plan(which, in_base, in_layout, out_base, out_layout, fftw::forward);
	}
	template<class InPtr, class In, class OutPtr, class Out>
	static auto backward(std::array<bool, In::rank_v> which, InPtr in_base, In in_layout, OutPtr out_base, Out out_layout) {
		return plan(which, in_base, in_layout, out_base, out_layout, fftw::backward);
	}

	template<class I, class O>
	void execute(I* in, O* out) const {
		static_assert(sizeof(in->imag()) == sizeof(double));
		static_assert(sizeof(out->imag()) == sizeof(double));

		static_assert(sizeof(*in) == sizeof(fftw_complex));
		static_assert(sizeof(*out) == sizeof(fftw_complex));

		::fftw_execute_dft(
			const_cast<fftw_plan>(impl_.get()),  // NOLINT(cppcoreguidelines-pro-type-const-cast) https://www.fftw.org/fftw3_doc/Thread-safety.html
			const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(in)),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) : to interface with legacy fftw
			reinterpret_cast<fftw_complex*>(out)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) : to interface with legacy fftw
		);
	}

	// template<class I, class O> void execute(I&& in, O&& out) const { execute_dft(std::forward<I>(in), std::forward<O>(out)); }
	// friend void                     execute(plan const& self) { self.execute(); }

	auto operator=(plan&&) -> plan&      = delete;
	auto operator=(plan const&) -> plan& = delete;

	[[nodiscard]] auto cost() const -> double { return fftw_cost(const_cast<fftw_plan>(impl_.get())); }  // NOLINT(cppcoreguidelines-pro-type-const-cast)
	[[nodiscard]] auto flops() const {
		struct ret_t {
			double add = 0.0;
			double mul = 0.0;
			double fma = 0.0;
			//  explicit operator double() const{return add + mul + 2*fma;}
		} ret{};
		fftw_flops(const_cast<fftw_plan>(impl_.get()), &ret.add, &ret.mul, &ret.fma);  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		return ret;
	}

#if HAVE_FFTW3_THREADS
 public:
	static void make_thread_safe() {
		fftw_make_planner_thread_safe();  // needs linking to -lfftw3_threads, requires FFTW-3.3.6 or greater
		is_thread_safe_ = true;
	}
	static int with_nthreads(int n) {
		fftw_plan_with_nthreads(n);
		nthreads_ = n;
		return n;
	}
	static int with_nthreads() {
		int n = std::thread::hardware_concurrency();
		return with_nthreads(n ? n : 2);
	}
	static bool is_thread_safe() { return is_thread_safe_; }
	static bool nthreads() { return nthreads_; }

 private:
	static bool is_thread_safe_;
	static int  nthreads_;
	static bool initialized_threads_;
#else
	static constexpr auto is_thread_safe() -> bool { return false; }
	static constexpr auto nthreads() -> bool { return true; }
	static constexpr auto with_nthreads() -> int { return 1; }
#endif
};

template<class In, class Out>
auto environment::make_plan_forward(std::array<bool, +In::rank_v> which, In const& in, Out&& out) {
	return plan::forward(which, in, std::forward<Out>(out));
}

template<class In, class Out>
auto environment::make_plan_backward(std::array<bool, +In::rank_v> which, In const& in, Out&& out) {
	return plan::backward(which, in, std::forward<Out>(out));
}

#if HAVE_FFTW3_THREADS
bool plan::is_thread_safe_ = (plan::make_thread_safe(), true);
int  plan::nthreads_       = (initialize_threads(), with_nthreads());
#endif

// enum strategy: decltype(FFTW_ESTIMATE){ estimate = FFTW_ESTIMATE, measure = FFTW_MEASURE };

// template<class In, class Out, dimensionality_type = In::rank_v>
// auto dft(In const& in, Out&& out, int dir)
// ->decltype(fftw::plan{in.layout(), out.layout(), dir}.execute(in.base(), out.base()), std::forward<Out>(out)) {
//  return fftw::plan{in.layout(), out.layout(), dir}.execute(in.base(), out.base()), std::forward<Out>(out); }

using std::decay_t;

template<class In, class Out, std::size_t D = In::rank_v>
auto dft(std::array<bool, +D> which, In const& in, Out&& out, sign dir)
	-> decltype(plan{which, in.base(), in.layout(), out.base(), out.layout(), dir}.execute(in.base(), out.base()), std::forward<Out>(out)) {
	return plan{which, in.base(), in.layout(), out.base(), out.layout(), dir}.execute(in.base(), out.base()), std::forward<Out>(out);
}

// template<typename T, dimensionality_type D, class... Args>
// auto rotate(multi::array<T, D, Args...>& inout) -> decltype(auto) {
//  multi::array_ref<T, D, typename multi::array<T, D, Args...>::element_ptr> before(data_elements(inout), extensions(inout));
//  inout.reshape(extensions(rotated(before)));
//  fftw::dft(before, inout, fftw::none);
//  return inout;
// }

template<typename In, multi::dimensionality_type D = std::decay_t<In>::rank::value,
         std::enable_if_t<std::is_assignable_v<decltype(*std::declval<In&&>().base()), typename std::decay_t<In>::element>, int> = 0>
auto dft(std::array<bool, +D> which, In&& in, sign dir)
	-> decltype(dft(which, in, in, dir), std::forward<In>(in)) {
	return dft(which, in, in, dir), std::forward<In>(in);
}

// template<typename In, std::size_t D = In::rank::value, class R = typename In::decay_type>
// void dft(std::array<bool, +D> which, In const& in) = delete;

// template<dimensionality_type Rank /*not deduced*/, typename In, class R = typename In::decay_type>
// [[nodiscard]]  // ("when second argument is const")
// auto
// dft(In const& in, sign dir) -> R {
//  static_assert(Rank <= In::rank_v, "!");
//  return dft<Rank>(in, R(extensions(in), get_allocator(in)), dir);
// }

// template<typename... A> auto dft_forward(A&&... array)
//  -> decltype(fftw::dft(std::forward<A>(array)..., fftw::forward)) {
//  return fftw::dft(std::forward<A>(array)..., fftw::forward);
// }

// template<typename BoolArray, typename A>
// [[nodiscard]]  // ("when input argument is read only")
// auto
// dft_forward(BoolArray which, A const& array)
//  -> decltype(fftw::dft(which, array, fftw::forward)) {
//  return fftw::dft(which, array, fftw::forward);
// }

// template<class A, multi::dimensionality_type D = A::rank::value>
// [[nodiscard]]  // ("when input argument is read only")
// auto
// dft_forward(std::array<bool, +D> which, A const& array)
//  -> decltype(fftw::dft(which, array, fftw::forward)) {
//  return fftw::dft(which, array, fftw::forward);
// }

template<class A, class O, multi::dimensionality_type D = A::rank_v>
auto dft_forward(std::array<bool, +D> which, A const& in, O&& out)
	-> decltype(fftw::dft(which, in, std::forward<O>(out), fftw::forward)) {
	return fftw::dft(which, in, std::forward<O>(out), fftw::forward);
}

template<class A, class O, multi::dimensionality_type D = A::rank::value>
auto dft_backward(std::array<bool, +D> which, A const& in, O&& out)
	-> decltype(fftw::dft(which, in, std::forward<O>(out), fftw::backward)) {
	return fftw::dft(which, in, std::forward<O>(out), fftw::backward);
}

template<typename... A> auto dft_backward(A&&... args)
	-> decltype(dft(std::forward<A>(args)..., fftw::backward)) {
	return dft(std::forward<A>(args)..., fftw::backward);
}

#if 0
template<typename In, class R=typename std::decay_t<In>::decay_type>
auto move(In&& in) {
	if(in.is_compact()) {
		multi::array_ref<typename In::element, In::rank_v, typename In::element_ptr> Ref(
			in.base(), extensions(in)
		);
		copy(in, Ref);
		return R(
			multi::array_ref<typename In::element, In::rank_v, std::move_iterator<typename In::element_ptr>>(std::make_move_iterator(in.mbase()), ((in.mbase()=0), extensions(Ref)))
		);
	}
	return copy(std::forward<In>(in));
}
#endif

template<class T, boost::multi::dimensionality_type D>
using static_array = ::boost::multi::static_array<T, D, fftw::allocator<T>>;

template<class T, multi::dimensionality_type D>
using array = ::boost::multi::array<T, D, fftw::allocator<T>>;

template<typename T, dimensionality_type D, class P, class R = typename multi::fftw::array<T, D>>
auto copy(multi::subarray<T, D, multi::move_ptr<T, P>>&& array) -> R {
	if(array.is_compact()) {
		return fftw::copy(
			       array.template static_array_cast<T, T*>(),
			       multi::array_ref<T, D, T*>(array.base().base(), array.extensions())
		)
			.template static_array_cast<T, multi::move_ptr<T>>();
	}
	return fftw::copy(array.template static_array_cast<T, P>());
}

template<class Array>
auto transpose(Array& array)
	-> decltype(fftw::copy(transposed(array), array.reshape(extensions(layout(array).transpose())))) {
	multi::array_ref<typename Array::element, Array::rank::value, typename Array::element_ptr> const ref(array.base(), extensions(array));
	return fftw::copy(ref.transposed(), array.reshape(layout(array).transpose().extensions()));
}

#if 0
// TODO(correaa) investigate why this doesn't work as expected
template<class Array>
auto rotate(Array& a)
->decltype(fftw::copy(rotated(a), a.reshape(extensions(layout(a).transpose())))){
	multi::array_ref<typename Array::element, Array::dimensionality, typename Array::element_ptr> r(a.base(), extensions(a));
	auto&& ro = r.rotated();
	return fftw::copy(ro, a.reshape(layout(a).rotate().extensions()));
}
#endif

}  // end namespace fftw
}  // end namespace boost::multi

namespace boost::multi::fftw {

template<class MDIterator>
class fft_iterator {
	MDIterator                                      base_;
	std::array<fftw::sign, MDIterator::rank::value> which_ = {};

 public:
	using iterator_type = MDIterator;

	using difference_type = typename std::iterator_traits<MDIterator>::difference_type;
	using value_type      = typename std::iterator_traits<MDIterator>::value_type;
	using pointer         = void*;
	class reference {
		typename MDIterator::reference::extensions_type x_;
		explicit reference(typename MDIterator::reference const& ref) : x_{ref.extensions()} {}
		friend class fft_iterator;

	 public:
		using extensions_type = typename MDIterator::reference::extensions_type;
		auto extensions() const -> extensions_type { return x_; }
	};

	using iterator_category = std::random_access_iterator_tag;  // using iterator_category = std::input_iterator_tag;

	explicit fft_iterator(iterator_type base, std::array<fftw::sign, MDIterator::rank_v> which) noexcept : base_{std::move(base)}, which_{which} {}

	friend auto operator-(fft_iterator const& self, fft_iterator const& other) -> difference_type {
		return self.base_ - other.base_;
	}

	// template<class ItOut>
	// friend auto copy(fft_iterator first, fft_iterator last, ItOut d_first) {
	//  assert(first.which_ == last.which_);
	//  fftw::dft(
	//      first.which_,
	//      multi::ref(first.base_, last.base_),
	//      multi::ref(d_first, d_first + (last.base_ - first.base_))
	//  );

	//  return d_first + (last.base_ - first.base_);
	// }
	// template<class ItOut>
	// friend auto uninitialized_copy(fft_iterator first, fft_iterator last, ItOut d_first) {
	//  return copy(first, last, d_first);
	// }

	auto operator*() const { return reference{*base_}; }
};

// template<class Origin, class Array>
// class fft_range {
//  Origin origin_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
//  Array  ref_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
//  using which_type = std::array<fftw::sign, std::decay_t<Array>::rank::value>;
//  which_type which_;

//  public:
//  using iterator_type = typename std::decay_t<Array>::const_iterator;

//  using size_type = typename std::decay_t<Array>::size_type;
//  using iterator  = fft_iterator<iterator_type>;

//  using decay_type = typename std::decay_t<Array>::decay_type;

//  explicit fft_range(Origin&& origin, Array&& in, which_type which)
//  : origin_{std::forward<Origin>(origin)}, ref_{std::forward<Array>(in)}, which_{which} {}

//  operator decay_type() && {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
//      if constexpr(std::is_same_v<Origin, decay_type&&>) {
//          decay_type the_ret{std::forward<Origin>(origin_)};
//          the_ret.reshape(this->extensions());

//          fftw::dft(
//              which_,
//              ref_,
//              the_ret
//          );

//          return the_ret;
//      } else {
//          return decay_type{this->begin(), this->end()};
//      }
//  }

//  auto operator+() const& { return static_cast<decay_type>(*this); }
//  auto operator+() && { return static_cast<decay_type>(std::move(*this)); }

//  auto begin() const { return iterator{ref_.begin(), which_}; }
//  auto end() const { return iterator{ref_.end(), which_}; }

//  auto size() const { return ref_.size(); }
//  auto extensions() const { return ref_.extensions(); }
//  auto num_elements() const { return ref_.num_elements(); }

//  auto base() const { return ref_.base(); }

//  auto rotated() const {
//      auto new_which = which_;
//      std::rotate(new_which.begin(), new_which.begin() + 1, new_which.end());
//      return fft_range<Origin, std::decay_t<decltype(ref_.rotated())>>{origin_, ref_.rotated(), new_which};
//  }
//  auto unrotated() const {
//      auto new_which = which_;
//      std::rotate(new_which.rbegin(), new_which.rbegin() + 1, new_which.rend());
//      return fft_range<Origin, std::decay_t<decltype(ref_.unrotated())>>{origin_, ref_.unrotated(), new_which};
//  }
//  auto transposed() const {
//      auto new_which = which_;
//      std::swap(std::get<0>(new_which), std::get<1>(new_which));
//      return fft_range<Origin, std::decay_t<decltype(ref_.transposed())>>{std::forward<Origin>(origin_), ref_.transposed(), new_which};
//  }

//  template<class... FBNs>
//  auto operator()(FBNs... fbns) const {
//      static_assert(sizeof...(fbns) <= std::decay_t<Array>::rank::value, "too many arguments");
//      auto                                    new_which = which_;
//      std::array<fftw::sign, sizeof...(fbns)> fbna{fbns...};
//      std::transform(fbna.begin(), fbna.end(), new_which.begin(), new_which.begin(),
//                     [](auto fbn, auto nw) {
//                         if(fbn == fftw::none) {
//                             return nw;
//                         }
//                         assert(nw == fftw::none);
//                         return fbn;
//                     }
//      );
//      return fft_range<Origin, std::decay_t<decltype(ref_())>>{std::forward<Origin>(origin_), ref_(), new_which};
//  }
// };

// template<class Array>
// auto ref(Array&& in) {
//  return fft_range<Array&&, Array&&>{
//      std::forward<Array>(in),
//      std::forward<Array>(in),
//      {}};
// }

// template<class Array> auto move(Array& in) { return fftw::ref(std::move(in)); }

// template<dimensionality_type ND = 1, class Array>
// auto fft(Array&& in) {
//  std::array<fftw::sign, std::decay_t<Array>::rank::value> which{};
//  std::fill_n(which.begin(), ND, fftw::forward);
//  return fft_range<Array&&, Array&&>{std::forward<Array>(in), std::forward<Array>(in), which};
// }

// template<dimensionality_type ND = 1, class Array>
// auto ifft(Array&& in) {
//  std::array<fftw::sign, Array::rank::value> which{};
//  std::fill_n(which.begin(), ND, fftw::backward);
//  return fft_range{in, which};
// }

}  // end namespace boost::multi::fftw

#endif
