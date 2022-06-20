#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0x$OXX `pkg-config --cflags --libs fftw3 cuda-11.0` -lboost_timer -lboost_unit_test_framework&&$0x$OXX&&rm $0x$OXX;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP

#include "../adaptors/../array.hpp"
#include "../adaptors/../config/NODISCARD.hpp"


#include "../detail/tuple_zip.hpp"

#include<algorithm> // sort
#include<complex>
#include<numeric> // accumulate

#if HAVE_FFTW3_THREADS
#include <thread>
#endif

#include<fftw3.h> // external fftw3 library

namespace boost{
namespace multi{

namespace fftw{
//	template<class T> auto alignment_of(T* p){return ::fftw_alignment_of((double*)p);}
#if  __cpp_lib_as_const >= 201510
using std::as_const;
#else
template<class T> constexpr std::add_const_t<T>& as_const(T& t) noexcept{return t;}
#endif

struct flags{
	using underlying_type = decltype(FFTW_PRESERVE_INPUT); // NOLINT(hicpp-signed-bitwise) : macro definition in external library
private:
	underlying_type underlying_;
public:
	constexpr explicit flags(underlying_type underlying) : underlying_{underlying}{}
	constexpr explicit operator underlying_type() const{return underlying_;}
	friend constexpr auto operator|(flags f1, flags f2){return flags{f1.underlying_ | f2.underlying_};}
};


constexpr flags estimate      {FFTW_ESTIMATE      }; // NOLINT(hicpp-signed-bitwise) : defined in an external lib 1U << 6
constexpr flags preserve_input{FFTW_PRESERVE_INPUT}; // NOLINT(hicpp-signed-bitwise) : defined in an external lib 1U << 4

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

template<class T, class Tpl> constexpr auto to_array(Tpl const& t){
	return detail::to_array_impl<T>(t, std::make_index_sequence<std::tuple_size<Tpl>{}>{});
}

#if 0
#if(__cpp_if_constexpr>=201606)
//https://stackoverflow.com/a/35110453/225186
template<class T> constexpr auto _constx(T&&t) -> std::remove_reference_t<T>{return t;}
#define logic_assert(C, M) \
	if constexpr(noexcept(_constx(C))) static_assert((C), M); else assert((C)&&(M));
#else
#define logic_assert(ConditioN, MessagE) assert(ConditioN && MessagE);
#endif
#endif

template<typename It1, class It2, std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0
>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign, fftw::flags flags)
->decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double>*>*/(base(d_first))), fftw_plan{}){ // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code

	static_assert( sizeof(*base(  first)) == sizeof(real(*base(  first))) + sizeof(imag(*base(  first))) and sizeof(*base(  first)) == sizeof(fftw_complex), 
		"input  must have complex pod layout" );
	static_assert( sizeof(*base(d_first)) == sizeof(real(*base(d_first))) + sizeof(imag(*base(d_first))) and sizeof(*base(d_first)) == sizeof(fftw_complex), 
		"output must have complex pod layout");

	assert(strides(*first) == strides(*last));
	assert(sizes(*first)==sizes(*d_first));

	auto const ssn_tuple = multi::detail::tuple_zip(strides(*first  ), strides(*d_first), sizes(*first));
	auto ssn = std::apply([](auto... e){return std::array<std::tuple<int, int, int>, sizeof...(e)>{
		std::make_tuple(static_cast<int>(std::get<0>(e)), static_cast<int>(std::get<1>(e)), static_cast<int>(std::get<2>(e)))...
	};}, ssn_tuple);
	std::sort(ssn.begin(), ssn.end(), std::greater<>{});

	auto const istrides = [&](){
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> istrides{};
		std::transform(ssn.begin(), ssn.end(), istrides.begin(), [](auto e){return std::get<0>(e);});
		return istrides;
	}();

	auto const ostrides = [&](){
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> ostrides{};
		std::transform(ssn.begin(), ssn.end(), ostrides.begin(), [](auto e){return std::get<1>(e);});
		return ostrides;
	}();
	assert( std::is_sorted(ostrides.begin(), ostrides.end(), std::greater<>{}) ); // otherwise ordering is incompatible

	auto const ion      = [&](){
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value> ion     {};
		std::transform(ssn.begin(), ssn.end(), ion     .begin(), [](auto e){return std::get<2>(e);});
		return ion;
	}();

	auto const inembed = [&](){
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value + 1> inembed{};
		std::adjacent_difference(
			istrides.rbegin(), istrides.rend(), inembed.rbegin(), [](auto a, auto b){assert(b != 0 and a%b == 0); return a/b;}
		);
		return inembed;
	}();

	auto const onembed = [&](){
		std::array<int, std::decay_t<decltype(*It1{})>::rank::value + 1> onembed{};
		std::adjacent_difference(
			ostrides.rbegin(), ostrides.rend(), onembed.rbegin(), [](auto a, auto b){assert(b != 0 and a%b == 0); return a/b;}
		);
		return onembed;
	}();

	auto ret = ::fftw_plan_many_dft(
		/*int rank*/ ion.size(),
		/*const int* n*/ ion.data(),
		/*int howmany*/ last - first,
		/*fftw_complex * in */ reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(static_cast<std::complex<double> const*>(base(first)))),   // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) input data
		/*const int *inembed*/ inembed.data(),
		/*int istride*/ istrides.back(),
		/*int idist*/ stride(first),
		/*fftw_complex * out */ reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(d_first))), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) adapt types
		/*const int *onembed*/ onembed.data(),
		/*int ostride*/ ostrides.back(),
		/*int odist*/ stride(d_first),
		/*int*/ sign, /*unsigned*/ static_cast<unsigned>(flags)
	);
	assert(ret); // if you get null here it could be because your library doesn't support this fftw call mode
	return ret;
}

template<typename It1, class It2, std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0
>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign)
->decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double>*>*/(base(d_first))), fftw_plan{}){ // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
	return fftw_plan_many_dft(first, last, d_first, sign, fftw::estimate);
}

template<
	class In, class Out, dimensionality_type D = std::decay_t<In>::rank_v,
	class=std::enable_if_t<D==std::decay_t<Out>::rank_v>,
	class=decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(base(std::declval<Out&>())))) // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
>
auto fftw_plan_dft(std::array<bool, +D> which, In&& in, Out&& out, int sign, fftw::flags flags) -> fftw_plan{
	static_assert( sizeof(*base(in )) == sizeof((*base(in )).real()) + sizeof((*base(in)).imag()) and sizeof(*base(in)) == sizeof(fftw_complex), 
		"input must have complex pod layout" );
	static_assert( sizeof(*base(out)) == sizeof((*base(out)).real()) + sizeof((*base(in)).imag()) and sizeof(*base(out)) == sizeof(fftw_complex), 
		"output must have complex pod layout" );

	assert(in.sizes() == out.sizes());

	auto const sizes_tuple   = in.sizes();
	auto const istride_tuple = in.strides();
	auto const ostride_tuple = out.strides();

	auto which_iodims = std::apply([](auto... e){
		return std::array<std::pair<bool, fftw_iodim64>, sizeof...(e)>{
			std::pair<bool, fftw_iodim64>{
				std::get<0>(e),
				fftw_iodim64{std::get<1>(e), std::get<2>(e), std::get<3>(e)}
			}...
		};
	}, boost::multi::detail::tuple_zip(which, sizes_tuple, istride_tuple, ostride_tuple));
	auto p = std::stable_partition(which_iodims.begin(), which_iodims.end(), [](auto e){return std::get<0>(e);});

	std::array<fftw_iodim64, D> dims{};
	auto const dims_end         = std::transform(which_iodims.begin(), p,         dims.begin(), [](auto e){return e.second;});

	std::array<fftw_iodim64, D> howmany_dims{};
	auto const howmany_dims_end = std::transform(p, which_iodims.end()  , howmany_dims.begin(), [](auto e){return e.second;});

//	auto const ion      = std::apply([](auto...e){return std::array{e...};}, in .sizes  ());
//	auto const istrides = std::apply([](auto...e){return std::array{e...};}, in .strides());
//	auto const ostrides = std::apply([](auto...e){return std::array{e...};}, out.strides());

//	std::array<fftw_iodim64, D> dims{};
//	auto l_dims = dims.begin();

//	std::array<fftw_iodim64, D> howmany{};
//	auto l_howmany = howmany.begin();

//	std::for_each(
//		which.begin(), which.end(), 
//		[&, i = 0](auto e) mutable{
//			if(e){
//				*l_dims    = {ion[i], istrides[i], ostrides[i]};
//				++l_dims;
//			}else{
//				*l_howmany = {ion[i], istrides[i], ostrides[i]};
//				++l_howmany;
//			}
//			++i;
//		}
//	);
//	for(int i=0; i != D; ++i){
//		if(which[i]){
//			*l_dims    = {ion[i], istrides[i], ostrides[i]};
//			++l_dims;
//		}else{
//			*l_howmany = {ion[i], istrides[i], ostrides[i]};
//			++l_howmany;
//		}
//	}

	assert( in .base() );
	assert( out.base() );

	assert( in.extensions() == out.extensions() );

	assert( (sign == -1) or (sign == +1) );

	fftw_plan ret = fftw_plan_guru64_dft(
		/*int rank*/ dims_end - dims.begin(), //p - which_iodims.begin(), //l_dims - dims.begin(),
		/*const fftw_iodim64 *dims*/ dims.data(),
		/*int howmany_rank*/ howmany_dims_end - howmany_dims.begin(), //which_iodims.endl_howmany - howmany.begin(),
		/*const fftw_iodim *howmany_dims*/ howmany_dims.data(), //howmany.data(), //nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(/*static_cast<std::complex<double> const *>*/(in.base()))), // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast) FFTW is taken as non-const while it is really not touched
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(out.base())), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
		sign, static_cast<unsigned>(flags) // | FFTW_ESTIMATE
	);

	assert(ret &&"fftw lib returned a null plan, if you are using MKL check the limitations of their fftw interface"); 
	//https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/appendix-d-fftw-interface-to-intel-math-kernel-library/fftw3-interface-to-intel-math-kernel-library/using-fftw3-wrappers.html
	return ret;
}

template<
	class In, class Out, dimensionality_type D = std::decay_t<In>::rank_v,
	class=std::enable_if_t<D==std::decay_t<Out>::rank_v>,
	class=decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(base(std::declval<Out&>())))) // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
>
auto fftw_plan_dft(std::array<bool, +D> which, In&& in, Out&& out, int sign) -> fftw_plan{
	return fftw_plan_dft(which, std::forward<In>(in), std::forward<Out>(out), sign, fftw::estimate);
}

template<class In, class Out, dimensionality_type D = In::rank_v, typename = decltype(reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(base(std::declval<Out&>()))))> // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
auto fftw_plan_dft(In const& in, Out&& out, int s, fftw::flags flags){
	static_assert( D == std::decay_t<Out>::rank_v , "!");
	using multi::sizes; using multi::strides; assert(sizes(in) == sizes(out));

	assert( in.sizes() == out.sizes() );

	auto const dims = std::apply([](auto... e){
		return std::array<fftw_iodim64, sizeof...(e)>{
			fftw_iodim64{std::get<0>(e), std::get<1>(e), std::get<2>(e)}
			...
		};
	}, boost::multi::detail::tuple_zip(in.sizes(), in.strides(), out.strides()));

	auto ret = fftw_plan_guru64_dft(
		/*int rank*/ s?D:0,
		/*const fftw_iodim64 *dims*/ dims.data(),
		/*int howmany_rank*/ 0,
		/*const fftw_iodim *howmany_dims*/ nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(base(in)))), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) : interact with legacy code
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(base(out))), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
		s, static_cast<unsigned>(flags)
	);
	assert(ret);
	return ret;
}

template<class In, class Out, dimensionality_type D = In::rank_v, typename = decltype(reinterpret_cast<fftw_complex*>(multi::implicit_cast<std::complex<double>*>(base(std::declval<Out&>()))))> // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : interact with legacy code
auto fftw_plan_dft(In const& in, Out&& out, int s){
	return fftw_plan_dft(in, out, s, fftw::estimate);
}


namespace fftw{

#if HAVE_FFTW3_THREADS
inline void initialize_threads(){int good = fftw_init_threads(); assert(good); (void)good;}
#else
inline void initialize_threads(){}
#endif

inline void cleanup(){fftw_cleanup();}

struct environment{
	environment() = default;
	environment(environment const&) = delete;
	environment(environment&&) = delete;
	auto operator=(environment const&) = delete;
	auto operator=(environment&&) = delete;
	~environment(){fftw_cleanup();}
};

class plan{
	plan() : impl_{nullptr, &fftw_destroy_plan}{}
	std::unique_ptr<std::remove_pointer_t<fftw_plan>, decltype(&fftw_destroy_plan)> impl_;
public:
	plan(plan const&) = delete;//default;
	plan(plan&&) = default;
	~plan() = default;
	template<typename... As, 
		typename = decltype(fftw_plan_dft(std::declval<As&&>()...))
	>
	explicit plan(As&&... as) : impl_{fftw_plan_dft(std::forward<As>(as)...), &fftw_destroy_plan}{
		assert(impl_);
	}
	template<typename... As>
	static auto many(As&&... as)
	->std::decay_t<decltype(fftw_plan_many_dft(std::forward<As>(as)...) , std::declval<plan>())>
	{
		plan r; r.impl_.reset(fftw_plan_many_dft(std::forward<As>(as)...)); return r; // this produces a compilation error in icc++17
	}

private:
	void execute() const{fftw_execute(impl_.get());} //TODO(correaa): remove const
	template<class I, class O>
	void execute_dft(I&& i, O&& o) const{
		::fftw_execute_dft(impl_.get(), const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(base(i)))), reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(o)))); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) : to interface with legacy fftw
	}
	template<class I, class O> void execute(I&& i, O&& o) const{execute_dft(std::forward<I>(i), std::forward<O>(o));}
	friend void execute(plan const& p){p.execute();}

public:
	auto operator=(plan&&) -> plan& = default;
	auto operator=(plan const&) -> plan& = delete;

	template<class I, class O> 
	void operator()(I&& i, O&& o) const{execute(std::forward<I>(i), std::forward<O>(o));}
	void operator()()             const{execute();} // http://www.fftw.org/fftw3_doc/Thread-safety.html#Thread-safety

	[[nodiscard]] auto cost() const -> double{return fftw_cost(impl_.get());}
	[[nodiscard]] auto flops() const{
		struct ret_t{
			double add = 0.;
			double mul = 0.;
			double fma = 0.;
		//	explicit operator double() const{return add + mul + 2*fma;}
		} ret{};
		fftw_flops(impl_.get(), &ret.add, &ret.mul, &ret.fma);
		return ret;
	}

	//std::string string_print() const{
	//	return std::unique_ptr<char>{fftw_sprint_plan(impl_.get())}.get();
	//}
	//friend std::ostream& operator<<(std::ostream& os, plan const& p){return os<<p.string_print()<<'\n';}
#if HAVE_FFTW3_THREADS
public:
	static void make_thread_safe(){
		fftw_make_planner_thread_safe(); // needs linking to -lfftw3_threads, requires FFTW-3.3.6 or greater
		is_thread_safe_ = true;
	}
	static int with_nthreads(int n){fftw_plan_with_nthreads(n); nthreads_ = n; return n;}
	static int with_nthreads(){
		int n=std::thread::hardware_concurrency(); return with_nthreads(n?n:2);
	}
	static bool is_thread_safe(){return is_thread_safe_;}
	static bool nthreads(){return nthreads_;}
private:
	static bool is_thread_safe_;
	static int nthreads_;
	static bool initialized_threads_;
#else
	static constexpr auto is_thread_safe() -> bool{return false;}
	static constexpr auto nthreads() -> bool{return true;}
	static constexpr auto with_nthreads() -> int{return 1;}
#endif
};

#if HAVE_FFTW3_THREADS
bool plan::is_thread_safe_ = (plan::make_thread_safe(), true);
int plan::nthreads_ = (initialize_threads(), with_nthreads());
#endif

using sign = int;
constexpr sign forward = FFTW_FORWARD;
constexpr sign none = 0;
constexpr sign backward = FFTW_BACKWARD;

static_assert( forward != none and none != backward and backward != forward, "!");

//enum strategy: decltype(FFTW_ESTIMATE){ estimate = FFTW_ESTIMATE, measure = FFTW_MEASURE };

template<class In, class Out>
auto dft(In const& i, Out&& o, int s)
->decltype(fftw::plan{i, o, s}(), std::forward<Out>(o)){
	return fftw::plan{i, o, s}(), std::forward<Out>(o);}

using std::decay_t;

template<class In, class Out, std::size_t D=In::rank_v>
auto dft(std::array<bool, +D> which, In const& i, Out&& o, sign s)
->decltype(plan{which, i, o, s}(), std::forward<Out>(o)){
	return plan{which, i, o, s}(), std::forward<Out>(o);}

template<typename In, class Out, dimensionality_type D=In::rank_v, dimensionality_type=std::decay_t<Out>::rank_v>
auto dft(std::array<sign, +D> w, In const& i, Out&& o){
//	std::array<bool, D> non;

	std::array<bool, D> fwd{};
	std::transform(begin(w), end(w), begin(fwd), [](auto e){return e==FFTW_FORWARD;});
	dft(fwd, i, o, fftw::forward);

	std::array<bool, D> bwd{};
	std::transform(begin(w), end(w), begin(bwd), [](auto e){return e==FFTW_BACKWARD;}); 
	if(std::accumulate(begin(bwd), end(bwd), false)){dft(bwd, o, o, FFTW_BACKWARD);}

	return std::forward<Out>(o);
}

template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, int sign)
->decltype(plan::many(first, last, d_first, sign)(), d_first + (last - first)){
	return plan::many(first, last, d_first, sign)(), d_first + (last - first);}

template<typename In, class R=typename In::decay_type>
NODISCARD("when first argument is const")
auto dft(In const& i, sign s)
->std::decay_t<decltype(dft(i, R(extensions(i), get_allocator(i)), s))>{
	return dft(i, R(extensions(i), get_allocator(i)), s);}

template<typename T, dimensionality_type D, class... Args>
auto rotate(multi::array<T, D, Args...>& i) -> decltype(auto){
	multi::array_ref<T, D, typename multi::array<T, D, Args...>::element_ptr> before(data_elements(i), extensions(i));
	i.reshape(extensions(rotated(before) ));
	fftw::dft(before, i, fftw::none);
	return i;
}

template<typename In, dimensionality_type D = In::rank_v, class R=typename In::decay_type>
NODISCARD("when first argument is const")
auto dft(std::array<bool, +D> which, In const& i, sign s)
->std::decay_t<decltype(fftw::dft(which, i, R(extensions(i), get_allocator(i)), s))>{
	return fftw::dft(which, i, R(extensions(i), get_allocator(i)), s);}

template<typename In, multi::dimensionality_type D = std::decay_t<In>::rank_v>
auto dft(std::array<bool, +D> which, In&& i, sign s)
->decltype(dft(which, i, i, s), std::forward<In>(i)){
	return dft(which, i, i, s), std::forward<In>(i);}

template<typename In, std::size_t D = In::rank_v, class R=typename In::decay_type>
void dft(std::array<bool, +D> which, In const& i) = delete;

template<dimensionality_type Rank /*not deduced*/, typename In, class R=typename In::decay_type>
NODISCARD("when second argument is const")
auto dft(In const& i, sign s) -> R{
	static_assert( Rank <= In::rank_v, "!" );
	return dft<Rank>(i, R(extensions(i), get_allocator(i)), s);
}

template<typename... A> auto            dft_forward(A&&... a)
->decltype(fftw::dft(std::forward<A>(a)..., fftw::forward)){
	return fftw::dft(std::forward<A>(a)..., fftw::forward);}

template<typename BoolArray, typename A>
NODISCARD("when input argument is read only")
auto dft_forward(BoolArray which, A const& a)
->decltype(fftw::dft(which, a, fftw::forward)){
	return fftw::dft(which, a, fftw::forward);}

template<class A, multi::dimensionality_type D = A::rank_v>
NODISCARD("when input argument is read only")
auto dft_forward(std::array<bool, +D> which, A const& a)
->decltype(fftw::dft(which, a, fftw::forward)){
	return fftw::dft(which, a, fftw::forward);}

template<class A, class O, multi::dimensionality_type D = A::rank_v>
auto dft_forward(std::array<bool, +D> which, A const& a, O&& o)
->decltype(fftw::dft(which, a, std::forward<O>(o), fftw::forward)){
	return fftw::dft(which, a, std::forward<O>(o), fftw::forward);}

template<typename A>
NODISCARD("when input argument is read only") auto dft_forward(A const& a)
->decltype(fftw::dft(a, fftw::forward)){
	return fftw::dft(a, fftw::forward);}

template<typename... A> auto            dft_backward(A&&... a)
->decltype(dft(std::forward<A>(a)..., fftw::backward)){
	return dft(std::forward<A>(a)..., fftw::backward);}

template<class In> auto dft_inplace(In&& i, sign s) -> In&&{
	fftw::plan{i, i, static_cast<int>(s)}();//(i, i); 
	return std::forward<In>(i);
}

template<class In, class Out, dimensionality_type D = In::rank_v>
auto copy(In const& i, Out&& o)
->decltype(dft(std::array<bool, D>{}, i, std::forward<Out>(o), fftw::forward)){
	return dft(std::array<bool, D>{}, i, std::forward<Out>(o), fftw::forward);}

template<typename In, class R=typename In::decay_type>
NODISCARD("when argument is const")
auto copy(In const& i) -> R
{//->decltype(copy(i, R(extensions(i), get_allocator(i))), R()){
	return copy(i, R(extensions(i), get_allocator(i)));}

template<typename In, class R=typename std::decay_t<In>::decay_type>
auto move(In&& in){
	if(in.is_compact()){
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

template<typename T, dimensionality_type D, class P, class R=typename multi::array<T, D>>
auto copy(multi::basic_array<T, D, multi::move_ptr<T, P>>&& a) -> R{
	if(a.is_compact()){
		return 
			fftw::copy(
				a.template static_array_cast<T, T*>(), 
				multi::array_ref<T, D, T*>(a.base().base(), a.extensions())
			).template static_array_cast<T, multi::move_ptr<T>>()
		;
	}
	return fftw::copy(a.template static_array_cast<T, P>());
}

template<class Array>
auto transpose(Array& a)
->decltype(fftw::copy(transposed(a), a.reshape(extensions(layout(a).transpose())))){
	multi::array_ref<typename Array::element, Array::rank_v, typename Array::element_ptr> r(a.base(), extensions(a));
	return fftw::copy(r.transposed(), a.reshape(layout(a).transpose().extensions()));
}


#if 0
// TODO investigate why this doesn't work as expected
template<class Array>
auto rotate(Array& a)
->decltype(fftw::copy(rotated(a), a.reshape(extensions(layout(a).transpose())))){
	multi::array_ref<typename Array::element, Array::dimensionality, typename Array::element_ptr> r(a.base(), extensions(a));
	auto&& ro = r.rotated();
	return fftw::copy(ro, a.reshape(layout(a).rotate().extensions()));
}
#endif

} // end namespace fftw
} // end namespace multi
} // end namespace boost

#endif

