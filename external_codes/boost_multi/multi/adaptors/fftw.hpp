#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && time c++ -O3 -g -rdynamic -std=c++17 -Wall -Wextra -Wpedantic -D_TEST_MULTI_ADAPTORS_FFTW $0.cpp catch_main.o -o $0.x -lfftw3 && $0.x $@ && rm $0.x $0.cpp; exit
#endif
#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP
//  (C) Copyright Alfredo A. Correa 2018

#include<fftw3.h>

#include "../../multi/utility.hpp"
#include "../../multi/array.hpp"

#include<cmath>
#include<complex>
#include<memory>

namespace boost{
namespace multi{

namespace fftw{
	template<class T>
	auto alignment_of(T* p){return fftw_alignment_of((double*)p);}
}

namespace{
	using std::complex;
}

template<typename Size>
auto fftw_plan_dft_1d(
	Size N, 
	complex<double> const* in, complex<double>* out, int sign, 
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
	complex<double>* in, std::complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert( fftw::alignment_of(in) == fftw::alignment_of(out) );
	return ::fftw_plan_dft_1d(N, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_2d(
	Size N1, Size N2, 
	complex<double> const* in, complex<double>* out, int sign, 
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
	complex<double>* in, complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_2d(N1, N2, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_3d(
	Size N1, Size N2, Size N3, 
	complex<double>* in, complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename Size>
auto fftw_plan_dft_3d(
	Size N1, Size N2, Size N3, 
	complex<double> const* in, complex<double>* out, int sign, 
	unsigned flags = FFTW_ESTIMATE
){
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags | FFTW_PRESERVE_INPUT);
}

template<typename Rank>
auto fftw_plan_dft(
	Rank r, int* ns, 
	complex<double>* in, std::complex<double>* out, 
	int sign, unsigned flags = FFTW_ESTIMATE
){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft(r, ns, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename RankType>
auto fftw_plan_dft(
	RankType r, int* ns, 
	complex<double> const* in, complex<double>* out, 
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

//////////////////////////////////////////////////////////////////////////////

#if 0
template<class InIt, typename Size, class OutIt>
auto fftw_plan_dft_1d_n(InIt it, Size n, OutIt dest, decltype(FFTW_FORWARD) sign, decltype(FFTW_ESTIMATE) flags = FFTW_ESTIMATE){
	return fftw_plan_dft_1d(n, data(it), data(dest), sign, flags);
}
template<class InIt, class OutIt>
auto fftw_plan_dft_1d(InIt first, InIt last, OutIt dest, decltype(FFTW_FORWARD) sign, decltype(FFTW_ESTIMATE) flags = FFTW_ESTIMATE){
	using std::distance;
	return fftw_plan_dft_1d_n(first, distance(first, last), dest, sign, flags);
}
#endif

///////////////////////////////////////////////////////////////////////////////

template<typename In, typename Out>
auto fftw_plan_dft_1d(
	In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE
){
	static_assert(in.dimensionality == 1, "!"); assert(size(in) == size(out));
	return multi::fftw_plan_dft_1d(size(in), data_elements(in), data_elements(out), sign, flags);
}

template<class In, class Out>
auto fftw_plan_dft_2d(
	In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE
){
	static_assert(in.dimensionality == 2, "!"); assert(in.sizes() == out.sizes());
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
	return multi::fftw_plan_dft_3d(
		sizes(in)[0], sizes(in)[1], sizes(in)[2],
		data(in), data(out),
		sign, flags
	);
}

/*
namespace detail {
template<class T, class Tuple, std::size_t... I>
constexpr std::array<std::remove_cv_t<T>, std::tuple_size<Tuple>{}> 
to_array_impl(Tuple&& t, std::index_sequence<I...>){
	return { static_cast<std::remove_cv_t<T>>(std::get<I>(std::forward<Tuple>(t)))... };
}
}
*/

template<class T, class Tuple> constexpr auto to_array(Tuple&& t){
    return detail::to_array_impl<T>(std::forward<Tuple>(t), std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}

template<class ArrayIn, class ArrayOut>
auto fftw_plan_dft(
	ArrayIn&& in, ArrayOut&& out, 
	int sign, unsigned flags = FFTW_ESTIMATE
){
	using multi::dimensionality;
	static_assert(dimensionality(in) == dimensionality(out), "!");
	using multi::sizes;
	assert(sizes(in) == sizes(out));
	auto ns = to_array<int>(sizes(in));
	return fftw_plan_dft(
		dimensionality(in), ns.data(), 
		origin(std::forward<ArrayIn>(in)), origin(std::forward<ArrayOut>(out)), 
		sign, flags
	);
}


#if 0
template<class ItIn, class ItOut>
fftw_plan_many_dft(ItIn first, ItIn last, ItOut d_first, int sign, unsigned flags = FFTW_ESTIMATE){
	assert(first->dimensionality == last->dimensionality);
	assert(first->dimensionality == d_first->dimensionality);
	assert(first->extensions() == last->extensions());
	auto ns = to_array<int>(sizes(*first));
//	for(int i = 0; i != first->dimensionality; ++i) 
	return fftw_plan_many_dft(
		first->dimensionality, ns.data(), std::distance(first, last),
		first.data(), nullptr, 
		
	)
}
	ItIn->dimensionality, 
fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
                             fftw_complex *in, const int *inembed,
                             int istride, int idist,
                             fftw_complex *out, const int *onembed,
                             int ostride, int odist,
                             int sign, unsigned flags);
#endif

///////////////////////////////////////////////////////////////////////////////

void fftw_execute_dft(
	fftw_plan p, std::complex<double>* in, std::complex<double>* out){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out) );
	::fftw_execute_dft(p, (fftw_complex*)in, (fftw_complex*)out);
}

template<class ArrayIn, class ArrayOut>
auto fftw_execute_dft(fftw_plan p, ArrayIn&& in, ArrayOut&& out)
->decltype(multi::fftw_execute_dft(p, data_elements(in), data_elements(out))){
	return multi::fftw_execute_dft(p, data_elements(in), data_elements(out));}

namespace fftw{

struct plan{
	std::unique_ptr<std::pointer_traits<fftw_plan>::element_type, decltype(&fftw_destroy_plan)> impl_;
	template<typename... Args>
	plan(Args&&... args) : impl_{multi::fftw_plan_dft(std::forward<Args>(args)...), &fftw_destroy_plan}{}
	void execute() const{fftw_execute(impl_.get());}
	template<class I, class O>
	void execute_dft(I&& i, O&& o) const{fftw_execute_dft(std::forward<I>(i), std::forward<O>(o));}
	template<class I, class O>
	void execute(I&& i, O&& o) const{execute_dft(std::forward<I>(i), std::forward<O>(o));}
	void operator()() const{execute();}
	template<typename I, typename O>
	void operator()(I&& i, O&& o) const{return execute(std::forward<I>(i), std::forward<O>(o));}
	friend void execute(plan const& p){p.execute();}
	double cost() const{return fftw_cost(impl_.get());}
	auto flops() const{
		struct{double add; double mul; double fma;} r;
		fftw_flops(impl_.get(), &r.add, &r.mul, &r.fma);
		return r;
	}
};

enum sign : int {forward = FFTW_FORWARD, backward = FFTW_BACKWARD };
enum strategy : unsigned { estimate = FFTW_ESTIMATE, measure = FFTW_MEASURE };

template<typename I, class O>
auto dft(I&& i, O&& o, sign s, strategy st = fftw::estimate){// unsigned flags = FFTW_ESTIMATE){
	execute(fftw::plan{std::forward<I>(i), std::forward<O>(o), (int)s, (unsigned)st | FFTW_PRESERVE_INPUT});//.execute();
//	fftw_plan p = multi::fftw_plan_dft(in, out, sign, flags); //	fftw_execute(p); //	fftw_destroy_plan(p);
}

template<typename I, typename O = multi::array<typename std::decay_t<I>::element_type, std::decay_t<I>::dimensionality>>
O dft(I const& i, sign s, strategy st = fftw::estimate){
	O ret(extensions(i));
	dft(i, ret, s, st);
	return ret;
}

template<typename I>
I& dft_inplace(I& i, sign s, strategy st = fftw::estimate){
	execute(fftw::plan{i, i, (int)s, (unsigned)st | FFTW_PRESERVE_INPUT});
	return i;
}
//template<typename I>
//I& dft(I& i, sign s, strategy st = fftw::estimate){return dft_inplace(i, s, st);}

template<typename I, typename = std::enable_if_t<std::is_rvalue_reference<I&&>{}> >
I&& dft(I&& i, sign s, strategy st = fftw::estimate){
	execute(fftw::plan{i, i, (int)s, (unsigned)st | FFTW_PRESERVE_INPUT});
	return i;
}

}

}}

#if _TEST_MULTI_ADAPTORS_FFTW

#include "../adaptors/fftw/allocator.hpp"
#include<iostream>
#include "../array.hpp"
#include<complex>
#include<numeric>

//#define CATCH_CONFIG_MAIN
#include<catch2/catch.hpp>

namespace{
	using namespace Catch::literals;

	using std::cout;
	namespace multi = boost::multi;
	namespace fftw = multi::fftw;

	using complex = std::complex<double>;

	template<class M>
	auto power(M const& m){
		auto sum_norm = [](auto& a, auto& b){return a + std::norm(b);};
		using multi::num_elements; using multi::data_elements; using std::accumulate;
		return accumulate(data_elements(m), data_elements(m) + num_elements(m), double{}, sum_norm);
	}

	constexpr int N = 16;
}

TEST_CASE("fftw 1D power", "[report]"){
	multi::array<complex, 1> in({N}, 0.); assert( size(in) == N );
	std::iota(begin(in), end(in), 1.);
	multi::array<complex, 1> out(extensions(in));
	auto p = multi::fftw_plan_dft(in, out, FFTW_FORWARD, FFTW_PRESERVE_INPUT);
	fftw_execute(p); 
	fftw_destroy_plan(p);
	REQUIRE( (power(in) - power(out)/num_elements(out)) == (0_a).margin(1e-17) );
}

TEST_CASE("fftw 1D allocator power", "[report]"){
	using multi::fftw::allocator;
	multi::array<complex, 1, allocator<>> in({16}, 0.); std::iota(begin(in), end(in), 1.);
	assert( size(in) == N );
	multi::array<complex, 1, allocator<>> out(extensions(in));
	auto p = multi::fftw_plan_dft(in, out, FFTW_FORWARD, FFTW_PRESERVE_INPUT);
	fftw_execute(p);
	fftw_destroy_plan(p);
	REQUIRE( (power(in) - power(out)/num_elements(out)) == (0_a).margin(1e-12) );
}

TEST_CASE("fftw 2D power", "[report]"){
	multi::array<complex, 2> in({N, N});
	std::iota(data(in), data(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	auto p = multi::fftw_plan_dft(in, out, FFTW_FORWARD, FFTW_PRESERVE_INPUT);
	fftw_execute(p); fftw_destroy_plan(p);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-12) );
}

TEST_CASE("fftw 2D power plan", "[report]"){
	multi::array<complex, 2> in({16, 16});
	std::iota(data(in), data(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const p{in, out, FFTW_FORWARD, FFTW_PRESERVE_INPUT};
	p(); //execute(p); //p.execute();
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-8) );
}

TEST_CASE("fftw 2D power dft", "[report]"){
	multi::array<complex, 2> in({16, 16}); std::iota(data(in), data(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft(in, out, multi::fftw::forward);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-8) );
}

TEST_CASE("fftw 2D power dft out", "[report]"){
	multi::array<complex, 2> in({16, 16}); std::iota(data(in), data(in) + num_elements(in), 1.2);
	auto out = multi::fftw::dft(in, multi::fftw::forward, multi::fftw::estimate);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-8) );
}

TEST_CASE("fftw 2D power dft out default", "[report]"){
	multi::array<complex, 2> in({16, 16}); std::iota(data(in), data(in) + num_elements(in), 1.2);
	auto out = multi::fftw::dft(in, fftw::forward);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-8) );
}

TEST_CASE("fftw 2D c-array power", "[report]"){
	int const N = 16;
	complex in[N][N];
	using multi::data_elements;	using multi::num_elements;
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	complex out[N][N];
	auto p = multi::fftw_plan_dft(in, out, FFTW_FORWARD | FFTW_PRESERVE_INPUT);
	fftw_execute(p); fftw_destroy_plan(p);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-8) );
}

TEST_CASE("fftw 3D power", "[report]"){
	multi::array<complex, 3> in({4, 4, 4}); std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);
	multi::array<complex, 3> out = fftw::dft(in, fftw::forward);
	REQUIRE( power(in) - power(out)/num_elements(out) == (0_a).margin(1e-10) );
}

TEST_CASE("fftw 3D power in place", "[report]"){
	multi::array<complex, 3> io({4, 4, 4}); std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);
	auto powerin = power(io);
	fftw::dft_inplace(io, fftw::forward);
	REQUIRE( powerin - power(io)/num_elements(io) == (0_a).margin(1e-10) );
}

#if 0
{
	auto in3 = [](){
		multi::array<complex, 3> in3({N, N, N});
		std::iota(in3.data(), in3.data() + in3.num_elements(), 1.2);
		return in3;
	}();
	multi::array<complex, 3> out3(extensions(in3));
	auto p = 
//		multi::fftw_plan_dft_3d(in3, out3, FFTW_FORWARD);
		multi::fftw_plan_dft(in3, out3, FFTW_FORWARD| FFTW_PRESERVE_INPUT)
	;
	fftw_execute(p);
	assert( power_diff(in3, out3) < 1e-3 );	
	fftw_execute_dft(p, in3, out3);//, out3);
	fftw_destroy_plan(p);
	assert( power_diff(in3, out3) < 1e-3 );
}
{
	auto in4 = [](){
		multi::array<complex, 4> in4({5, 6, 7, 8});
		std::iota(in4.data(), in4.data() + in4.num_elements(), 10.2);
		return in4;
	}();
	multi::array<complex, 4> out4(extensions(in4));
	auto p = multi::fftw_plan_dft(in4, out4, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fftw_execute(p);
	fftw_destroy_plan(p);
	assert( power_diff(in4, out4) < 1e-3 );
}
#endif
#endif
#endif

