#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ `#-DNDEBUG` -std=c++17 -Wall -Wextra -I$HOME/prj -D_TEST_MULTI_ADAPTORS_FFTW $0x.cpp -o $0x.x -lfftw3 && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP

#include<memory> // allocator
#include<cassert>
#include<complex>
#include <fftw3.h>
#include<new>

#include "../utility.hpp"

namespace boost{
namespace multi{

template<class T>
struct fftw_allocator : std::allocator<T>{
	template< class U > struct rebind { typedef fftw_allocator<U> other; };
	T* allocate(typename fftw_allocator::size_type n){
		T* ret = static_cast<T*>(fftw_malloc(sizeof(T)*n));
		if(ret == nullptr) throw std::bad_alloc{};
		return ret;
	}
	void deallocate(T* p, typename fftw_allocator::size_type){fftw_free(p);}
};

template<class T>
struct fftw_uallocator : fftw_allocator<T>{
	template< class U > struct rebind { typedef fftw_uallocator<U> other; };
	template<class U> void construct(U*) noexcept{
		static_assert(
			std::is_trivially_destructible<T>{}, 
			"uallocator cannot be used with non trivial types"
		);
	}
};


//////////////////////////////////////////////////////////////////////////////

template<typename Size>
auto fftw_plan_dft_1d(Size N, std::complex<double> const* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE){
	assert( flags & FFTW_ESTIMATE );
	assert( flags & FFTW_PRESERVE_INPUT );
#ifndef NDEBUG
	auto check = in[N/3];
#endif
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out) );
	auto ret=::fftw_plan_dft_1d(N, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
	assert(check == in[N/3]); // check that const data has not been overwritten
	return ret;
}
template<typename Size>
auto fftw_plan_dft_1d(Size N, std::complex<double>* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_1d(N, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_2d(Size N1, Size N2, std::complex<double> const* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE |FFTW_PRESERVE_INPUT){
	assert( flags & FFTW_ESTIMATE );
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_2d(N1, N2, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename Size>
auto fftw_plan_dft_2d(Size N1, Size N2, std::complex<double>* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_2d(N1, N2, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Size>
auto fftw_plan_dft_3d(Size N1, Size N2, Size N3, std::complex<double>* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename Size>
auto fftw_plan_dft_3d(Size N1, Size N2, Size N3, std::complex<double> const* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT){
	assert( flags & FFTW_ESTIMATE );
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft_3d(N1, N2, N3, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

template<typename Rank>
auto fftw_plan_dft(Rank r, int* ns, std::complex<double>* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft(r, ns, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}
template<typename Rank>
auto fftw_plan_dft(Rank r, int* ns, std::complex<double> const* in, std::complex<double>* out, int sign, unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT){
	assert( flags & FFTW_ESTIMATE );
	assert( flags & FFTW_PRESERVE_INPUT );
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out));
	return ::fftw_plan_dft(r, ns, (fftw_complex*)in, (fftw_complex*)out, sign, flags);
}

//////////////////////////////////////////////////////////////////////////////

template<class InIt, typename Size, class OutIt>
auto fftw_plan_dft_1d_n(InIt it, Size n, OutIt dest, decltype(FFTW_FORWARD) sign, decltype(FFTW_ESTIMATE) flags = FFTW_ESTIMATE){
	return fftw_plan_dft_1d(n, data(it), data(dest), sign, flags);
}
template<class InIt, class OutIt>
auto fftw_plan_dft_1d(InIt first, InIt last, OutIt dest, decltype(FFTW_FORWARD) sign, decltype(FFTW_ESTIMATE) flags = FFTW_ESTIMATE){
	using std::distance;
	return fftw_plan_dft_1d_n(first, distance(first, last), dest, sign, flags);
}

///////////////////////////////////////////////////////////////////////////////

template<class ArrayIn, class ArrayOut>
auto fftw_plan_dft_1d(ArrayIn&& in, ArrayOut&& out, decltype(FFTW_FORWARD) sign, decltype(FFTW_MEASURE) flags = FFTW_ESTIMATE)
->decltype(boost::multi::fftw_plan_dft_1d(size(in), data(in), data(out), sign, flags)){
	static_assert(in.dimensionality == 1, "!");
	assert(size(in) == size(out));
	return boost::multi::fftw_plan_dft_1d(size(in), data(in), data(out), sign, flags);}

template<class ArrayIn, class ArrayOut>
auto fftw_plan_dft_2d(ArrayIn&& in, ArrayOut&& out, decltype(FFTW_FORWARD) sign, decltype(FFTW_MEASURE) flags = FFTW_ESTIMATE){
	static_assert(in.dimensionality == 2, "!");
	assert(in.sizes() == out.sizes());
	if(std::is_const<decltype(out)>{}) assert(flags & FFTW_MEASURE);
	return fftw_plan_dft_2d(
		sizes(in)[0], sizes(in)[1],
		data(in), data(out),
		sign, flags
	);
}

template<class ArrayIn, class ArrayOut>
auto fftw_plan_dft_3d(ArrayIn&& in, ArrayOut&& out, decltype(FFTW_FORWARD) sign, decltype(FFTW_MEASURE) flags = FFTW_ESTIMATE){
	static_assert(in.dimensionality == 3, "!");
	assert(in.sizes() == out.sizes());
	if(std::is_const<decltype(out)>{}) assert(flags & FFTW_MEASURE);
	return fftw_plan_dft_3d(
		sizes(in)[0], sizes(in)[1], sizes(in)[2],
		data(in), data(out),
		sign, flags
	);
}

namespace detail {
template<class T, class Tuple, std::size_t... I>
constexpr std::array<std::remove_cv_t<T>, std::tuple_size<Tuple>{}> 
to_array_impl(Tuple&& t, std::index_sequence<I...>){
	return { static_cast<std::remove_cv_t<T>>(std::get<I>(std::forward<Tuple>(t)))... };
}
}
 
template<class T, class Tuple> constexpr auto to_array(Tuple&& t){
    return detail::to_array_impl<T>(std::forward<Tuple>(t), std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}

template<class ArrayIn, class ArrayOut>
auto fftw_plan_dft(
	ArrayIn&& in, ArrayOut&& out, 
	int sign, unsigned flags = FFTW_ESTIMATE
){
	using boost::multi::dimensionality;
	static_assert(dimensionality(in) == dimensionality(out), "!");
//	static_assert(in.dimensionality == out.dimensionality, "!");
	using boost::multi::sizes;
	assert(sizes(in) == sizes(out));
	static_assert(dimensionality(out) > 0, "!");
	static_assert(dimensionality(in) > 0, "!");
	std::array<int, dimensionality(in)> ns = to_array<int>(sizes(in));
//	for(auto i = 0; i != dimensionality(in); ++i) 
//		ns[i] = sizes(in)[i];
	return fftw_plan_dft(dimensionality(in), ns.data(), origin(std::forward<ArrayIn>(in)), origin(std::forward<ArrayOut>(out)), sign, flags);
}

///////////////////////////////////////////////////////////////////////////////

void fftw_execute_dft(
	fftw_plan p, std::complex<double>* in, std::complex<double>* out){
	assert(fftw_alignment_of((double*)in) == fftw_alignment_of((double*)out) );
	::fftw_execute_dft(p, (fftw_complex*)in, (fftw_complex*)out);
}

template<class ArrayIn, class ArrayOut>
auto fftw_execute_dft(fftw_plan p, ArrayIn&& in, ArrayOut&& out)
->decltype(boost::multi::fftw_execute_dft(p, data(in), data(out))){
	boost::multi::fftw_execute_dft(p, data(in), data(out));
}

}}

#if _TEST_MULTI_ADAPTORS_FFTW

#include<iostream>
#include "../array.hpp"
#include<complex>

using std::cout;
namespace multi = boost::multi;
using complex = std::complex<double>;

template<class M> 
auto power(M const& m){
	auto sum_norm = [](auto& a, auto& b){return a + std::norm(b);};
	return accumulate(data(m), data(m) + num_elements(m), complex{}, sum_norm);
}
template<class M1, class M2> 
auto power_diff(M1 const& m1, M2 const& m2){
	using std::abs;
	return abs(power(m1) - power(m2)/double(num_elements(m2)));
}

constexpr int N = 16;

int main(){
	std::allocator_traits<typename multi::fftw_uallocator<double>>::rebind_alloc<complex> aa{};
	complex* p = aa.allocate(10);
	complex* q = new complex[10];//aa.allocate(10);
	assert(fftw_alignment_of((double*)p) == fftw_alignment_of((double*)q));	
{
	auto const in = [](){
		multi::array<complex, 1, multi::fftw_allocator<complex>> in({N});
		std::iota(begin(in), end(in), 1.2);
		return in;
	}();
	multi::array<complex, 1> out(in.extensions());
	auto p = 
//		multi::fftw_plan_dft_1d(in, out, FFTW_FORWARD)
		multi::fftw_plan_dft(in, out, FFTW_FORWARD | FFTW_PRESERVE_INPUT)
	;
	fftw_execute(p);
	fftw_destroy_plan(p);
	assert( power_diff(in, out) < 1e-10 );
}
{
	auto const in2 = [](){
		multi::array<complex, 2, multi::fftw_allocator<complex>> in2({N, N});
		std::iota(in2.data(), in2.data() + in2.num_elements(), 1.2);
		return in2;
	}();
	multi::array<complex, 2, multi::fftw_uallocator<complex>> out2(extensions(in2));
	auto p = 
//		multi::fftw_plan_dft_2d(in2, out2, FFTW_FORWARD)
		multi::fftw_plan_dft(in2, out2, FFTW_FORWARD | FFTW_PRESERVE_INPUT)
	;
	fftw_execute(p); fftw_destroy_plan(p);
	assert( power_diff(in2, out2) < 1e-4 );
}
{
	complex in2[32][32];
	complex out2[32][32];
	auto p = multi::fftw_plan_dft(in2, out2, FFTW_FORWARD | FFTW_PRESERVE_INPUT);
	fftw_execute(p); fftw_destroy_plan(p);
	assert( power_diff(in2, out2) < 1e-4 );
}
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
}

#endif
#endif

