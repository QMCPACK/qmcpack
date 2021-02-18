#ifdef COMPILATION// -*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
$CXX $0 -o $0x -lcudart -lcufft `pkg-config --libs fftw3` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_ADAPTORS_CUFFTW_HPP
#define MULTI_ADAPTORS_CUFFTW_HPP

#include "../adaptors/../utility.hpp"
#include "../adaptors/../array.hpp"
#include "../adaptors/../config/NODISCARD.hpp"

#include "../adaptors/cuda.hpp"

#include<numeric>

#include<tuple> // std::apply
#include<array>

#include<vector>

#include "../complex.hpp"

//#include<execution>
#include<future>
#include<cufft.h>

namespace boost{
namespace multi{
namespace memory{
namespace cuda{

#if 0
template<class T1, class T1const, class T2, class T2const>
auto copy(
	array_iterator<T1, 1, managed::ptr<T1const>> first, 
	array_iterator<T1, 1, managed::ptr<T1const>> last, 
	array_iterator<T2, 1, managed::ptr<T2const>> d_first
){
	assert(first.stride() == last.stride());
	auto s = cudaMemcpy2D(raw_pointer_cast(d_first.data()), d_first.stride()*sizeof(T2), raw_pointer_cast(first.data()), first.stride()*sizeof(T2), sizeof(T2), last - first, cudaMemcpyDefault);
	switch(s){
		case cudaSuccess: break;
		case cudaErrorInvalidValue: assert(0);
		case cudaErrorInvalidPitchValue: assert(0);
		case cudaErrorInvalidDevicePointer: assert(0);
		case cudaErrorInvalidMemcpyDirection: assert(0);
		default: assert(0); // unknown error
	}
	return d_first + (last - first);
}
#endif

}}}}

namespace boost{
namespace multi{
namespace cufft{

class sign{
	int impl_;
public:
	sign() = default;
	constexpr sign(int i) : impl_{i}{}
	constexpr operator int() const{return impl_;}
};

constexpr sign forward{CUFFT_FORWARD};
constexpr sign none{0};
constexpr sign backward{CUFFT_INVERSE};

static_assert(forward != none and none != backward and backward != forward, "!");

class plan{
	using complex_type = cufftDoubleComplex;
	complex_type const* idata_     = nullptr;
	complex_type*       odata_     = nullptr;
	int                 direction_ = 0;
	cufftHandle h_;
	plan() = default;
	plan(plan const&) = delete;
	plan(plan&& other) : 
		idata_{std::exchange(other.idata_, nullptr)}, 
		odata_{std::exchange(other.odata_, nullptr)}, 
		direction_{std::exchange(other.direction_, 0)}, 
		h_{std::exchange(other.h_, {})}
	{} // needed in <=C++14 for return
	void ExecZ2Z(complex_type const* idata, complex_type* odata, int direction) const{
		++tl_execute_count;
	//	assert(idata_ and odata_); 
	//	assert(direction_!=0);
		cufftResult r = ::cufftExecZ2Z(h_, const_cast<complex_type*>(idata), odata, direction); 
		switch(r){
			case CUFFT_SUCCESS        : break;// "cuFFT successfully executed the FFT plan."
			case CUFFT_INVALID_PLAN   : throw std::runtime_error{"The plan parameter is not a valid handle."};
		//	case CUFFT_ALLOC_FAILED   : throw std::runtime_error{"CUFFT failed to allocate GPU memory."};
		//	case CUFFT_INVALID_TYPE   : throw std::runtime_error{"The user requests an unsupported type."};
			case CUFFT_INVALID_VALUE  : throw std::runtime_error{"At least one of the parameters idata, odata, and direction is not valid."};
			case CUFFT_INTERNAL_ERROR : throw std::runtime_error{"Used for all internal driver errors."};
			case CUFFT_EXEC_FAILED    : throw std::runtime_error{"CUFFT failed to execute an FFT on the GPU."};
			case CUFFT_SETUP_FAILED   : throw std::runtime_error{"The cuFFT library failed to initialize."};
		//	case CUFFT_INVALID_SIZE   : throw std::runtime_error{"The user specifies an unsupported FFT size."};
		//	case CUFFT_UNALIGNED_DATA : throw std::runtime_error{"Unaligned data."};
		//	case CUFFT_INCOMPLETE_PARAMETER_LIST: throw std::runtime_error{"Incomplete parameter list."};
		//	case CUFFT_INVALID_DEVICE : throw std::runtime_error{"Invalid device."};
		//	case CUFFT_PARSE_ERROR    : throw std::runtime_error{"Parse error."};
		//	case CUFFT_NO_WORKSPACE   : throw std::runtime_error{"No workspace."};
		//	case CUFFT_NOT_IMPLEMENTED: throw std::runtime_error{"Not implemented."};
		//	case CUFFT_LICENSE_ERROR  : throw std::runtime_error{"License error."};
		//	case CUFFT_NOT_SUPPORTED  : throw std::runtime_error{"CUFFT_NOT_SUPPORTED"};
			default                   : throw std::runtime_error{"cufftExecZ2Z unknown error"};
		}
	//	if(cudaDeviceSynchronize() != cudaSuccess) throw std::runtime_error{"Cuda error: Failed to synchronize"};
	}
	void swap(plan& other){ 
		using std::swap;
		swap(idata_, other.idata_);
		swap(odata_, other.odata_);
		swap(direction_, other.direction_);
		swap(h_, other.h_);
	}
public:
	thread_local static int tl_execute_count;
	plan& operator=(plan other){swap(other); return *this;}
	void operator()() const{ExecZ2Z(idata_, odata_, direction_);}
	template<class I, class O>
	O&& execute_dft(I&& i, O&& o, int direction) const{
		ExecZ2Z(
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(i))),
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(o))),
			direction
		);
		return std::forward<O>(o);
	}
	template<class I, class O>
	void execute_dft(I&& i, O&& o) const{execute_dft(std::forward<I>(i), std::forward<O>(o), direction_);}
	~plan(){if(h_) cufftDestroy(h_);}
	using size_type = int;
	using ssize_type = int;

	template<class I, class O, //std::enable_if_t<(I::dimensionality < 4), int> =0,
		dimensionality_type D = I::dimensionality,  
		typename = decltype(raw_pointer_cast(base(std::declval<I const&>())), reinterpret_cast<complex_type*      >(raw_pointer_cast(base(std::declval<O&>()))))
	>
	plan(I const& i, O&& o, sign s) : 
		idata_{                          reinterpret_cast<complex_type const*>(raw_pointer_cast(base(i))) },
		odata_{const_cast<complex_type*>(reinterpret_cast<complex_type*      >(raw_pointer_cast(base(o))))},
		direction_{s}
	{
		assert( I::dimensionality < 4 );
		assert( CUFFT_FORWARD == s or CUFFT_INVERSE == s or s == 0 );
		assert( sizes(i) == sizes(o) );

//		using std::experimental::apply;// using std::experimental::make_array;
		auto ion      = std::apply([](auto... t){return std::array< size_type, D>{static_cast< size_type>(t)...};}, sizes  (i));
		auto istrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(i));
		auto ostrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(o));

		std::array<std::tuple<int, int, int>, I::dimensionality> ssn;
		for(std::size_t i = 0; i != ssn.size(); ++i) ssn[i] = std::make_tuple(istrides[i], ostrides[i], ion[i]);
		std::sort(ssn.begin(), ssn.end(), std::greater<>{});

		for(std::size_t i = 0; i != ssn.size(); ++i){
			istrides[i] = std::get<0>(ssn[i]);
			ostrides[i] = std::get<1>(ssn[i]);
			ion[i]      = std::get<2>(ssn[i]);
		}// = std::tuple<int, int, int>(istrides[i], ostrides[i], ion[i]);

		int istride = istrides.back();
		auto inembed = istrides; inembed.fill(0);
		int ostride = ostrides.back();
		auto onembed = ostrides; onembed.fill(0);	
		for(std::size_t i = 1; i != onembed.size(); ++i){
			assert(ostrides[i-1] >= ostrides[i]); // otherwise ordering is incompatible
			assert(ostrides[i-1]%ostrides[i]==0);
			onembed[i]=ostrides[i-1]/ostrides[i]; //	assert( onembed[i] <= ion[i] );
			assert(istrides[i-1]%istrides[i]==0);
			inembed[i]=istrides[i-1]/istrides[i]; //	assert( inembed[i] <= ion[i] );
		}

		direction_ = s;
		idata_ =                           reinterpret_cast<complex_type const*>(raw_pointer_cast(base(i))) ;
		odata_ = const_cast<complex_type*>(reinterpret_cast<complex_type*      >(raw_pointer_cast(base(o))));

		switch(::cufftPlanMany(
			/*cufftHandle *plan*/ &h_, 
			/*int rank*/          ion.size(), 
			/*int *n*/            ion.data(), //	/*NX*/      last - first,
			/*int *inembed*/      inembed.data(),
			/*int istride*/       istride, 
			/*int idist*/         1, //stride(first), 
			/*int *onembed*/      onembed.data(), 
			/*int ostride*/       ostride, 
			/*int odist*/         1, //stride(d_first), 
			/*cufftType type*/    CUFFT_Z2Z, 
			/*int batch*/         1 //BATCH
		)){
			case CUFFT_SUCCESS        : break;// "cuFFT successfully executed the FFT plan."
			case CUFFT_ALLOC_FAILED   : throw std::runtime_error{"CUFFT failed to allocate GPU memory."};
			case CUFFT_INVALID_VALUE  : throw std::runtime_error{"At least one of the parameters idata, odata, and direction is not valid."};
			case CUFFT_INTERNAL_ERROR : throw std::runtime_error{"Used for all internal driver errors."};
			case CUFFT_SETUP_FAILED   : throw std::runtime_error{"The cuFFT library failed to initialize."};
			case CUFFT_INVALID_SIZE   : throw std::runtime_error{"The user specifies an unsupported FFT size."};
			default                   : throw std::runtime_error{"cufftPlanMany unknown error"};
		}
	}
#ifndef __INTEL_COMPILER
	template<class It1, class It2, dimensionality_type D = decltype(*It1{})::dimensionality>
	static auto many(It1 first, It1 last, It2 d_first, int sign = 0, unsigned = 0)
	->std::decay_t<decltype(const_cast<complex_type*>(reinterpret_cast<complex_type*>(raw_pointer_cast(base(d_first)))), std::declval<plan>())>
#else
	template<class It1, class It2, 
		dimensionality_type D = decltype(*It1{})::dimensionality, 
		typename TT = decltype(const_cast<complex_type*>(reinterpret_cast<complex_type*>(It2{}.base().raw_pointer_cast())))
	>
	static auto many(It1 first, It1 last, It2 d_first, int sign = 0, unsigned = 0)
#endif
	{
		assert( CUFFT_FORWARD == sign or CUFFT_INVERSE == sign or sign == 0 );
		assert(sizes(*first)==sizes(*d_first));
		auto ion      = std::apply([](auto... t){return std::array< size_type, D>{static_cast< size_type>(t)...};}, sizes  (*  first));

		assert(strides(*first) == strides(*last));
		auto istrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(*  first));
		auto ostrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(*d_first));

		std::array<std::tuple<int, int, int>, std::decay_t<decltype(*It1{})>::dimensionality> ssn;
		for(std::size_t i = 0; i != ssn.size(); ++i) ssn[i] = std::make_tuple(istrides[i], ostrides[i], ion[i]);
		std::sort(ssn.begin(), ssn.end(), std::greater<>{});

		for(std::size_t i = 0; i != ssn.size(); ++i){
			istrides[i] = std::get<0>(ssn[i]);
			ostrides[i] = std::get<1>(ssn[i]);
			ion[i]      = std::get<2>(ssn[i]);
		}

		int istride = istrides.back();
		auto inembed = istrides; inembed.fill(0);
		int ostride = ostrides.back();
		auto onembed = ostrides; onembed.fill(0);	
		for(std::size_t i = 1; i != onembed.size(); ++i){
			assert(ostrides[i-1] >= ostrides[i]); // otherwise ordering is incompatible
			assert(ostrides[i-1]%ostrides[i]==0);
			onembed[i]=ostrides[i-1]/ostrides[i]; //	assert( onembed[i] <= ion[i] );
			assert(istrides[i-1]%istrides[i]==0);
			inembed[i]=istrides[i-1]/istrides[i]; //	assert( inembed[i] <= ion[i] );
		}

		plan ret;
		ret.direction_ = sign;
		ret.idata_ =                           reinterpret_cast<complex_type const*>(  first.base().raw_pointer_cast()) ;
		ret.odata_ = const_cast<complex_type*>(reinterpret_cast<complex_type*      >(d_first.base().raw_pointer_cast()));

		switch(::cufftPlanMany(
			/*cufftHandle *plan*/ &ret.h_, 
			/*int rank*/          ion.size(), 
			/*int *n*/            ion.data(), //	/*NX*/      last - first,
			/*int *inembed*/      inembed.data(),
			/*int istride*/       istride, 
			/*int idist*/         stride(first), 
			/*int *onembed*/      onembed.data(), 
			/*int ostride*/       ostride, 
			/*int odist*/         stride(d_first), 
			/*cufftType type*/    CUFFT_Z2Z, 
			/*int batch*/         last - first //BATCH
		)){
			case CUFFT_SUCCESS        : break;// "cuFFT successfully executed the FFT plan."
		//	case CUFFT_INVALID_PLAN   : throw std::runtime_error{"The plan parameter is not a valid handle."};
			case CUFFT_ALLOC_FAILED   : throw std::runtime_error{"CUFFT failed to allocate GPU memory."};
		//	case CUFFT_INVALID_TYPE   : throw std::runtime_error{"The user requests an unsupported type."};
			case CUFFT_INVALID_VALUE  : throw std::runtime_error{"At least one of the parameters idata, odata, and direction is not valid."};
			case CUFFT_INTERNAL_ERROR : throw std::runtime_error{"Used for all internal driver errors."};
		//	case CUFFT_EXEC_FAILED    : throw std::runtime_error{"CUFFT failed to execute an FFT on the GPU."};
			case CUFFT_SETUP_FAILED   : throw std::runtime_error{"The cuFFT library failed to initialize."};
			case CUFFT_INVALID_SIZE   : throw std::runtime_error{"The user specifies an unsupported FFT size."};
		//	case CUFFT_UNALIGNED_DATA : throw std::runtime_error{"Unaligned data."};
		//	case CUFFT_INCOMPLETE_PARAMETER_LIST: throw std::runtime_error{"Incomplete parameter list."};
		//	case CUFFT_INVALID_DEVICE : throw std::runtime_error{"Invalid device."};
		//	case CUFFT_PARSE_ERROR    : throw std::runtime_error{"Parse error."};
		//	case CUFFT_NO_WORKSPACE   : throw std::runtime_error{"No workspace."};
		//	case CUFFT_NOT_IMPLEMENTED: throw std::runtime_error{"Not implemented."};
		//	case CUFFT_LICENSE_ERROR  : throw std::runtime_error{"License error."};
		//	case CUFFT_NOT_SUPPORTED  : throw std::runtime_error{"CUFFT_NOT_SUPPORTED"};
			default                   : throw std::runtime_error{"cufftPlanMany unknown error"};
		}
		return ret;
	}
};

thread_local int plan::tl_execute_count = 0;

template<typename In, class Out>
auto dft(In const& i, Out&& o, int s)
->decltype(cufft::plan{i, o, s}(), std::forward<Out>(o)){
	return cufft::plan{i, o, s}(), std::forward<Out>(o);}

template<typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
R dft(In const& i, int s){
	static_assert(std::is_trivially_default_constructible<typename In::element_type>{}, "!");
	R ret(extensions(i), get_allocator(i));
	cufft::dft(i, ret, s);
	return ret;
}

#ifndef __INTEL_COMPILER
template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, sign s)
->decltype(plan::many(first, last, d_first, s)(), d_first + (last - first)){
	return plan::many(first, last, d_first, s)(), d_first + (last - first);}
#else
template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, sign s)
->decltype(plan::many(first, last, d_first, s)(), d_first + (last - first)){
	return plan::many(first, last, d_first, s)(), d_first + (last - first);}
#endif

template<typename In, class Out,  std::size_t D = In::dimensionality, std::enable_if_t<(D==1), int> = 0>
Out&& dft(std::array<bool, D> which, In const& i, Out&& o, int s){
	if(which[0]) return cufft::dft(i, std::forward<Out>(o), s);
	else return std::forward<Out>(std::forward<Out>(o) = i);
}

template <class Array, std::size_t... Ns>
constexpr auto array_tail_impl(Array const& t, std::index_sequence<Ns...>){
	return std::array<typename Array::value_type, std::tuple_size<Array>{} - 1>{std::get<Ns + 1>(t)...};
}

template<class Array>
constexpr auto array_tail(Array const& t)
->decltype(array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>())){
	return array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>());}

template<typename In, class Out, std::size_t D = In::dimensionality, std::enable_if_t<(D>1), int> = 0> 
auto dft(std::array<bool, D> which, In const& i, Out&& o, int s)
->decltype(many_dft(i.begin(), i.end(), o.begin(), s),std::forward<Out>(o))
{
	assert(extension(i) == extension(o));
	auto ff = std::find(begin(which)+1, end(which), false);
	if(which[0] == true){
		if(ff==end(which)) cufft::dft(i, std::forward<Out>(o), s);
		else{
			auto n = ff - which.begin();
			std::rotate(begin(which), ff, end(which));
			dft(which, i<<n, o<<n, s);
		}
	}else if(which[0]==false){
		if(D==1 or std::none_of(begin(which)+1, end(which), [](auto e){return e;})){
			if(base(o) != base(i)) std::forward<Out>(o) = i;
			else if(o.layout() != i.layout()) std::forward<Out>(o) = +i;
		}
		else if(ff==end(which)) many_dft(i.begin(), i.end(), o.begin(), s);
		else{
			std::array<bool, D-1> tail = array_tail(which);
			if(which[1] == false and i.is_flattable() and o.is_flattable()) cufft::dft(tail, i.flatted(), o.flatted(), s);
			else{
				auto d_min = 0; auto n_min = size(i);
				for(auto d = 0; d != D - 1; ++d){
					if((size(i<<d) < n_min) and (tail[d]==false)){
						n_min = size(i<<d); 
						d_min = d;
					}
				}
				if( d_min!=0 ){
					std::rotate(which.begin(), which.begin()+d_min, which.end());
					dft(which, i<<d_min, o<<d_min, s);
				}else
				{
					if(base(i) == base(o) and i.layout() != o.layout()){
						auto const tmp = +i;
						for(auto idx : extension(i)) cufft::dft(tail, tmp[idx], o[idx], s);
					}else for(auto idx : extension(i)) cufft::dft(tail, i[idx], o[idx], s);
				}
			}
		}
	}
	return std::forward<Out>(o);
}

template<typename In,  std::size_t D = In::dimensionality>
NODISCARD("when passing a const argument")
auto dft(std::array<bool, D> which, In const& i, int sign)->std::decay_t<decltype(
dft(which, i, typename In::decay_type(extensions(i), get_allocator(i)), sign))>{return 
dft(which, i, typename In::decay_type(extensions(i), get_allocator(i)), sign);}

template<typename In,  std::size_t D = In::dimensionality>
auto dft(std::array<bool, D> which, In&& i, int sign)
->decltype(dft(which, i, i, sign), std::forward<In>(i)){
	return dft(which, i, i, sign), std::forward<In>(i);}	

//template<typename... A> auto            dft_forward(A&&... a)
//->decltype(cufft::dft(std::forward<A>(a)..., cufft::forward)){
//	return cufft::dft(std::forward<A>(a)..., cufft::forward);}

template<typename Array, typename A> NODISCARD("when passing a const argument")
auto dft_forward(Array arr, A const& a) 
->decltype(cufft::dft(arr, a, cufft::forward)){
	return cufft::dft(arr, a, cufft::forward);}

template<typename Array, dimensionality_type D> NODISCARD("when passing a const argument")
auto dft_forward(Array arr, multi::cuda::array<std::complex<double>, D>&& a) 
->decltype(cufft::dft(arr, a, cufft::forward), multi::cuda::array<std::complex<double>, D>{}){//assert(0);
	return cufft::dft(arr, a, cufft::forward), std::move(a);}

template<typename A> NODISCARD("when passing a const argument")
auto dft_forward(A const& a)
->decltype(cufft::dft(a, cufft::forward)){
	return cufft::dft(a, cufft::forward);}

template<typename... A> auto            dft_backward(A&&... a)
->decltype(cufft::dft(std::forward<A>(a)..., cufft::backward)){
	return cufft::dft(std::forward<A>(a)..., cufft::backward);}

template<typename Array, typename A> NODISCARD("when passing a const argument")
auto dft_backward(Array arr, A const& a) 
->decltype(cufft::dft(arr, a, cufft::backward)){
	return cufft::dft(arr, a, cufft::backward);}

template<typename A> NODISCARD("when passing a const argument")
auto dft_backward(A const& a)
->decltype(cufft::dft(a, cufft::backward)){
	return cufft::dft(a, cufft::backward);}


}

}}


#if not __INCLUDE_LEVEL__ // TEST BELOW

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuFFT adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

#include "../adaptors/cuda.hpp"
#include "../adaptors/fftw.hpp"
#include "../adaptors/cufft.hpp"

//#include "../adaptors/fft.hpp"

#include<complex>
#include<thrust/complex.h>
#include "../complex.hpp"

#include<cuda_runtime.h> // cudaDeviceSynchronize

#include<iostream>

namespace multi = boost::multi;
using complex = std::complex<double>;
namespace utf = boost::unit_test;


template <class T>
__attribute__((always_inline)) inline void DoNotOptimize(const T &value) {
	asm volatile("" : "+m"(const_cast<T &>(value)));
}

struct watch : private std::chrono::high_resolution_clock{
	std::string label_; time_point  start_;
	watch(std::string label ="") : label_{label}, start_{}{
		cudaDeviceSynchronize();
		start_ = now();
	}
	~watch(){
		cudaDeviceSynchronize();
		auto const count = std::chrono::duration<double>(now() - start_).count();
		std::cerr<< label_<<": "<< count <<" sec"<<std::endl;
	}
};

constexpr complex I{0, 1};

#if 1

BOOST_AUTO_TEST_CASE(cufft_2D, *boost::unit_test::tolerance(0.0001)){

	multi::array<complex, 2> const in_cpu = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> fw_cpu(extensions(in_cpu));
	multi::fftw::dft(in_cpu, fw_cpu, multi::fftw::forward);

	multi::cuda::array<complex, 2> const in_gpu = in_cpu;
	multi::cuda::array<complex, 2> fw_gpu(extensions(in_gpu));
	multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);

	BOOST_TEST( imag(static_cast<complex>(fw_gpu[3][2]) - fw_cpu[3][2]) == 0. );

	auto fw2_gpu = multi::cufft::dft(in_gpu, multi::cufft::forward);
	BOOST_TEST( imag(static_cast<complex>(fw2_gpu[3][1]) - fw_cpu[3][1]) == 0. );

	multi::cuda::managed::array<complex, 2> const in_mng = in_cpu;
	multi::cuda::managed::array<complex, 2> fw_mng(extensions(in_gpu));
	multi::cufft::dft(in_mng, fw_mng, multi::cufft::forward);

	BOOST_TEST( imag(fw_mng[3][2] - fw_cpu[3][2]) == 0. );

	auto fw2_mng = multi::fftw::dft(in_mng, multi::fftw::forward);
	BOOST_TEST( imag(fw2_mng[3][1] - fw_cpu[3][1]) == 0. );

}

BOOST_AUTO_TEST_CASE(cufft_3D_timing, *boost::unit_test::tolerance(0.0001)){

	auto x = std::make_tuple(300, 300, 300);
	{
		multi::array<complex, 3> const in_cpu(x, 10.); 
		BOOST_ASSERT( in_cpu.num_elements()*sizeof(complex) < 2e9 );
		multi::array<complex, 3> fw_cpu(extensions(in_cpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t;  // 1.041691s wall, 1.030000s user + 0.000000s system = 1.030000s CPU (98.9%)
			multi::fftw::dft(in_cpu, fw_cpu, multi::fftw::forward);
			BOOST_TEST( fw_cpu[8][9][10] != 99. );
		}
	}
	{
		multi::cuda::array<complex, 3> const in_gpu(x, 10.); 
		multi::cuda::array<complex, 3> fw_gpu(extensions(in_gpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::fftw::forward);

			BOOST_TEST( static_cast<complex>(fw_gpu[8][9][10]) != 99. );
		}
	}
	{
		multi::cuda::managed::array<complex, 3> const in_gpu(x, 10.); 
		multi::cuda::managed::array<complex, 3> fw_gpu(extensions(in_gpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);
		//	BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);
		//	BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
	}
}

BOOST_AUTO_TEST_CASE(cufft_combinations, *utf::tolerance(0.00001)){

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::clog<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	multi::cuda::array<complex, 4> const in_gpu = in;
	multi::cuda::managed::array<complex, 4> const in_mng = in;

	using std::clog;
	for(auto c : std::vector<std::array<bool, 4>>{
		{false, true , true , true }, 
		{false, true , true , false}, 
		{true , false, false, false}, 
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	}){
		std::clog<<"case "; copy(begin(c), end(c), std::ostream_iterator<bool>{std::clog,", "}); std::clog<<std::endl;
		multi::array<complex, 4> out = in;
		multi::array<complex, 4> in_rw = in;
		[&, _ = watch{"cpu_opl "}]{
			multi::fftw::dft_forward(c, in, out);
		}();
		[&, _ = watch{"cpu_ipl "}]{
			multi::fftw::dft(c, in_rw, multi::fftw::forward);
			BOOST_TEST( abs( static_cast<multi::complex<double>>(in_rw[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );			
		}();
		{
			multi::array<complex, 4> in_rw2 = in;
			[&, _ = watch{"cpu_mov "}]{
				multi::array<complex, 4> const out_mov = multi::fftw::dft_forward(c, std::move(in_rw2));
			//	what(out_mov);
				BOOST_TEST( abs( static_cast<multi::complex<double>>(out_mov[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );			
				BOOST_REQUIRE( is_empty(in_rw2) );
				BOOST_REQUIRE( extensions(out_mov) == extensions(in) );
			}();
		}


		[&, _ = watch{"cpu_new "}]{
			auto const out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST( abs( static_cast<multi::complex<double>>(out_cpy[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );
		}();
		multi::cuda::array<complex, 4> out_gpu(extensions(in_gpu));
		[&, _ = watch{"gpu_opl "}]{
			multi::cufft::dft(c, in_gpu   , out_gpu, multi::cufft::forward);
			BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
		}();
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_ipl "}]{
				multi::cufft::dft(c, in_rw_gpu, multi::cufft::forward);
				BOOST_TEST( abs( static_cast<complex>(in_rw_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::cuda::array<complex, 4> const out_mov = multi::cufft::dft_forward(c, std::move(in_rw_gpu));
			//	BOOST_REQUIRE( in_rw_gpu.empty() );
			//	BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::cuda::array<complex, 4> out_mov = std::move(in_rw_gpu);
				multi::cufft::dft(c, out_mov, multi::cufft::forward);
			//	BOOST_REQUIRE( in_rw_gpu.empty() );
			//	BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		cudaDeviceSynchronize();
		[&, _ = watch{"gpu_new "}]{
			multi::cuda::array<complex, 4> const out_cpy = multi::cufft::dft(c, in_gpu, multi::cufft::forward);
		}();
		multi::cuda::managed::array<complex, 4> out_mng(extensions(in_mng));
		[&, _ = watch{"mng_cld "}]{
			multi::cufft::dft(c, in_mng, out_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
		[&, _ = watch{"mng_hot "}]{
			multi::cufft::dft(c, in_mng   , out_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
		[&, _ = watch{"mng_new "}]{
			auto const out_mng = multi::cufft::dft(c, in_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
	}
	std::clog<<std::endl;

}

BOOST_AUTO_TEST_CASE(cufft_many_3D, *utf::tolerance(0.00001) ){

	auto const in_cpu = []{
		multi::array<complex, 4> ret({45, 18, 32, 16});
		std::generate(
			ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();

	multi::cuda::array<complex, 4> const in = in_cpu;
	multi::cuda::array<complex, 4>       out(extensions(in));

#if 0
	multi::cufft::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), +1);

	multi::array<complex, 4> out_cpu(extensions(in));
	multi::fft::many_dft(begin(unrotated(in_cpu)), end(unrotated(in_cpu)), begin(unrotated(out_cpu)), +1);

	BOOST_TEST( imag( static_cast<complex>(out[5][4][3][2]) - out_cpu[5][4][3][2]) == 0. );
#endif
}

#if 0

BOOST_AUTO_TEST_CASE(cufft_4D, *utf::tolerance(0.00001)){
	auto const in = []{
		multi::array<complex, 3> ret({10, 10, 10});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();

	multi::array<complex, 3> out(extensions(in));
//	multi::fftw::dft({true, false, true}, in, out, multi::fftw::forward);
	multi::fft::many_dft(begin(in<<1), end(in<<1), begin(out<<1), multi::fftw::forward);

	multi::cuda::array<complex, 3> in_gpu = in;
	multi::cuda::array<complex, 3> out_gpu(extensions(in));

//	multi::cufft::dft({true, false, true}, in_gpu, out_gpu, multi::fft::forward);//multi::cufft::forward);	
	multi::cufft::many_dft(begin(in_gpu<<1), end(in_gpu<<1), begin(out_gpu<<1), multi::fftw::forward);
	BOOST_TEST( imag( static_cast<complex>(out_gpu[5][4][3]) - out[5][4][3]) == 0. );	
}



//BOOST_AUTO_TEST_CASE(cu

#if 0
BOOST_AUTO_TEST_CASE(cufft_4D){
	auto const in = []{
		multi::array<complex, 3> ret({10, 10, 10});
		ret[2][3][4] = 99.;
		return ret;
	}();
	multi::array<complex, 3> out(extensions(in));

	multi::fftw::dft({true, true, false}, in, out, multi::fftw::forward);
	
//	auto fwd = multi::fftw::dft({true, true, true, true}, in, out, multi::fftw::forward);
//	BOOST_REQUIRE(in[2][3][4][5] == 99.);
	std::cout << out[9][1][2] << std::endl;
	for(auto i = 0; i != out.num_elements(); ++i) std:cout << (out.data_elements()[i]) <<' ';

#if 0
	multi::cuda::array<complex, 3> in_gpu = in;//[]{
//		multi::cuda::array<complex, 4> ret({10, 10, 10, 10});
//		ret[2][3][4][5] = 99.;
//		return ret;
//	}();
	multi::cuda::array<complex, 3> out_gpu(extensions(in));
	multi::cufft::dft({true, true, false}, in_gpu, out_gpu, multi::cufft::forward);

	std::cout << out_gpu[5][4][3].operator complex() << std::endl;

//	multi::cufft::dft({true, true, true, true}, in_gpu, out_gpu, multi::cufft::forward);
//	multi::cufft::dft({true, true, true, true}, in_gpu, out_gpu, multi::cufft::forward);
#endif

}
#endif
#endif
#endif

#endif
#endif

