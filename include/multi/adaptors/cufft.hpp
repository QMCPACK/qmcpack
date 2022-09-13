// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_CUFFTW_HPP
#define MULTI_ADAPTORS_CUFFTW_HPP

#include "../config/MARK.hpp"

#include "../adaptors/../utility.hpp"
#include "../adaptors/../array.hpp"
#include "../adaptors/../config/NODISCARD.hpp"

#include "../adaptors/cuda.hpp"

#include<numeric>

#include<tuple> // std::apply
#include<array>

#include<vector>

#include "../complex.hpp"

#include<cufft.h>

namespace boost{
namespace multi{
namespace cufft{

class sign {
	int impl_;

 public:
	sign() = default;
	constexpr sign(int i) : impl_{i} {}
	constexpr operator int() const {return impl_;}
};

constexpr sign forward{CUFFT_FORWARD};
constexpr sign none{0};
constexpr sign backward{CUFFT_INVERSE};

static_assert(forward != none and none != backward and backward != forward, "!");

class plan {
	using complex_type = cufftDoubleComplex;
	complex_type const* idata_     = nullptr;
	complex_type*       odata_     = nullptr;
	int                 direction_ = 0;
	cufftHandle h_;

public:
	template<multi::dimensionality_type D>
	static
	std::map<std::tuple<std::array<int, D>, std::array<int, D>, int, int, std::array<int, D>, int, int, int>, cufftHandle>&
	cache() {
		static std::map<std::tuple<std::array<int, D>, std::array<int, D>, int, int, std::array<int, D>, int, int, int>, cufftHandle> cache_;
		return cache_;
	}

 private:
	plan() = default;
	plan(plan const&) = delete;
	plan(plan&& other)
	: idata_{std::exchange(other.idata_, nullptr)}
	, odata_{std::exchange(other.odata_, nullptr)}
	, direction_{std::exchange(other.direction_, 0)}
	, h_{std::exchange(other.h_, {})}
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
		cudaDeviceSynchronize();
	}
	void swap(plan& other) {
		using std::swap;
		swap(idata_, other.idata_);
		swap(odata_, other.odata_);
		swap(direction_, other.direction_);
		swap(h_, other.h_);
	}

 public:
	thread_local static int tl_execute_count;
	plan& operator=(plan other) {swap(other); return *this;}
	void operator()() const {ExecZ2Z(idata_, odata_, direction_);}
	template<class I, class O>
	O&& execute_dft(I&& i, O&& o, int direction) const {
		ExecZ2Z(
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(i))),
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(o))),
			direction
		);
		return std::forward<O>(o);
	}
	template<class I, class O>
	void execute_dft(I&& i, O&& o) const{execute_dft(std::forward<I>(i), std::forward<O>(o), direction_);}
	~plan() {
		MULTI_MARK_SCOPE("cufft plan dtor");
	//	if(h_) cufftDestroy(h_);
	}
	using size_type = int;
	using ssize_type = int;

	template<class I, class O, //std::enable_if_t<(I::dimensionality < 4), int> =0,
		dimensionality_type D = I::dimensionality,  
		typename = decltype(raw_pointer_cast(base(std::declval<I const&>())), reinterpret_cast<complex_type*      >(raw_pointer_cast(base(std::declval<O&>()))))
	>
	plan(I const& i, O&& o, sign s)
	: idata_{                          reinterpret_cast<complex_type const*>(raw_pointer_cast(base(i))) }
	, odata_{const_cast<complex_type*>(reinterpret_cast<complex_type*      >(raw_pointer_cast(base(o))))}
	, direction_{s}
	{
		MULTI_MARK_SCOPE("cufft plan ctor");

		assert( I::dimensionality < 4 );
		assert( CUFFT_FORWARD == s or CUFFT_INVERSE == s or s == 0 );
		assert( sizes(i) == sizes(o) );

		auto ion      = std::apply([](auto... t){return std::array< size_type, D>{static_cast< size_type>(t)...};}, sizes  (i));
		auto istrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(i));
		auto ostrides = std::apply([](auto... t){return std::array<ssize_type, D>{static_cast<ssize_type>(t)...};}, strides(o));

		std::array<std::tuple<int, int, int>, I::dimensionality> ssn;
		for(std::size_t i = 0; i != ssn.size(); ++i) {ssn[i] = std::make_tuple(istrides[i], ostrides[i], ion[i]);}
		std::sort(ssn.begin(), ssn.end(), std::greater<>{});

		for(std::size_t i = 0; i != ssn.size(); ++i) {
			istrides[i] = std::get<0>(ssn[i]);
			ostrides[i] = std::get<1>(ssn[i]);
			ion[i]      = std::get<2>(ssn[i]);
		}// = std::tuple<int, int, int>(istrides[i], ostrides[i], ion[i]);

		int istride = istrides.back();
		auto inembed = istrides; inembed.fill(0);
		int ostride = ostrides.back();
		auto onembed = ostrides; onembed.fill(0);
		for(std::size_t i = 1; i != onembed.size(); ++i) {
			assert(ostrides[i-1] >= ostrides[i]); // otherwise ordering is incompatible
			assert(ostrides[i-1]%ostrides[i]==0);
			onembed[i]=ostrides[i-1]/ostrides[i]; //	assert( onembed[i] <= ion[i] );
			assert(istrides[i-1]%istrides[i]==0);
			inembed[i]=istrides[i-1]/istrides[i]; //	assert( inembed[i] <= ion[i] );
		}

		direction_ = s;
		idata_ =                           reinterpret_cast<complex_type const*>(raw_pointer_cast(base(i))) ;
		odata_ = const_cast<complex_type*>(reinterpret_cast<complex_type*      >(raw_pointer_cast(base(o))));

		auto it = cache<D>().find(std::make_tuple(ion, inembed, istride, 1, onembed, ostride, 1, 1));
		if(it != cache<D>().end()) {
			h_ = it->second;
		}else{
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
			)) {
				case CUFFT_SUCCESS        : break;// "cuFFT successfully executed the FFT plan."
				case CUFFT_ALLOC_FAILED   : throw std::runtime_error{"CUFFT failed to allocate GPU memory."};
				case CUFFT_INVALID_VALUE  : throw std::runtime_error{"At least one of the parameters idata, odata, and direction is not valid."};
				case CUFFT_INTERNAL_ERROR : throw std::runtime_error{"Used for all internal driver errors."};
				case CUFFT_SETUP_FAILED   : throw std::runtime_error{"The cuFFT library failed to initialize."};
				case CUFFT_INVALID_SIZE   : throw std::runtime_error{"The user specifies an unsupported FFT size."};
				default                   : throw std::runtime_error{"cufftPlanMany unknown error"};
			}
			cache<D>().insert(std::make_pair(std::make_tuple(ion, inembed, istride, 1, onembed, ostride, 1, 1), h_));
		}
		assert(h_);
	}
#ifndef __INTEL_COMPILER
	template<class It1, class It2, dimensionality_type D = decltype(*It1{})::dimensionality>
	static auto many(It1 first, It1 last, It2 d_first, int sign = 0, unsigned = 0)
	->std::decay_t<decltype(const_cast<complex_type*>(reinterpret_cast<complex_type*>(raw_pointer_cast(base(d_first)))), std::declval<plan>())>
#else
	template<class It1, class It2, dimensionality_type D = decltype(*It1{})::dimensionality, 
		typename TT = decltype(const_cast<complex_type*>(reinterpret_cast<complex_type*>(raw_pointer_cast(It2{}.base())))
	>
	static auto many(It1 first, It1 last, It2 d_first, int sign = 0, unsigned = 0)
#endif
	{
		MULTI_MARK_SCOPE("cufft plan many factory");

		assert( CUFFT_FORWARD == sign or CUFFT_INVERSE == sign or sign == 0 );
		assert( sizes(*first) == sizes(*d_first) );

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
		for(std::size_t i = 1; i != onembed.size(); ++i) {
			assert(ostrides[i-1] >= ostrides[i]); // otherwise ordering is incompatible
			assert(ostrides[i-1]%ostrides[i]==0);
			onembed[i]=ostrides[i-1]/ostrides[i]; //	assert( onembed[i] <= ion[i] );
			assert(istrides[i-1]%istrides[i]==0);
			inembed[i]=istrides[i-1]/istrides[i]; //	assert( inembed[i] <= ion[i] );
		}

		plan ret;
		ret.direction_ = sign;
		ret.idata_ =                           reinterpret_cast<complex_type const*>(raw_pointer_cast(  first.base())) ;
		ret.odata_ = const_cast<complex_type*>(reinterpret_cast<complex_type*      >(raw_pointer_cast(d_first.base())));

		auto it = cache<D>().find(std::make_tuple(ion, inembed, istride, stride(first), onembed, ostride, stride(d_first), last - first));
		if(it != cache<D>().end()) {
			ret.h_ = it->second;
		} else {
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
			)) {
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
				default                   : throw std::logic_error{"cufftPlanMany unknown error"};
			}
			cache<D>().insert(std::make_pair(std::make_tuple(ion, inembed, istride, stride(first), onembed, ostride, stride(d_first), last - first), ret.h_));
		}
		assert(ret.h_);
		return ret;
	}
};

thread_local int plan::tl_execute_count = 0;

template<typename In, class Out>
auto dft(In const& i, Out&& o, int s)
->decltype(cufft::plan{i, o, s}(), std::forward<Out>(o)) {
	return cufft::plan{i, o, s}(), std::forward<Out>(o); }

template<typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
R dft(In const& i, int s) {
	static_assert(std::is_trivially_default_constructible<typename In::element_type>{}, "!");
	R ret(extensions(i), get_allocator(i));
	cufft::dft(i, ret, s);
	if(cudaDeviceSynchronize() != cudaSuccess) throw std::runtime_error{"Cuda error: Failed to synchronize"};
	return ret;
}

#ifndef __INTEL_COMPILER
template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, sign s)
->decltype(plan::many(first, last, d_first, s)(), d_first + (last - first)) {
	return plan::many(first, last, d_first, s)(), d_first + (last - first); }
#else
template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, sign s)
->decltype(plan::many(first, last, d_first, s)(), d_first + (last - first)) {
	return plan::many(first, last, d_first, s)(), d_first + (last - first); }
#endif

template<typename In, class Out,  std::size_t D = In::dimensionality, std::enable_if_t<(D==1), int> = 0>
Out&& dft(std::array<bool, D> which, In const& i, Out&& o, int s) {
	if(which[0]) return cufft::dft(i, std::forward<Out>(o), s);
	else return std::forward<Out>(std::forward<Out>(o) = i);
}

template <class Array, std::size_t... Ns>
constexpr auto array_tail_impl(Array const& t, std::index_sequence<Ns...>) {
	return std::array<typename Array::value_type, std::tuple_size<Array>{} - 1>{std::get<Ns + 1>(t)...};
}

template<class Array>
constexpr auto array_tail(Array const& t)
->decltype(array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>())) {
	return array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>()); }

template<typename In, class Out, std::size_t D = In::dimensionality, std::enable_if_t<(D>1), int> = 0>
auto dft(std::array<bool, D> which, In const& i, Out&& o, int s)
->decltype(many_dft(i.begin(), i.end(), o.begin(), s),std::forward<Out>(o))
{
	assert(extension(i) == extension(o));
	auto ff = std::find(begin(which)+1, end(which), false);
	if(which[0] == true) {
		if(ff==end(which)) {cufft::dft(i, std::forward<Out>(o), s);}
		else {
			auto const n = ff - which.begin();
			std::rotate(begin(which), ff, end(which));
			// TODO(correaa) : make this more elegant
			switch(n) {
				case 0: dft(which, i                              , o                              , s); break;
				case 1: dft(which, i.rotated()                    , o.rotated()                    , s); break;
				case 2: dft(which, i.rotated().rotated()          , o.rotated().rotated()          , s); break;
				case 3: dft(which, i.rotated().rotated().rotated(), o.rotated().rotated().rotated(), s); break;
				default: assert(0);
			}
		}
	} else if(which[0]==false) {
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
				for(auto d = 0; d != D - 1; ++d) {
					switch(d) {
						case 0: if( (size(i                              ) < n_min) and (tail[d] == false)) {n_min = size(i                              ); d_min = d;} break;
						case 1: if( (size(i.rotated()                    ) < n_min) and (tail[d] == false)) {n_min = size(i.rotated()                    ); d_min = d;} break;
						case 2: if( (size(i.rotated().rotated()          ) < n_min) and (tail[d] == false)) {n_min = size(i.rotated().rotated()          ); d_min = d;} break;
						case 3: if( (size(i.rotated().rotated().rotated()) < n_min) and (tail[d] == false)) {n_min = size(i.rotated().rotated().rotated()); d_min = d;} break;
						default: assert(0);
					}
				//  if((size(i<<d) < n_min) and (tail[d]==false)) {n_min = size(i<<d); d_min = d;}
				}
				if( d_min!=0 ) {
					std::rotate(which.begin(), which.begin()+d_min, which.end());
					switch(d_min) {
						case 0: dft(which, i, o, s); break;
						case 1: dft(which, i.rotated()                    , o.rotated()                    , s); break;
						case 2: dft(which, i.rotated().rotated()          , o.rotated().rotated()          , s); break;
						case 3: dft(which, i.rotated().rotated().rotated(), o.rotated().rotated().rotated(), s); break;
						default: assert(0);
					}
				//  dft(which, i<<d_min, o<<d_min, s);
				} else {
					if(base(i) == base(o) and i.layout() != o.layout()){
						auto const tmp = +i;
						for(auto idx : extension(i)) cufft::dft(tail, tmp[idx], o[idx], s);
					}else for(auto idx : extension(i)){
						MULTI_MARK_SCOPE("cufft inner loop");
						cufft::dft(tail, i[idx], o[idx], s);
					}
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
#endif
