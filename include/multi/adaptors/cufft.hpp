// Copyright 2020-2024 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_CUFFTW_HPP
#define MULTI_ADAPTORS_CUFFTW_HPP

#include "../adaptors/../utility.hpp"
#include "../adaptors/../array.hpp"
#include "../adaptors/../config/NODISCARD.hpp"

// #include "../adaptors/cuda.hpp"

#include<tuple>
#include<array>

// #include "../complex.hpp"

#include<thrust/memory.h>  // for raw_pointer_cast

#if not defined(__HIP_ROCclr__)
#include <cufft.h>
#include <cufftXt.h>
#endif

namespace boost{
namespace multi{
namespace cufft{

// cuFFT API errors
static char const* _cudaGetErrorEnum(cufftResult error) {
    switch (error) {
        case CUFFT_SUCCESS:        return "CUFFT_SUCCESS";

        case CUFFT_ALLOC_FAILED:   return "CUFFT_ALLOC_FAILED";
		case CUFFT_EXEC_FAILED:    return "CUFFT_EXEC_FAILED";
		case CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
		case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
		case CUFFT_INVALID_DEVICE: return "CUFFT_INVALID_DEVICE";
        case CUFFT_INVALID_PLAN:   return "CUFFT_INVALID_PLAN";
		case CUFFT_INVALID_SIZE:   return "CUFFT_INVALID_SIZE";
		case CUFFT_INVALID_TYPE:   return "CUFFT_INVALID_TYPE";
		case CUFFT_INVALID_VALUE:  return "CUFFT_INVALID_VALUE";
		case CUFFT_NO_WORKSPACE:   return "CUFFT_NO_WORKSPACE";
		case CUFFT_NOT_IMPLEMENTED:return "CUFFT_NOT_IMPLEMENTED";
		case CUFFT_NOT_SUPPORTED : return "CUFFT_NOT_SUPPORTED";
		case CUFFT_PARSE_ERROR:    return "CUFFT_PARSE_ERROR";
        case CUFFT_SETUP_FAILED:   return "CUFFT_SETUP_FAILED";
		case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
		// case CUFFT_LICENSE_ERROR:  return "CUFFT_LICENSE_ERROR";
    }
    return "<unknown>";
}

#define cufftSafeCall(err) implcufftSafeCall(err, __FILE__, __LINE__)
inline void implcufftSafeCall(cufftResult err, const char *file, const int line) {
	if( CUFFT_SUCCESS != err) {
		std::cerr <<"CUFFT error in file "<< __FILE__ <<", line "<< __LINE__ <<"\nerror "<< err <<": "<<_cudaGetErrorEnum(err)<<"\n";
		//fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n", __FILE__, __LINE__, err, 
        //                        _cudaGetErrorEnum(err));
		cudaDeviceReset()==cudaSuccess?void():assert(0);
		assert(0);
	}
}

class sign {
	int impl_ = 0;

 public:
	sign() = default;
	constexpr sign(int i) : impl_{i} {}
	constexpr operator int() const {return impl_;}
};

constexpr sign forward{CUFFT_FORWARD};
constexpr sign none{0};
constexpr sign backward{CUFFT_INVERSE};
// constexpr sign backward{CUFFT_BACKWARD};

static_assert(forward != none and none != backward and backward != forward, "!");

template<dimensionality_type DD = -1, class Alloc = void*>
struct plan {
	Alloc alloc_;
	::size_t workSize_ = 0;
	void* workArea_;

	using complex_type = cufftDoubleComplex;
	cufftHandle h_;  // TODO(correaa) put this in a unique_ptr
	std::array<std::pair<bool, fftw_iodim64>, DD + 1> which_iodims_{};
	int first_howmany_;

public:
	using allocator_type = Alloc;

	plan(plan&& other) noexcept :
		h_{std::exchange(other.h_, {})},
		which_iodims_{std::exchange(other.which_iodims_, {})},
		first_howmany_{std::exchange(other.first_howmany_, {})},
		workSize_{std::exchange(other.workSize_, {})},
		workArea_{std::exchange(other.workArea_, {})},
		alloc_{std::move(other.alloc_)}
	{}

    template<
		class ILayout, class OLayout, dimensionality_type D = std::decay_t<ILayout>::rank::value,
		class=std::enable_if_t<D == std::decay_t<OLayout>::rank::value>
	>
	plan(std::array<bool, +D> which, ILayout const& in, OLayout const& out, allocator_type const& alloc = {}) : alloc_{alloc} {

		assert(in.sizes() == out.sizes());

		auto const sizes_tuple   = in.sizes();
		auto const istride_tuple = in.strides();
		auto const ostride_tuple = out.strides();

		using boost::multi::detail::get;
		auto which_iodims = std::apply([](auto... elems) {
			return std::array<std::pair<bool, fftw_iodim64>, sizeof...(elems) + 1>{  // TODO(correaa) added one element to avoid problem with gcc 13 static analysis (out-of-bounds)
				std::pair<bool, fftw_iodim64>{
					get<0>(elems),
					fftw_iodim64{get<1>(elems), get<2>(elems), get<3>(elems)}
				}...,
				std::pair<bool, fftw_iodim64>{}
			};
		}, boost::multi::detail::tuple_zip(which, sizes_tuple, istride_tuple, ostride_tuple));

		std::stable_sort(which_iodims.begin(), which_iodims.end() - 1, [](auto const& a, auto const& b){return get<1>(a).is > get<1>(b).is;});

		auto const part = std::stable_partition(which_iodims.begin(), which_iodims.end() - 1, [](auto elem) {return std::get<0>(elem);});

		std::array<fftw_iodim64, D> dims{};
		auto const dims_end         = std::transform(which_iodims.begin(), part,         dims.begin(), [](auto elem) {return elem.second;});

		std::array<fftw_iodim64, D> howmany_dims{};
		auto const howmany_dims_end = std::transform(part, which_iodims.end() -1, howmany_dims.begin(), [](auto elem) {return elem.second;});

		which_iodims_ = which_iodims;
		first_howmany_ = part - which_iodims.begin();

		////////////////////////////////////////////////////////////////////////

		std::array<int, D> istrides{};
		std::array<int, D> ostrides{};
		std::array<int, D> ion{};

		auto const istrides_end = std::transform(dims.begin(), dims_end, istrides.begin(), [](auto elem) {return elem.is;});
		auto const ostrides_end = std::transform(dims.begin(), dims_end, ostrides.begin(), [](auto elem) {return elem.os;});
		auto const ion_end      = std::transform(dims.begin(), dims_end, ion.begin(),      [](auto elem) {return elem.n;});

		int istride = *(istrides_end -1);
		auto inembed = istrides; inembed.fill(0);
		int ostride = *(ostrides_end -1);
		auto onembed = ostrides; onembed.fill(0);

		for(std::size_t i = 1; i != ion_end - ion.begin(); ++i) {
			assert(ostrides[i-1] >= ostrides[i]);
			assert(ostrides[i-1]%ostrides[i]==0);
			onembed[i]=ostrides[i-1]/ostrides[i];
			assert(istrides[i-1]%istrides[i]==0);
			inembed[i]=istrides[i-1]/istrides[i];
		}

		if(dims_end == dims.begin()) {throw std::runtime_error{"no ffts in any dimension is not supported"};}

		while(first_howmany_ < D - 1) {
			int nelems = 1;
			// for(int i = D - 1; i != first_howmany_ + 1; --i) {
			//  nelems *= which_iodims_[i].second.n;
			//  if(
			//      which_iodims_[i - 1].second.is == nelems and
			//      which_iodims_[i - 1].second.os == nelems
			//  ) {
			//      for(int j = i - 1; j != first_howmany_; --j) {
			//          which_iodims_[j].second.n *= which_iodims_[j + 1].second.n;
			//      }
			//  }

			// }
			for(int i = first_howmany_ + 1; i != D; ++i) {nelems *= which_iodims_[i].second.n;}
			if(
				which_iodims_[first_howmany_].second.is == nelems and
				which_iodims_[first_howmany_].second.os == nelems
			) {
				which_iodims_[first_howmany_ + 1].second.n *= which_iodims_[first_howmany_].second.n;
				++first_howmany_;
			} else {
				break;
			}
		}

		if(first_howmany_ == D) {
			if constexpr(std::is_same_v<Alloc, void*>) {
				cufftSafeCall(::cufftPlanMany(
					/*cufftHandle *plan*/ &h_,
					/*int rank*/          dims_end - dims.begin(),
					/*int *n*/            ion.data(),
					/*int *inembed*/      inembed.data(),
					/*int istride*/       istride,
					/*int idist*/         1, //stride(first),
					/*int *onembed*/      onembed.data(),
					/*int ostride*/       ostride,
					/*int odist*/         1, //stride(d_first),
					/*cufftType type*/    CUFFT_Z2Z,
					/*int batch*/         1 //BATCH
				));
			} else {
				cufftSafeCall(cufftCreate(&h_));
				cufftSafeCall(cufftSetAutoAllocation(h_, false));
				cufftSafeCall(cufftMakePlanMany(
					/*cufftHandle *plan*/ h_,
					/*int rank*/          dims_end - dims.begin(),
					/*int *n*/            ion.data(),
					/*int *inembed*/      inembed.data(),
					/*int istride*/       istride,
					/*int idist*/         1, //stride(first),
					/*int *onembed*/      onembed.data(),
					/*int ostride*/       ostride,
					/*int odist*/         1, //stride(d_first),
					/*cufftType type*/    CUFFT_Z2Z,
					/*int batch*/         1, //BATCH
					/*size_t **/          &workSize_
				));
				cufftSafeCall(cufftGetSize(h_, &workSize_));
				workArea_ = ::thrust::raw_pointer_cast(alloc_.allocate(workSize_)); static_assert(sizeof(Alloc) == 1000);
				// auto s = cudaMalloc(&workArea_, workSize_);
				// if(s != cudaSuccess) {throw std::runtime_error{"L212"};}
				cufftSafeCall(cufftSetWorkArea(h_, workArea_));
			}
			if(not h_) {throw std::runtime_error{"cufftPlanMany null"};}
			return;
		}

		std::sort(which_iodims_.begin() + first_howmany_, which_iodims_.begin() + D, [](auto const& a, auto const& b){return get<1>(a).n > get<1>(b).n;});

		if(first_howmany_ <= D - 1) {
			if constexpr(std::is_same_v<Alloc, void*>) {  // NOLINT(bugprone-branch-clone) workaround bug in DeepSource
				cufftSafeCall(::cufftPlanMany(
					/*cufftHandle *plan*/ &h_,
					/*int rank*/          dims_end - dims.begin(),
					/*int *n*/            ion.data(),
					/*int *inembed*/      inembed.data(),
					/*int istride*/       istride,
					/*int idist*/         which_iodims_[first_howmany_].second.is,
					/*int *onembed*/      onembed.data(),
					/*int ostride*/       ostride,
					/*int odist*/         which_iodims_[first_howmany_].second.os,
					/*cufftType type*/    CUFFT_Z2Z,
					/*int batch*/         which_iodims_[first_howmany_].second.n
				));
			} else {
				cufftSafeCall(cufftCreate(&h_));
				cufftSafeCall(cufftSetAutoAllocation(h_, false));
				cufftSafeCall(cufftMakePlanMany(
					/*cufftHandle *plan*/ h_,
					/*int rank*/          dims_end - dims.begin(),
					/*int *n*/            ion.data(),
					/*int *inembed*/      inembed.data(),
					/*int istride*/       istride,
					/*int idist*/         which_iodims_[first_howmany_].second.is,
					/*int *onembed*/      onembed.data(),
					/*int ostride*/       ostride,
					/*int odist*/         which_iodims_[first_howmany_].second.os,
					/*cufftType type*/    CUFFT_Z2Z,
					/*int batch*/         which_iodims_[first_howmany_].second.n,
					/*size_t **/          &workSize_
				));
				cufftSafeCall(cufftGetSize(h_, &workSize_));
				workArea_ = ::thrust::raw_pointer_cast(alloc_.allocate(workSize_));
				cufftSafeCall(cufftSetWorkArea(h_, workArea_));
			}
			if(not h_) {throw std::runtime_error{"cufftPlanMany null"};}
			++first_howmany_;
			return;
		}
		// throw std::runtime_error{"cufft not implemented yet"};
	}

 private:
	plan() = default;
	plan(plan const&) = delete;
	void ExecZ2Z(complex_type const* idata, complex_type* odata, int direction) const{
		cufftSafeCall(::cufftExecZ2Z(h_, const_cast<complex_type*>(idata), odata, direction));
		// cudaDeviceSynchronize();
	}

 public:
	template<class IPtr, class OPtr>
	void execute(IPtr idata, OPtr odata, int direction) {  // TODO(correaa) make const
		if(first_howmany_ == DD) {
			ExecZ2Z((complex_type const*)::thrust::raw_pointer_cast(idata), (complex_type*)::thrust::raw_pointer_cast(odata), direction);
			return;
		}
		if(first_howmany_ == DD - 1) {
			if( which_iodims_[first_howmany_].first) {throw std::runtime_error{"logic error"};}
			for(int i = 0; i != which_iodims_[first_howmany_].second.n; ++i) {
				::cufftExecZ2Z(
					h_,
					const_cast<complex_type*>((complex_type const*)::thrust::raw_pointer_cast(idata + i*which_iodims_[first_howmany_].second.is)),
					                          (complex_type      *)::thrust::raw_pointer_cast(odata + i*which_iodims_[first_howmany_].second.os) ,
					direction
				);
			}
			return;
		}
		if(first_howmany_ == DD - 2) {
			if( which_iodims_[first_howmany_ + 0].first) {throw std::runtime_error{"logic error0"};}
			if( which_iodims_[first_howmany_ + 1].first) {throw std::runtime_error{"logic error1"};}
			if(idata == odata) {throw std::runtime_error{"complicated inplace 2"};}
			for(int i = 0; i != which_iodims_[first_howmany_].second.n; ++i) {
				for(int j = 0; j != which_iodims_[first_howmany_ + 1].second.n; ++j) {
					::cufftExecZ2Z(
						h_,
						const_cast<complex_type*>((complex_type const*)::thrust::raw_pointer_cast(idata + i*which_iodims_[first_howmany_].second.is + j*which_iodims_[first_howmany_ + 1].second.is)),
												  (complex_type      *)::thrust::raw_pointer_cast(odata + i*which_iodims_[first_howmany_].second.os + j*which_iodims_[first_howmany_ + 1].second.os) ,
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
		ExecZ2Z((complex_type const*)::thrust::raw_pointer_cast(idata), (complex_type*)::thrust::raw_pointer_cast(odata), direction);
	}
	template<class I, class O>
	O&& execute_dft(I&& i, O&& o, int direction) const {
		ExecZ2Z(
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(i))),
			const_cast<complex_type*>(reinterpret_cast<complex_type const*>(base(o))),
			direction
		);
		return std::forward<O>(o);
	}

	~plan() {
		if constexpr(not std::is_same_v<Alloc, void*>) {
			if(workSize_ > 0) {alloc_.deallocate(typename std::allocator_traits<Alloc>::pointer((char*)workArea_), workSize_);}
		}
		if(h_) {cufftSafeCall(cufftDestroy(h_));}
	}
	using size_type = int;
	using ssize_type = int;
};

template<dimensionality_type D, class Alloc = void*>
struct cached_plan {
	typename std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc> >::iterator it;

	cached_plan(cached_plan const&) = delete;
	cached_plan(cached_plan&&) = delete;

	cached_plan(std::array<bool, D> which, boost::multi::layout_t<D, boost::multi::size_type> in, boost::multi::layout_t<D, boost::multi::size_type> out, Alloc const& alloc = {}) {
		static thread_local std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc> >& LEAKY_cache = *new std::map<std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>, plan<D, Alloc> >;
		it = LEAKY_cache.find(std::tuple<std::array<bool, D>, multi::layout_t<D>, multi::layout_t<D>>{which, in, out});
		if(it == LEAKY_cache.end()) {it = LEAKY_cache.insert(std::make_pair(std::make_tuple(which, in, out), plan<D, Alloc>(which, in, out, alloc))).first;}
	}
	template<class IPtr, class OPtr>
	void execute(IPtr idata, OPtr odata, int direction) {
		// assert(it != LEAKY_cache.end());
		it->second.execute(idata, odata, direction);
	}
};

template<typename In, class Out, dimensionality_type D = In::rank::value, std::enable_if_t<not multi::has_get_allocator<In>::value, int> =0>
auto dft(std::array<bool, +D> which, In const& i, Out&& o, int s)
->decltype(cufft::cached_plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), s), std::forward<Out>(o)) {
	return cufft::cached_plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), s), std::forward<Out>(o); }

template<typename In, class Out, dimensionality_type D = In::rank::value, std::enable_if_t<    multi::has_get_allocator<In>::value, int> =0>
auto dft(std::array<bool, +D> which, In const& i, Out&& o, int s)
->decltype(cufft::cached_plan<D /*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/ >{which, i.layout(), o.layout()/*, i.get_allocator()*/}.execute(i.base(), o.base(), s), std::forward<Out>(o)) {
	return cufft::cached_plan<D /*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/ >{which, i.layout(), o.layout()/*, i.get_allocator()*/}.execute(i.base(), o.base(), s), std::forward<Out>(o); }

template<typename In, class Out, dimensionality_type D = In::rank::value>//, std::enable_if_t<not multi::has_get_allocator<In>::value, int> =0>
auto dft_forward(std::array<bool, +D> which, In const& i, Out&& o)  -> Out&& {
//->decltype(cufft::plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::forward), std::forward<Out>(o)) {
	return cufft::cached_plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::forward), std::forward<Out>(o); }

// template<typename In, class Out, dimensionality_type D = In::rank::value, class = typename In::allocator_type, std::enable_if_t<    multi::has_get_allocator<In>::value, int> =0>
// auto dft_forward(std::array<bool, +D> which, In const& i, Out&& o) -> Out&& {
// //->decltype(cufft::plan<D, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char> >{which, i.layout(), o.layout(), i.get_allocator()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o)) {
//  return cufft::cached_plan<D/*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/>{which, i.layout(), o.layout()/*, i.get_allocator()*/}.execute(i.base(), o.base(), cufft::forward), std::forward<Out>(o); }

template<typename In, class Out, dimensionality_type D = In::rank::value>//, std::enable_if_t<not multi::has_get_allocator<In>::value, int> =0>
auto dft_backward(std::array<bool, +D> which, In const& i, Out&& o) -> Out&& {
//->decltype(cufft::plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o)) {
	return cufft::cached_plan<D>{which, i.layout(), o.layout()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o); }

// template<typename In, class Out, dimensionality_type D = In::rank::value, class = typename In::allocator_type, std::enable_if_t<    multi::has_get_allocator<In>::value, int> =0>
// auto dft_backward(std::array<bool, +D> which, In const& i, Out&& o) -> Out&& {
// //->decltype(cufft::plan<D, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char> >{which, i.layout(), o.layout(), i.get_allocator()}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o)) {
//  return cufft::cached_plan<D/*, typename std::allocator_traits<typename In::allocator_type>::rebind_alloc<char>*/>{which, i.layout(), o.layout()/*, i.get_allocator()*/}.execute(i.base(), o.base(), cufft::backward), std::forward<Out>(o); }

template<typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
R dft(In const& i, int s) {
	static_assert(std::is_trivially_default_constructible<typename In::element_type>{});
	R ret(extensions(i), get_allocator(i));
	cufft::dft(i, ret, s);
	// if(cudaDeviceSynchronize() != cudaSuccess) throw std::runtime_error{"Cuda error: Failed to synchronize"};
	return ret;
}

template <class Array, std::size_t... Ns>
constexpr auto array_tail_impl(Array const& t, std::index_sequence<Ns...>) {
	return std::array<typename Array::value_type, std::tuple_size<Array>{} - 1>{std::get<Ns + 1>(t)...};
}

template<class Array>
constexpr auto array_tail(Array const& t)
->decltype(array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>())) {
	return array_tail_impl(t, std::make_index_sequence<std::tuple_size<Array>{} - 1>()); }

// template<typename In, class Out, std::size_t D = In::dimensionality, std::enable_if_t<(D>1), int> = 0>
// auto dft_forward(std::array<bool, +D> which, In const& i, Out&& o)
// ->decltype(dft(which, i, std::forward<Out>(o), cufft::forward)) {
//  return dft(which, i, std::forward<Out>(o), cufft::forward); }

// template<typename In, class Out, std::size_t D = In::dimensionality, std::enable_if_t<(D>1), int> = 0>
// auto dft_backward(std::array<bool, +D> which, In const& i, Out&& o)
// ->decltype(dft(which, i, std::forward<Out>(o), cufft::backward)) {
//  return dft(which, i, std::forward<Out>(o), cufft::backward); }

template<typename In,  std::size_t D = In::dimensionality>
NODISCARD("when passing a const argument")
auto dft(std::array<bool, D> which, In const& i, int sign)->std::decay_t<decltype(
dft(which, i, typename In::decay_type(extensions(i), get_allocator(i)), sign))>{return
dft(which, i, typename In::decay_type(extensions(i), get_allocator(i)), sign);}

template<typename In,  std::size_t D = In::dimensionality>
auto dft(std::array<bool, D> which, In&& i, int sign)
->decltype(dft(which, i, i, sign), std::forward<In>(i)){
	return dft(which, i, i, sign), std::forward<In>(i);}

template<typename Array, typename A> NODISCARD("when passing a const argument")
auto dft_forward(Array arr, A const& a) 
->decltype(cufft::dft(arr, a, cufft::forward)){
	return cufft::dft(arr, a, cufft::forward);}

// template<typename Array, dimensionality_type D> NODISCARD("when passing a const argument")
// auto dft_forward(Array arr, multi::cuda::array<std::complex<double>, D>&& a) 
// ->decltype(cufft::dft(arr, a, cufft::forward), multi::cuda::array<std::complex<double>, D>{}){//assert(0);
//  return cufft::dft(arr, a, cufft::forward), std::move(a);}

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
