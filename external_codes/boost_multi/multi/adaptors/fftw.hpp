#ifdef COMPILATION_INSTRUCTIONS
$CXX $0 -o $0x -lcudart `pkg-config --libs fftw3` -lboost_timer -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2019

#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP

#include<fftw3.h> // external fftw3 library
	
#include "../../multi/utility.hpp"
#include "../../multi/array.hpp"

#include "../../multi/config/NODISCARD.hpp"

#include<cmath>
#include<complex>
#include<memory>
#include<numeric> // accumulate

#if HAVE_FFTW3_THREADS
#include <thread>
#endif

#include<experimental/tuple> // experimental::apply

#include<utility>
#include<type_traits>

namespace boost{
namespace multi{

namespace fftw{
//	template<class T> auto alignment_of(T* p){return ::fftw_alignment_of((double*)p);}
#if  __cpp_lib_as_const >= 201510
using std::as_const;
#else
template<class T> constexpr std::add_const_t<T>& as_const(T& t) noexcept{return t;}
#endif

}

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

#if(__cpp_if_constexpr>=201606)
//https://stackoverflow.com/a/35110453/225186
template<class T> constexpr std::remove_reference_t<T> _constx(T&&t){return t;}
#define logic_assert(C, M) \
	if constexpr(noexcept(_constx(C))) static_assert((C), M); else assert((C)&& M);
#else
#define logic_assert(ConditioN, MessagE) assert(ConditioN && MessagE);
#endif

template<typename It1, class It2, std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign, unsigned flags = FFTW_ESTIMATE)
->decltype(reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(d_first))), fftw_plan{}){
	assert(sizes(*first)==sizes(*d_first));
	auto ion      = to_array<int>(sizes(*first));

	assert(strides(*first) == strides(*last));
	auto istrides = to_array<int>(strides(*first));
	auto ostrides = to_array<int>(strides(*d_first));

//	auto inelemss = to_array<int>(first->nelemss());
//	auto onelemss = to_array<int>(d_first->nelemss());

	std::array<std::tuple<int, int, int>, std::decay_t<decltype(*It1{})>::dimensionality> ssn;
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

	return ::fftw_plan_many_dft(
		/*int rank*/ ion.size(), 
		/*const int* n*/ ion.data(),
		/*int howmany*/ last - first,
		/*fftw_complex * in */ reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(static_cast<std::complex<double> const*>(base(first)))), 
		/*const int *inembed*/ inembed.data(),
		/*int*/ istride, 
		/*int idist*/ stride(first),
		/*fftw_complex * out */ reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(d_first))),
		/*const int *onembed*/ onembed.data(),
		/*int*/ ostride, 
		/*int odist*/ stride(d_first),
		/*int*/ sign, /*unsigned*/ flags
	);
}

template<class In, class Out, std::size_t D = std::decay_t<In>::dimensionality,
typename = std::enable_if_t<D == std::decay_t<Out>::dimensionality>,
typename = decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(base(std::declval<Out&>()))))
>
fftw_plan fftw_plan_dft(std::decay_t<std::array<bool, D>> which, In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE){
	using multi::sizes; using multi::strides; assert(sizes(in) == sizes(out));
	auto ion      = to_array<ptrdiff_t>(sizes(in));
	auto istrides = to_array<ptrdiff_t>(strides(in));
	auto ostrides = to_array<ptrdiff_t>(strides(out));
	std::array<fftw_iodim64, D> dims   ; auto l_dims = dims.begin();
	std::array<fftw_iodim64, D> howmany; auto l_howmany = howmany.begin();
	for(int i = 0; i != D; ++i) 
		(which[i]?*l_dims++:*l_howmany++) = fftw_iodim64{ion[i], istrides[i], ostrides[i]};
	return fftw_plan_guru64_dft(
		/*int rank*/ sign?(l_dims - dims.begin()):0, 
		/*const fftw_iodim64 *dims*/ dims.data(), 
		/*int howmany_rank*/ l_howmany - howmany.begin(),
		/*const fftw_iodim *howmany_dims*/ howmany.data(), //nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const *>(base(in)))), 
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(base(out))),
		sign, flags// | FFTW_ESTIMATE
	);
}

template<class To, class From, std::enable_if_t<std::is_convertible<From, To>{},int> =0>
To implicit_cast(From&& f){return static_cast<To>(f);}

template<class In, class Out, dimensionality_type D = In::dimensionality, typename = decltype(reinterpret_cast<fftw_complex*>(implicit_cast<std::complex<double>*>(base(std::declval<Out&>()))))>
auto fftw_plan_dft(In const& in, Out&& out, int s, unsigned flags = FFTW_ESTIMATE){
	static_assert( D == std::decay_t<Out>::dimensionality , "!");
	using multi::sizes; using multi::strides; assert(sizes(in) == sizes(out));
	auto 
		ion      = to_array<ptrdiff_t>(sizes(in)), 
		istrides = to_array<ptrdiff_t>(strides(in)),
		ostrides = to_array<ptrdiff_t>(strides(out))
	;
	std::array<fftw_iodim64, D> dims;
	for(int i=0; i!=D; ++i) dims[i] = {ion[i], istrides[i], ostrides[i]};
	return fftw_plan_guru64_dft(
		/*int rank*/ s?D:0,
		/*const fftw_iodim64 *dims*/ dims.data(),
		/*int howmany_rank*/ 0,
		/*const fftw_iodim *howmany_dims*/ nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(base(in)))), 
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(implicit_cast<std::complex<double>*>(base(out))),
		s, flags
	);
}

//std::complex<double> const* base(std::complex<double> const& c){return &c;}

namespace fftw{

#if HAVE_FFTW3_THREADS
void initialize_threads(){int good = fftw_init_threads(); assert(good);}
#else
void initialize_threads(){}
#endif

class plan{
	plan() : impl_{nullptr, &fftw_destroy_plan}{}
	std::unique_ptr<std::remove_pointer_t<fftw_plan>, decltype(&fftw_destroy_plan)> impl_;
public:
	plan(plan const&) = delete;//default;
	plan(plan&&) = default;
	template<typename... As, 
		typename = decltype(fftw_plan_dft(std::declval<As&&>()...))
	> plan(As&&... as) : impl_{fftw_plan_dft(std::forward<As>(as)...), &fftw_destroy_plan}{
		assert(impl_);
	}
	template<typename... As>
	static auto many(As&&... as)
	->std::decay_t<decltype(fftw_plan_many_dft(std::forward<As>(as)...) , std::declval<plan>())>
	{
		plan r; r.impl_.reset(fftw_plan_many_dft(std::forward<As>(as)...)); return r;
	}
private:
	void execute() const{fftw_execute(impl_.get());}
	template<class I, class O>
	void execute_dft(I&& i, O&& o) const{
		::fftw_execute_dft(impl_.get(), const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(base(i)))), reinterpret_cast<fftw_complex*>(static_cast<std::complex<double>*>(base(o))));
	}
	template<class I, class O> void execute(I&& i, O&& o) const{execute_dft(std::forward<I>(i), std::forward<O>(o));}
	friend void execute(plan const& p){p.execute();}
public:
	plan& operator=(plan&&) = default;
	plan& operator=(plan const&) = delete;//default;
	void operator()() const{execute();} // http://www.fftw.org/fftw3_doc/Thread-safety.html#Thread-safety
	template<class I, class O> void operator()(I&& i, O&& o) const{return execute(std::forward<I>(i), std::forward<O>(o));}
	double cost() const{return fftw_cost(impl_.get());}
	auto flops() const{
		struct{double add; double mul; double fma; operator double() const{return add + mul + 2*fma;}} r;
		fftw_flops(impl_.get(), &r.add, &r.mul, &r.fma);
		return r;
	}
#if HAVE_FFTW3_THREADS
public:
	static void make_thread_safe(){
		fftw_make_planner_thread_safe(); // needs linking to -lfftw3_threads, requires FFTW-3.3.6 or greater
		is_thread_safe_ = true;
	}
	static int with_nthreads(int n){fftw_plan_with_nthreads(n); nthreads_ = n; return n;}
	static int with_nthreads(){
		int n=std::thread::hardware_concurrency(); return with_nthreads(n?n:2);
	//	return with_nthreads(std::thread::hardware_concurrency()?:2);
	}
	static bool is_thread_safe(){return is_thread_safe_;}
	static bool nthreads(){return nthreads_;}
private:
	static bool is_thread_safe_;
	static int nthreads_;
	static bool initialized_threads_;
#else
	static constexpr bool is_thread_safe(){return false;}
	static constexpr bool nthreads(){return 1;}
	static constexpr int with_nthreads(){return 1;}
#endif
};

#if HAVE_FFTW3_THREADS
bool plan::is_thread_safe_ = (plan::make_thread_safe(), true);
int plan::nthreads_ = (initialize_threads(), with_nthreads());
#endif

//enum sign: decltype(FFTW_FORWARD){forward = FFTW_FORWARD, none = 0, backward = FFTW_BACKWARD };
using sign = int;
constexpr sign forward = FFTW_FORWARD;
constexpr sign none = 0;
constexpr sign backward = FFTW_BACKWARD; 

static_assert( forward != none and none != backward and backward != forward, "!");

enum strategy: decltype(FFTW_ESTIMATE){ estimate = FFTW_ESTIMATE, measure = FFTW_MEASURE };


template<class In, class Out>
auto dft(In const& i, Out&& o, int s)
->decltype(fftw::plan{i, o, s}(), std::forward<Out>(o)){
	return fftw::plan{i, o, s}(), std::forward<Out>(o);}

template<class In, class Out, std::size_t = std::decay_t<In>::dimensionality>
Out&& transpose(In const& i, Out&& o){
	return dft(i, std::forward<Out>(o), fftw::none);
}

//template<class In, class Out>//, std::size_t D = >
//decltype(auto) dft(std::array<bool, In::dimensionality> which, In const& i, Out&& o, sign s)
//->decltype(plan{which, i, o, s}(), std::forward<Out>(o)){
//{	return plan{which, i, o, s}(), std::forward<Out>(o);}

//template<class In, class Out, std::size_t W>//class Array = std::array<bool, In::dimensionality> >//, std::size_t D = >
//decltype(auto) dft(std::array<bool, W> which, In const& i, Out&& o, sign s)
//->decltype(plan{which, i, o, s}(), std::forward<Out>(o)){
//{	return plan{which, i, o, s}(), std::forward<Out>(o);}

using std::decay_t;

template<class In, class Out, std::size_t D = In::dimensionality>
auto dft(std::array<bool, D> which, In const& i, Out&& o, sign s)
->decltype(plan{which, i, o, s}(), std::forward<Out>(o)){
	return plan{which, i, o, s}(), std::forward<Out>(o);}


/*
template<dimensionality_type R, class In, class Out, std::size_t D = std::decay_t<In>::dimensionality>
Out&& dft(In const& i, Out&& o, sign s){
	static_assert( R <= D , "dimension of transpformation cannot be larger than total dimension" );
	std::array<bool, D> which; std::fill(std::fill_n(begin(which), R, false), end(which), true);
	plan{which, i, o, s}();//(i, std::forward<Out>(o)); 
	return std::forward<Out>(o);
}
*/

template<typename In, class Out, std::size_t D = In::dimensionality, std::size_t = std::decay_t<Out>::dimensionality>
auto dft(std::array<sign, D> w, In const& i, Out&& o){
	std::array<bool, D> fwd, /*non,*/ bwd;

	std::transform(begin(w), end(w), begin(fwd), [](auto e){return e==FFTW_FORWARD;});
	dft(fwd, i, o, fftw::forward);

//	std::transform(begin(w), end(w), begin(non), [](auto e){return e==sign::none;});
	std::transform(begin(w), end(w), begin(bwd), [](auto e){return e==FFTW_BACKWARD;}); 
	if(std::accumulate(begin(bwd), end(bwd), false)) dft(bwd, o, o, FFTW_BACKWARD);

	return std::forward<Out>(o);
}

template<typename It1, typename It2>
auto many_dft(It1 first, It1 last, It2 d_first, int sign)
->decltype(plan::many(first, last, d_first, sign)(), d_first + (last - first)){
	return plan::many(first, last, d_first, sign)(), d_first + (last - first);}

template<typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
auto dft(In const& i, sign s)
->std::decay_t<decltype(dft(i, R(extensions(i), get_allocator(i)), s))>{
	return dft(i, R(extensions(i), get_allocator(i)), s);}

template<typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
R transpose(In const& i){
	return transpose(i, R(extensions(i), get_allocator(i)));
}

template<typename T, dimensionality_type D, class... Args>
decltype(auto) rotate(multi::array<T, D, Args...>& i, int = 1){
	multi::array_ref<T, D, typename multi::array<T, D, Args...>::element_ptr> before(data_elements(i), extensions(i));
//	std::cout << "1. "<< size(i) <<' '<< size(rotated(i)) << std::endl;
	i.reshape(extensions(rotated(before) ));
//	auto x = extensions(i);
//	std::cout << "2. "<< size(i) <<' '<< size(rotated(i)) << std::endl;
	fftw::dft(before, i, fftw::none);
//	std::cout << "3. "<< size(i) <<' '<< size(rotated(i)) << std::endl;
	return i;
//	assert( extensions(i) == x );
//	return i;
}

template<typename In, std::size_t D = In::dimensionality, typename R = multi::array<typename In::element_type, D, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when first argument is const")
auto dft(std::array<bool, D> which, In const& i, sign s)
->std::decay_t<decltype(fftw::dft(which, i, R(extensions(i), get_allocator(i)), s))>{
	return fftw::dft(which, i, R(extensions(i), get_allocator(i)), s);}

template<typename In, std::size_t D = std::decay_t<In>::dimensionality>
auto dft(std::array<bool, std::decay_t<In>::dimensionality> which, In&& i, sign s)
->decltype(dft(which, i, i, s), std::forward<In>(i)){
	return dft(which, i, i, s), std::forward<In>(i);}

/*
template<typename In, std::size_t D = In::dimensionality, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when second argument is const")
R dft(std::array<sign, D> which, In const& i){
	return dft(which, i, R(extensions(i), get_allocator(i)));
}*/

template<typename In, std::size_t D = In::dimensionality, typename R = multi::array<typename In::element_type, D, decltype(get_allocator(std::declval<In>()))>>
void dft(std::array<bool, D> which, In const& i) = delete;

template<dimensionality_type Rank, typename In, typename R = multi::array<typename In::element_type, In::dimensionality, decltype(get_allocator(std::declval<In>()))>>
NODISCARD("when second argument is const")
R dft(In const& i, sign s){
	return dft<Rank>(i, R(extensions(i), get_allocator(i)), s);
}

/*
template<typename T, dimensionality_type D, class... As, typename R = multi::array<T, D, As...>>//typename std::decay_t<In>::element_type, std::decay_t<In>::dimensionality>>
NODISCARD("when first argument can be destroyed")
R dft(multi::array<T, D, As...>&& i, sign s){
//	R ret(extensions(i), get_allocator(i));
//	plan{i, ret, s, static_cast<unsigned>(fftw::estimate) | FFTW_DESTROY_INPUT}();//(i, ret); // to do destroy input for move iterators
	return R{std::move(dft(i, s))};
}
*/

template<typename T> decltype(auto) dft(std::initializer_list<T> il, sign s){return dft(multi::array<T, 1>(il), s);}
template<typename T> decltype(auto) dft(std::initializer_list<std::initializer_list<T>> il, sign s){return dft(multi::array<T, 2>(il), s);}

template<typename... A> auto            dft_forward(A&&... a)
->decltype(fftw::dft(std::forward<A>(a)..., fftw::forward)){
	return fftw::dft(std::forward<A>(a)..., fftw::forward);}

template<typename Array, typename A>
NODISCARD("when input argument is read only")
auto dft_forward(Array which, A const& a)
->decltype(fftw::dft(which, a, fftw::forward)){
	return fftw::dft(which, a, fftw::forward);}

template<typename Array, typename A> 
NODISCARD("when input argument is read only")
auto dft_forward(A const& a)
->decltype(fftw::dft(a, fftw::forward)){
	return fftw::dft(a, fftw::forward);}

template<typename... A> auto            dft_backward(A&&... a)
->decltype(dft(std::forward<A>(a)..., fftw::backward)){
	return dft(std::forward<A>(a)..., fftw::backward);}

template<typename T, typename... As> decltype(auto) dft_forward(As&... as, std::initializer_list<T> il){return dft_forward(std::forward<As>(as)..., multi::array<T, 1>(il));}
template<typename T, typename... As> decltype(auto) dft_forward(As&... as, std::initializer_list<std::initializer_list<T>> il){return dft_forward(std::forward<As>(as)..., multi::array<T, 2>(il));}

template<typename T, typename... As> decltype(auto) dft_backward(As&... as, std::initializer_list<T> il){return dft_backward(std::forward<As>(as)..., multi::array<T, 1>(il));}
template<typename T, typename... As> decltype(auto) dft_backward(As&... as, std::initializer_list<std::initializer_list<T>> il){return dft_backward(std::forward<As>(as)..., multi::array<T, 2>(il));}

template<class In> In&& dft_inplace(In&& i, sign s){
	fftw::plan{i, i, (int)s}();//(i, i); 
	return std::forward<In>(i);
}

}

namespace fft{
	using fftw::many_dft;
	using fftw::dft;
	using fftw::dft_forward;
	using fftw::dft_backward;

	static constexpr int forward = fftw::forward;//FFTW_FORWARD;
	static constexpr int none = 0;
	static constexpr int backward = fftw::backward;//FFTW_BACKWARD;

	static_assert( forward != none and none != backward and backward != forward, "!");
}

}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_FFTW

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

//#include "../adaptors/fftw/allocator.hpp"
#include<iostream>
#include "../array.hpp"
#include<complex>
#include<numeric>

#include<experimental/array>
#include<experimental/tuple>

#include<random>

#include "../adaptors/cuda.hpp"

#include<chrono>

#include "../../multi/complex.hpp"

namespace{

	using std::cout;
	namespace multi = boost::multi;
	namespace fftw = multi::fftw;

	using complex = std::complex<double>;
	complex const I{0, 1};

	template<class M>
	auto power(M const& m){
		auto sum_norm = [](auto& a, auto& b){return a + std::norm(b);};
		using multi::num_elements; using multi::data_elements; using std::accumulate;
		return accumulate(data_elements(m), data_elements(m) + num_elements(m), double{}, sum_norm);
	}

	constexpr int N = 16;
}

struct watch : private std::chrono::high_resolution_clock{
	std::string label_; time_point  start_;
	watch(std::string label ="") : label_{label}, start_{now()}{}
	~watch(){
		std::cerr<< label_<<": "<< std::chrono::duration<double>(now() - start_).count() <<" sec"<<std::endl;
	}
};

template<class T> struct randomizer{
	template<class M> void operator()(M&& m) const{for(auto&& e:m) operator()(e);}
	void operator()(T& e) const{
		static std::random_device r; static std::mt19937 g{r()}; static std::normal_distribution<T> d;
		e = d(g);
	}
};

template<class T> struct randomizer<std::complex<T>>{
	template<class M> void operator()(M&& m) const{for(auto&& e:m) operator()(e);}
	void operator()(std::complex<T>& e) const{
		static std::random_device r; static std::mt19937 g{r()}; static std::normal_distribution<T> d;
		e = std::complex<T>(d(g), d(g));
	}
};

BOOST_AUTO_TEST_CASE(fftw_2D_identity_2, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> fwd(extensions(in));
//	multi::fftw::dft({false, false}, in, fwd, multi::fftw::forward);
//	multi::fftw::dft<2>(in, fwd, multi::fftw::forward);
//	multi::fftw::dft({multi::fftw::none, multi::fftw::none}, in, fwd);
	multi::fftw::dft(in, fwd, fftw::none);
//	multi::fftw::transpose(in, fwd);
//	BOOST_REQUIRE( fwd == in );
}

BOOST_AUTO_TEST_CASE(fftw_1D){
	multi::array<complex, 1> in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};
	auto fwd = multi::fftw::dft(in, fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_TEST( size(fwd) == size(in) );
//	auto fwd = multi::fftw::dft({multi::fftw::forward}, in); // Fourier[in, FourierParameters -> {1, -1}]
//	auto fwd = multi::fftw::dft<0>(in, multi::fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]

	BOOST_REQUIRE(fwd[2] == -2. - 2.*I);
	BOOST_REQUIRE( in[1] == 2. + 3.*I );

	auto bwd = multi::fftw::dft(in, FFTW_BACKWARD); // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE(bwd[2] == -2. - 2.*I);
}
#if 1
/*
BOOST_AUTO_TEST_CASE(fftw_1D_cuda){
	multi::cuda::managed::array<complex, 1> in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};
	auto fwd = multi::fftw::dft(in, multi::fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
//	auto fwd = multi::fftw::dft(in, multi::fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
//	auto fwd = multi::fftw::dft({multi::fftw::forward}, in); // Fourier[in, FourierParameters -> {1, -1}]
//	auto fwd = multi::fftw::dft<0>(in, multi::fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE(fwd[2] == -2. - 2.*I);
	BOOST_REQUIRE( in[1] == 2. + 3.*I );

	auto bwd = multi::fftw::dft(in, multi::fftw::backward); // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE(bwd[2] == -2. - 2.*I);
}
*/

BOOST_AUTO_TEST_CASE(fftw_2D_identity, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	auto fwd = multi::fftw::dft(in, fftw::none);
//	auto fwd = multi::fftw::dft<0>(in, multi::fftw::none);
//	BOOST_REQUIRE( fwd == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
//	using multi::fftw::forward;
	auto fwd = multi::fftw::dft(in, fftw::forward);
//	auto fwd = multi::fftw::dft<0>(in, forward);
//	auto fwd = multi::fftw::dft({forward, forward}, in);
//	auto fwd = dft({true, true}, in, forward);

	BOOST_TEST( real(fwd[3][1]) == -19.0455 ); // Fourier[in, FourierParameters -> {1, -1}][[4]][[2]]
	BOOST_TEST( imag(fwd[3][1]) == - 2.22717 );

	multi::array<complex, 1> const in0 = { 1. + 2.*I, 9. - 1.*I, 2. + 4.*I};
	using multi::fftw::dft_forward;

	BOOST_REQUIRE( dft_forward(in[0]) == dft_forward(in0) );
//	BOOST_REQUIRE( dft_forward(in[3]) == dft_forward({3.-1.*I, 8.+7.*I, 2.+1.*I}) );
//	BOOST_REQUIRE( dft_forward(rotated(in)[0]) == dft_forward({1.+2.*I, 3.+3.*I, 4. + 1.*I,  3. - 1.*I, 31. - 1.*I}) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_rotated, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
//	using multi::fftw::forward;
	auto fwd = multi::fftw::dft(in, fftw::forward);
//	auto fwd = multi::fftw::dft<0>(in, forward);
//	auto fwd = dft({true, true}, in, forward);
//	auto fwd = multi::fftw::dft({forward, forward}, in);

	using multi::fftw::dft_forward;
//	BOOST_REQUIRE( dft_forward(rotated(in)[0]) == dft_forward({1.+2.*I, 3.+3.*I, 4. + 1.*I,  3. - 1.*I, 31. - 1.*I}) );
//	BOOST_REQUIRE( dft_forward(rotated(in)) == rotated(fwd) );// rotated(fwd) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_many, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> out(extensions(in));
//	multi::fftw::dft<1>(in, out, multi::fftw::forward);
//	dft({false, true}, in, out, multi::fftw::forward);
	multi::fftw::dft({fftw::none, fftw::forward}, in, out);

	using multi::fftw::dft_forward;
	BOOST_REQUIRE( dft_forward(in[0]) == out[0] );

	multi::fftw::dft({false, true}, rotated(in), rotated(out), fftw::forward);
	BOOST_REQUIRE( dft_forward(rotated(in)[0]) == rotated(out)[0] );

	multi::fftw::dft({false, false}, rotated(in), rotated(out), fftw::forward);
	BOOST_REQUIRE( in == out );

	multi::fftw::many_dft(begin(in), end(in), begin(out), fftw::forward);
	using multi::fftw::dft_forward;
	BOOST_REQUIRE( dft_forward(in[0]) == out[0] );
}

BOOST_AUTO_TEST_CASE(fftw_many1_from_2){
	multi::array<complex, 2> in({3, 10}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({3, 10});
	fftw::dft({false, true}, in, out, fftw::forward);

	multi::array<complex, 2> out2({3, 10});
	for(int i = 0; i!=size(in); ++i)
		fftw::dft(in[i], out2[i], fftw::forward);

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_3){
	multi::array<complex, 3> in({3, 5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 3> out({3, 5, 6});
	fftw::dft({false, true, true}, in, out, FFTW_FORWARD);

	multi::array<complex, 3> out2({3, 5, 6});
	for(int i = 0; i!=size(in); ++i)
		fftw::dft(in[i], out2[i], FFTW_FORWARD);

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_2){
	multi::array<complex, 2> in({5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({5, 6});
	fftw::dft({true, true}, in, out, FFTW_FORWARD);

	multi::array<complex, 2> out2({5, 6});
	fftw::dft(in, out2, FFTW_FORWARD);
	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_3D){
	multi::array<complex, 3> in({10, 10, 10});
	in[2][3][4] = 99.;
	auto fwd = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_4D){
	multi::array<complex, 4> const in = []{
		multi::array<complex, 4> in({10, 10, 10, 10}); in[2][3][4][5] = 99.; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, true}, in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4][5] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_4D_many){

	auto const in = []{multi::array<complex, 4> in({97, 95, 101, 10}); in[2][3][4][5] = 99.; return in;}();
	auto fwd = multi::fftw::dft({true, true, true, false}, in, fftw::forward);
	BOOST_REQUIRE( in[2][3][4][5] == 99. );

	multi::array<complex, 4> out(extensions(in));
	multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
	BOOST_REQUIRE( fwd == out );

//	multi::array<complex, 4> out2({10, 97, 95, 101});
//	multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(out2), multi::fftw::forward);
//	BOOST_REQUIRE( fwd == rotated(out2) );

}

BOOST_AUTO_TEST_CASE(cufft_many_2D){
	auto const in = []{
		multi::array<complex, 3> ret({10, 10, 10});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	multi::array<complex, 3> out(extensions(in));
	multi::fftw::many_dft((in<<1).begin(), (in<<1).end(), (out<<1).begin(), multi::fftw::forward);

	multi::array<complex, 3> out2(extensions(in));	
	multi::fftw::dft({true, false, true}, in, out2, multi::fftw::forward);
	BOOST_REQUIRE( out == out2 );
}

BOOST_AUTO_TEST_CASE(fftw_5D){
	multi::array<complex, 5> in({4, 5, 6, 7, 8});
	in[2][3][4][5][6] = 99.;
	auto fwd = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4][5][6] == 99.);

	BOOST_REQUIRE( std::get<2>(sizes(in)) == 6 );
	auto sizes_as_int = std::experimental::apply(
		[](auto... n){
			auto safe = [](auto i){assert(i<=std::numeric_limits<int>::max()); return static_cast<int>(i);};
			return std::array<int, sizeof...(n)>{safe(n)...};
		}, 
		sizes(in)
	);
	BOOST_REQUIRE( sizes_as_int[2] == 6 );
}

BOOST_AUTO_TEST_CASE(fftw_1D_power){
	multi::array<complex, 1> in(N, 0.); assert( size(in) == N );
	std::iota(begin(in), end(in), 1.);
	multi::array<complex, 1> out(extensions(in));
	static_assert(dimensionality(in)==dimensionality(out), "!");
	auto p = multi::fftw_plan_dft(in, out, fftw::forward, FFTW_PRESERVE_INPUT);
	fftw_execute(p); 
	fftw_destroy_plan(p);
	BOOST_REQUIRE( (power(in) - power(out)/num_elements(out)) < 1e-17 );
}

/*
BOOST_AUTO_TEST_CASE(fftw_1D_allocator_power){
	using multi::fftw::allocator;
	multi::array<complex, 1, allocator<complex>> in(16, 0.); std::iota(begin(in), end(in), 1.);
	assert( size(in) == N );
	multi::array<complex, 1, allocator<complex>> out(extensions(in));
	auto p = multi::fftw_plan_dft(in, out, fftw::forward, FFTW_PRESERVE_INPUT);
	fftw_execute(p);
	fftw_destroy_plan(p);
	BOOST_REQUIRE( (power(in) - power(out)/num_elements(out)) < 1e-12 );
}
*/

BOOST_AUTO_TEST_CASE(fftw_2D_power){
	multi::array<complex, 2> in({N, N});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	auto p = multi::fftw_plan_dft(in, out, fftw::forward, FFTW_PRESERVE_INPUT);
	fftw_execute(p); fftw_destroy_plan(p);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-12 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan){
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const p{in, out, fftw::forward, FFTW_PRESERVE_INPUT};
	p(); //execute(p); //p.execute();
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft){
	multi::array<complex, 2> in({16, 16}); std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft(in, out, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out){
	multi::array<complex, 2> in({16, 16}); std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out_default){
	multi::array<complex, 2> in({16, 16}); std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

/*
BOOST_AUTO_TEST_CASE(fftw_2D_carray_power){
	int const N = 16;
	complex in[N][N];
	using multi::data_elements;	using multi::num_elements;
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);
	complex out[N][N];
	auto p = multi::fftw_plan_dft(in, out, fftw::forward | FFTW_PRESERVE_INPUT);
	fftw_execute(p); fftw_destroy_plan(p);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}
*/

BOOST_AUTO_TEST_CASE(fftw_3D_power){
	multi::array<complex, 3> in({4, 4, 4}); std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);
	multi::array<complex, 3> out = fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( std::abs(power(in) - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place){
	multi::array<complex, 3> io({4, 4, 4}); std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);
	auto powerin = power(io);
	fftw::dft_inplace(io, fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place_over_ref_inplace){
	multi::array<complex, 3> io({4, 4, 4}); std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);
	auto powerin = power(io);
//	fftw::dft_inplace(multi::array_ref<complex, 3>(io.data(), io.extensions()), fftw::forward);
	fftw::dft_inplace(multi::array_ref<complex, 3>(data_elements(io), extensions(io)), fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_ref){
	multi::array<complex, 3> in({4, 4, 4}); std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);
	multi::array<complex, 3> out({4, 4, 4});
	multi::array_ref<complex, 3>(data_elements(out), extensions(out)) = fftw::dft(multi::array_cref<complex, 3>(data_elements(in), extensions(in)), fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_temporary){
	double powerin;
	auto f = [&](){
		multi::array<complex, 3> in({4, 4, 4}); 
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);
		powerin = power(in);
		return in;
	};
	auto out = fftw::dft(f(), fftw::forward);
	BOOST_REQUIRE( std::abs(powerin - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_transposition_square_inplace){
	multi::array<complex, 2> in = {
		{11., 12.},
		{21., 22.}
	};
	BOOST_REQUIRE( in[1][0] == 21. );

	multi::fftw::transpose(in, rotated(in));
//	BOOST_REQUIRE( in[1][0] == 12. );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_benchmark){
	multi::array<complex, 3> in({100, 100, 10});
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);
	complex sum = 1.2;
	{//boost::timer::auto_cpu_timer t;
		for(int i = 0; i != 100; ++i){
			multi::array<complex, 3> out = fftw::dft({true, true, true}, in, fftw::forward);
			sum += out[0][0][0];
		}
	}
	std::cout << sum << std::endl;
}

namespace utf = boost::unit_test::framework;

BOOST_AUTO_TEST_CASE(fft_combinations, *boost::unit_test::tolerance(0.00001)){

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	std::vector<std::array<bool, 4>> cases = {
		{false, true , true , true }, 
		{false, true , true , false}, 
		{true , false, false, false}, 
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	using std::cout;
	for(auto c : cases){
		cout<<"case "; copy(begin(c), end(c), std::ostream_iterator<bool>{cout,", "}); cout<<"\n";
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"cpu_oplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in, out);
		}
		{
			multi::fftw::plan p(c, in, out, fftw::forward);
			boost::timer::auto_cpu_timer t{"cpu_oplac planned %ws wall, CPU (%p%)\n"};
			p();
		}
		{
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_iplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in_rw);
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			multi::fftw::plan p(c, in_rw, in_rw, fftw::forward);
			boost::timer::auto_cpu_timer t{"cpu_iplac planned %ws wall, CPU (%p%)\n"};
			p();
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			multi::fftw::plan p(c, in_rw, in_rw, fftw::forward);// | FFTW_MEASURE);
			boost::timer::auto_cpu_timer t{"cpu_iplac planned measured %ws wall, CPU (%p%)\n"};
			p();
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			boost::timer::auto_cpu_timer t{"cpu_alloc %ws wall, CPU (%p%)\n"}; 
			auto out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST( abs( out_cpy[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_move %ws wall, CPU (%p%)\n"}; 
			auto out_cpy = multi::fftw::dft_forward(c, std::move(in_rw));
			BOOST_REQUIRE( in_rw.empty() );
			BOOST_TEST( abs( out_cpy[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark, *boost::unit_test::disabled() ){
	auto x = multi::array<complex, 4>::extensions_type({64, 128, 128, 128});
	multi::array<complex, 4> in(x);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_REQUIRE( in[0][0][0][0] == 1.2 );
	std::array<bool, 4> c = {false, true, true, true};
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}]{
		fftw::dft(c, in, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}]{
		fftw::dft(c, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_REQUIRE( in0000 != 1.2 );


	multi::array<complex, 4> out(x);
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}]{
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}]{
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	BOOST_REQUIRE( in0000 == in[0][0][0][0] );
}

#endif
#endif
#endif

