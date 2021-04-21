#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS $0 -o $0x$OXX `pkg-config --cflags --libs fftw3 cuda-11.0` -lboost_timer -lboost_unit_test_framework&&$0x$OXX&&rm $0x$OXX;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MULTI_ADAPTORS_FFTW_HPP
#define MULTI_ADAPTORS_FFTW_HPP
	
#include "../adaptors/../array.hpp"
#include "../adaptors/../config/NODISCARD.hpp"

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

template<typename It1, class It2, std::enable_if_t<std::is_pointer<decltype(base(It2{}))>{} or std::is_convertible<decltype(base(It2{})), std::complex<double>*>{}, int> = 0
>
auto fftw_plan_many_dft(It1 first, It1 last, It2 d_first, int sign, unsigned flags = FFTW_ESTIMATE)
->decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double>*>*/(base(d_first))), fftw_plan{}){

	static_assert( sizeof(*base(  first)) == sizeof(real(*base(  first))) + sizeof(imag(*base(  first))) and sizeof(*base(  first)) == sizeof(fftw_complex), 
		"input  must have complex pod layout" );
	static_assert( sizeof(*base(d_first)) == sizeof(real(*base(d_first))) + sizeof(imag(*base(d_first))) and sizeof(*base(d_first)) == sizeof(fftw_complex), 
		"output must have complex pod layout");

	assert(sizes(*first)==sizes(*d_first));
	auto ion      = to_array<int>(sizes(*first));

	assert(strides(*first) == strides(*last));
	auto istrides = to_array<int>(strides(*first));
	auto ostrides = to_array<int>(strides(*d_first));

	std::array<std::array<int, 3>, std::decay_t<decltype(*It1{})>::rank::value> ssn;
	for(std::size_t i = 0; i != ssn.size(); ++i) ssn[i] = {istrides[i], ostrides[i], ion[i]};
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

	auto ret = ::fftw_plan_many_dft(
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
	assert(ret);
	return ret;
}

template<
	class In, class Out, dimensionality_type D = std::decay_t<In>::dimensionality,
	class=std::enable_if_t<D==std::decay_t<Out>::dimensionality>,
	class=decltype(reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(base(std::declval<Out&>()))))
>
fftw_plan fftw_plan_dft(std::array<bool, +D> which, In&& in, Out&& out, int sign, unsigned flags = FFTW_ESTIMATE){
	static_assert( sizeof(*base(in )) == sizeof((*base(in )).real()) + sizeof((*base(in)).imag()) and sizeof(*base(in)) == sizeof(fftw_complex), 
		"input must have complex pod layout" );
	static_assert( sizeof(*base(out)) == sizeof((*base(out)).real()) + sizeof((*base(in)).imag()) and sizeof(*base(out)) == sizeof(fftw_complex), 
		"output must have complex pod layout" );

	using multi::sizes;
	assert(sizes(in) == sizes(out));

	using multi::strides;
	auto ion      = to_array<ptrdiff_t>(in.sizes());
	auto istrides = to_array<ptrdiff_t>(in.strides());
	auto ostrides = to_array<ptrdiff_t>(out.strides());

	std::array<fftw_iodim64, D> dims   ; 
	auto l_dims = dims.begin();

	std::array<fftw_iodim64, D> howmany; 
	auto l_howmany = howmany.begin();

	for(int i=0; i!=D; ++i) *(which[i]?l_dims:l_howmany)++ = {ion[i], istrides[i], ostrides[i]};

	assert( D == l_dims - dims.begin() + l_howmany - howmany.begin() );
	assert(in.base()); assert(out.base()); assert( in.extensions() == out.extensions() ); 
	assert( (sign == -1) or (sign == +1) );
	fftw_plan ret = fftw_plan_guru64_dft(
		/*int rank*/ l_dims - dims.begin(),
		/*const fftw_iodim64 *dims*/ dims.data(),
		/*int howmany_rank*/ l_howmany - howmany.begin(),
		/*const fftw_iodim *howmany_dims*/ howmany.data(), //nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(/*static_cast<std::complex<double> const *>*/(in.base()))), 
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(/*static_cast<std::complex<double> *>*/(out.base())),
		sign, flags// | FFTW_ESTIMATE
	);
	assert(ret &&"fftw lib returned a null plan, if you are using MKL check the limitations of their fftw interface"); 
	//https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/appendix-d-fftw-interface-to-intel-math-kernel-library/fftw3-interface-to-intel-math-kernel-library/using-fftw3-wrappers.html
	return ret;
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
	auto ret = fftw_plan_guru64_dft(
		/*int rank*/ s?D:0,
		/*const fftw_iodim64 *dims*/ dims.data(),
		/*int howmany_rank*/ 0,
		/*const fftw_iodim *howmany_dims*/ nullptr, //howmany_dims.data(), //;//nullptr,
		/*fftw_complex *in*/ const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(static_cast<std::complex<double> const*>(base(in)))), 
		/*fftw_complex *out*/ reinterpret_cast<fftw_complex*>(implicit_cast<std::complex<double>*>(base(out))),
		s, flags
	);
	assert(ret);
	return ret;
}

namespace fftw{

#if HAVE_FFTW3_THREADS
void initialize_threads(){int good = fftw_init_threads(); assert(good); (void)good;}
#else
void initialize_threads(){}
#endif

void cleanup(){fftw_cleanup();}

struct environment{
	~environment(){cleanup();}
};

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
		plan r; r.impl_.reset(fftw_plan_many_dft(std::forward<As>(as)...)); return r; // this produces a compilation error in icc++17
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
	static constexpr bool is_thread_safe(){return false;}
	static constexpr bool nthreads(){return 1;}
	static constexpr int with_nthreads(){return 1;}
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

enum strategy: decltype(FFTW_ESTIMATE){ estimate = FFTW_ESTIMATE, measure = FFTW_MEASURE };

template<class In, class Out>
auto dft(In const& i, Out&& o, int s)
->decltype(fftw::plan{i, o, s}(), std::forward<Out>(o)){
	return fftw::plan{i, o, s}(), std::forward<Out>(o);}

using std::decay_t;

template<class In, class Out, std::size_t D=In::dimensionality>
auto dft(std::array<bool, +D> which, In const& i, Out&& o, sign s)
->decltype(plan{which, i, o, s}(), std::forward<Out>(o)){
	return plan{which, i, o, s}(), std::forward<Out>(o);}

template<typename In, class Out, dimensionality_type D=In::dimensionality, dimensionality_type=std::decay_t<Out>::dimensionality>
auto dft(std::array<sign, +D> w, In const& i, Out&& o){
	std::array<bool, D> fwd, /*non,*/ bwd;

	std::transform(begin(w), end(w), begin(fwd), [](auto e){return e==FFTW_FORWARD;});
	dft(fwd, i, o, fftw::forward);

	std::transform(begin(w), end(w), begin(bwd), [](auto e){return e==FFTW_BACKWARD;}); 
	if(std::accumulate(begin(bwd), end(bwd), false)) dft(bwd, o, o, FFTW_BACKWARD);

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
decltype(auto) rotate(multi::array<T, D, Args...>& i, int = 1){
	multi::array_ref<T, D, typename multi::array<T, D, Args...>::element_ptr> before(data_elements(i), extensions(i));
	i.reshape(extensions(rotated(before) ));
	fftw::dft(before, i, fftw::none);
	return i;
}

template<typename In, dimensionality_type D = In::dimensionality, class R=typename In::decay_type>
NODISCARD("when first argument is const")
auto dft(std::array<bool, +D> which, In const& i, sign s)
->std::decay_t<decltype(fftw::dft(which, i, R(extensions(i), get_allocator(i)), s))>{
	return fftw::dft(which, i, R(extensions(i), get_allocator(i)), s);}

template<typename In, multi::dimensionality_type D = std::decay_t<In>::dimensionality>
auto dft(std::array<bool, +D> which, In&& i, sign s)
->decltype(dft(which, i, i, s), std::forward<In>(i)){
	return dft(which, i, i, s), std::forward<In>(i);}

template<typename In, std::size_t D = In::dimensionality, class R=typename In::decay_type>
void dft(std::array<bool, +D> which, In const& i) = delete;

template<dimensionality_type Rank /*not deduced*/, typename In, class R=typename In::decay_type>
NODISCARD("when second argument is const")
R dft(In const& i, sign s){
	static_assert( Rank <= In::dimensionality, "!" );
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

template<class A, multi::dimensionality_type D = A::dimensionality>
NODISCARD("when input argument is read only")
auto dft_forward(std::array<bool, +D> which, A const& a)
->decltype(fftw::dft(which, a, fftw::forward)){
	return fftw::dft(which, a, fftw::forward);}

template<class A, class O, multi::dimensionality_type D = A::dimensionality>
auto dft_forward(std::array<bool, +D> which, A const& a, O&& o)
->decltype(fftw::dft(which, a, std::forward<O>(o), fftw::forward)){
	return fftw::dft(which, a, std::forward<O>(o), fftw::forward);}

template<typename A>
NODISCARD("when input argument is read only")
auto dft_forward(A const& a)
->decltype(fftw::dft(a, fftw::forward)){
	return fftw::dft(a, fftw::forward);}

template<typename... A> auto            dft_backward(A&&... a)
->decltype(dft(std::forward<A>(a)..., fftw::backward)){
	return dft(std::forward<A>(a)..., fftw::backward);}

template<class In> In&& dft_inplace(In&& i, sign s){
	fftw::plan{i, i, (int)s}();//(i, i); 
	return std::forward<In>(i);
}

template<class In, class Out, dimensionality_type D = In::dimensionality>
auto copy(In const& i, Out&& o)
->decltype(dft(std::array<bool, D>{}, i, std::forward<Out>(o), fftw::forward)){
	return dft(std::array<bool, D>{}, i, std::forward<Out>(o), fftw::forward);}

template<typename In, class R=typename In::decay_type>
NODISCARD("when argument is const")
R copy(In const& i)
{//->decltype(copy(i, R(extensions(i), get_allocator(i))), R()){
	return copy(i, R(extensions(i), get_allocator(i)));}
	
template<typename In, class R=typename std::decay_t<In>::decay_type>
auto move(In&& in){
	if(in.is_compact()){
		multi::array_ref<typename In::element, In::dimensionality, typename In::element_ptr> ref(
			in.base(), extensions(in)
		);
		copy(in, ref);
		return R(
			multi::array_ref<typename In::element, In::dimensionality_type, std::move_iterator<typename In::element_ptr>>(std::make_move_iterator(in.mbase()), ((in.mbase()=0), extensions(ref)))
		);
	}else return copy(std::forward<In>(in));
}

template<typename T, dimensionality_type D, class P, class R=typename multi::array<T, D>>
R copy(multi::basic_array<T, D, multi::move_ptr<T, P>>&& a){
	if(a.is_compact()){
		return 
			fftw::copy(
				a.template static_array_cast<T, T*>(), 
				multi::array_ref<T, D, T*>(a.base().base(), a.extensions())
			).template static_array_cast<T, multi::move_ptr<T>>()
		;
	}else return fftw::copy(a.template static_array_cast<T, P>());
}

template<class Array>
auto transpose(Array& a)
->decltype(fftw::copy(transposed(a), a.reshape(extensions(layout(a).transpose())))){
	multi::array_ref<typename Array::element, Array::dimensionality, typename Array::element_ptr> r(a.base(), extensions(a));
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

}}}

////////////////////////////////////////////////////////////////////////////////

#if not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../adaptors/../complex.hpp"

#include<chrono>
#include<random>

#include<thrust/complex.h>

namespace{

	namespace multi = boost::multi;
	namespace fftw = multi::fftw;

	using complex = std::complex<double>; MAYBE_UNUSED complex const I{0, 1};

	template<class M> auto power(M const& m)->decltype(std::norm(m)){return std::norm(m);}

	template<class M, DELETE((M::dimensionality < 1))> double power(M const& m){return accumulate(begin(m), end(m), 0., [](auto const& a, auto const& b){return a + power(b);});}

	struct sum_power{
		template<class A, class B> auto operator()(A const& a, B const& b) const{return a+power(b);}
	};

	MAYBE_UNUSED constexpr int N = 16;
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

struct fftw_fixture : fftw::environment{
	void setup(){} 
	void teardown(){}//fftw_cleanup();}
};

BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

BOOST_AUTO_TEST_CASE(fftw_3D){
	using complex = std::complex<double>; //TODO make it work with thrust
	multi::array<complex, 3> in({10, 10, 10});
	in[2][3][4] = 99.;
	auto fwd = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_1D_const){
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto fwd = multi::fftw::dft(in, fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2. - 2.*I  );
	BOOST_REQUIRE( in[1]  == +2. + 3.*I  );

	auto bwd = multi::fftw::dft(in, fftw::forward); // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE( bwd[2] == -2. - 2.*I  );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity_2, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft({false, false}, in, out, fftw::forward); // out = in;
	BOOST_REQUIRE( power(in) == power(out) );
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	auto fwd = multi::fftw::dft({}, in, fftw::forward);
	BOOST_REQUIRE( fwd == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	
	namespace fftw = multi::fftw;
	auto fwd = fftw::dft_forward(in);
	BOOST_TEST_REQUIRE( fwd[3][1].real() == -19.0455  ); // Fourier[in, FourierParameters -> {1, -1}][[4]][[2]]
	BOOST_TEST_REQUIRE( fwd[3][1].imag() == - 2.22717 );

	multi::array<complex, 1> const in0 = {1. + 2.*I, 9. - 1.*I, 2. + 4.*I};

	auto b = multi::fftw::dft_forward(in0);
	auto a = multi::fftw::dft_forward(in[0]);
	BOOST_REQUIRE( fftw::dft_forward(in[0]) == fftw::dft_forward(in0) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_rotated, *boost::unit_test::tolerance(0.0001)){
	using multi::array;
	array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	using multi::fftw::dft_forward;
	auto fwd = dft_forward(in);
	BOOST_REQUIRE(
		dft_forward(rotated(in)[0])
			== dft_forward(array<complex, 1>{1.+2.*I, 3.+3.*I, 4. + 1.*I,  3. - 1.*I, 31. - 1.*I})
	);
	BOOST_REQUIRE( dft_forward(rotated(in)) == rotated(fwd) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_many, *boost::unit_test::tolerance(0.0001)){
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> out(extensions(in));

	using multi::fftw::dft_forward;

	multi::fftw::dft({fftw::none, fftw::forward}, in, out);
	BOOST_REQUIRE( dft_forward(in[0]) == out[0] );

	multi::fftw::dft({false, true}, rotated(in), rotated(out), fftw::forward);
	BOOST_REQUIRE( dft_forward(rotated(in)[0]) == rotated(out)[0] );

	multi::fftw::dft_forward({false, false}, rotated(in), rotated(out));
	BOOST_REQUIRE( in == out );

	multi::fftw::many_dft(begin(in), end(in), begin(out), fftw::forward);
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

BOOST_AUTO_TEST_CASE(fftw_4D){
	multi::array<complex, 4> const in = []{
		multi::array<complex, 4> in({10, 10, 10, 10}); in[2][3][4][5] = 99.; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, true}, in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4][5] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_4D_many){

	auto const in = []{
		multi::array<complex, 4> in({97, 95, 101, 10}, 0.); 
		in[2][3][4][5] = 99.; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, false}, in, fftw::forward);
	BOOST_REQUIRE( in[2][3][4][5] == 99. );

	multi::array<complex, 4> out(extensions(in));
	multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
	BOOST_REQUIRE( out == fwd );

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

	multi::fftw::copy(in, rotated(in));
	BOOST_TEST( in[0][1].real() == 21. );
	BOOST_TEST( in[0][1].imag() ==  0. );
}

BOOST_AUTO_TEST_CASE(fftw_4D_inq_poisson){

	multi::array<complex, 4> const in = []{
		multi::array<complex, 4> in({50, 100, 137, 1}); 
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);
		return in;
	}();
	
	multi::array<complex, 4> out(extensions(in));
	multi::fftw::dft({0, 1, 1, 0}, in, out);

	BOOST_TEST( power(in) == power(out)/std::get<1>(sizes(out))/std::get<2>(sizes(out)) , boost::test_tools::tolerance(1e-10) );

}


#endif
#endif

