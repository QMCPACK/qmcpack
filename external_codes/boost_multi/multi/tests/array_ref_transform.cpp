#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif

#include "../array.hpp"

#include<complex>
#include<iostream>
#include<boost/iterator/transform_iterator.hpp> //might need -DBOOST_RESULT_OF_USE_DECLTYPE

using std::cout; using std::cerr;
namespace multi = boost::multi;

auto neg = [](auto const& x){return -x;};

//auto inverse(decltype(neg)){return [](auto&& x){return -x;};}
//auto bconj = [](auto const& x){using std::conj; return conj(x);};
//auto inverse(decltype(bconj)){return [](auto const& x){using std::conj; return conj(x);};}

template<class It, class F> class involuter;

#if __has_cpp_attribute(no_unique_address)>=201803
#define CPP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
#define CPP_NO_UNIQUE_ADDRESS
#endif

template<class T, std::enable_if_t<std::is_empty<T>{}, int> =0> 
T default_construct(){return *(T*)(&default_construct<T>);}

template<class T, std::enable_if_t<not std::is_empty<T>{}, int> =0> 
T default_construct(){return {};}

template<class Ref, class Involution>
class involuted{
protected:
	Ref r_;
	CPP_NO_UNIQUE_ADDRESS Involution f_;
public:
	using decay_type =std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	explicit involuted(Ref r) : involuted(
		std::forward<Ref>(r), 
		default_construct<Involution>()
//		[&]()->Involution{if constexpr(std::is_empty<Involution>{}) return*(Involution*)this;else return {};}()
	){
	///	static_assert( std::is_trivially_default_constructible<std::decay_t<Involution>>{}, "!" );
	//	static_assert( std::is_empty<Involution>{}, "!" );
	}
	explicit involuted(Ref r, Involution f) : r_{std::forward<Ref>(r)}, f_{f}{}
	involuted& operator=(involuted const& other)=delete;//{r_ = other.r_; return *this;}
public:
	involuted(involuted const&) = delete;
#if __cplusplus >= 201703L
	involuted(involuted&&) = delete;
#else
	involuted(involuted&&) = default;
#endif
	operator decay_type() const&{return f_(r_);}
	decltype(auto) operator&()&&{return involuter<decltype(&std::declval<Ref>()), Involution>{&r_, f_};}
	template<class DecayType>
	auto operator=(DecayType&& other)&&
	->decltype(r_=f_(std::forward<DecayType>(other)), *this){
		return r_=f_(std::forward<DecayType>(other)), *this;}
	template<class DecayType>
	auto operator=(DecayType&& other)&
	->decltype(r_=f_(std::forward<DecayType>(other)), *this){
		return r_=f_(std::forward<DecayType>(other)), *this;}
	template<class OtherRef>
	auto operator=(involuted<OtherRef, Involution> const& o)&
	->decltype(r_=f_==o.f_?std::forward<decltype(o.r_)>(o.r_):f_(o), *this){
		return r_=f_==o.f_?std::forward<decltype(o.r_)>(o.r_):f_(o), *this;}
	template<class DecayType>
	auto operator==(DecayType&& other) const&
	->decltype(this->operator decay_type()==other){
		return this->operator decay_type()==other;}
	template<class Any> friend auto operator<<(Any&& a, involuted const& self)->decltype(a << std::declval<decay_type>()){return a << self.operator decay_type();}
};

#if __cpp_deduction_guides
template<class T, class F> involuted(T&&, F)->involuted<T const, F>;
//template<class T, class F> involuted(T&, F)->involuted<T&, F>;
//template<class T, class F> involuted(T const&, F)->involuted<T const&, F>;
#endif

template<class It, class F>
class involuter : public std::iterator_traits<It>{
	It it_;
	CPP_NO_UNIQUE_ADDRESS F f_;
public:
	involuter(It it) : involuter(
		std::move(it), 
		default_construct<F>()
//		[&]{if constexpr(sizeof(F)<=1) return *(F*)(this); else return F{};}()
	){}
	involuter(It it, F f) : it_{std::move(it)}, f_{std::move(f)}{}
	using reference = involuted<typename std::iterator_traits<It>::reference, F>;
	auto operator*() const{return reference{*it_, f_};}
	involuter& operator+=(typename involuter::difference_type n){it_+=n; return *this;}
	involuter operator+(typename involuter::difference_type n) const{return {it_+n, f_};}
	decltype(auto) operator->() const{
		return involuter<typename std::iterator_traits<It>::pointer, F>{&*it_, f_};
	}
};

#undef CPP_NO_UNIQUE_ADDRESS

template<class ComplexRef> using negated = involuted<ComplexRef, decltype(neg)>;
template<class ComplexIt> using negater = involuter<ComplexIt, decltype(neg)>;

auto const conj = [](std::complex<double> const& a){return std::conj(std::forward<decltype(a)>(a));};

template<class ComplexRef> struct conjd : involuted<ComplexRef, decltype(conj)>{
	conjd(ComplexRef r) : involuted<ComplexRef, decltype(conj)>(r){}
	decltype(auto) real() const{return this->r_.real();}
	decltype(auto) imag() const{return negated<decltype(this->r_.imag())>(this->r_.imag());}//-this->r_.imag();}//negated<std::decay_t<decltype(this->r_.imag())>>(this->r_.imag());} 
	friend decltype(auto) real(conjd const& self){using std::real; return real(static_cast<typename conjd::decay_type>(self));}
	friend decltype(auto) imag(conjd const& self){using std::imag; return imag(static_cast<typename conjd::decay_type>(self));}
};

#if __cpp_deduction_guides
template<class T> conjd(T&&)->conjd<T>;
#endif

template<class Complex> using conjr = involuter<Complex, decltype(conj)>;

#if 0
template<class It, class F, class InvF = decltype(inverse(std::declval<F>()))>
struct bitransformer{
	It it_;
	F f_;
	InvF invf_;
	bitransformer(It it) : bitransformer(
		std::move(it), 
		[]{if constexpr(sizeof(F)==sizeof(std::true_type)) return *static_cast<F*>(nullptr); else return F{};}()
	){}
	bitransformer(It it, F f) : bitransformer(it, f, inverse(f)){}
	bitransformer(It it, F f, InvF invf) : it_{std::move(it)}, f_{f}, invf_{invf}{}
	using difference_type = typename std::iterator_traits<It>::difference_type;
	using value_type = typename std::iterator_traits<It>::value_type;
	using pointer  = typename std::iterator_traits<It>::pointer;
	class reference{
		typename std::iterator_traits<It>::reference r_;
		F f_;
		InvF invf_;
		reference(typename std::iterator_traits<It>::reference r, F f, InvF invf) : r_{r}, f_{std::move(f)}, invf_{std::move(invf)}{}
		reference(reference const&) = delete;
		reference(reference&&) = default;
	public:
		operator typename std::iterator_traits<It>::value_type() &&{return f_(r_);}
		template<class T, typename = decltype(*(std::declval<It>()) = std::declval<T>())> 
		reference&& operator=(T&& t)&&{r_ = invf_(std::forward<T>(t)); return std::move(*this);}
		friend struct bitransformer;
	};
	using iterator_category = typename std::iterator_traits<It>::iterator_category;
	reference operator*() const{return reference{*it_, f_, invf_};}
	bitransformer operator+(std::ptrdiff_t n) const{return {it_ + n, f_, invf_};}
};
#endif

template<class P = std::complex<double>*>
struct indirect_real : std::iterator_traits<typename std::pointer_traits<P>::element_type::value_type*>{
	P impl_;
	indirect_real(P const& p) : impl_{p}{}
    auto operator+(std::ptrdiff_t n) const{return indirect_real{impl_ + n};}
	double& operator*() const{return reinterpret_cast<double(&)[2]>(*impl_)[0];}
};

struct A{
	A(A const&)=delete;
	A(A&&)=delete;
};

template<class T> void what(T&&);

int main(){

	{
		using complex = std::complex<double>;
		complex c{1., 2.};
		auto&& z = conjd<complex&>{c};
		cout<< z << std::endl;
		auto pz = &z;
		cout << real(*pz) <<' '<< imag(*pz) << std::endl;
		cout << pz->real() << std::endl;
		cout << pz->imag() << std::endl;
	}
	{
		double a = 5;
		double& b = a; assert( b == 5 );
		auto&& c = involuted<double&, decltype(neg)>(a, neg);
		assert( c == -5. );
		c = 10.; assert( c == 10. );
		assert( a = -10. );
		
	//	double&& aa = 5.;
	//	auto neg_fun = [](auto const& x){return -x;};
	//	auto m5 = involuted<double const&, decltype(neg_fun)>(5., neg_fun);
	//	m5 = 4.;
	//	std::cout << m5 << std::endl;
	//	assert( m5 == -5. );
	//	c = m5;
	}

#if 0
	double const d2D[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	double const md2D[4][5] {
		{ -0,  -1,  -2,  -3,  -4}, 
		{ -5,  -6,  -7,  -8,  -9}, 
		{-10, -11, -12, -13, -14}, 
		{-15, -16, -17, -18, -19}
	};
#endif
//#if __cpp_deduction_guides
//	multi::array_ref d2DA({4, 5}, boost::transform_iterator{&d2D[0][0], neg});
//	multi::array_ref d2DB({4, 5}, &md2D[0][0]);
//#else
//	auto d2DA = multi::make_array_ref(
//		boost::make_transform_iterator(&d2D[0][0], neg), 
//		{4, 5}
//	);
//	auto d2DB = multi::make_array_ref(&md2D[0][0], {4, 5});
//#endif
//	d2DA[0][0] = 4.;
//	assert( d2DA == d2DB );

{
	multi::array<double, 1> A = { 0,  1,  2,  3,  4};
	auto&& A_ref = multi::static_array_cast<double, double const*>(A);
	assert( A_ref[2] == A[2] );
//	assert( mA_ref == mA );
//	assert( mA_ref[1][1] == mA[1][1] );
//	assert( mA_ref[1] == mA[1] );
//	assert( mA_ref == mA ); 
}

{
	multi::array<double, 1> A = { 0,  1,  2,  3,  4};
	multi::array<double, 1> mA = { -0,  -1,  -2,  -3, -4};
	auto&& mA_ref = multi::static_array_cast<double, negater<double*>>(A);
	assert( mA_ref[2] == mA[2] );
//	assert( mA_ref == mA );
//	assert( mA_ref[1][1] == mA[1][1] );
//	assert( mA_ref[1] == mA[1] );
//	assert( mA_ref == mA ); 
}
{
	multi::array<double, 2> A = {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	multi::array<double, 2> mA = {
		{ -0,  -1,  -2,  -3, -4}, 
		{ -5,  -6,  -7,  -8, -9}, 
		{-10, -11, -12, -13, -14}, 
		{-15, -16, -17, -18, -19}
	};
	auto&& mA_ref = multi::static_array_cast<double, negater<double*>>(A);
	assert( mA_ref[1][1] == mA[1][1] );
//	assert( mA_ref[1] == mA[1] );
//	assert( mA_ref == mA ); 
}

{
#if __cpp_deduction_guides
	double Z[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	auto d2DC = multi::make_array_ref(involuter<double*, decltype(neg)>{&Z[0][0], neg}, {4, 5});
//	multi::array_ref d2DC{bitransformer<decltype(neg), decltype(&Z[0][0])>{&Z[0][0], neg}, {4, 5}};
	cout<< d2DC[1][1] <<'\n';
	d2DC[1][1] = -66;
	cout<< d2DC[1][1] <<'\n';
	assert( Z[1][1] == 66 );
#endif

	{
		using complex = std::complex<double>;
		multi::array<complex, 2> d2D = {
			{ {0, 3}, {1, 9}, {2, 4},  3,  4}, 
			{  5    , {6, 3}, {7, 5},  8,  9}, 
			{ {1, 4}, {9, 1}, 12    , 13, 14}, 
			{  15   ,  16   , 17    , 18, 19}
		};

		using multi::reinterpret_array_cast;
		auto&& d2Dreal = reinterpret_array_cast<double>(d2D);
		assert( d2Dreal[2][1] == 9. );
		d2Dreal[2][1] = 12.;
		assert( d2D[2][1] == complex(12., 1.) ); 

		auto&& d2DrealT = reinterpret_array_cast<double>(rotated(d2D));
		assert( d2DrealT[2][1] == 7. );

		multi::array<double, 2> d2real_copy = reinterpret_array_cast<double>(d2D);//d2Dreal;

		using multi::static_array_cast;
		auto&& d2Dreal2 = static_array_cast<double, indirect_real<>>(d2D);
		assert( d2Dreal2[2][1] == 12. );

		struct indirect_imag{
			std::complex<double>* underlying; using element_type = double;
			indirect_imag(std::complex<double>* underlying) : underlying{underlying}{}
			indirect_imag operator+(std::ptrdiff_t n) const{return {underlying + n};}
			double& operator*() const{return reinterpret_cast<double(&)[2]>(*underlying)[1];}
			operator double*() const{return &(*(*this));}
			using difference_type = std::ptrdiff_t;
			using value_type = double;
			using pointer = double*;
			using reference = double&;
			using iterator_category = std::random_access_iterator_tag;
		};
		auto&& d2imag2 = static_array_cast<double, indirect_imag>(d2D);
		assert( d2imag2[2][1] == 1. );
		double* p = d2imag2.base();
		assert( *p == 3 );

	}
	{
		using complex = std::complex<double>;
		constexpr auto const I = complex{0., 1.};
		multi::array<complex, 2> A = {
			{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
			{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
		};
	//	[[maybe_unused]] 
		auto Aconj = multi::static_array_cast<complex, conjr<complex*>>(A);
		assert( Aconj[1][2] == conj(A[1][2]) );
	}
}

}

