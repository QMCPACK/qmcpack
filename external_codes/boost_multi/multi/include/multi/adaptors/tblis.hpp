#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX -std=c++17 -O0 -I/home/correaa/include/tblis -I/home/correaa/tblis/src/external/tci -L/home/correaa/lib -Wl,-rpath=/home/correaa/lib -ltblis $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2021

#include "tblis/tblis.h"

#include "../array.hpp"

#include<tuple>

namespace boost::multi::tblis{

struct index{
	char v_;
	explicit operator char() const{return v_;}
};

namespace indices{
	constexpr index a{'a'}, b{'b'}, c{'c'}, d{'d'}, e{'e'}, f{'f'}, g{'g'}, h{'h'}, i{'i'}, j{'j'}, k{'k'}, l{'l'}, m{'m'}, n{'n'}, o{'o'}, p{'p'}, q{'q'}, r{'r'}, s{'s'}, t{'t'}, u{'u'}, v{'v'}, w{'w'}, x{'x'}, y{'y'}, z{'z'};

namespace abc{
	using indices::a;
	using indices::b;
	using indices::c;
	using indices::d;
	using indices::e;
	using indices::f;
	using indices::h;
}

namespace ijk{
	using indices::i;
	using indices::j;
	using indices::k;
	using indices::l;
	using indices::m;
	using indices::n;
	using indices::o;
	using indices::p;
	using indices::q;
	using indices::r;
	using indices::s;
	using indices::t;
}

namespace uvw{
	using indices::u;
	using indices::v;
	using indices::w;
	using indices::x;
	using indices::y;
	using indices::z;
}

namespace xyz{
	using indices::x;
	using indices::y;
	using indices::z;
}

namespace xyzt{
	using indices::x;
	using indices::y;
	using indices::z;
	using indices::t;
}

namespace greek{
	constexpr index alpha{'1'}, beta{'2'}, gamma{'3'}, delta{'4'}, epsilon{'5'}, zeta{'6'}, mu{'7'}, nu{'8'}, sigma{'9'};
#if defined(__clang__)
	constexpr index α = alpha, β = beta, γ = gamma, δ = delta, ε = epsilon, ζ = zeta, μ = mu, ν = nu, σ = sigma;
#endif
}

#if defined(__clang__)
namespace αβγδ{
	using indices::greek::α;
	using indices::greek::β;
	using indices::greek::γ;
	using indices::greek::δ;
	using indices::greek::ε;
	using indices::greek::ζ;
}

namespace μνσ{
	using indices::greek::μ;
	using indices::greek::ν; 
	using indices::greek::σ;
	using indices::greek::α;
	using indices::greek::β;
	using indices::greek::γ;
}
#endif


}








template<class T> auto init_matrix = std::enable_if_t<sizeof(T*)==0>{};
template<> auto init_matrix<float               > = ::tblis::tblis_init_matrix_s;
template<> auto init_matrix<double              > = ::tblis::tblis_init_matrix_d;
template<> auto init_matrix<std::complex<float >> = ::tblis::tblis_init_matrix_c;
template<> auto init_matrix<std::complex<double>> = ::tblis::tblis_init_matrix_z;

template<class T> auto init_tensor = std::enable_if_t<sizeof(T*)==0>{};
template<> auto init_tensor<float               > = ::tblis::tblis_init_tensor_s;
template<> auto init_tensor<double              > = ::tblis::tblis_init_tensor_d;
template<> auto init_tensor<std::complex<float >> = ::tblis::tblis_init_tensor_c;
template<> auto init_tensor<std::complex<double>> = ::tblis::tblis_init_tensor_z;

template<class Element, multi::dimensionality_type D>
struct indexed_tensor;

template<class Element, multi::dimensionality_type D>
struct tensor : ::tblis::tblis_tensor{

	std::array<::tblis::len_type   , D> lens_;
	std::array<::tblis::stride_type, D> strides_;
	template<class A, std::enable_if_t<not std::is_base_of<tensor, std::decay_t<A>>{}, int> =0>
	explicit tensor(A&& a) : 
		lens_   (std::apply([](auto... s){return std::array<::tblis::len_type   , D>{s...};}, sizes  (a))),
		strides_(std::apply([](auto... s){return std::array<::tblis::stride_type, D>{s...};}, strides(a))){
		tblis::init_tensor<std::decay_t<Element>>(this, D, lens_.data(), const_cast<std::decay_t<Element>*>(base(a)), strides_.data());
	}
	tensor(tensor const&) = delete;
	tensor(tensor&& other) : lens_{other.lens_}, strides_{other.strides_}{
		tblis::init_tensor<std::decay_t<Element>>(this, D, lens_.data(), const_cast<std::decay_t<Element>*>(other.data()), strides_.data());
	}
	using dimensionality_type = multi::dimensionality_type;
	static constexpr dimensionality_type dimensionality(){return D;}
	indexed_tensor<Element, D> operator[](std::string_view indices)&&;
	Element* data() const{return static_cast<Element*>(::tblis::tblis_tensor::data);}
	template<class... Rest>
	auto operator()(Rest...)&&;
};

template<class A, class P = decltype(std::declval<A&&>().base())> tensor(A&&)->tensor<typename std::pointer_traits<P>::element_type, std::decay_t<A>::dimensionality>;

template<class Array, class... Indices>
auto paren(Array&& arr, Indices&&... indices){
	return tblis::tensor(std::forward<Array>(arr))[std::string{char(indices)...}];
}
//->decltype(tensor<typename std::pointer_traits<decltype(arr.base())>::element_type, std::decay_t<Array>::dimensionality>(std::forward<Array>(arr))(indices...)){
//	return tensor<typename std::pointer_traits<decltype(arr.base())>::element_type, std::decay_t<Array>::dimensionality>(std::forward<Array>(arr))(indices...);}

template<class Element, multi::dimensionality_type D>
struct indexed_tensor{
	tensor<Element, D> tensor_;
	std::string indices_;
	indexed_tensor(tensor<Element, D>&& t, std::string indices) : tensor_(std::move(t)), indices_{std::move(indices)}{}
	indexed_tensor(indexed_tensor&& other) = default;
	tensor<Element, D>& tensor_part()&{return tensor_;}
	std::string indices() const{return indices_;}
};

template<class Element, multi::dimensionality_type D>
indexed_tensor<Element, D> tensor<Element, D>::operator[](std::string_view indices)&&{
	return indexed_tensor<Element, D>{std::move(*this), std::string{indices}};
}

template<class Element, multi::dimensionality_type D>
template<class... Rest>
auto tensor<Element, D>::operator()(Rest... rest)&&{
	return indexed_tensor<Element, D>{std::move(*this), std::string{char(rest)...}};
}

template<class TensorA, class TensorB, class TensorC>
auto mult(TensorA&& a, std::string a_indices, TensorB&& b, std::string b_indices, TensorC&& c, std::string c_indices)
->decltype(tblis_tensor_mult(NULL, NULL, &a, a_indices.data(), &b, b_indices.data(), &c, c_indices.data())){
	assert( std::string::size_type(a.dimensionality()) == a_indices.size() );
	assert( std::string::size_type(b.dimensionality()) == b_indices.size() );
	assert( std::string::size_type(c.dimensionality()) == c_indices.size() );
	tblis_tensor_mult(NULL, NULL, &a, a_indices.data(), &b, b_indices.data(), &c, c_indices.data());
}

template<class ITensorA, class ITensorB, class ITensorC>
auto mult(ITensorA&& aijk, ITensorB&& bijk, ITensorC&& cijk)
->decltype(mult(aijk.tensor_part(), aijk.indices(), bijk.tensor_part(), bijk.indices(), cijk.tensor_part(), cijk.indices())){
	return mult(aijk.tensor_part(), aijk.indices(), bijk.tensor_part(), bijk.indices(), cijk.tensor_part(), cijk.indices());}

template<class Element>
struct matrix : ::tblis::tblis_matrix{
public:
	template<class A, std::enable_if_t<not std::is_base_of<matrix<Element>, std::decay_t<A>>{}, int> =0>
	matrix(A&& a){
		init_matrix<Element>(this, 
			std::get<0>(a.sizes()), std::get<1>(a.sizes()), const_cast<double*>(a.base()), 
			std::get<0>(a.strides()), std::get<1>(a.strides())
		);
	}
//	template<class EE> matrix(matrix<EE> const& other) : ::tblis::tblis_matrix
	matrix(matrix const&) = delete;
	matrix(matrix&&) = default;
};

template<class A, class P = typename std::decay_t<A>::element_ptr> matrix(A&&)->matrix<typename std::pointer_traits<P>::element_type>;

template<class AElement, class BElement, class CElement>
void mult(matrix<AElement> const& A, matrix<BElement> const& B, matrix<CElement>&& C){
	::tblis::tblis_matrix_mult(nullptr, nullptr, &A, &B, &C);
}

template<class A, class B, class C>
auto mult(A const& a, B const& b, C&& c)
->decltype(mult(matrix(a), matrix(b), matrix(c))){
	return mult(matrix(a), matrix(b), matrix(c));}

}


