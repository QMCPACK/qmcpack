#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
${CXX:-c++} -std=c++17 $CXXFLAGS -I../include $0 -o $0.$X&&$0.$X&&rm $0.$X;exit
#endif
// Copyright 2018-2022 Alfredo A. Correa

#include "./multi/array.hpp"

#include <thrust/iterator/transform_iterator.h>
#include <thrust/functional.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/functional/hash.hpp>
#include <boost/utility/result_of.hpp>

#include<complex>

namespace multi = boost::multi;

constexpr auto conj = [](auto const& c) -> auto const {return std::conj(c);};

template<class T> struct conjr : boost::transform_iterator<decltype(conj), T*> {
	template<class... As> conjr(As const&... as) : boost::transform_iterator<decltype(conj), T*>{as...} {}
};

template<class Array2D, class Complex = typename Array2D::element_type>
auto hermitized(Array2D const& arr) {
	return arr
		.transposed() // lazily tranposes the array
		.template static_array_cast<Complex, conjr<Complex>>(conj)  // lazy conjugate elements
	;
}

int main(){
{
    using namespace std::complex_literals;
	multi::array A = {
		{ 1. + 2.i,  3. +  4.i},
		{ 8. + 9.i, 10. + 11.i}
	};

	auto const& Ah = hermitized(A);

	assert( Ah[1][0] == std::conj(A[0][1]) );
}

{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x) {return x + 1;};
	std::vector<double> v(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7. );
}
{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x) {return x + 1;};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7. );
}
{
	multi::array<double, 1> v(10);
	auto r = extension(v);
	auto f = [](auto x) {return x * 2;};
	v.assign(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 2. );
}
{
	multi::array<double, 1> v(10);
	multi::array<double, 1> r = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
	auto f = [](auto x) {return x * 2;};
	v.assign(
		boost::make_transform_iterator(r.base()           , f),
		boost::make_transform_iterator(r.base() + r.size(), f)
	);
	assert( v[1] == 4. );
}
{
	auto r = multi::make_extension_t(10l);
	auto f = [](auto x){
		std::size_t seed = 1234;
	//	boost::hash_combine(seed, );
		seed ^= boost::hash<multi::index>{}(x) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		return static_cast<double>(seed)/static_cast<double>(std::numeric_limits<std::size_t>::max());
	};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);

	std::size_t seed = 12349l;
	//	boost::hash_combine(seed, );
//	seed ^= boost::hash<std::size_t>{}(13) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	boost::hash_combine(seed, 13);

	assert( v.size() == r.size() );
	assert( v[1] >= 0. );
	assert( v[1] < 1.  );
	assert( std::all_of(begin(v), end(v), [](auto x){
		return x >= 0. and x < 1.;
	}) );
}

//struct conj_t : thrust::unary_function<std::complex<double>, std::complex<double>> {
//	std::complex<double> operator()(std::complex<double> const& e) const {return std::conj(e);}
//} conj;

{
    using namespace std::complex_literals;
	multi::array<std::complex<double>, 1> A = { 1. + 2.i,  3. +  4.i,  5. +  7.i};

	auto const conj = [](auto e) {return std::conj(e);};

	std::vector<std::complex<double>> v(thrust::make_transform_iterator(A.elements().begin(), conj), thrust::make_transform_iterator(A.elements().end(), conj));
	std::cout << v[1] << std::endl;
	assert( v[1] == 3. - 4.i );

//	using ittc = std::iterator_traits<conjr<std::complex<double>>>::iterator_category;

//	using rt = boost::result_of<const std::complex<double> &(std::complex<double> &)>::type;
//	using conjugater = decltype(boost::make_transform_iterator(A.data_elements(), conj));
//	conjugater ll(A.data_elements(), conj);
//	auto conjA = A.template static_array_cast<std::complex<double>, thrust::transform_iterator<conj_t, std::complex<double>*> >();
//	A.static_array_cast<std::complex<double>, transformer<std::complex<double>, decltype(conj)> >(conj);


//		return {this->layout(), P2{this->base(), std::forward<Args>(args)...}};
//	boost::multi::basic_array<std::complex<double>, 1, std::decay_t<decltype(thrust::make_transform_iterator(A.base(), conj))>> bb{A.layout(), thrust::make_transform_iterator(A.base(), conj)};
}

}
