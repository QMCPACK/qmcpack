#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
//  © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi member cast"
#ifdef BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#else
#include<boost/test/included/unit_test.hpp>
#endif

#include "../array.hpp"

#include "../complex.hpp"

namespace multi = boost::multi;

namespace boost{
namespace multi{

namespace{
	struct priority_0{}; struct priority_1 : priority_0{};
	template<class Array, typename E = typename std::decay_t<Array>::element, typename R = decltype(std::real(E{}))>
	decltype(auto) real_(Array&& a, priority_0){return a;}
	template<class Array, typename E = typename std::decay_t<Array>::element, typename R = decltype(std::real(E{})), typename I = decltype(E{}.imag())>
	decltype(auto) real_(Array&& a, priority_1){
		struct C{R real; I imag;}; static_assert(sizeof(E) == sizeof(C), "!");
		return member_array_cast<R>(reinterpret_array_cast<C>(a), &C::real);
	}
}

template<class Array> decltype(auto) Real(Array&& a){return real_(a, priority_1{});}

template<class Array, typename E = typename std::decay_t<Array>::element, typename R = decltype(std::real(E{})), typename I = decltype(E{}.imag())>
decltype(auto) Imag(Array&& a){
	struct C{R real; I imag;}; static_assert(sizeof(E) == sizeof(C), "!");
	return member_array_cast<I>(reinterpret_array_cast<C>(a), &C::imag);
}

}}

BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos){
{

	using v3d = std::array<double, 3>;

	struct particles_SoA{
		multi::array<double,2> masses; 
		multi::array<v3d,2> positions;
		struct reference{double& mass; v3d& position ;};
		reference operator()(int i, int j){return {masses[i][j], positions[i][j]};}
	};

	struct particle{
		double mass;
		v3d position alignas(2*sizeof(double));  // __attribute__((aligned(2*sizeof(double))))
		particle() = default;
		particle(double m, v3d v) : mass{m}, position{v}{}
		particle(particle const&) = default;
		particle(particles_SoA::reference&& r) : mass{r.mass}, position{r.position}{}
	};

	multi::array<particle, 2> AoS({2, 2}); 
	AoS[1][1] = particle{99., {1.,2.}};

	auto&& masses = multi::member_array_cast<double>(AoS, &particle::mass);
	BOOST_REQUIRE( size(masses) == 2 );
	BOOST_REQUIRE( masses[1][1] == 99. );

	multi::array<double, 2> masses_copy = masses;
	particles_SoA SoA = {
		multi::member_array_cast<double>(AoS, &particle::mass), 
		multi::member_array_cast<v3d>(AoS, &particle::position)
	};
	BOOST_REQUIRE(SoA(1, 1).mass == 99. );

	particle p11 = SoA(1, 1); 
	BOOST_REQUIRE(p11.mass == 99. );
}
{
	struct alignas(sizeof(std::string)) employee{
		std::string name;
		short salary;
		std::size_t age;
		employee(std::string name, short salary, std::size_t age) : name{name}, salary{salary}, age{age}{}
	};

	multi::array<employee, 2> d2D = {
		{ {"Al"  , 1430, 35}, {"Bob"  , 3212, 34} }, 
		{ {"Carl", 1589, 32}, {"David", 2300, 38} }
	};

	using multi::member_array_cast;
	auto&& d2D_names = member_array_cast<std::string>(d2D, &employee::name);

	BOOST_REQUIRE( size(d2D_names) == size(d2D) ); 
	BOOST_REQUIRE( d2D_names[1][1] == "David" );

	multi::static_array<std::string, 2> d2D_names_copy{d2D_names};
	BOOST_REQUIRE( d2D_names == d2D_names_copy );
	BOOST_REQUIRE( base(d2D_names) != base(d2D_names_copy) );

}
}

BOOST_AUTO_TEST_CASE(member_array_cast_complex){

	using complex = std::complex<double>;
	multi::array<complex, 2> A = {
		{ {1.,2.}, {3.,4.} },
		{ {22.,33.}, {5.,9.} }
	};
	struct Complex{
		double real;
		double imag;
	};
	{
		struct complex{double real; double imag;};
		auto&& Areal = multi::member_array_cast<double>(A, &complex::real);
		auto&& Aimag = multi::member_array_cast<double>(A, &complex::imag);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
		auto&& Areal = multi::member_array_cast<double>(A, &multi::complex<double>::real);
		auto&& Aimag = multi::member_array_cast<double>(A, &multi::complex<double>::imag);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
		auto&& Acast = multi::reinterpret_array_cast<Complex>(A);
		auto&& Areal = multi::member_array_cast<double>(Acast, &Complex::real);
		auto&& Aimag = multi::member_array_cast<double>(Acast, &Complex::imag);
		assert( Areal[1][0] == 22. and std::get<1>(strides(Areal)) == 2 );
		assert( Aimag[1][0] == 33. and std::get<1>(strides(Aimag)) == 2 );
	}
	{
		auto&& Areal = multi::member_array_cast<double>(multi::reinterpret_array_cast<Complex>(A), &Complex::real);
		auto&& Aimag = multi::member_array_cast<double>(multi::reinterpret_array_cast<Complex>(A), &Complex::imag);
		assert( Areal[1][0] == 22. and std::get<1>(strides(Areal)) == 2 );
		assert( Aimag[1][0] == 33. and std::get<1>(strides(Aimag)) == 2 );
		Areal[1][0] = 55.;
	}
	{
		auto&& Areal = Real(A);
		auto&& Aimag = Imag(A);
		auto Areal_copy = decay(Real(A));
		assert( Areal[1][0] == 55. and std::get<1>(strides(Areal)) == 2 );
		assert( Aimag[1][0] == 33. and std::get<1>(strides(Aimag)) == 2 );
		Areal[1][0] = 888.;
	}
	{
		multi::array<double, 2> A = {
			{  1., 3.},
			{ 22., 5.}
		};
		auto&& Areal = Real(A);
		assert( Areal[1][1] == 5. );
		Areal[1][1] = 55.;
	}
	{
		multi::array<double, 2> const A = {
			{  1., 3.},
			{ 22., 5.}
		};
		auto&& Areal = Real(A);
		assert( Areal[1][1] == 5. );
	}
}

