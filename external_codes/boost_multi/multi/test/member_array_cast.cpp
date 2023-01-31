// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi member cast"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

using v3d = std::array<double, 3>;

BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos) {
// some members might need explicit padding to work well with member_cast
struct particle{
	double mass;
	v3d position alignas(2*sizeof(double));  //  __attribute__((aligned(2*sizeof(double))))
};

class particles_soa {
	multi::array<double, 2> masses_;
	multi::array<v3d, 2> positions_;

 public:
	// NOLINTNEXTLINE(runtime/explicit)
	particles_soa(multi::array<particle, 2> const& AoS)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : particle_soa can represent a particles' AoS
	: masses_   {AoS.member_cast<double>(&particle::mass    )}
	, positions_{AoS.member_cast<v3d   >(&particle::position)} {}

	struct reference {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		double& mass;   // NOLINT(misc-non-private-member-variables-in-classes): exposed by design
		v3d& position;  // NOLINT(misc-non-private-member-variables-in-classes): exposed by design
		operator particle() const {return {mass, position};} // NOLINT(google-explicit-constructor, hicpp-explicit-conversions): allow equal assignment
		auto operator+() const {return operator particle();}

		reference(double& mss, v3d& pos) : mass{mss}, position{pos} {}  // NOLINT(google-runtime-references)

	 private: // NOLINT(whitespace/indent) : bug in cpplint 1.5.5
		friend class particles_soa;
		reference(reference const&) = default;
	//  reference(reference&&) = default;

	 public:  // NOLINT(whitespace/indent) : bug in cpplint 1.5.5
	//  ~reference() noexcept = default;  // lints cppcoreguidelines-special-member-functions,hicpp-special-member-functions
	//  #endif

		// NOLINTNEXTLINE(cert-oop54-cpp, fuchsia-trailing-return): simulate reference
		auto operator=(reference const& other) -> reference& {
			std::tie(mass, position) = std::tie(other.mass, other.position);
			return *this;
		}
		// NOLINTNEXTLINE(fuchsia-trailing-return): simulate reference
		auto operator=(reference&& other) noexcept -> reference& {operator=(other); return *this;}

		auto operator==(reference const& other) const {return std::tie(mass, position) == std::tie(other.mass, other.position);}
		auto operator!=(reference const& other) const {return std::tie(mass, position) != std::tie(other.mass, other.position);}
	};

	auto operator()(int eye, int jay){return reference{masses_[eye][jay], positions_[eye][jay]};}
};

	multi::array<particle, 2> AoS({2, 2}, particle{});
	AoS[1][1] = particle{99., v3d{{1., 2.}} };

	auto&& masses = AoS.member_cast<double>(&particle::mass);
	BOOST_REQUIRE( size(masses) == 2 );
	BOOST_REQUIRE( masses[1][1] == 99. );

	multi::array<double, 2> masses_copy = masses;
	BOOST_REQUIRE( &masses_copy[1][1] != &masses[1][1] );

	particles_soa SoA{AoS};

	BOOST_REQUIRE(SoA(1, 1).mass == 99. );

	particle p11 = SoA(1, 1);
	BOOST_REQUIRE(p11.mass == 99. );

	auto autop11 = +SoA(1, 1);
	BOOST_REQUIRE(autop11.mass == 99. );

	SoA(1, 1).mass = 88;
	BOOST_REQUIRE(SoA(1, 1).mass == 88. );

	SoA(1, 1) = SoA(0, 0);
	BOOST_REQUIRE(SoA(1, 1).mass == SoA(0, 0).mass );
	BOOST_REQUIRE(SoA(1, 1) == SoA(0, 0) );
	BOOST_REQUIRE(not (SoA(1, 1) != SoA(0, 0)) );
}

struct alignas(32) employee {
	std::string name;
	int16_t salary;
	std::size_t age;
//	private: //	char padding_[9];// std::array<char, 9> padding_; // use alignment or padding to allow member_cast
};

BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos_employee) {
	multi::array<employee, 1> d1D = { {"Al"  , 1430, 35}, {"Bob"  , 3212, 34} };
	auto&& d1D_names = d1D.member_cast<std::string>(&employee::name);
	BOOST_REQUIRE( size(d1D_names) == size(d1D) );
	BOOST_REQUIRE(  d1D_names[1] ==  d1D[1].name );
	BOOST_REQUIRE( &d1D_names[1] == &d1D[1].name );

	multi::array<employee, 2> d2D = {
		{ {"Al"  , 1430, 35}, {"Bob"  , 3212, 34} },
		{ {"Carl", 1589, 32}, {"David", 2300, 38} }
	};
	BOOST_REQUIRE( d2D[0][0].name == "Al" );
	BOOST_REQUIRE( d2D[0][0].salary == 1430 );
	BOOST_REQUIRE( d2D[0][0].age == 35 );

	auto&& d2D_names = d2D.member_cast<std::string>(&employee::name);
	BOOST_REQUIRE( size(d2D_names) == size(d2D) );
	BOOST_REQUIRE( d2D_names[1][1] == "David" );

	multi::array<std::string, 2> d2D_names_copy{d2D_names};
	BOOST_REQUIRE( d2D_names == d2D_names_copy );
	BOOST_REQUIRE( base(d2D_names) != base(d2D_names_copy) );
}

