// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4324)  // Explicit padding required
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos) {
	using v3d = std::array<double, 3>;

	// some members might need explicit padding to work well with member_cast
	struct particle {
		int mass;
		v3d    position alignas(2 * sizeof(double));  // __attribute__((aligned(2*sizeof(double))))
	};

	class particles_soa {
		multi::array<int, 2> masses_;
		multi::array<v3d, 2>    positions_;

	 public:  // NOLINT(whitespace/indent) nested class
		// NOLINTNEXTLINE(runtime/explicit)
		explicit particles_soa(multi::array<particle, 2> const& AoS)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : particle_soa can represent a particles' AoS
		: masses_{AoS.member_cast<int>(&particle::mass)}, positions_{AoS.member_cast<v3d>(&particle::position)} {}

		struct reference {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions) // NOSONAR
			int& mass;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members) exposed by design
			v3d&    position;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members) exposed by design

			operator particle() const { return {mass, position}; }  // NOLINT(google-explicit-constructor, hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) allow direct assignment
			auto operator+() const { return operator particle(); }

			reference(int& mss, v3d& pos) : mass{mss}, position{pos} {}  // NOLINT(google-runtime-references)
			// unused: explicit reference(particle& other) : reference{other.mass, other.position} {}

		 private:  // NOLINT(whitespace/indent) nested class
			friend class particles_soa;

		 public:  // NOLINT(whitespace/indent) nested class
			auto operator=(reference const& other) && -> reference& {
				if(this == std::addressof(other)) {
					return *this;
				}
				std::tie(mass, position) = std::tie(other.mass, other.position);
				return *this;
			}

			auto operator==(reference const& other) const { return std::tie(mass, position) == std::tie(other.mass, other.position); }
			auto operator!=(reference const& other) const { return std::tie(mass, position) != std::tie(other.mass, other.position); }
		};

		auto operator()(int eye, int jay) { return reference{masses_[eye][jay], positions_[eye][jay]}; }
	};

	multi::array<particle, 2> AoS({2, 2}, particle{});
	AoS[1][1] = particle{99, v3d{{1.0, 2.0}}};

	auto&& masses = AoS.member_cast<int>(&particle::mass);
	BOOST_REQUIRE(size(masses) == 2);
	BOOST_REQUIRE(masses[1][1] == 99 );

	multi::array<int, 2> masses_copy = masses;
	BOOST_REQUIRE(&masses_copy[1][1] != &masses[1][1]);

	particles_soa SoA{AoS};

	BOOST_REQUIRE( SoA(1, 1).mass == 99 );

	particle const p11 = SoA(1, 1);
	BOOST_REQUIRE(p11.mass == 99 );

	auto autop11 = +SoA(1, 1);
	BOOST_REQUIRE(autop11.mass == 99 );

	SoA(1, 1).mass = 88;
	BOOST_REQUIRE( SoA(1, 1).mass == 88 );

	SoA(1, 1) = SoA(0, 0);
	BOOST_REQUIRE( SoA(1, 1).mass == SoA(0, 0).mass);
	BOOST_REQUIRE( SoA(1, 1) == SoA(0, 0));
	BOOST_REQUIRE( ! (SoA(1, 1) != SoA(0, 0)));
}

struct employee_dummy {
	std::string name;
	// NOLINTNEXTLINE(runtime/int)
	short       salary;  // NOLINT(google-runtime-int)
	std::size_t age;
};

struct employee {
	std::string name;
	// NOLINTNEXTLINE(runtime/int)
	short       salary;  // NOLINT(google-runtime-int)
	std::size_t age;
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	char padding_[(((offsetof(employee_dummy, age) + sizeof(age)) / sizeof(std::string) + 1) * sizeof(std::string) - (offsetof(employee_dummy, age) + sizeof(age)))] = {};
};

// TODO(correaa) this doesn't work with NVCC (triggered by adl fill)
#if !(defined(__NVCC__) || defined(__HIPCC__))
BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos_employee) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	multi::array<employee, 1> d1D = {
		{ "Al"s, 1430, 35},
		{"Bob"s, 3212, 34},
	};

	auto&& d1D_names = d1D.member_cast<std::string>(&employee::name);
	BOOST_REQUIRE(size(d1D_names) == size(d1D));
	BOOST_REQUIRE(d1D_names[1] == d1D[1].name);
	BOOST_REQUIRE(&d1D_names[1] == &d1D[1].name);

	multi::array<employee, 2> d2D = {
		{  {"Al"s, 1430, 35},   {"Bob"s, 3212, 34}},
		{{"Carl"s, 1589, 32}, {"David"s, 2300, 38}},
	};
	BOOST_REQUIRE(d2D[0][0].name == "Al");
	BOOST_REQUIRE(d2D[0][0].salary == 1430);
	BOOST_REQUIRE(d2D[0][0].age == 35);

	auto&& d2D_names = d2D.member_cast<std::string>(&employee::name);
	BOOST_REQUIRE(size(d2D_names) == size(d2D));
	BOOST_REQUIRE(d2D_names[1][1] == "David");

#if !(defined(__clang__) && defined(__CUDACC__))
#if !defined(__circle_build__) || (__circle_build__ > 200 )
	multi::array<std::string, 2> d2D_names_copy_members = d2D.element_transformed(&employee::name);
	BOOST_REQUIRE(d2D_names_copy_members[1][1] == "David");
	BOOST_REQUIRE(d2D_names_copy_members       == d2D_names);
#endif
#endif

#if !(defined(__clang__) && defined(__CUDACC__))
	multi::array<std::string, 2> d2D_names_copy{d2D_names};
	BOOST_REQUIRE(d2D_names == d2D_names_copy);
	BOOST_REQUIRE(base(d2D_names) != base(d2D_names_copy));
#endif
}
#endif

#if !defined(__circle_build__) || (__circle_build__ > 200 )
BOOST_AUTO_TEST_CASE(element_transformed_from_member) {
	struct record {
		int    id;
		double data;
	};

	multi::array<record, 2> const recs = {
		{{1, 1.1}, {2, 2.2}},
		{{3, 3.3}, {4, 4.4}},
	};

	// multi::array<int, 2> ids = recs.element_transformed(std::mem_fn(& A::id));
	multi::array<int, 2> ids = recs.element_transformed(&record::id);

	BOOST_REQUIRE( ids[1][1] == 4 );
	BOOST_REQUIRE( ids == recs.member_cast<int>(&record::id) );

	// recs.element_transformed(std::mem_fn(& A::id) )[1][1] = 5;  // not assignable, ok
	// BOOST_REQUIRE( recs[1][1].id == 5 );
}
#endif

// TODO(correaa) this doesn't work with NVCC (triggered by adl fill)
#if !(defined(__NVCC__) || defined(__HIPCC__))
BOOST_AUTO_TEST_CASE(element_transformed_from_member_no_amp) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	multi::array<employee, 2> d2D = {
		{  {"Al"s, 1430, 35},   {"Bob"s, 3212, 34}},
		{{"Carl"s, 1589, 32}, {"David"s, 2300, 38}},
	};

	// multi::array<std::size_t, 2> d2D_ages_copy =
	d2D.element_transformed(std::mem_fn(&employee::age));
	BOOST_REQUIRE( d2D.element_transformed(std::mem_fn(&employee::age)) == d2D.element_transformed(&employee::age) );
}
#endif
