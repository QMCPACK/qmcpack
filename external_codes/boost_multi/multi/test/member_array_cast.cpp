// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, transform_ptr, static_array

// IWYU pragma: no_include <algorithm>                        // for equal  // bug in iwyu 14.0.6? with GNU stdlib
// IWYU pragma: no_include <utility>                          // for addressof  // bug in iwyu 14.0.6? with GNU stdlib
#include <array>     // for array, operator==
#include <cstddef>   // for offsetof, size_t
#include <functional>  // for mem_fn  // IWYU pragma: keep
#include <iterator>  // for size
#include <memory>    // for addressof  // IWYU pragma: keep
#include <string>    // for operator""s, allocator, char_traits
#include <tuple>     // for tie, operator==, tuple

#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Winvalid-offsetof"  // Explicit padding, for particle example
#elif defined(_MSC_VER)
	#pragma warning(push)
	#pragma warning(disable : 4324)  // Explicit padding, for particle example
#endif

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos) {
		using v3d = std::array<double, 3>;

		// some members might need explicit padding to work well with member_cast
		struct particle {
			int mass;
			v3d position alignas(2 * sizeof(double));  // __attribute__((aligned(2*sizeof(double))))
		};

		class particles_soa {
			multi::array<int, 2> masses_;
			multi::array<v3d, 2> positions_;

		 public:  // NOLINT(whitespace/indent) nested class
			// NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : particle_soa can represent a particles' AoS
			explicit particles_soa(multi::array<particle, 2> const& AoS)  // NOLINTNEXTLINE(runtime/explicit)
			: masses_{AoS.member_cast<int>(&particle::mass)}, positions_{AoS.member_cast<v3d>(&particle::position)} {}

			// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)  // NOSONAR
			struct reference {
				int& mass;      // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members) exposed by design
				v3d& position;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members) exposed by design

				// NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions)
				operator particle() const { return {mass, position}; }  // NOSONAR(cpp:S1709) allow direct assignment
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
		BOOST_TEST(size(masses) == 2);
		BOOST_TEST(masses[1][1] == 99 );

		multi::array<int, 2> masses_copy{masses};
		BOOST_TEST(&masses_copy[1][1] != &masses[1][1]);

		particles_soa SoA{AoS};

		BOOST_TEST( SoA(1, 1).mass == 99 );

		particle const p11 = SoA(1, 1);
		BOOST_TEST(p11.mass == 99 );

		auto autop11 = +SoA(1, 1);
		BOOST_TEST(autop11.mass == 99 );

		SoA(1, 1).mass = 88;
		BOOST_TEST( SoA(1, 1).mass == 88 );

		SoA(1, 1) = SoA(0, 0);
		BOOST_TEST( SoA(1, 1).mass == SoA(0, 0).mass);
		BOOST_TEST( SoA(1, 1) == SoA(0, 0));
		BOOST_TEST( !(SoA(1, 1) != SoA(0, 0)));
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

		// clang-format off
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	char padding_[
		((offsetof(employee_dummy, age) + sizeof(age)) / sizeof(std::string) + 1) * sizeof(std::string)
		- (offsetof(employee_dummy, age) + sizeof(age))
	] = {};
		// clang-format on
	};

// TODO(correaa) this doesn't work with NVCC (triggered by adl fill)
#if !(defined(__NVCC__) || defined(__HIPCC__))
	BOOST_AUTO_TEST_CASE(member_array_cast_soa_aos_employee) {
		using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

		// NOLINTBEGIN(misc-include-cleaner) bug in clang-tidy 18
		multi::array<employee, 1> d1D = {
			{ "Al"s, 1430, 35},
			{"Bob"s, 3212, 34},
		};
		// NOLINTEND(misc-include-cleaner) bug in clang-tidy 18

		auto&& d1D_names = d1D.member_cast<std::string>(&employee::name);
		BOOST_TEST(size(d1D_names) == size(d1D));
		BOOST_TEST(d1D_names[1] == d1D[1].name);
		BOOST_TEST(&d1D_names[1] == &d1D[1].name);

		multi::array<employee, 2> d2D = {
			{  {"Al"s, 1430, 35},   {"Bob"s, 3212, 34}},
			{{"Carl"s, 1589, 32}, {"David"s, 2300, 38}},
		};
		BOOST_TEST(d2D[0][0].name == "Al");
		BOOST_TEST(d2D[0][0].salary == 1430);
		BOOST_TEST(d2D[0][0].age == 35);

		auto&& d2D_names = d2D.member_cast<std::string>(&employee::name);
		BOOST_TEST(size(d2D_names) == size(d2D));
		BOOST_TEST(d2D_names[1][1] == "David");

	#if !(defined(__clang__) && defined(__CUDACC__))
		multi::array<std::string, 2> d2D_names_copy_members{d2D.element_transformed(&employee::name)};
		BOOST_TEST(d2D_names_copy_members[1][1] == "David");
		BOOST_TEST(d2D_names_copy_members       == d2D_names);

		multi::array<std::string, 2> d2D_names_copy{d2D_names};
		BOOST_TEST( d2D_names == d2D_names_copy);
		BOOST_TEST( d2D_names.base() != d2D_names_copy.base() );
	#endif
	}
#endif

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
		multi::array<int, 2> ids{recs.element_transformed(&record::id)};

		BOOST_TEST( ids[1][1] == 4 );
		BOOST_TEST( ids == recs.member_cast<int>(&record::id) );

		// recs.element_transformed(std::mem_fn(& A::id) )[1][1] = 5;  // not assignable, ok
		// BOOST_TEST( recs[1][1].id == 5 );
	}

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
		BOOST_TEST( d2D.element_transformed(std::mem_fn(&employee::age)) == d2D.element_transformed(&employee::age) );
	}
#endif

	return boost::report_errors();
}

#if defined(__clang__)
	#pragma clang diagnostic pop
#elif defined(_MSC_VER)
	#pragma warning(pop)
#endif
