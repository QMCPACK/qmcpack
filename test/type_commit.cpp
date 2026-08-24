// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/type.hpp>

#include <mpi3/detail/tuple_offset.hpp>

#include <boost/core/lightweight_test.hpp>

#include <tuple>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	{
		using Tuple = std::tuple<int, double, int, char, float, long double>;
		Tuple tup;
		auto offset4 = mpi3::detail::element_offset<4, Tuple>();
		BOOST_TEST( reinterpret_cast<char*>(&tup) + offset4 == reinterpret_cast<char*>(&std::get<4>(tup)) );  // NOLINT(cert-dcl03-c,hicpp-static-assert,misc-static-assert,cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-bounds-pointer-arithmetic) for some compiler this is not a constexpr
	}

	mpi3::type const t = mpi3::int_[100];  // mpi3::type::int_.contiguous(100)

	return boost::report_errors();
}
