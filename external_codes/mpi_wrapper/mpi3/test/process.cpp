// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>  // NOLINT(misc-include-cleaner) used indirectly

#include <algorithm>
#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <string>
#include <vector>

// nontrivial nonpod class
struct B {                          // NOLINT(readability-identifier-naming)
	std::string name_ = "unnamed";  // NOLINT(misc-non-private-member-variables-in-classes)
	std::size_t n_    = 0;          // NOLINT(misc-non-private-member-variables-in-classes)
	double*     data  = nullptr;    // NOLINT(misc-non-private-member-variables-in-classes)

	B() = default;
	explicit B(std::size_t n) : n_{n}, data{new double[n]} { std::fill_n(data, n, 0.0); }
	B(B const& other) : name_{other.name_}, n_{other.n_}, data{new double[other.n_]} {}
	B(B&&) = delete;

	B& operator=(B&&) = default;
	B& operator=(B const& other) {
		if(this == &other) {
			return *this;
		}
		name_ = other.name_;
		n_    = other.n_;
		delete[] data;
		data = new double[other.n_];  // NOLINT(cppcoreguidelines-owning-memory)
		std::copy_n(other.data, other.n_, data);
		return *this;
	}

	~B() { delete[] data; }
};

// nonintrusive serialization
template<class Archive>
static void save(Archive& ar, B const& b, unsigned int const /*version*/) {  // NOLINT(misc-use-anonymous-namespace) needed by boost serialization
	ar << b.name_ << b.n_ << boost::serialization::make_array(b.data, b.n_);
}

template<class Archive>
static void load(Archive& ar, B& b, unsigned int const /*version*/) {  // NOLINT(misc-use-anonymous-namespace) needed by boost serialization
	ar >> b.name_ >> b.n_;
	delete[] b.data;            // NOLINT(cppcoreguidelines-owning-memory)
	b.data = new double[b.n_];  // NOLINT(cppcoreguidelines-owning-memory)
	ar >> boost::serialization::make_array(b.data, b.n_);
}
BOOST_SERIALIZATION_SPLIT_FREE(B)  // cppcheck-suppress unknownMacro

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		int const a = 5;
		world[1] << a;
		break;
	}
	case 1: {
		int a = -1;
		world[0] >> a;  // specific source (any tag)
		BOOST_TEST(a == 5);
		break;
	}
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		int const a = 7;
		world[1] << a;
		break;
	}
	case 1: {
		int a = -1;
		world >> a;  // any source (any tag)
		BOOST_TEST(a == 7);
		break;
	}
	default: {}
	}

	int b = world.rank();
	world[1] & b;  // broadcast (from rank 1)
	BOOST_TEST(b == 1);

	// if(world.root()) {
	//  B b1(4);
	//  b1.data[2] = 4.5;
	//  world[1] << b1;
	// } else {
	//  B b2;
	//  world[0] >> b2;
	//  BOOST_TEST(b2.data[2] == 4.5);
	// }

	{
		switch(world.rank()) {
		case 0: {
			int const a = 7;
			world[1] << a;
			break;
		}
		case 1: {
			int a = -1;
			world >> a;  // any source (any tag)
			BOOST_TEST(a == 7);
			break;
		}
		default: {
		}
		}

		switch(world.rank()) {
		case 0: {
			world[1] << true;
			break;
		}
		case 1: {
			bool bo = false;
			world >> bo;
			BOOST_TEST(bo == true);
			break;
		}
		default: {}
		}

		{
			std::vector<double> const v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

			switch(world.rank()) {
			case 0: {
				world[1] << v;
				break;
			}
			case 1: {
				std::vector<double> w;
				world >> w;
				BOOST_TEST(v == w);
			}
			default: {
			}
			}
		}
		{
			std::vector<bool> const v_send = {false, true, false};
			switch(world.rank()) {
			case 0: {
				world[1] << v_send;
				break;
			}
			case 1: {
				std::vector<bool> v_recv;
				world >> v_recv;
				BOOST_TEST(v_recv == v_send);
			}
			default: {
			}
			}
		}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
