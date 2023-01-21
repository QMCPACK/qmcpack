// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>
#include <mpi3/process.hpp>

#include<boost/serialization/vector.hpp>

// nontrivial nonpod class
struct B {  // NOLINT(readability-identifier-naming)
	std::string name_ = "unnamed";  // NOLINT(misc-non-private-member-variables-in-classes)
	std::size_t n_    = 0;  // NOLINT(misc-non-private-member-variables-in-classes)
	double*     data  = nullptr;  // NOLINT(misc-non-private-member-variables-in-classes)

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
void save(Archive& ar, B const& b, unsigned int const /*version*/) {
	ar << b.name_ << b.n_ << boost::serialization::make_array(b.data, b.n_);
}
template<class Archive>
void load(Archive& ar, B& b, unsigned int const /*version*/) {
	ar >> b.name_ >> b.n_;
	delete[] b.data;  // NOLINT(cppcoreguidelines-owning-memory)
	b.data = new double[b.n_];  // NOLINT(cppcoreguidelines-owning-memory)
	ar >> boost::serialization::make_array(b.data, b.n_);
}
BOOST_SERIALIZATION_SPLIT_FREE(B)  // cppcheck-suppress unknownMacro

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {

	assert(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		int a = 5;
		world[1] << a;
		break;
	}
	case 1: {
		int a = -1;
		world[0] >> a;  // specific source (any tag)
		assert(a == 5);
		break;
	}
	}

	switch (world.rank()) {
	case 0: {
		int a = 7;
		world[1] << a;
		break;
	}
	case 1: {
		int a = -1;
		world >> a;  // any source (any tag)
		assert(a == 7);
		break;
	}
	}

	int b = world.rank();
	world[1] & b;  // broadcast (from rank 1)
	assert(b == 1);

	// if(world.root()) {
	// 	B b1(4);
	// 	b1.data[2] = 4.5;
	// 	world[1] << b1;
	// } else {
	// 	B b2;
	// 	world[0] >> b2;
	// 	assert(b2.data[2] == 4.5);
	// }

	{
		switch(world.rank()) {
		case 0: {
			int a = 7;
			world[1] << a;
			break;
		}
		case 1: {
			int a = -1;
			world >> a;  // any source (any tag)
			assert(a == 7);
			break;
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
			assert(bo == true);
			break;
		}
		}

		{
			std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

			switch(world.rank()) {
			case 0: {
				world[1] << v;
				break;
			}
			case 1: {
				std::vector<double> w;
				world >> w;
				assert(v == w);
			}
			}
		}

		{
			std::vector<bool> v = {false, true, false};
			switch(world.rank()) {
			case 0: {
				world[1] << v;
				break;
			}
			case 1: {
				std::vector<bool> v;
				world >> v;
				assert((v == std::vector<bool>{false, true, false}));
			}
			}
		}
	}
	return 0;
} catch(...) { return 1; }
