// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>
//#include <mpi3/detail/package_archive.hpp"

#include <boost/serialization/utility.hpp>  // serialize std::pair
#include <boost/serialization/vector.hpp>

#include <set>

namespace mpi3 = boost::mpi3;

struct A {  // NOLINT(readability-identifier-naming) example name
 private:
	std::string               name_ = "unnamed";
	std::size_t               n_    = 0;
	std::unique_ptr<double[]> data_;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) notation

 public:
	A() = default;
	explicit A(std::size_t n) : n_{n}, data_{new double[n_]} {}
	A(A const& other) : name_{other.name_}, n_{other.n_}, data_{new double[n_]} {}
	A(A&&) = delete;
	~A()   = default;

	// auto operator=(A&&) = delete;
	auto operator=(A const& other) -> A& {
		if(this == &other) {
			return *this;
		}
		name_ = other.name_;
		n_    = other.n_;
		data_.reset(new double[other.n_]);  // NOLINT(cppcoreguidelines-owning-memory)
		std::copy_n(other.data_.get(), n_, data_.get());
		return *this;
	}
	auto operator=(A&&) -> A& = default;

	auto operator[](std::ptrdiff_t i) -> double& { return data_.get()[i]; }
	// intrusive serialization
	template<class Archive>
	void save(Archive& ar, unsigned int const /*version*/) const {
		ar << name_ << n_ << boost::serialization::make_array(data_.get(), n_);
	}
	template<class Archive>
	void load(Archive& ar, unsigned int const /*version*/) {
		ar >> name_ >> n_;
		data_.reset(new double[n_]);  // NOLINT(cppcoreguidelines-owning-memory)
		ar >> boost::serialization::make_array(data_.get(), n_);
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

struct B {  // NOLINT(readability-identifier-naming) example name
	std::string               name_ = "unnamed";  // NOLINT(misc-non-private-member-variables-in-classes) exposed for serialization
	std::size_t               n_    = 0;  // NOLINT(misc-non-private-member-variables-in-classes)
	std::unique_ptr<double[]> data;  // NOLINT(misc-non-private-member-variables-in-classes, cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

	B() = default;
	explicit B(std::size_t n) : n_{n}, data{std::make_unique<double[]>(n_)} {}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	B(B const& other) : name_{other.name_}, n_{other.n_}, data{std::make_unique<double[]>(n_)} {}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	B(B&&) = delete;

	auto operator=(B const& other) -> B& {
		if(this == &other) {
			return *this;
		}
		name_ = other.name_;
		n_    = other.n_;
		data  = std::make_unique<double[]>(other.n_);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		std::copy_n(other.data.get(), n_, data.get());
		return *this;
	}
	auto operator=(B&& other) -> B& = default;

	auto operator[](std::ptrdiff_t i) const -> double& { return data.get()[i]; }
	~B() = default;
};

// nonintrusive serialization
template<class Archive>
void save(Archive& ar, B const& b, unsigned int const /*version*/) {
	ar << b.name_ << b.n_ << boost::serialization::make_array(b.data.get(), b.n_);
}
template<class Archive>
void load(Archive& ar, B& b, unsigned int const /*version*/) {  // NOLINT(google-runtime-references): serialization protocol
	ar >> b.name_ >> b.n_;
	b.data = std::make_unique<double[]>(b.n_);  // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	ar >> boost::serialization::make_array(b.data.get(), b.n_);
}
BOOST_SERIALIZATION_SPLIT_FREE(B)

template<> struct mpi3::datatype<
	std::pair<std::complex<double>, std::complex<double>>
> : struct_<
	std::complex<double>,
	std::complex<double>
> {};

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {

	assert(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		std::vector<std::vector<double>> buffer(10, std::vector<double>(20));
		buffer[4][5] = 6.1;
		world.send(buffer.begin(), buffer.end(), 1, 123);
		break;
	};
	case 1: {
		std::vector<std::vector<double>> buffer(10, std::vector<double>(20));
		world.receive(buffer.begin(), buffer.end(), 0, 123);
		assert(buffer[4][5] == 6.1);
		break;
	};
	}
	switch(world.rank()) {
	case 0: {
		std::vector<double> buffer(10);
		iota(begin(buffer), end(buffer), 0.0);
		world.send(begin(buffer), end(buffer), 1);
		break;
	};
	case 1: {
		std::vector<double> v(10);
		auto                it = world.receive_n(begin(v), 10, 0);
		assert(it == end(v) and v[3] == 3.0);
		break;
	};
	}
	switch(world.rank()) {
	case 0: {
		std::map<int, std::vector<double>> m;
		m[2] = std::vector<double>(2);
		m[5] = std::vector<double>(5);
		world.send(begin(m), end(m), 1, 123);
		break;
	};
	case 1: {
		std::vector<std::pair<int, std::vector<double>>> v(2);
		world.receive(begin(v), end(v), 0, 123);
		assert((v[1] == std::pair<int, std::vector<double>>{5, std::vector<double>(5)}));
		break;
	};
	}

	switch(world.rank()) {
	case 0: {
		std::vector<A> v(5, A(3));
		v[2][2] = 3.14;
		world.send(begin(v), end(v), 1, 123);
		break;
	};
	case 1: {
		std::vector<A> v(5);
		world.receive(begin(v), end(v), 0, 123);
		assert(v[2][2] == 3.14);
		break;
	};
	}

	switch(world.rank()) {
	case 0: {
		std::vector<B> v(5, B(3));
		v[2][2] = 3.14;
		world.send(begin(v), end(v), 1, 123);
		break;
	};
	case 1: {
		std::vector<B> v(5);
		world.receive(begin(v), end(v), 0, 123);
		assert(v[2][2] == 3.14);
		break;
	};
	}

	switch(world.rank()) {
	case 0: {
		std::vector<std::pair<std::complex<double>, std::complex<double>> > v(5);
		v[2] = std::make_pair(std::complex<double>{3.14, 6.28}, std::complex<double>{4.0, 5.0});
		world.send(begin(v), end(v), 1);
		break;
	};
	case 1: {
		std::vector<std::pair<std::complex<double>, std::complex<double>> > v(5);
		world.receive(begin(v), end(v), 0);
		assert( v[2] == std::make_pair(std::complex<double>{3.14, 6.28}, std::complex<double>{4.0, 5.0}) );
		break;
	};
	}

	return 0;
} catch(...) {
	return 1;
}
