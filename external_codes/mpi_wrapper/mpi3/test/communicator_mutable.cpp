// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <list>
#include <utility>
#include <vector>

namespace mpi3 = boost::mpi3;

struct projector {
	projector() = default;

	explicit projector(mpi3::communicator comm) : n_{5}, comm_{std::move(comm)} {}

	projector(projector     &&) = default;
#ifdef _MSC_VER
	projector(projector const& other) : n_{other.n_}, comm_{other.comm_} {}  // mutable member: const& gives non-const access. Necessary for MSVC
	auto operator=(projector const& other) -> projector& { n_ = other.n_; comm_ = std::move(mpi3::communicator{other.comm_}); return *this; }  // NOLINT(clang-diagnostic-deprecated-declarations)
#else
	projector(projector const&) = default;
	auto operator=(projector const&) -> projector& = default;  // NOLINT(clang-diagnostic-deprecated-declarations) TODO(correaa) deprecate copy assigment
#endif
	auto operator=(projector     &&) -> projector& = default;  // NOLINT(clang-diagnostic-deprecated-declarations) TODO(correaa) deprecate move assigment
//  auto operator=(projector      &) -> projector& = default;

	friend auto operator==(projector const& a, projector const& b) {return a.n_ == b.n_;} //  a.comm_ == b.comm_;}
	friend auto operator!=(projector const& a, projector const& b) {return a.n_ != b.n_;} //  a.comm_ == b.comm_;}

	decltype(auto) get_comm() const {return comm_;}
	auto get_n() const -> int{return n_;}
	~projector() = default;

 private:
	int n_ = 0;
	mutable mpi3::communicator comm_;
};

struct projector2 {
	projector2() = default;

	explicit projector2(mpi3::communicator comm) : n_{5}, comm_{std::move(comm)} {}

	projector2(projector2     &&) = default;
#ifdef _MSC_VER
	projector2(projector2 const& other) : n_{other.n_}, comm_{other.comm_} {}  // mutable member: const& gives non-const access. Necessary for MSVC
#else
	projector2(projector2 const&) = default;
#endif

	auto operator=(projector2 const& other) -> projector2& {
		if(this == &other) {return *this;}
		BOOST_TEST(comm_.compare(other.comm_) == mpi3::detail::equality::congruent);
		n_ = other.n_;
		return *this;
	}
	auto operator=(projector2     &&) -> projector2& = default;  // NOLINT(clang-diagnostic-deprecated-declarations) TODO(correaa) deprecate move assigment
//  auto operator=(projector2      &) -> projector2& = default;

	friend auto operator==(projector2 const& a, projector2 const& b) {return a.n_ == b.n_;} //  a.comm_ == b.comm_;}
	friend auto operator!=(projector2 const& a, projector2 const& b) {return a.n_ != b.n_;} //  a.comm_ == b.comm_;}

	decltype(auto) get_comm() const { return comm_; }
	auto get_n() const -> int { return n_; }
	~projector2() = default;

 private:
	int n_ = 0;
	mutable mpi3::communicator comm_;
};

struct projector3 {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	projector3() = default;

	explicit projector3(mpi3::communicator comm) : n_{5}, comm_{std::move(comm)} {}

	projector3(projector3     &&) = default;
#ifdef _MSC_VER
	projector3(projector3 const& other) : n_{other.n_}, comm_{other.comm_} {}  // mutable member: const& gives non-const access. Necessary for MSVC
#else
	projector3(projector3 const&) = default;
#endif

	auto operator=(projector3 other) noexcept -> projector3& {
		swap(other);
		return *this;
	}
//  auto operator=(projector3     &&) -> projector3& = default;
//  auto operator=(projector3      &) -> projector3& = default;

	friend auto operator==(projector3 const& a, projector3 const& b) {return a.n_ == b.n_;} //  a.comm_ == b.comm_;}
	friend auto operator!=(projector3 const& a, projector3 const& b) {return a.n_ != b.n_;} //  a.comm_ == b.comm_;}

	decltype(auto) get_comm() const { return comm_; }
	auto get_n() const -> int { return n_; }
	~projector3() = default;

 private:
	int n_ = 0;
	mutable mpi3::communicator comm_;
	void swap(projector3& other) noexcept {
		std::swap(n_, other.n_);
		std::swap(comm_, other.comm_);
	}
};

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

//  {
//      projector const p{world};
//      projector p2;
//      p2 = p;
//  }
	{
		std::list<mpi3::communicator> v;
		v.emplace_back(world);
		v.emplace_back(world);
	}
//  { // doesn't compile, communicator is not copiable
//      std::vector<mpi3::communicator> v = {world, world};
//      v.emplace_back(world);
//      v.emplace_back(world);
//  }
	{ // but this works because the member is mutable
		std::vector<projector> v = {projector{world}, projector{world}};
		v.emplace_back(world);
		v.emplace_back(world);
		v.emplace_back(world);
		v[1] = v[0];
		BOOST_TEST( v[1] == v[0] );
		v[1] = std::move(v[0]);
		BOOST_TEST( v[0].get_comm().is_empty() );
		v[1] = static_cast<projector const&>(v[2]);
		BOOST_TEST( v[2].get_comm() == world );
	}
	{ // but this works because the member is mutable
		std::vector<projector2> v = {projector2{world}, projector2{world}};
		v.emplace_back(world);
		v.emplace_back(world);
		v.emplace_back(world);
		v[1] = v[0];
		BOOST_TEST( v[1] == v[0] );
		v[1] = std::move(v[0]);
		BOOST_TEST( v[0].get_comm().is_empty() );
		v[1] = static_cast<projector2 const&>(v[2]);
		BOOST_TEST( v[2].get_comm() == world );
	}
	{ // but this works because the member is mutable
		std::vector<projector3> v = {projector3{world}, projector3{world}};
		v.emplace_back(world);
		v.emplace_back(world);
		v.emplace_back(world);
		v[1] = v[0];
		BOOST_TEST( v[1] == v[0] );
		v[1] = std::move(v[0]);
		BOOST_TEST( v[0].get_comm().is_empty() );
		v[1] = static_cast<projector3 const&>(v[2]);
		BOOST_TEST( v[2].get_comm() == world );
	}
	{
		mpi3::communicator const comm;
		BOOST_TEST( ! comm );
		BOOST_TEST( comm.is_empty() );
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
