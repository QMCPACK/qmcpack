// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/request.hpp>

#include <array>
#include <boost/core/lightweight_test.hpp>
#include <tuple>
#include <vector>

namespace mpi3 = boost::mpi3;

namespace {

template<class CV, class CR>
void halo_sayan(mpi3::communicator& comm, CV& cv, CR& cr) {
	std::vector<mpi3::request> rs;

	if (comm.rank() != 0) {
		rs.push_back(comm.ireceive(cr[0].begin(), cr[0].end(), comm.rank() - 1));
		rs.push_back(comm.isend   (cv[1].begin(), cv[1].end(), comm.rank() - 1));
	}

	if (comm.rank() != comm.size()-1) {
		rs.push_back(comm.ireceive(cr[1].begin(), cr[1].end(), comm.rank() + 1));
		rs.push_back(comm.isend   (cv[0].begin(), cv[0].end(), comm.rank() + 1));
	}

	wait_all(rs.begin(), rs.end());
}

template<class CV, class CR>
void halo_npos(mpi3::communicator& comm, CV& cv, CR& cr) {
	std::vector<mpi3::request> rs;

	rs.push_back(comm.ireceive(cr[0].begin(), cr[0].end(), (comm.rank()!=0             )?comm.rank() - 1:mpi3::communicator::nproc));
	rs.push_back(comm.isend   (cv[1].begin(), cv[1].end(), (comm.rank()!=0             )?comm.rank() - 1:mpi3::communicator::nproc));

	rs.push_back(comm.ireceive(cr[1].begin(), cr[1].end(), (comm.rank()!= comm.size()-1)?comm.rank() + 1:mpi3::communicator::nproc));
	rs.push_back(comm.isend   (cv[0].begin(), cv[0].end(), (comm.rank()!= comm.size()-1)?comm.rank() + 1:mpi3::communicator::nproc));

	wait_all(rs.begin(), rs.end());
}

template<class CV, class CR>
void halo_init(mpi3::communicator& comm, CV& cv, CR& cr) {
	auto rs = std::tuple{
		comm.ireceive(cr[0].begin(), cr[0].end(), comm.prev()),
		comm.isend   (cv[1].begin(), cv[1].end(), comm.prev()),
		comm.ireceive(cr[1].begin(), cr[1].end(), comm.next()),
		comm.isend   (cv[0].begin(), cv[0].end(), comm.next())
	};

	using std::get;
	wait(std::move(get<0>(rs)), std::move(get<1>(rs)), std::move(get<2>(rs)), std::move(get<3>(rs)));
}

template<class CV, class CR>
void halo_wait(mpi3::communicator& comm, CV& cv, CR& cr) {
	wait(
		comm.ireceive(cr[0].begin(), cr[0].end(), comm.prev()),
		comm.isend   (cv[1].begin(), cv[1].end(), comm.prev()),
		comm.ireceive(cr[1].begin(), cr[1].end(), comm.next()),
		comm.isend   (cv[0].begin(), cv[0].end(), comm.next())
	);
}

template<class CV, class CR>
void halo_temporaries(mpi3::communicator& comm, CV& cv, CR& cr) {
	(
		comm.ireceive(cr[0].begin(), cr[0].end(), comm.prev()),
		comm.isend   (cv[1].begin(), cv[1].end(), comm.prev()),
		comm.ireceive(cr[1].begin(), cr[1].end(), comm.next()),
		comm.isend   (cv[0].begin(), cv[0].end(), comm.next())
	);
}

}  // end namespace

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	std::vector<std::array<int, 4>> cv(2); cv[0][0] = world.rank()  ; cv[1][0] = -world.rank();
	std::vector<std::array<int, 4>> cr(2); cr[0][0] = world.rank()*2; cr[1][0] = -world.rank()*3;

	halo_sayan(world, cv, cr);

	{
		std::vector<std::array<int, 4>> _cv(2);  _cv[0][0] = world.rank()  ; _cv[1][0] = -world.rank();
		std::vector<std::array<int, 4>> _cr(2);  _cr[0][0] = world.rank()*2; _cr[1][0] = -world.rank()*3;

		halo_npos(world, _cv, _cr);

		BOOST_TEST( cv == _cv );
		BOOST_TEST( cr == _cr );
	}
	{
		std::vector<std::array<int, 4>> _cv(2);  _cv[0][0] = world.rank()  ; _cv[1][0] = -world.rank();
		std::vector<std::array<int, 4>> _cr(2);  _cr[0][0] = world.rank()*2; _cr[1][0] = -world.rank()*3;

		halo_init(world, _cv, _cr);

		BOOST_TEST( cv == _cv );
		BOOST_TEST( cr == _cr );
	}
	{
		std::vector<std::array<int, 4>> _cv(2);  _cv[0][0] = world.rank()  ; _cv[1][0] = -world.rank();
		std::vector<std::array<int, 4>> _cr(2);  _cr[0][0] = world.rank()*2; _cr[1][0] = -world.rank()*3;

		halo_wait(world, _cv, _cr);

		BOOST_TEST( cv == _cv );
		BOOST_TEST( cr == _cr );
	}
	{
		std::vector<std::array<int, 4>> _cv(2);  _cv[0][0] = world.rank()  ; _cv[1][0] = -world.rank();
		std::vector<std::array<int, 4>> _cr(2);  _cr[0][0] = world.rank()*2; _cr[1][0] = -world.rank()*3;

		halo_temporaries(world, _cv, _cr);

		BOOST_TEST( cv == _cv );
		BOOST_TEST( cr == _cr );
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
