// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <list>
#include <numeric>  // for std::iota

namespace mpi3 = boost::mpi3;

auto main() -> int try {
	mpi3::environment env;
	mpi3::communicator world = env.world();

	{
		std::vector<double> v(static_cast<std::size_t>(world.size()));
		iota(begin(v), end(v), 0);
		std::vector<double> w(1);
		world.scatter(begin(v), end(v), begin(w), 0);
		assert( w[0] == world.rank() );
	}
	{
		std::vector<double> v(world.root()?static_cast<std::size_t>(world.size()):0);
		iota(begin(v), end(v), 0);
		std::vector<double> w(1);
	//	auto e = 
		world.scatter_from(begin(w), end(w), begin(v), 0);
	//	assert( e == end(v) );
		assert( w[0] == world.rank() );
	}
	{
		std::vector<double> v(world.root()?static_cast<std::size_t>(world.size()):0);
		iota(begin(v), end(v), 0);
		double w = -1;
	//	auto e = 
		world.scatter_value_from(w, begin(v), 0);
	//	assert( e == end(v) );
		assert( w == world.rank() );
	}
	{
		std::vector<double> v(world.root()?static_cast<std::size_t>(world.size()):0);
		iota(begin(v), end(v), 0);
		double w = world.scatter(begin(v), end(v), 0);
		assert( w == world.rank() );
	}
	{
		std::vector<double> v(world.root()?static_cast<std::size_t>(world.size()):0);
		iota(begin(v), end(v), 0);
		double w = v / world;
		assert( w == world.rank() );
	}
	{
		std::list<double> l(world.root()?static_cast<std::size_t>(world.size()):0);
		iota(begin(l), end(l), 0);
		double w = l / world;
		assert( w == world.rank() );
	}
#if 0
	{
		std::vector<double> v = {1, 2, 2, 3, 3, 3}; if(!world.root()) v.clear();
		std::vector<double> w(world.rank() + 1);
		auto e = world.scatterv_from(begin(w), end(w), begin(v), 0);
		assert( end(v) == e );
		switch(world.rank()){
			case 0: assert((w==std::vector<double>{1}    )); break;
			case 1: assert((w==std::vector<double>{2,2}  )); break;
			case 2: assert((w==std::vector<double>{3,3,3})); break;
		}
	}
	{
		if(auto duo = (world < 2)){
			assert( duo.size() == 2 );
		//	std::vector<std::vector<double>> vs = { {1}, {2, 2} };
			std::vector<double> v = {1, 2, 2};
			std::vector<int> ns = {1, 2}; if(!duo.root()) ns.clear();
			std::vector<std::vector<double>::iterator> its = {v.begin(), v.begin()+1}; if(!world.root()) its.clear();
			std::vector<double> w(duo.rank()+1);
			duo.scatterv_n(begin(its), begin(ns), begin(w));
			switch(duo.rank()){
				case 0: assert(( w == std::vector<double>{1} )); break;
				case 1: std::cerr <<"w2:" << w[0] << ", " << w[1] << std::endl; break;
			}
		}
	}
	{
		if(auto duo = (world < 2)){
			assert( duo.size() == 2 );
			std::vector<std::vector<double>> vs = { {1}, {2, 2} }; if(!duo.root()) vs.clear();
			std::vector<double> w(duo.rank()+1);
			duo.scatterv(begin(vs), begin(w));
		//	switch(duo.rank()){
		//		case 0: assert(( w == std::vector<double>{1} )); break;
		//		case 1: std::cerr <<"w2:" << w[0] << ", " << w[1] << std::endl; break;
		//	}
		}
	}
#endif
	return 0;
} catch(...) { return 1; }
