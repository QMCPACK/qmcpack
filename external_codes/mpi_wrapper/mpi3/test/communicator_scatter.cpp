// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>

#include <cstddef>
#include <list>
#include <numeric>  // for std::iota
#include <vector>

namespace mpi3 = boost::mpi3;

auto main() -> int try {
	mpi3::environment env;

	mpi3::communicator world = env.world();

	{
		std::vector<double> v(static_cast<std::size_t>(world.size()));
		iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		std::vector<double> w(1);
		world.scatter(begin(v), end(v), begin(w), 0);
		BOOST_TEST( w[0] == world.rank() );
	}
	{
		std::vector<double> v(world.root() ? static_cast<std::size_t>(world.size()) : 0);
		iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		std::vector<double> w(1);
		//  auto e =
		world.scatter_from(begin(w), end(w), begin(v), 0);
		//  BOOST_TEST( e == end(v) );
		BOOST_TEST( w[0] == world.rank() );
	}
	{
		std::vector<double> v(world.root() ? static_cast<std::size_t>(world.size()) : 0);
		iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		double w = -1;
		//  auto e =
		world.scatter_value_from(w, begin(v), 0);
		//  BOOST_TEST( e == end(v) );
		BOOST_TEST( w == world.rank() );
	}
	{
		std::vector<double> v(world.root() ? static_cast<std::size_t>(world.size()) : 0);
		iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		double const w = world.scatter(begin(v), end(v), 0);
		BOOST_TEST( w == world.rank() );
	}
	{
		std::vector<double> v(world.root() ? static_cast<std::size_t>(world.size()) : 0);
		iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		double const w = v / world;
		BOOST_TEST( w == world.rank() );
	}
	{
		std::list<double> l(world.root() ? static_cast<std::size_t>(world.size()) : 0);
		iota(begin(l), end(l), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		double const w = l / world;
		BOOST_TEST( w == world.rank() );
	}
// #if 0
//  {
//      std::vector<double> v = {1, 2, 2, 3, 3, 3}; if(!world.root()) v.clear();
//      std::vector<double> w(world.rank() + 1);
//      auto e = world.scatterv_from(begin(w), end(w), begin(v), 0);
//      BOOST_TEST( end(v) == e );
//      switch(world.rank()){
//          case 0: BOOST_TEST((w==std::vector<double>{1}    )); break;
//          case 1: BOOST_TEST((w==std::vector<double>{2,2}  )); break;
//          case 2: BOOST_TEST((w==std::vector<double>{3,3,3})); break;
//      }
//  }
//  {
//      if(auto duo = (world < 2)){
//          BOOST_TEST( duo.size() == 2 );
//      //  std::vector<std::vector<double>> vs = { {1}, {2, 2} };
//          std::vector<double> v = {1, 2, 2};
//          std::vector<int> ns = {1, 2}; if(!duo.root()) ns.clear();
//          std::vector<std::vector<double>::iterator> its = {v.begin(), v.begin()+1}; if(!world.root()) its.clear();
//          std::vector<double> w(duo.rank()+1);
//          duo.scatterv_n(begin(its), begin(ns), begin(w));
//          switch(duo.rank()){
//              case 0: BOOST_TEST(( w == std::vector<double>{1} )); break;
//              case 1: std::cerr <<"w2:" << w[0] << ", " << w[1] << std::endl; break;
//          }
//      }
//  }
//  {
//      if(auto duo = (world < 2)){
//          BOOST_TEST( duo.size() == 2 );
//          std::vector<std::vector<double>> vs = { {1}, {2, 2} }; if(!duo.root()) vs.clear();
//          std::vector<double> w(duo.rank()+1);
//          duo.scatterv(begin(vs), begin(w));
//      //  switch(duo.rank()){
//      //      case 0: BOOST_TEST(( w == std::vector<double>{1} )); break;
//      //      case 1: std::cerr <<"w2:" << w[0] << ", " << w[1] << std::endl; break;
//      //  }
//      }
//  }
// #endif

	return boost::report_errors();
} catch(...) {
	return 1;
}
