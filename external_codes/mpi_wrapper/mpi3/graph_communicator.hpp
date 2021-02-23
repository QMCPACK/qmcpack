#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_GRAPH_COMMUNICATOR $0x.cpp -o $0x.x && time mpirun -n 5 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_GRAPH_COMMUNICATOR_HPP
#define BOOST_MPI3_GRAPH_COMMUNICATOR_HPP

#include "../mpi3/communicator.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <iterator>
#include <boost/graph/graphviz.hpp>

namespace boost{
namespace mpi3{

struct graph_communicator : communicator{
	using graph = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>;

	std::vector<int> neighbors(int rank) const{
		int nneighbors;
		MPI_Graph_neighbors_count(impl_, rank, &nneighbors);
		std::vector<int> ret(nneighbors);
		MPI_Graph_neighbors(impl_, rank, ret.size(), ret.data());
		return ret;
	}

	graph topology() const{
		int nnodes;
		int nedges;
		MPI_Graphdims_get(impl_, &nnodes, &nedges);
		std::vector<int> index(nnodes);
		std::vector<int> edges(nedges);
		MPI_Graph_get(impl_, nnodes, nedges, index.data(), edges.data());
		auto ite = std::begin(edges);
		int source = 0;
		graph ret;
	//	std::vector<decltype(add_vertex(g))> v_ptr;
	//	for(int i = 0; i != nnodes; ++i){
	//		v_ptr.push_back(add_vertex(g));
	//	}
	//	std::cout << " index: ";
	//	for(auto e: index) std::cout << e << " ";
	//	std::cout << '\n';

	//	std::cout << " edges: ";
	//	for(auto e: edges) std::cout << e << " ";
	//	std::cout << '\n';
	//	return ret;
		int j = 0;
		for(int i : index){
			for(; j != i; ++j){
				add_edge(source, edges[j], ret);
			}
			++source;
		}
		return ret;
	}
};

/*
template<class Graph>
graph_communicator communicator::make_graph(Graph const& g) const{
	graph_communicator ret;
	std::vector<int> indx;
	std::vector<int> edges;
	std::pair<typename Graph::vertex_iterator, typename Graph::vertex_iterator> vs = boost::vertices(g);
	for(auto vit = vs.first; vit != vs.second; ++vit){
		auto neighbors = boost::adjacent_vertices(*vit, g);
		for(auto nit = neighbors.first; nit != neighbors.second; ++nit) edges.push_back(*nit);
		indx.push_back(edges.size()); // http://stackoverflow.com/a/32833024/225186
	}
	MPI_Graph_create(impl_, num_vertices(g), indx.data(), edges.data(), true, &ret.impl_);
	return ret;
}*/

/*
template<class Graph>
int communicator::graph_rank(Graph const& g) const{
	std::vector<int> indx;
	std::vector<int> edges;
	std::pair<typename Graph::vertex_iterator, typename Graph::vertex_iterator> vs = boost::vertices(g);
	for(auto vit = vs.first; vit != vs.second; ++vit){
		auto neighbors = boost::adjacent_vertices(*vit, g);
		for(auto nit = neighbors.first; nit != neighbors.second; ++nit) edges.push_back(*nit);
		indx.push_back(edges.size()); // http://stackoverflow.com/a/32833024/225186
	}
	int ret;
	MPI_Graph_map(impl_, num_vertices(g), indx.data(), edges.data(), &ret);
	return ret;
}*/

}}

#ifdef _TEST_BOOST_MPI3_GRAPH_COMMUNICATOR

#include<iostream>
#include "../mpi3/main.hpp"
#include "../mpi3/version.hpp"
#include "../mpi3/ostream.hpp"

#include<boost/graph/adjacency_list.hpp>
#include<boost/graph/random.hpp>
#include<random>

#include<boost/graph/graphviz.hpp>
#include <boost/graph/isomorphism.hpp>

int boost::mpi3::main(int argc, char* argv[], boost::mpi3::communicator world){
	if(world.size() != 5) return 1;

	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> graph;
	graph g;

	std::default_random_engine rd;
	boost::generate_random_graph(g, 5, 8, rd);

	if(world.rank()==0){
		write_graphviz(std::cout, g);
		std::cout << "---\n";
	}

//	auto world_graph = world.make_graph(g);

/*	auto outg = world_graph.topology();

	if(world.rank()==0){
		write_graphviz(std::cout, outg);
		std::cout << "---\n";
	}

	assert(boost::isomorphism(outg, g));
*/

	return 0;

}

#endif
#endif


