#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -O3 -Wall -std=c++17`#-Wextra` -Wfatal-errors  -I${HOME}/prj -D_TEST_BOOST_MULTI_ORDERING $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_ORDERING_HPP
#define BOOST_MULTI_ORDERING_HPP

#include<cstddef> // ptrdiff_t
#include<utility> // index_sequence

namespace boost{
namespace multi{

using dimensionality_type = std::ptrdiff_t;

//template<dimensionality_type D>
//using fortran_ordering = std::make_index_sequence<D>;

template<std::size_t C, std::size_t N, std::size_t... Is>
struct c_ordering_aux : c_ordering_aux<C - 1, N, N - C, Is...>{};

//template<std::size_t N, std::size_t... Is>
//struct c_ordering_aux<0, N, Is...>{using type = std::index_sequence<N, Is...>;};

//template<dimensionality_type D>
//using c_ordering = typename c_ordering_aux<D-1, D-1>::type;

}}


#if _TEST_BOOST_MULTI_ORDERING

#include<cassert>
#include<iostream>

namespace bm = boost::multi;

struct node{
	int val;
	node const& next;
};

struct Node{
	int val;
	Node const* const nextp;
};

int main(){

	node s{3, {4, {5, s}}};

	assert(&(s.next.next.next) == &s);

	assert(s.val == 3);
	assert(s.next.val == 4);
	assert(s.next.next.val == 5);


	Node tt{3, &tt};
//	Node t{3, &t};
//	Node tt{3, &Node{4, &t}}
//	Node t1;
//	Node t2;
//	t1 = {3, &t2};
//	t2 = {4, &t1};
//	Node [t1, t2] = std::pair{Node{3, &t2}, Node{4, &t1}};
}

#endif
#endif

