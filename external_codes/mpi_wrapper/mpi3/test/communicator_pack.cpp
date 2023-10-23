#if COMPILATION_INSTRUCTIONS
mpic++ -std=c++14 -O3 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"

#include "../../mpi3/detail/package.hpp"

#include<list>
#include<vector>

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	{
		mpi3::detail::package p(world);

		int i = 123;
		std::vector<int> c(10); iota(begin(c), end(c), 0);

		p.pack_n(&i, 1);
		std::istringstream iss("0 1 2 3 4 5 6 7 8 9");
		p.pack(std::istream_iterator<int>(iss), std::istream_iterator<int>());
		int i2;
		std::list<int> c2(10);
		p.unpack_n(&i2, 1);
		p.unpack(begin(c2), end(c2));
		assert(i2 == i);
		assert(equal(begin(c2), end(c2), begin(c)));
	}
	{
		assert(world.size() > 1);

		mpi3::detail::package p(world);

		int i;
		std::vector<char> c(100);

		switch(world.rank()){
			case 0:
				i = 123;
				iota(c.begin(), c.end(), 0);
				p.pack_n(&i, 1);
				p.pack(begin(c), end(c));
				p.send(1);
				break;
			case 1:
				p.receive(0);
				p.unpack_n(&i, 1);
				p.unpack(begin(c), end(c));
				assert(i == 123);
				std::vector<char> c2(100);
				iota(begin(c2), end(c2), 0);
				assert(c == c2);
				break;
		}

	}
	return 0;
}

