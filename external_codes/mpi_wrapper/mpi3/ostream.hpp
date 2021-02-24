#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_MPI3_OSTREAM $0x.cpp -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_OSTREAM_HPP
#define MPI3_OSTREAM_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

#include<boost/icl/split_interval_map.hpp>

#include<iostream>
#include<sstream>

namespace boost{
namespace mpi3{

struct ostream : public std::ostream{ // http://stackoverflow.com/a/2212940/225186
	class streambuf : public std::stringbuf{
		communicator& comm_;
		std::ostream& output;
		std::string msg_;
		public:
		streambuf(communicator& comm, std::ostream& strm = std::cout) : 
			comm_(comm), output(strm)
		{}
		virtual int sync(){
			// following code can be improved by a custom reduce operation
			if(comm_.root()){
				boost::icl::interval_map<int, std::string> messages;
				messages.insert(std::make_pair(0, str()));
				for(int i = 1; i != comm_.size(); ++i){
					match m = comm_.matched_probe(i);
					msg_.resize(m.count<char>());
					m.receive(msg_.begin());
					messages.insert(std::make_pair(i, msg_));
				}
				for(auto& m : messages){
					output << comm_.name();
					if((int)size(m.first) < (int)comm_.size()){
						if(size(m.first) == 1) output<<"["<< lower(m.first) <<"]";
						else output<<"["<< m.first <<"]";
					}
					output<<"\t: "<< m.second;
				}
			}else comm_.send_n(str().begin(), str().size(), 0);
			str("");
			output.flush();
			comm_.barrier();
			return 0;
		}
	};
	streambuf buffer;
public:
	ostream(communicator& comm, std::ostream& os = std::cout) : 
		std::ostream(&buffer), buffer(comm, os)
	{}
	~ostream(){flush();}
};

}}

#ifdef _TEST_MPI3_OSTREAM

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::ostream wout(world);
	wout << "hello" << std::endl;
	wout << "hello, I am rank " << world.rank() << " in " << world.name() << std::endl;
	wout << "hello, my rank/2 is " << world.rank()/2 << std::endl;
	wout << (not world.root()?"not root":"") << std::endl;
	mpi3::communicator firsttwo = (world < 2);
	if(firsttwo) firsttwo.name("firsttwo");

	if(firsttwo){
		mpi3::ostream fout(firsttwo);
		fout << "hola, I am rank " << firsttwo.rank() << " in " << firsttwo.name() << " and also rank " << world.rank() << " in " << world.name() << std::endl;
	}

	return 0;
}

/* output: 
world: hello
world[0]: hello, I am rank 0 in world
world[1]: hello, I am rank 1 in world
world[2]: hello, I am rank 2 in world
world[3]: hello, I am rank 3 in world
world[0-1]: hello, my rank/2 is 0
world[2-3]: hello, my rank/2 is 1
firsttwo[0]: hola, I am rank 0 in firsttwo
firsttwo[1]: hola, I am rank 1 in firsttwo
firsttwo[2]: hola, I am rank 2 in firsttwo
*/

#endif
#endif

