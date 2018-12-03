#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_MPI3_MESSAGE $0x.cpp -o $0x.x && time mpirun -np 8 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_MESSAGE_HPP
#define MPI3_MESSAGE_HPP

#include "../mpi3/detail/iterator_traits.hpp"
#include "../mpi3/detail/value_traits.hpp"

#include<mpi.h>

namespace boost{
namespace mpi3{

class message{
public:
	MPI_Message impl_;
	MPI_Message operator&() const{return impl_;}
	template<class It, typename Size>
	auto receive_n(
		It it, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count
	){
		int s = MPI_Mrecv(detail::data(it), count, detail::basic_datatype<typename std::iterator_traits<It>::value_type>{}, &impl_, MPI_STATUS_IGNORE);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot message receive");
	}
	template<class It, typename Size>
	auto receive_n(It it, Size count){
		return receive_n(
			it, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count
		);
	}
	template<class It>
	auto receive(It first, It last, detail::random_access_iterator_tag){
		return receive_n(first, std::distance(first, last));
	}
	template<class It>
	auto receive(It first, It last){
		return receive(
			first, last,
				detail::iterator_category_t<It>{}
		);
	}
//protected:
};

}}

#ifdef _TEST_MPI3_MESSAGE

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
}

#endif
#endif

