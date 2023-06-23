/* -*- indent-tabs-mode: t -*- */

#ifndef MPI3_MESSAGE_HPP
#define MPI3_MESSAGE_HPP

#include "../mpi3/detail/iterator_traits.hpp"
#include "../mpi3/detail/value_traits.hpp"

#include<mpi.h>

namespace boost {
namespace mpi3 {

class message {
 public:
	MPI_Message impl_;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa)

	MPI_Message operator&() const {return impl_;}  // NOLINT(google-runtime-operator) design

	template<class It, typename Size>
	auto receive_n(
		It it,
		/**/ detail::contiguous_iterator_tag /*contiguous*/,
		/**/ detail::basic_tag /*basic*/,
		Size count
	) {
		MPI_(Mrecv)(
			detail::data(it), static_cast<int>(count), detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},  // TODO(correaa) use safe cast
			&impl_, MPI_STATUS_IGNORE
		);
	}
	template<class It, typename Size>
	auto receive_n(It it, Size count) {
		return receive_n(
			it,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count
		);
	}
	template<class It>
	auto receive(
		It first, It last, 
			detail::random_access_iterator_tag /*random_access*/
	) {
		return receive_n(first, std::distance(first, last));
	}
	template<class It>
	auto receive(It first, It last) {
		return receive(
			first, last,
				detail::iterator_category_t<It>{}
		);
	}
};

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_MESSAGE

//#include "../mpi3/main.hpp"

//namespace mpi3 = boost::mpi3;

//int mpi3::main(int, char*[], mpi3::communicator world){
//}

//#endif
#endif

