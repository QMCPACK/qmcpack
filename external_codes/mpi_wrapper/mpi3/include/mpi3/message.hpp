// Copyright 2018-2025 Alfredo A. Correa

#ifndef BMPI3_MESSAGE_HPP
#define BMPI3_MESSAGE_HPP

#include <mpi3/detail/iterator_traits.hpp>
#include <mpi3/detail/value_traits.hpp>
#include <mpi3/vector.hpp>

#include <mpi3/detail/mpi_impl.h>

#include <iterator>

namespace boost {
namespace mpi3 {

#ifndef EXAMPI
class message {
 public:
	MPI_Message impl_;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa)
	MPI_Message operator&() const { return impl_; }  // NOLINT(google-runtime-operator) design

	template<class It, typename Size>
	auto receive_n(
		It it,
		/**/ detail::contiguous_iterator_tag /*contiguous*/,
		/**/ detail::basic_tag /*basic*/,
		Size count
	) -> It {
		// Avoid dereferencing iterator when count==0 (MSVC debug asserts on end iterators)
		typename std::iterator_traits<It>::pointer ptr = (count != 0) ? detail::data(it) : nullptr;
		MPI_(Mrecv)(
			ptr, static_cast<int>(count), detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},  // TODO(correaa) use safe cast
#define OMPI_SKIP_MPICXX 1  // NOLINT(cppcoreguidelines-macro-usage) // https://github.com/open-mpi/ompi/issues/5157
			&impl_, MPI_STATUS_IGNORE  // cppcheck-suppress[cstyleCast,intToPointerCast]
		);
		return std::next(it, static_cast<typename std::iterator_traits<It>::difference_type>(count));
	}

	template<class It, typename Size>
	auto receive_n(
		It it,
		/**/ boost::mpi3::detail::forward_iterator_tag /*contiguous*/,
		/**/ detail::basic_tag /*basic*/,
		Size count
	) -> It {
		boost::mpi3::uvector<typename std::iterator_traits<It>::value_type> rc(static_cast<typename boost::mpi3::uvector<typename std::iterator_traits<It>::value_type>::size_type>(static_cast<std::size_t>(count)));
		receive_n(rc.begin(), count);
		return std::move(rc.begin(), rc.end(), it);
	}

	template<class It, typename Size>
	auto receive_n(It it, Size count) -> It {
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
	) -> void {
		receive_n(first, std::distance(first, last));
	}
	template<class It>
	auto receive(It first, It last) -> void {
		receive(
			first, last,
			detail::iterator_category_t<It>{}
		);
	}
};
#endif

}  // end namespace mpi3
}  // end namespace boost
#endif
