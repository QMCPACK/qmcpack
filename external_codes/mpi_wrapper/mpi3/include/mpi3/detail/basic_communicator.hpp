// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef MPI3_DETAIL_BASIC_COMMUNICATOR_HPP
#define MPI3_DETAIL_BASIC_COMMUNICATOR_HPP

#include "../../mpi3/vector.hpp"

#include "../../mpi3/detail/buffer.hpp"
#include "../../mpi3/detail/iterator_traits.hpp"
#include "../../mpi3/detail/value_traits.hpp"

#include "../../mpi3/match.hpp"
#include<mpi.h>

#include<algorithm>

#include<cassert>

namespace boost {
namespace mpi3 {
namespace detail {

class basic_communicator{
 protected:
	using impl_t = MPI_Comm;
	impl_t impl_ = MPI_COMM_NULL;  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes,misc-non-private-member-variables-in-classes) : TODO(correaa) make private

	explicit basic_communicator(MPI_Comm impl) noexcept : impl_(impl){}

 public:
	basic_communicator() noexcept = default; //: impl_(MPI_COMM_NULL){}
	~basic_communicator() = default;

	auto operator=(basic_communicator const&) -> basic_communicator& = delete;
	auto operator=(basic_communicator&&)      -> basic_communicator& = delete;

	basic_communicator(basic_communicator const&) = delete; // communicators are not copyable, only duplicable, if you know what you are doing, use `mutable`
	basic_communicator(basic_communicator& other) {
		if(MPI_COMM_NULL != other.impl_) {MPI_(Comm_dup)(other.impl_, &impl_);}
	}
	basic_communicator(basic_communicator&& other) noexcept :
		impl_{std::exchange(other.impl_, MPI_COMM_NULL)} {}

	// [[deprecated]] void swap(basic_communicator& o) noexcept {std::swap(impl_, o.impl_);}
	// [[deprecated]] friend void swap(basic_communicator& self, basic_communicator& other) noexcept {self.swap(other);}

	template<class T>
	int pack_size(int count, detail::basic_tag /*tag*/) {
		int size = -1;
		MPI_(Pack_size)(count, detail::basic_datatype<T>{}, impl_, &size);
		return size;
	}
	template<class T>
	int pack_size(int count) {
		return pack_size<T>(count, detail::value_category_t<T>{});
	}
	template<class It, typename Size>
	auto pack_n(
		It first,
			detail::contiguous_iterator_tag /*contiguous*/,
			detail::basic_tag /*basic*/,
		Size count,
		uvector<detail::packed>& p, int pos  // NOLINT(misc-unused-parameters) bug in clang-tidy 12
	) {
		using value_type = typename std::iterator_traits<It>::value_type;
		MPI_(Pack)(
			detail::data(first), static_cast<int>(count), detail::basic_datatype<value_type>{},
			p.data(), static_cast<int>(p.size()),  // TODO(correaa) use safe cast
			&pos, impl_
		);
		return pos;
	}
	template<class It, typename Size>
	auto pack_n(
		It first, Size count, 
		uvector<detail::packed>& b, int pos
	){
		return pack_n(
			first, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			b, pos
		);
	}
	template<class It, typename Size>
	auto pack_n(It first, Size count, uvector<detail::packed>& p) {
		assert(p.size() < std::numeric_limits<int>::max());
		int const pos = static_cast<int>(p.size());
		p.resize(p.size() + static_cast<std::size_t>(pack_size<typename std::iterator_traits<It>::value_type>(static_cast<int>(count))));
		return pack_n(first, count, p, pos);
	}
	template<class It>
	auto pack(
		It first, It second,
		detail::random_access_iterator_tag /*random_access*/,
		detail::value_unspecified_tag /*value_unspecified*/,
		uvector<detail::packed>& b
	) {
		return pack_n(first, std::distance(first, second), b);
	}

	template<class It>
	auto pack(
		It first, It last,
		/**/ detail::input_iterator_tag /*input*/,
		uvector<detail::packed>& b
	) {
		std::for_each(first, last, [this, &b](auto& e) { pack_n(std::addressof(e), 1, b); });
		return b.size();
	}

	template<class It>
	auto pack(It first, It second, uvector<detail::packed>& b){
		return pack(
			first, second,
				typename detail::iterator_category<It>::type{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			b
		);
	}
	template<class It, typename Size>
	auto unpack_n(
		uvector<detail::packed>& b, int pos,
		It first,
		/**/ detail::contiguous_iterator_tag /*contiguous*/,
		/**/ detail::basic_tag /*basic*/,
		Size count
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		MPI_(Unpack)(
			b.data(), static_cast<int>(b.size()), &pos,
			detail::data(first), static_cast<int>(count),  // TODO(correaa) use safe cast
			detail::basic_datatype<value_type>{},
			impl_
		);
		return pos;
	}
	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last,
			detail::random_access_iterator_tag /*random_access*/
	) {
		return unpack_n(b, pos, first, std::distance(first, last));
	}

	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last,
			detail::forward_iterator_tag /*forward*/
	) {
		std::for_each(first, last, [&b, &pos, this](auto& e) {pos = unpack_n(b, pos, std::addressof(e), 1);});
		// while(first != last){
		// 	pos = unpack_n(b, pos, std::addressof(*first), 1);
		// 	++first;
		// }
		return pos;
	}
	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last
	) {
		return unpack(
			b, pos,
			first, last,
			/**/ detail::iterator_category_t<It>{}
		);
	}
	template<class It>
	auto unpack(detail::buffer& b, It first, It last){
		return b.pos = unpack(b, b.pos, first, last);
	}
	template<class It, typename Size>
	auto unpack_n(uvector<detail::packed>& b, int pos, It first, Size count){
		return unpack_n(
			b, pos,
			first,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename  std::iterator_traits<It>::value_type>{},
			count
		);
	}
	template<class It, typename Size>
	auto unpack_n(detail::buffer& b, It first, Size count) {
	//	assert(0);
		b.pos = unpack_n(b, b.pos, first, count);
		return b.pos;
	}
	template<class It, typename Size>
	auto send_n(
		It first,
			detail::contiguous_iterator_tag /*contiguous*/,
			detail::basic_tag /*basic*/,
		Size count,
		int dest, int tag
	) {
		MPI_(Send)(
			detail::data(first), static_cast<int>(count),  // TODO(correaa) use safe cast
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, tag, 
			impl_
		);
	}
	template<class It, typename Size>
	auto send_n(It first, Size count, int dest, int tag = 0) {
		return send_n(
			first, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			dest, tag
		);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::random_access_iterator_tag /*random_access*/,
			detail::value_unspecified_tag /*value_unspecified*/,
		int dest, int tag
	) {
		return send_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::input_iterator_tag /*input*/,
			detail::basic_tag /*basic*/,
		int dest, int tag
	) {
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		return send_n(buffer.data(), buffer.size(), dest, tag);
	}
	template<class It>
	auto send(It first, It last, int dest, int tag = 0) {
		return send(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);
	}
	auto send(uvector<detail::packed> const& p, int dest, int tag = 0) {
		return send_n(p.data(), p.size(), dest, tag);
	}
	match matched_probe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		match m;
		MPI_(Mprobe)(source, tag, impl_, &m.message::impl_, &m.status::impl_);
		return m;
	}
	template<class It, typename Size>
	auto receive_n(
		It dest,
		/**/ detail::forward_iterator_tag /*forward*/,
		/**/ detail::basic_tag /*basic*/,
		Size n,
		int source, int tag
	) {
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(n);
		receive_n(buffer.data(), buffer.size(), source, tag);
		return std::copy_n(buffer.begin(), n, dest);
	}
	auto receive(uvector<detail::packed>& b, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		match m = matched_probe(source, tag);
		auto const count = static_cast<std::size_t>(m.count<detail::packed>());
		auto const size = static_cast<std::ptrdiff_t>(b.size());
		b.resize(b.size() + count);
		return m.receive_n(std::next(b.data(), size), count);
	}
	auto receive(detail::buffer& b, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		return receive(static_cast<uvector<detail::packed>&>(b), source, tag);
	}
	template<class It, typename Size, typename... Meta>
	auto send_receive_replace_n(
		It first,
			detail::forward_iterator_tag /*forward*/,
			detail::basic_tag /*basic*/,
		Size size, Meta... meta
	) {
		using value_type = typename std::iterator_traits<It>::value_type;
		mpi3::uvector<value_type> buffer(size);
		std::copy_n(first, buffer.size(), buffer.begin());
		send_receive_replace_n(buffer.begin(), buffer.size(), meta...);
		return std::copy_n(buffer.begin(), buffer.size(), first);
	}

	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::contiguous_iterator_tag /*contiguous*/,
			detail::basic_tag /*basic*/,
		Size size,
		int dest, int source,
		int sendtag, int recvtag
	) {
		using value_type = typename std::iterator_traits<It>::value_type;
		status ret;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) delayed initialization
		int s = MPI_Sendrecv_replace(
			detail::data(first), size, detail::basic_datatype<value_type>{}, 
			dest, sendtag, source, recvtag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) {throw std::runtime_error("cannot send_receive");}
		return first + size;
	}

	template<class It, class Size>
	auto send_receive_replace_n(
		It first, Size size,
		int dest, int source, // = MPI_ANY_SOURCE,
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	) {
		return send_receive_replace_n(
			first, 
				detail::iterator_category_t<It>{}, 
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			size,
			dest, source, sendtag, recvtag
		);
	}
	template<class It, class Size>
	auto send_receive_n(
		It first, Size size,
		int dest, int source, // = MPI_ANY_SOURCE,
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	) {
		return send_receive_replace_n(first, size, dest, source, sendtag, recvtag);
	}
};

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost
#endif
