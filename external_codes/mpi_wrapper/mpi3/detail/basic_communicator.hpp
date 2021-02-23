#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
mpic++ -x c++ $0 -o $0x&&mpirun -n 1 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MPI3_DETAIL_BASIC_COMMUNICATOR_HPP
#define MPI3_DETAIL_BASIC_COMMUNICATOR_HPP
#define OMPI_SKIP_MPICXX 1 // workaround for https://github.com/open-mpi/ompi/issues/5157

#include "../../mpi3/vector.hpp"

#include "../../mpi3/detail/buffer.hpp"
#include "../../mpi3/detail/iterator_traits.hpp"
#include "../../mpi3/detail/value_traits.hpp"

#include "../../mpi3/match.hpp"
#include<mpi.h>

#include<algorithm>

#include<cassert>

namespace boost{
namespace mpi3{
namespace detail{

class basic_communicator{// : public detail::caller<communicator, MPI_Comm>{
protected:
	using impl_t = MPI_Comm;
	impl_t impl_ = MPI_COMM_NULL;
	basic_communicator(MPI_Comm impl) noexcept : impl_(impl){}
public:
	basic_communicator() noexcept = default; //: impl_(MPI_COMM_NULL){}
	[[deprecated("communicators are not values, they are not copiable, only duplicable; they cannot be elements of std::vector")]] 
	basic_communicator(basic_communicator const& other){// = delete;/*{
		if(MPI_COMM_NULL != other.impl_){
			int s = MPI_Comm_dup(other.impl_, &impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot duplicate communicator");
		}
	}
	basic_communicator(basic_communicator& other){
		if(MPI_COMM_NULL != other.impl_){
			int s = MPI_Comm_dup(other.impl_, &impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot duplicate communicator");
		}
	}
	basic_communicator(basic_communicator&& other) noexcept : 
		impl_{std::exchange(other.impl_, MPI_COMM_NULL)}
	{}
	void swap(basic_communicator& other) noexcept{
		std::swap(impl_, other.impl_);
	}
	template<class T>
	int pack_size(int count, detail::basic_tag){
		int size = -1;
		int s = MPI_Pack_size(count, detail::basic_datatype<T>{}, impl_, &size);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot pack size");
		return size;
	}
	template<class T>
	int pack_size(int count){
		return pack_size<T>(count, detail::value_category_t<T>{});
	}
	template<class It, typename Size>
	auto pack_n(
		It first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		uvector<detail::packed>& p, int pos
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		int s = MPI_Pack(
			detail::data(first), count, 
			detail::basic_datatype<value_type>{}, 
			p.data(), p.size(), 
			&pos, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot pack_n");
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
	auto pack_n(It first, Size count, uvector<detail::packed>& p){
		int pos = p.size();
		p.resize(p.size() + pack_size<typename std::iterator_traits<It>::value_type>(count));
		return pack_n(first, count, p, pos);
	}
	template<class It>
	auto pack(
		It first, It second, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		uvector<detail::packed>& b
	){
		return pack_n(first, std::distance(first, second), b);
	}
	template<class It>
	auto pack(
		It first, It second, 
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		uvector<detail::packed>& b
	){
		while(first != second){
			pack_n(std::addressof(*first), 1, b);
			++first;
		}
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
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		int s = MPI_Unpack(
			b.data(), b.size(), &pos, 
			detail::data(first), count, 
			detail::basic_datatype<value_type>{}, 
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot unpack_n");
		return pos;
	}
	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last,
			detail::random_access_iterator_tag
	){
		return unpack_n(b, pos, first, std::distance(first, last));
	}
	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last, 
			detail::forward_iterator_tag
	){
		while(first != last){
			pos = unpack_n(b, pos, std::addressof(*first), 1);
			++first;
		}
		return pos;
	}
	template<class It>
	auto unpack(
		uvector<detail::packed>& b, int pos,
		It first, It last
	){
		return unpack(
			b, pos, 
			first, last, 
				detail::iterator_category_t<It>{}
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
	auto unpack_n(detail::buffer& b, It first, Size count){
	//	assert(0);
		b.pos = unpack_n(b, b.pos, first, count); 
		return b.pos;
	}
	template<class It, typename Size>
	auto send_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int dest, int tag
	){
		MPI_Send(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, tag, 
			impl_
		);
	}
	template<class It, typename Size>
	auto send_n(It first, Size count, int dest, int tag = 0){
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
			detail::random_access_iterator_tag, 
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return send_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::input_iterator_tag, 
			detail::basic_tag,
		int dest, int tag
	){
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		return send_n(buffer.data(), buffer.size(), dest, tag);
	}
	template<class It>
	auto send(It first, It last, int dest, int tag = 0){
		return send(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);
	}
	auto send(uvector<detail::packed> const& p, int dest, int tag = 0){
		return send_n(p.data(), p.size(), dest, tag);
	}
	match matched_probe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		match m;
		int s = MPI_Mprobe(source, tag, impl_, &m.message::impl_, &m.status::impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot mprobe");
		return m;
	}
	template<class It, typename Size>
	auto receive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		Size n, 
		int source, int tag
	){
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(n);
		receive_n(buffer.data(), buffer.size(), source, tag);
		return std::copy_n(buffer.begin(), n, dest);
	}
	auto receive(uvector<detail::packed>& b, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		match m = matched_probe(source, tag);
		auto count = m.count<detail::packed>();
		auto size = b.size();
		b.resize(b.size() + count);
		return m.receive_n(b.data() + size, count);
	}
	auto receive(detail::buffer& b, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(static_cast<uvector<detail::packed>&>(b), source, tag);
	}
	template<class It, typename Size, typename... Meta>
	auto send_receive_replace_n(
		It first, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		Size size, Meta... meta
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		mpi3::uvector<value_type> buffer(size);
		std::copy_n(first, buffer.size(), buffer.begin());
		send_receive_replace_n(buffer.begin(), buffer.size(), meta...);
		return std::copy_n(buffer.begin(), buffer.size(), first);
	}
	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size size,
		int dest, int source, 
		int sendtag, int recvtag
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		status ret;
		int s = MPI_Sendrecv_replace(
			detail::data(first), size, detail::basic_datatype<value_type>{}, 
			dest, sendtag, source, recvtag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send_receive");
		return first + size;
	}
	template<class It, class Size>
	auto send_receive_replace_n(
		It first, Size size, 
		int dest, int source, // = MPI_ANY_SOURCE, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
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
	){
		return send_receive_replace_n(first, size, dest, source, sendtag, recvtag);
	}
};

}}}

#if not __INCLUDE_LEVEL__ // def _TEST_MPI3_DETAIL_BASIC_COMMUNICATOR

#include "../../mpi3/version.hpp"

int main(int argc, char* argv[]){

    int rank, nprocs;

    MPI_Init(&argc,&argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    MPI_Finalize();
    return 0; 
}

#endif
#endif


