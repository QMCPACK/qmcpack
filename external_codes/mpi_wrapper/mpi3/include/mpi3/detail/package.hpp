#ifndef MPI3_DETAIL_PACKAGE_HPP
#define MPI3_DETAIL_PACKAGE_HPP

// TODO(correaa) move from detail to root

#include "../../mpi3/vector.hpp"

#include "../../mpi3/detail/basic_communicator.hpp"
#include "../../mpi3/detail/buffer.hpp"
#include "../../mpi3/detail/iterator.hpp"
#include "../../mpi3/detail/iterator_traits.hpp"
#include "../../mpi3/detail/value_traits.hpp"

namespace boost {
namespace mpi3 {

class communicator;

namespace detail{

struct package : buffer {
 private:
	basic_communicator& bcomm_;

 public:
	explicit package(communicator& comm, buffer::size_type n = 0)
	: buffer{n}, bcomm_{reinterpret_cast<basic_communicator&>(comm)} {  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa) break cyclic dependency of classes
		reserve(n);
	}
	package(package const&) = delete;
	package(package     &&) = delete;

	package& operator=(package const&) = delete;
	package& operator=(package     &&) = delete;

	~package() noexcept = default;

	template<class It, typename Size>
	void pack_n(It first, Size count){
		bcomm_.pack_n(first, count, static_cast<buffer&>(*this));
	}
	template<class It>
	auto pack(It first, It last){
		bcomm_.pack(first, last, static_cast<buffer&>(*this));
	}
	template<class It, class Size>
	void unpack_n(It first, Size count){
		bcomm_.unpack_n(static_cast<buffer&>(*this), first, count);
	}
	template<class It>
	void unpack(It first, It last){
		bcomm_.unpack(static_cast<buffer&>(*this), first, last);
	}
	explicit operator bool() const {return pos < static_cast<std::ptrdiff_t>(size());}

	template<class T>
	package& operator>>(T& t){
		unpack_n(std::addressof(t), 1);
		return *this;
	}
	auto send(int dest, int tag = 0) {
		return bcomm_.send(static_cast<buffer&>(*this), dest, tag);
	}
	auto receive(int dest, int tag = 0) {
		return bcomm_.receive(static_cast<buffer&>(*this), dest, tag);
	}

//	package const& send(int dest, int tag = 0) const;
/*	package const& send(int dest, int tag = 0) const{
		comm_.send_packed_n(buffer_.data(), in_, dest, tag);
		return *this;
	}*/
//	package& receive(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);
/*	package& receive(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		MPI_Status status;
		MPI_Message msg;
		int count = -1;
		MPI_Mprobe(source, tag, comm_.impl_, &msg, &status);
		MPI_Get_count(&status, MPI_PACKED, &count);
		buffer_.resize(count);
		MPI_Mrecv(buffer_.data(), count, MPI_PACKED, &msg, MPI_STATUS_IGNORE);
	//	int n = comm_.probe(source, tag).count<char>();
	//	buffer_.resize(n);
	//	comm_.receive_packed_n(buffer_.data(), n, source, tag);
		return *this;
	}*/
/*	package& broadcast(int root = 0){ // see https://www.researchgate.net/publication/228737912_Dynamically-Sized_Messages_in_MPI-3
		comm_.broadcast_value(in_, root);
		buffer_.resize(in_);
		comm_.broadcast_n(buffer_.data(), in_, root);
		return *this;
	}*/
//	package& gather(int root = 0){
//		
//	}

//	template<class T> int size(int n = 1) const{return comm_.pack_size<T>(n);}

};

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_PACKAGE

//#include "../../mpi3/communicator.hpp"
//#include "../../mpi3/main.hpp"

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int argc, char* argv[], mpi3::communicator world){

//	if(world.rank() == 0){
//		char buf[1000];
//		int i = 12;
//		int j = 13;
//		auto end = world.pack_n(&i, 1, buf);
//		     end = world.pack_n(&j, 1, end);
//		world.send_packed(buf, end, 1); //world.send_packed_n(buff, std::distance(buff, end), 1); //world.send_packed_n(buff, end - buff, 1);
//		world.send_packed_n(buf, 1000, 2);
//	}else if(world.rank() == 1){
//		std::vector<int> v(2);
//		world.receive(v.begin(), v.end(), 0);
//		assert(v[0] == 12);
//		assert(v[1] == 13);
//	}else if(world.rank() == 2){
//		char buf[1000];
//		world.receive_packed_n(buf, 1000, 0);
//		int i = -1;
//		int j = -1;
//		auto end = world.unpack_n(&i, 1, buf);
//		     end = world.unpack_n(&j, 1, end);
//		assert(i == 12);
//		assert(j == 13);
//	}
//	world.barrier();
////	return 0;

//	if(world.rank() == 0){
//		mpi3::package p(world);
//		int i = 12;
//		int j = 13;
//		(p << i << j).send(1).send(2);
//	//	p.send(1);
//	//	p.send(2);
//	}else if(world.rank() == 1){
//		std::vector<int> v(2, -1);
//		world.receive(v.begin(), v.end(), 0);
//		assert(v[0] = 12);
//		assert(v[1] == 13);
//	}else if(world.rank() == 2){
//		mpi3::package p(world);
//		int i = -1;
//		int j = -1;
//		p.receive(0) 
//			>> i 
//			>> j
//		;
//		assert(i == 12);
//		assert(j == 13);
//	}
//}

//#endif
#endif
