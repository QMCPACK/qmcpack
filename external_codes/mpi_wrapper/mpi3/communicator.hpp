#if COMPILATION // -*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
mpic++ -x c++ $0 -o $0x&&mpirun -n 1 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MPI3_COMMUNICATOR_HPP
#define MPI3_COMMUNICATOR_HPP

#include "../mpi3/communication_mode.hpp"
#include "../mpi3/generalized_request.hpp"
#include "../mpi3/handle.hpp"
#include "../mpi3/info.hpp"
#include "../mpi3/message.hpp"
#include "../mpi3/operation.hpp"
#include "../mpi3/port.hpp"
#include "../mpi3/request.hpp"
#include "../mpi3/status.hpp"
#include "../mpi3/type.hpp"
#include "../mpi3/error.hpp"
#include "../mpi3/group.hpp"
//#include "../mpi3/window.hpp"

#include "../mpi3/detail/equality.hpp"
#include "../mpi3/detail/basic_communicator.hpp"
#include "../mpi3/detail/buffer.hpp"
#include "../mpi3/detail/datatype.hpp"
#include "../mpi3/detail/iterator.hpp"
#include "../mpi3/detail/value_traits.hpp"
//#include "../mpi3/detail/strided.hpp"
#include "../mpi3/detail/package.hpp"

#include "../mpi3/config/NODISCARD.hpp"

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include <boost/optional.hpp> // TODO replace by std::optional (c++17)
//#include <boost/static_assert.hpp>
//#include <boost/type_traits/is_array.hpp> 

#define BOOST_PACKAGE_ARCHIVE_SOURCE

#include <boost/archive/detail/common_oarchive.hpp>
#include <boost/archive/detail/common_iarchive.hpp>
#include <boost/archive/basic_streambuf_locale_saver.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/archive/detail/auto_link_archive.hpp>
//#include <boost/archive/detail/abi_prefix.hpp> // must be the last header

#include <boost/serialization/item_version_type.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

#include <boost/mpl/placeholders.hpp>

#include<boost/any.hpp> // or <any> in C++17

// use this to avoid need for linking -lserialization
#ifdef _MAKE_BOOST_SERIALIZATION_HEADER_ONLY
//#include <boost/archive/detail/decl.hpp>
#if BOOST_VERSION > 106000 && BOOST_VERSION < 106600
#include "../mpi3/serialization_hack/singleton.cpp"
#endif
#if BOOST_VERSION < 105900
#define BOOST_ARCHIVE_DECL
#define BOOST_SERIALIZATION_DECL
#endif
#include "../mpi3/serialization_hack/archive_exception.cpp"
#include "../mpi3/serialization_hack/extended_type_info.cpp"
#include "../mpi3/serialization_hack/extended_type_info_typeid.cpp"
#include "../mpi3/serialization_hack/basic_archive.cpp"
#include "../mpi3/serialization_hack/basic_iarchive.cpp"
#include "../mpi3/serialization_hack/basic_oarchive.cpp"
#include "../mpi3/serialization_hack/basic_iserializer.cpp"
#include "../mpi3/serialization_hack/basic_oserializer.cpp"
#endif

#include "../mpi3/package_archive.hpp"

#include<cassert>
#include<iostream>
#include<iterator> // iterator_traits
#include<limits>
#include<map>
#include<numeric> // std::accumulate
#include<string>
#include<thread>
#include<type_traits>
#include<vector>
#include<type_traits> // is_same

namespace boost{
namespace mpi3{

inline bool check_mpi_(enum error e){
	if(e!=mpi3::error::success) return true;
	throw std::system_error{e, "cannot call function"};
	return false;
}

#define SAFE_MPI_(F) check_mpi_(MPI_##F)

#if !defined(OPEN_MPI) || (OMPI_MAJOR_VERSION < 2)
#define OMPI_COMM_TYPE_NODE     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_HWTHREAD MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CORE     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L1CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L2CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L3CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_SOCKET   MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_NUMA     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_BOARD    MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_HOST     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CU       MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CLUSTER  MPI_COMM_TYPE_SHARED
#endif

// https://www.open-mpi.org/doc/v4.0/man3/MPI_Comm_split_type.3.php#toc8
enum class communicator_type : int{
	shared    = MPI_COMM_TYPE_SHARED   ,/*synomym*/ node = OMPI_COMM_TYPE_NODE,
	hw_thread = OMPI_COMM_TYPE_HWTHREAD,
	core      = OMPI_COMM_TYPE_CORE    ,
	l1_cache  = OMPI_COMM_TYPE_L1CACHE ,
	l2_cache  = OMPI_COMM_TYPE_L2CACHE ,
	l3_cache  = OMPI_COMM_TYPE_L3CACHE ,
	socket    = OMPI_COMM_TYPE_SOCKET  ,
	numa      = OMPI_COMM_TYPE_NUMA    ,
	board     = OMPI_COMM_TYPE_BOARD   ,
	host      = OMPI_COMM_TYPE_HOST    ,
	cu        = OMPI_COMM_TYPE_CU      ,/*synomym*/ cpu = OMPI_COMM_TYPE_CU   ,
	cluster   = OMPI_COMM_TYPE_CLUSTER 
};

enum constant{
	undefined    = MPI_UNDEFINED ,
	process_null = MPI_PROC_NULL ,
	any_source   = MPI_ANY_SOURCE
};

enum key{ // for attributes
	tag_ub             = MPI_TAG_UB, 
	host               = MPI_HOST,
	io                 = MPI_IO,
	wtime_is_global    = MPI_WTIME_IS_GLOBAL,
	application_number = MPI_APPNUM,
	universe_size      = MPI_UNIVERSE_SIZE,
	last_used_code     = MPI_LASTUSEDCODE
};

template<int N = 10> struct overload_priority : overload_priority<N-1>{
//	using overload_priority<N-1>::overload_priority;
};
template<> struct overload_priority<0>{};

class environment;
class group;

using boost::optional;

struct error_handler;

template<class T>
struct shm_pointer;

template<class T>
struct pointer;

using address = MPI_Aint;
using intptr_t = MPI_Aint;
using size_t = MPI_Aint;
using ptrdiff_t = std::make_signed_t<size_t>;
struct request;

struct send_request;
struct receive_request;

struct FILE;

struct process;

struct ostream;

struct message_header{
	int tag;
	
};

struct graph_communicator;
struct shared_communicator; // intracommunicator

using boost::any;
using boost::any_cast;

//class communicator_ptr{};

template<class T> struct window;

class communicator : protected detail::basic_communicator{
	friend struct detail::package;
	friend struct window<void>;
protected:
	bool is_null() const{return MPI_COMM_NULL == impl_;}
	friend class mpi3::environment;
	detail::equality compare(communicator const& other) const{
		detail::equality ret;// = boost::mpi3::detail::unequal;
		MPI_(Comm_compare)(impl_, other.impl_, reinterpret_cast<int*>(&ret));
		return ret;
	}
public:
	communicator(communicator const& o, group const& g);
	communicator(group const& g);
	communicator(group const& g, int tag);
	explicit operator group() const;

	communicator duplicate(){
		communicator ret;
		MPI_(Comm_dup)(impl_, &ret.impl_);
		return ret;
	}

	template<class T = void*>
	struct keyval{
		static int delete_fn_(MPI_Comm /*comm*/, int /*keyval*/, void *attr_val, void */*extra_state*/){
			delete (T*)attr_val; attr_val = nullptr;
			return MPI_SUCCESS;
		}
		static int copy_fn_(
			MPI_Comm /*oldcomm*/, int /*keyval*/,
			void * /*extra_state*/, void *attribute_val_in,
			void *attribute_val_out, int *flag
		){
			*(void**)attribute_val_out = (void*)new T{*((T const*)attribute_val_in)};
			assert(flag); *flag = 1;
			return MPI_SUCCESS;
		}
		using mapped_type = T;
		int impl_;
		keyval(){
			MPI_Comm_create_keyval(copy_fn_, delete_fn_, &impl_, (void *)0);
		}
		keyval(keyval const&) = delete;
		~keyval(){MPI_Comm_free_keyval(&impl_);}
	};
	using detail::basic_communicator::send_receive_n;
	using detail::basic_communicator::matched_probe;
	template<class It, typename Size>
	auto send_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int dest, int tag
	){
		auto e = static_cast<enum error>(
			MPI_Send(
				detail::data(first), count, 
				detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
				dest, tag, impl_
			)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot send"};
	}
	template<class It, typename Size>
	auto isend_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int dest, int tag
	){
		mpi3::request r;
		int s = MPI_Isend(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, tag, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send");
		return r;
	}
	template<class It, typename Size>
	void send_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count,
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(count--) poa << *first++;
		send_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}
	template<class It, typename Size>
	auto isend_n(It first, Size count, int dest, int tag = 0){
		return isend_n(
			first, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			dest, tag
		);
	}
	template<class It, typename Size>
	auto isend_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(count--) poa << *first++;
		return isend_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}
	template<class T, class = decltype(T::dimensionality)> static std::true_type  has_dimensionality_aux(T const&);
	                                                       static std::false_type has_dimensionality_aux(...);
	template<class T> struct has_dimensionality : decltype(has_dimensionality_aux(T{})){};

	template<class It, typename Size, class = std::enable_if_t<(not has_dimensionality<It>{})> >
	void send_n(It first, Size count, int dest, int tag = 0){
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
			detail::contiguous_iterator_tag,
			detail::basic_tag,
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
		mpi3::vector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		return send_n(buffer.begin(), buffer.size(), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(first!=last) poa << *first++;
		send_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}

	template<class MultiIt>
	auto send_n(MultiIt first, typename MultiIt::difference_type count, int dest, int tag = 0)
	->decltype( MPI_Send (mpi3::base(first), count, mpi3::type{first}, dest, tag, impl_), first + count){
		return MPI_(Send)(mpi3::base(first), count, mpi3::type{first}, dest, tag, impl_), first + count;}

	template<class MultiIt>
	auto receive_n(MultiIt first, typename MultiIt::difference_type count, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG)
	->decltype( MPI_Recv (mpi3::base(first), count, mpi3::type{first}, source, tag, impl_, MPI_STATUS_IGNORE), first + count){
		return MPI_(Recv)(mpi3::base(first), count, mpi3::type{first}, source, tag, impl_, MPI_STATUS_IGNORE), first + count;}

#if 0
	template<
		class MultiIt, class Size = typename MultiIt::difference_type, class Stride = typename MultiIt::stride_type,
			std::size_t D = MultiIt::dimensionality, std::enable_if_t<(D>=2), int> = 0,
			typename Element = typename MultiIt::element, typename DataType = detail::basic_datatype<Element>,
			std::enable_if_t<detail::is_basic<Element>{}, int> = 0
	>
	MultiIt receive_n(MultiIt first, decltype(Size{}) count, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		status sta;
		auto e = static_cast<enum error>(
			MPI_Recv(
				mpi3::base(first), count, mpi3::type{first}.commit(),
				source, tag, impl_, &sta.impl_
			)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot receive"};
		return first + count;
	}
#endif

	template<class It>
	auto isend(
		It first, It last,
			detail::random_access_iterator_tag, 
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return isend_n(first, std::distance(first, last), dest, tag);
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
	template<class It>
	auto isend(It first, It last, int dest, int tag = 0){
		return isend(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);
	}
	using detail::basic_communicator::basic_communicator;
	communicator(communicator const&) = default;
	communicator(communicator&) = default;
	communicator(communicator&&) = default;
	communicator() = default;
/*	communicator& operator=(communicator const& other){
		communicator tmp(other);
		swap(tmp);
		return *this;
	}*/
//	communicator& operator=(communicator&& other){
//		communicator tmp(std::move(other));
//		swap(tmp);
//		return *this;
//	}
	communicator& operator=(communicator other){return swap(other), *this;}
	bool operator==(communicator const& other) const{
		return &*this==&other or compare(other)==detail::equality::congruent;
	//	auto eq = compare(other);
	//	return (eq == equality::identical) or (eq == equality::congruent);
	}
	bool operator!=(communicator const& other) const{return not(*this==other);}
	explicit operator bool() const{return not is_null();}
//	impl_t operator&() const{return impl_;}
	auto get() const{return impl_;}
	~communicator(){
		if(impl_ != MPI_COMM_WORLD and impl_ != MPI_COMM_NULL and impl_ != MPI_COMM_SELF){
			MPI_Comm_disconnect(&impl_); //this will wait for communications to finish communications, <s>if it gets to this point is probably an error anyway</s> <-- not true, it is necessary to synchronize the flow
		//	MPI_Comm_free(&impl_);
		}
	}
	int size() const{
		if(is_null()) return 0;//throw std::runtime_error("size called on null communicator");
		int size = -1; 
		int s = MPI_Comm_size(impl_, &size);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot get size"); 
		return size;
	}
	bool empty() const{return size() == 0;}
	void abort(int errorcode = 0) const{MPI_Abort(impl_, errorcode);}
	bool is_intercommunicator() const{
		int flag = false;
		MPI_Comm_test_inter(impl_, &flag);
		return flag;
	}
	communicator split(int color, int key) const{
		communicator ret;
		int s = MPI_Comm_split(impl_, color, key, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot split communicator");
		if(ret) ret.name(name() + std::to_string(color));// + std::to_string(key));
		if(ret) ret.attribute("color") = color;
		return ret;
	}
	communicator split(int color = MPI_UNDEFINED) const{
		return split(color, rank());
	}
	communicator keep(bool cond) const{return split(cond?0:mpi3::undefined);}

	shared_communicator split_shared(int key = 0) const;
	shared_communicator split_shared(communicator_type t, int key = 0) const;
	int remote_size() const{
		int ret = -1;
		int s = MPI_Comm_remote_size(impl_, &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot remote size");
		return ret;
	}
	communicator reversed() const{
		return split(0, size() - rank());
	}
	int cartesian_map(std::vector<int> const& dims, std::vector<int> const& periods) const{
		assert( dims.size() == periods.size() );
		int ret;
		int s = MPI_Cart_map(impl_, dims.size(), dims.data(), periods.data(), &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot map"};
		return ret;
	}
	int cartesian_map(std::vector<int> const& dimensions) const{
		return cartesian_map(dimensions, std::vector<int>(dimensions.size(), 0));
	}
	pointer<void> malloc(MPI_Aint size) const;
	template<class T = void> void deallocate_shared(pointer<T> p);
	template<class T = void> void deallocate(pointer<T>& p, MPI_Aint size = 0);
	void free(pointer<void>& p) const;

	bool similar(communicator const& o) const{return compare(o)==detail::equality::similar;}
	template<class Vector>//, typename = typename std::enable_if<std::is_same<decltype(Vector{}.data()), int*>{}>::type>
	communicator subcomm(Vector const& v) const{
		MPI_Group old_g;
		MPI_Comm_group(impl_, &old_g);
		MPI_Group new_g;
		MPI_Group_incl(old_g, v.size(), v.data(), &new_g);
		communicator ret; MPI_Comm_create(impl_, new_g, &ret.impl_);
		MPI_Group_free(&new_g);
		MPI_Group_free(&old_g);
		return ret;
	}
	communicator subcomm(std::initializer_list<int> l) const{
		return subcomm(std::vector<int>(l));
	}
	enum class topology{undefined = MPI_UNDEFINED, graph = MPI_GRAPH, cartesian = MPI_CART};
//	topology topo() const{return static_cast<topology>(call<&MPI_Topo_test>());}
	int rank() const{
		int rank = -1;
		if(impl_ == MPI_COMM_NULL) throw std::runtime_error("rank on null communicator");
		int s = MPI_Comm_rank(impl_, &rank);
		if(s != MPI_SUCCESS) MPI_Comm_call_errhandler(impl_, s);
		return rank;
	}
	int right() const{return (rank() + 1) % size();}
	int left() const{
		int left = rank() - 1;
		if(left < 0) left = size() - 1;
		return left;
	}
	int next(int n = 1) const{
		assert(rank() + n < size());
		return rank() + n;
	}
	int prev(int n = 1) const{
		assert(rank() - n > 0);
		return rank() - n;
	}
	communicator accept(port const& p, int root = 0) const{
		communicator ret;
		MPI_Comm_accept(p.name_.c_str(), MPI_INFO_NULL, root, impl_, &ret.impl_);
		return ret;
	}
	void barrier() const{MPI3_CALL(MPI_Barrier)(impl_);}
	communicator connect(port const& p, int root = 0) const{
		communicator ret;
		MPI_Comm_connect(p.name_.c_str(), MPI_INFO_NULL, root, impl_, &ret.impl_);
		return ret;
	}
	bool    root() const{return (not empty()) and (rank() == 0);}
	bool is_root() const{return root();}

	void set_error_handler(error_handler const& eh);
	error_handler get_error_handler() const;
//	template<typename T>
//	void set_attribute(keyval const& kv, int idx, T const& val);
//	void delete_attribute(keyval const& kv, int idx);
//	template<typename T>
//	void get_attribute(keyval const& kv, int idx, T& val);
//	template<typename T>
//	T const& get_attribute_as(keyval const& kv, int idx);
//	bool has_attribute(keyval const& kv, int idx);

	process operator[](int i);
protected:
	template<class T> void set_attribute(int kv_idx, T const& t){
		MPI_Comm_set_attr(impl_, kv_idx, new T{t});
	}
	inline void delete_attribute(int kv_idx){
		MPI_Comm_delete_attr(impl_, kv_idx);
	}
	void* get_attribute(int kvidx) const{
		void* v = nullptr; int flag = 0;
		MPI_Comm_get_attr(impl_, kvidx, &v, &flag);
		if(not flag){assert(!v); return nullptr;}
		return v;
	}
	bool has_attribute(int kvidx) const{
		void* v = nullptr; int flag = 0;
		MPI_Comm_get_attr(impl_, kvidx, &v, &flag);
		if(not flag) return false;
		return true;
	}
public:
	template<class T, class TT = T> void
	set_attribute(keyval<T> const& k, TT const& t = {}){set_attribute<T>(k.impl_, t);}
	template<class T>
	inline void delete_attribute(keyval<T> const& k){delete_attribute(k.impl_);}
	template<class T>
	T const& get_attribute(keyval<T> const& kv) const{return *(T*)get_attribute(kv.impl_);}
	template<class T> 
	T& get_attribute(keyval<T> const& kv){return *(T*)get_attribute(kv.impl_);}
	template<class T>
	bool has_attribute(keyval<T> const& kv) const{return has_attribute(kv.impl_);}
	template<class T>
	T& attribute(keyval<T> const& kv){
		if(not has_attribute(kv)) set_attribute(kv);
		return get_attribute(kv);
	}
	mpi3::any& attribute(std::string const& s);

	void call_error_handler(int errorcode){
		int status = MPI_Comm_call_errhandler(impl_, errorcode);
		if(status != MPI_SUCCESS) throw std::runtime_error{"cannot call error handler"};
	}
	void error(mpi3::error const& e){
		int s = MPI_Comm_call_errhandler(impl_, static_cast<int>(e));
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot call error handler"};
	}
	communicator divide_low(int n) const{
		return split(
			(rank() < size()/n*(n-size()%n))?
				rank()/(size()/n):
				n-size()%n + (rank() - (n-size()%n)*(size()/n))/((size()/n)+1)
		);
	}
	communicator divide_high(int n) const{
		int bat=size()/n; int residue = size()%n;
		int i = 0;
		for(int last = 0; ; i++){
			last += bat + ((i < residue)?1:0);
			if(rank() < last) break;
		}
		return split(i);
	}
	communicator operator/(int n) const{
		assert(n!=0);
		if(n > 0) return divide_high(n);
		return divide_low(n);
	}
	communicator operator%(int n) const{return split(rank()%n);}
	communicator divide_even(int n) const{
		return split(2*(rank()%n) > n?mpi3::undefined:rank()/n);
	}
	communicator operator/(double nn) const{return divide_even(nn);}
	communicator operator<(int n) const{return split((rank() < n)?0:MPI_UNDEFINED);}
	communicator operator<=(int n) const{return split((rank() <= n)?0:MPI_UNDEFINED);}
	communicator operator>(int n) const{return split((rank() > n)?0:MPI_UNDEFINED);}
	communicator operator>=(int n) const{return split((rank() >= n)?0:MPI_UNDEFINED);}
	communicator operator==(int n) const{return split((rank() == n)?0:MPI_UNDEFINED);}

	template<class T>
	void send_value(T const& t, int dest, int tag = 0){
		send(std::addressof(t), std::addressof(t) + 1, dest, tag);
	}
	template<class T>
	auto isend_value(T const& t, int dest, int tag = 0){
		return isend(std::addressof(t), std::addressof(t) + 1, dest, tag);
	}
	template<class T, std::size_t N>
	void send_value(T(&t)[N], int dest, int tag = 0){
		send(std::addressof(t[0]), std::addressof(t[N-1]) + 1, dest, tag);
	}
	template<class ContIt, class Size>
	request send_init_n(ContIt first, Size count, int dest, int tag = 0);
	template<class ContIt>
	request send_init(ContIt first, ContIt last, int dest, int tag = 0);

	template<class ContIt, class Size>
	receive_request receive_init_n(ContIt first, Size count, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);
	template<class ContIt>
	receive_request receive_init(ContIt first, ContIt last, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);

	template<typename InputIt, typename Size, 
		class value_type = typename std::iterator_traits<InputIt>::value_type, 
		class datatype = detail::basic_datatype<value_type>
	>
	auto pack_n_aux(std::random_access_iterator_tag, InputIt first, Size count, char* out, int max_size
	){
		assert(0);
		int position = 0;
		MPI_Pack(detail::data(first), count, datatype{}, out, max_size, &position, impl_);
		return out + position;
	}

	template<class InputIt, class Size, typename V = typename std::iterator_traits<InputIt>::value_type, class datatype = detail::basic_datatype<V>>
	auto unpack_from_n(InputIt first, Size count, char const* buffer){
		int position = 0;
		using detail::data;
	//	assert( count % pack_size<V>() == 0 );
		int s = MPI_Unpack(buffer, count*pack_size<V>(), &position, data(first), count, datatype{}, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot packaga unpack");
		return buffer + position;
	}
//	template<class InputIt, typename Size>
//	auto unpack_n(InputIt first, Size count, char const* buffer){
//		return unpack_from_n(first, count, buffer);
//	}
//	template<class V>
//	auto unpack_value(V& v, char const* buffer){
//		return unpack_n(&v, 1, buffer);
//	}
//	template<class V>
//	auto pack_value(V const& v, char const* buffer){
//		return pack_n(&v, 1, buffer);
//	}
	template<class It>
	auto pack(
		It first, It last, 
			detail::random_access_iterator_tag,
		uvector<detail::packed>& b, int pos
	){
		return pack_n(first, std::distance(first, last), b, pos);
	}
	template<class It>
	auto pack(It first, It last, uvector<detail::packed>& b, int pos){
		return pack(
			first, last,
				detail::iterator_category_t<It>{},
			b, pos
		);
	}
	template<class It, class Size>
	auto send_receive_replace_n(
		It first, Size size, 
		int dest, int source, // = MPI_ANY_SOURCE, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		return send_receive_replace_n(
			first, 
				detail::iterator_category_t<It>{}, 
				detail::value_category_t<value_type>{},
			size,
			dest, source, sendtag, recvtag
		);
	}
	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::random_access_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int source, int sendtag, int recvtag
	){
		status ret;
		int s = MPI_Sendrecv_replace(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{}, 
			dest, sendtag, source, recvtag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot sendrecv_replace"};
		return ret;
	}
	template<class It, typename Size>
	auto send_receive_n(
		It first, Size count, int dest, 
		It d_first, Size d_count, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive_n(
			first, count, dest, 
			d_first, d_count, source,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{}, 
			sendtag, recvtag
		);
	}
	template<class It, typename Size>
	auto send_receive_n(
		It first, Size count, int dest, 
		It d_first, Size d_count, int source,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int sendtag, int recvtag
	){
		status ret;
		int s = MPI_Sendrecv(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, sendtag, 
			detail::data(d_first), d_count,
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			source, recvtag,
			impl_,
			&ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot sendrecv"};
		return ret;
	}
	template<class It, typename Size, typename... Meta>
	auto send_receive_replace_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int dest, int source,
		int sendtag, int recvtag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		auto first_copy = first;
		while(count--) poa << *first_copy++;
		auto s = p.size();
		send_receive_replace_n(&s, 1, dest, source, sendtag, recvtag);
		detail::package p2(*this);
		p2.resize(s);
		send_receive_n(
			p.begin(), p.size(), dest, 
			p2.begin(), p2.size(), 
			source, sendtag, recvtag
		);
		package_iarchive pia(p2);
		while(p2) pia >> *first++;
		return first;
	}
	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int source, int sendtag, int recvtag
	){
		uvector<typename std::iterator_traits<It>::value_type> v(count);
		std::copy_n(first, count, v.begin());
		send_receive_replace_n(v.begin(), v.size(), dest, source, sendtag, recvtag);
		return std::copy_n(v.begin(), v.size(), first);
	}
	template<class It, class Size>
	auto send_receive_n(
		It first, Size size, 
		int dest, int source, // = MPI_ANY_SOURCE, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive_replace_n(
			first, size,
			dest, source, sendtag, recvtag
		);
	}
public:
	template<class It>
	auto send_receive(
		It first, It last, int dest, 
		It d_first, It d_last, int source, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int sendtag, int recvtag
	){
		return send_receive_n(
			first, std::distance(first, last), dest,
			d_first, std::distance(d_first, d_last), source,
			sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive(
		It first, It last, int dest, 
		It d_first, It d_last, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive(
			first, last, dest,
			d_first, d_last, source,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{}, 
			sendtag, recvtag
		);
	}
/*
		int size = comm_.pack_size<V>()*count;
		auto old_buffer_size = buffer_.size();
		buffer_.resize(old_buffer_size + size);
		auto curr = buffer_.data() + old_buffer_size;
		int position = -1;
		int outcount = this->pack_size<V>()*count;
		using detail::data;
//		auto end = comm_.pack_n(data(first), count, curr);
		MPI_Pack(detail::data(first), count, datatype{}, out, outcount, &position, impl_);
		end = out + position;
		assert(end == std::addressof(*buffer_.end()));
		in_ = buffer_.size();
	//	in_ += end - curr;
	}
*/

public:
	status probe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		status s;
		MPI_Probe(source, tag, impl_, &s.impl_);
		return s;
	}
	optional<status> iprobe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		status s;
		int flag;
		MPI_Iprobe(source, tag, impl_, &flag, &s.impl_);
		if(not flag) return optional<status>();
		return optional<status>(s);
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int dest, int source, int sendtag, int recvtag
	){
		return send_receive_replace_n(
			first, std::distance(first, last), 
			dest, source, sendtag, recvtag
		);		
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int dest, int source, int sendtag, int recvtag
	){
		mpi3::vector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		send_receive_replace_n(buffer.begin(), buffer.size(), dest, source, sendtag, recvtag);
		return std::copy(buffer.begin(), buffer.end(), first);
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int source, int sendtag, int recvtag
	){
		return send_receive_replace_n(
			first, std::distance(first, last), 
			dest, source, sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive_replace(It first, It last, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG){
		return send_receive_replace(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, source, sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive(It first, It last, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG){
		return send_receive_replace(first, last, dest, source, sendtag, recvtag);
	}
/*	template<class It, class... Args>
	auto send_receive_replace_category(std::random_access_iterator_tag const&, 
		It first, It last, Args&&... args){
		return send_receive_replace_n(first, std::distance(first, last), std::forward<Args>(args)...);
	}*/
/*
	template<class T1, class T2>
	auto send_receive_value(
		T1 const& t1, int dest, 
		T2&       t2, int source = MPI_ANY_SOURCE
	){
		return send_receive_n(&t1, 1, dest, &t2, 1, source);
	}
	template<
		class IOIterator, class Size, class... A, 
		class CategoryIO = typename std::iterator_traits<IOIterator>::iterator_category,
		class datatype = detail::basic_datatype<typename std::iterator_traits<IOIterator>::value_type>
	>
	auto send_receive_n_aux(
		std::random_access_iterator_tag, IOIterator IO, Size count, int dest, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_SOURCE
	) const{
		status ret;
		MPI_Sendrecv_replace(std::addressof(*IO), count, datatype{}, dest, sendtag, source, recvtag, impl_, &ret.impl_);
		return ret;
	}
*/
/*	template<
		class IOIterator, class Size, class... A, 
		class CategoryIO = typename std::iterator_traits<IOIterator>::iterator_category
	>
	auto send_receive_n(IOIterator IO, Size count, int dest, int source, A&&... a) const{
		return send_receive_n_aux(CategoryIO(), IO, count, dest, source, std::forward<A>(a)...);
	}*/

//	template<class T>
//	auto send_receive_value(T& t, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG) const{
//		return send_receive_n(&t, 1, dest, source, sendtag, recvtag);
//	}

	void send_packed_n(void const* begin, int n, int dest, int tag = 0){
		std::cout << "sending packet of size " << n << std::endl;
		MPI_Send(const_cast<void*>(begin), n, MPI_PACKED, dest, tag, impl_);
	}
//	void send_value(package const& p, int dest, int tag = 0);
//	template<class Char = char>
//	auto send_packed(Char const* begin, Char const* end, int dest, int tag = 0){
//		return send_packed_n(begin, (char*)end - (char*)begin, dest, tag);
//	}
	auto receive_packed_n(void* begin, int n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const{
		status ret;
		MPI_Recv(begin, n, MPI_PACKED, source, tag, impl_, &ret.impl_);
		return ret;
	}
	auto receive_packed(void* begin, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		MPI_Status status;
		MPI_Message msg;
		int count = -1;
		MPI_Mprobe(source, tag, impl_, &msg, &status);
		MPI_Get_count(&status, MPI_PACKED, &count);
		MPI_Mrecv(begin, count, MPI_PACKED, &msg, MPI_STATUS_IGNORE);
	//	auto n = probe(source, tag).count<char>();
	//	receive_packed_n(begin, n, source, tag);
		return (void*)((char*)begin + count);
	}
	template<class It, typename Size>
	auto receive_n(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		int source, int tag
	){
		status sta;
		auto e = static_cast<enum error>(
			MPI_Recv(
				detail::data(dest), count, 
				detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
				source, tag, impl_, &sta.impl_
			)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot receive"};
		return dest + count;
	}
	template<class It, typename Size>
	auto ireceive_n(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		int source, int tag
	){
		mpi3::request r;
		int s = MPI_Irecv(
			detail::data(dest), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			source, tag, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("receive_n");
		return r;
	//	return dest + count;
	}
	template<class It, typename Size>
	auto receive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}
/*	template<class It, typename Size>
	auto ireceive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}*/
/*	template<class It, typename Size>
	receive_request ireceive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}*/
	template<class It, typename Size
		, std::enable_if_t<not has_dimensionality<It>{}, int> =0// or (not detail::is_basic<typename std::iterator_traits<It>::value_type>{}), int> =0 // needed by intel commpiler
	>
	auto receive_n(It dest, Size n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive_n(
			dest, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			n,
			source, tag
		);
	}
	template<class It, typename Size>
	mpi3::request ireceive_n(
		It dest, Size n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG
	){
		return ireceive_n(
			dest,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			n,
			source, tag
		);
	}
	template<class It>
	auto receive(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		match m = matched_probe(source, tag);
		auto count = m.count<typename std::iterator_traits<It>::value_type>();
		m.receive_n(dest, count);
		return dest + count;
	}
	template<class It>
	auto receive(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(p) pia >> *dest++;
		return dest;
	}
	template<class It>
	auto receive(
		It dest, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
	//	match m = matched_probe(source, tag);
	//	auto count = m.count<typename std::iterator_traits<It>::value_type>();
	//	mpi3::uvector<typename std::iterator_traits<It>::value_type> v(count);
		return matched_probe(source, tag).receive_n(dest);//, count);
	}
	template<class It>
	auto receive(It dest, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(
			dest,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			source, tag
		);
	}
	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		return receive_n(d_first, std::distance(d_first, d_last), source, tag);
	}
	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		return receive_n(std::addressof(*d_first), std::distance(d_first, d_last), source, tag);
	//	return std::copy(buffer.begin(), buffer.end(), d_first);
	}

	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(std::distance(d_first, d_last));
		receive_n(buffer.begin(), buffer.size(), source, tag);
		return std::copy(buffer.begin(), buffer.end(), d_first);
	}
	class ir_req{
		boost::mpi3::status query(){
			boost::mpi3::status ret;
			ret.set_source(MPI_UNDEFINED);
			ret.set_tag(MPI_UNDEFINED);
			ret.set_cancelled();
			ret.set_elements<char>(0);
			return ret;
		}
		void free(){
			std::cout << "free" << std::endl;
		}
		void cancel(int complete){
			std::cout << "cancel " << complete << std::endl;
		}
	};
	template<class It>
	struct receive_args{ 
		communicator* commP; 
		It d_first;
	//	It d_last;
		int source;
		int tag; 
		MPI_Request* requestP; 
	};
	struct receive_state{
		int cancelled = 0;
		int source = MPI_UNDEFINED;
		int tag = MPI_UNDEFINED;
	};
	template<class It>
	inline static void* receive_thread(void* ptr){
		receive_args<It>* args = (receive_args<It>*)ptr;
		args->commP->receive(args->d_first, args->source, args->tag);//, /*args->d_last,*/ );
		MPI_Grequest_complete(*args->requestP);
		::free(ptr);
		return NULL;
	}
	inline static int query_fn(void* extra_state, MPI_Status *status){
		receive_state* rs = (receive_state*)extra_state;
		/* always send just one int */ 
		MPI_Status_set_elements(status, MPI_INT, 1);
		/* can never cancel so always true */ 
		MPI_Status_set_cancelled(status, rs->cancelled); 
		/* choose not to return a value for this */
		status->MPI_SOURCE = rs->source; 
		/* tag has not meaning for this generalized request */ 
		status->MPI_TAG = rs->tag; 
		/* this generalized request never fails */ 
		return MPI_SUCCESS; 
	}
	inline static int free_fn(void* extra_state){ 
		/* this generalized request does not need to do any freeing */ 
		/* as a result it never fails here */
		::free(extra_state);
		return MPI_SUCCESS; 
	}
	inline static int cancel_fn(void* /*extra_state*/, int complete) 
	{ 
		/* This generalized request does not support cancelling. 
		   Abort if not already done.  If done then treat as if cancel failed. */ 
		if(!complete){
			fprintf(stderr, "Cannot cancel generalized request - aborting program\n"); 
			MPI_Abort(MPI_COMM_WORLD, 99); 
		} 
		return MPI_SUCCESS;
	}
	template<class It>
	auto ireceive(It d_first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		// based on http://liinwww.ira.uka.de/courses/spprakt/mpi2-html-doc/node157.html
		mpi3::request ret; /*	receive_args<It>* args = (receive_args<It>*)::malloc(sizeof(receive_args<It>)); args->commP = this; args->d_first = d_first; //	args->d_last = d_last; args->source = source; args->tag = tag; args->requestP = &ret.impl_;*/
		receive_state* rs = (receive_state*)::malloc(sizeof(receive_state));
		rs->cancelled = 0;
		rs->source = source;
		rs->tag = tag;
		MPI_Grequest_start(query_fn, free_fn, cancel_fn, rs, &ret.impl_);//args->requestP);
		std::thread( //	static_cast<void*(*)(void*)>(receive_thread<It>), args
			[this, d_first, source, tag, &ret](){
				this->receive(d_first, source, tag); //	receive_args<It>* args = (receive_args<It>*)ptr; //	args->commP->receive(args->d_first, args->source, args->tag);//, /*args->d_last,*/ );
				MPI_Grequest_complete(ret.impl_); //	MPI_Grequest_complete(*args->requestP); //	::free(ptr);
			}
		).detach();	//	t.detach(); //	pthread_t thread; //	pthread_create(&thread, NULL, static_cast<void*(*)(void*)>(receive_thread<It>), args); //	pthread_detach(thread);
		return ret;		
	}
	template<class It>
	auto ireceive(
		It d_first, It d_last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		return ireceive_n(d_first, std::distance(d_first, d_last), source, tag);
	}
	template<class It>
	auto receive(It d_first, It d_last, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(
			d_first, d_last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			source, tag
		);
	}
	template<typename InputIt, typename Size, class V = typename std::iterator_traits<InputIt>::value_type>
	auto bsend_init_n(
		InputIt first, Size count, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int dest, int tag
	){
		request ret;
		int s = MPI_Bsend_init(
			std::addressof(*first), count, detail::basic_datatype<V>{},
			dest, tag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot bsend init"};
		return ret;
	}
	template<typename It, typename Size>
	auto bsend_init_n(It first, Size count, int dest, int tag = 0){
		return bsend_init_n(
			first, count,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);				
	}
	template<class It>
	auto bsend_init(
		It first, It last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return bsend_init_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto bsend_init(It first, It last, int dest, int tag = 0){
		bsend_init(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);				
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto bsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(buffered_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class InputIt, class V = typename std::iterator_traits<InputIt>::value_type>
	auto dynamic_receive(InputIt first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
	//	auto count = probe(source, tag).count<V>();
	//	return receive(first, first + count, source, tag);
		MPI_Status status;
	    MPI_Message msg;
    	int count = -1;
    	MPI_Mprobe(source, tag, impl_, &msg, &status);
    	MPI_Get_count(&status, detail::basic_datatype<V>{}, &count);
    //	int* buffer = (int*)malloc(sizeof(int) * count);
    	using detail::data;
    	MPI_Mrecv(data(first), count, detail::basic_datatype<V>{}, &msg, MPI_STATUS_IGNORE);
	}

/*	template<class InputIterator>
	auto receive(InputIterator first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return dynamic_receive(first, source, tag);
	}*/ 

	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto breceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(buffered_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto ibsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(buffered_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto ibreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(buffered_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}

	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto ssend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(synchronous_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto sreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(synchronous_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto issend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(synchronous_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class CommunicationMode, class It1, class Size>
	auto isend_n(
		CommunicationMode, 
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int tag
	){
		mpi3::request ret;
		CommunicationMode{}.ISend(
			std::addressof(*first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{}, dest, tag, impl_, &ret.impl_
		);
		return ret;
	}
	template<class It1, class Size>
	auto issend_n(It1 first, Size count, int dest, int tag = 0){
		return isend_n(
			synchronous_communication_mode{}, 
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			dest, tag
		);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto isreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(synchronous_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}

	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto rsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(ready_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto rreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(ready_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto irsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(ready_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto irreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(ready_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}
	template<class CommunicationMode, class BlockingMode, class InputIterator>
	auto send_category(CommunicationMode cm, BlockingMode bm, std::input_iterator_tag, 
		InputIterator first, InputIterator last, int dest, int tag
	);
	template<class CommunicationMode, class BlockingMode, class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto send(CommunicationMode cm, BlockingMode bm, InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send_category(cm, bm, category{}, It1, It2, dest, tag);
	}
	
private:
/*	{
		auto it = mpi3::type::registered.find(typeid(V)); assert(it != boost::mpi3::type::registered.end());
		int s = cm.Send(detail::data(first), n, it->second.impl_, dest, tag, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
	}*/
//	{
//		auto p = make_package();
//		assert(0);
#if 0
		auto it = boost::mpi3::type::registered.find(typeid(ValueType));
		if(it != boost::mpi3::type::registered.end()){
			int s = cm.Send(detail::data(first), n, it->second.impl_, dest, tag, impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else if(std::is_pod<ValueType>{}){
#ifdef BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION
			throw std::logic_error("cannot send unregistered type " + std::string(typeid(ValueType).name()) + " (pod comm disallowed)");
#endif
			// allow seamless receive of pod types
			int s = cm.Send(detail::data(first), n*sizeof(ValueType), MPI_BYTE, dest, tag, impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else throw std::logic_error("cannot send unregistered non pod type " + std::string(typeid(ValueType).name()));
#endif
//	}

//	template<class CommunicationMode, class ContiguousIterator, class Size, class V = typename std::iterator_traits<ContiguousIterator>::value_type>
//	void receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag);
/*	{
		status stts;
		auto it = mpi3::type::registered.find(typeid(V)); assert(it != boost::mpi3::type::registered.end());
		int s = cm.Recv(detail::data(first), n, it->second.impl_, dest, tag, impl_, &stts.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
	}*/
#if 0
{
		auto it = boost::mpi3::type::registered.find(typeid(ValueType));
		if(it != boost::mpi3::type::registered.end()){
			int s = cm.Recv(detail::data(first), n, it->second.impl_, dest, tag, impl_, MPI_STATUS_IGNORE);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else if(std::is_pod<ValueType>{}){
#ifdef BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION
			throw std::logic_error("cannot receive not registered type " + std::string(typeid(ValueType).name()) + " (pod comm disallowed)");
#endif
			// allow seamless receive of pod types
			int s = cm.Recv(detail::data(first), sizeof(ValueType)*n, MPI_BYTE, dest, tag, impl_, MPI_STATUS_IGNORE);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else throw std::logic_error("cannot receive not registered non pod type " + std::string(typeid(ValueType).name()));
	}
#endif

	template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = typename detail::basic_datatype<ValueType> >
	request send_n_randomaccess_contiguous_builtin(CommunicationMode cm, nonblocking_mode, std::true_type, ContiguousIterator first, Size n, int dest, int tag){
		request r;
		int s = cm.ISend(detail::data(first), n, datatype{}, dest, tag, impl_, &r.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot send random access iterators"};
		return r;
	}
	template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = typename detail::basic_datatype<ValueType> >
	request receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, nonblocking_mode, std::true_type, ContiguousIterator first, Size n, int dest, int tag){
		request r;
		int s = cm.IRecv(detail::data(first), n, datatype{}, dest, tag, impl_, &r.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		return r;
	}
private:
	template<class It1, typename Size, class It2>
	auto all_to_all_n(
		It1 first,
			detail::contiguous_iterator_tag,
			detail::basic_tag, 
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag
	){
		assert( count % size() == 0 );
		MPI_(Alltoall)( 
			detail::data(first), count/size(),
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), count/size(),
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			impl_
		);
		return d_first + count;
	}
	template<class It1, typename Size>
	auto all_to_all_n(
		It1 first,
			detail::contiguous_iterator_tag,
			detail::basic_tag, 
		Size count 
	)
	->decltype(MPI_(Alltoall)( 
			MPI_IN_PLACE, 0*count/size(),
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(first), count/size(),
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			impl_
		), first + count)
	{
		assert( count % size() == 0 );
		return MPI_(Alltoall)( 
			MPI_IN_PLACE, 0*count/size(),
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(first), count/size(),
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			impl_
		), first + count;}
public:
	template<class It1, typename Size, class It2>
	auto all_to_all_n(It1 first, Size count, It2 d_first){
		assert( count % size() == 0 );
		return 
			all_to_all_n(
				first, 
					detail::iterator_category_t<It1>{},
					detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
				count,
				d_first,
					detail::iterator_category_t<It2>{},
					detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}
			);
	}
	template<class It1, typename Size>
	auto all_to_all_n(It1 first, Size count)
	->decltype(
			all_to_all_n(
				first, 
					detail::iterator_category_t<It1>{},
					detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
				count
			))
	{
		assert( count % size() == 0 );
		return 
			all_to_all_n(
				first, 
					detail::iterator_category_t<It1>{},
					detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
				count
			);}
	template<class It1, class It2>
	auto all_to_all(It1 first, It2 d_first){return all_to_all_n(first, size(), d_first);}
	template<class It1> 
	auto all_to_all(It1 first)
	->decltype(all_to_all_n(first, size())){
		return all_to_all_n(first, size());}
public:
	template<class T>
	auto operator()(T&& t)
	->decltype(communicator::all_to_all(begin(std::forward<T>(t))), void()){
		assert(t.size() == size());
		auto e = communicator::all_to_all(begin(std::forward<T>(t)));
		using std::end;
		assert( e == end(t) );
		return std::forward<T>(t);
	}
	template<class T> 
	NODISCARD("do not ignore result when argument is const")
	auto operator()(T const& t)
	->decltype(communicator::all_to_all(t.begin(), std::declval<T>().begin()), T(communicator::size())){
		assert(t.size() == communicator::size());
		T ret(communicator::size()); 
		communicator::all_to_all(t.begin(), ret.begin());
		return ret;
	}
private:
/*	template<class BlockingMode, class InputIt, class OutputIt, class InputCategory = typename std::iterator_traits<OutputIt>::iterator_category, class OutputCategory = typename std::iterator_traits<OutputIt>::iterator_category >
	auto gather(BlockingMode bm, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, int root = 0){
		return gather_category(bm, InputCategory{}, OutputCategory{}, first, last, d_first, d_last, root);
	}
	template<class BlockingMode, class InputIt, class OutputIt>
	auto gather_category(BlockingMode bm, std::random_access_iterator_tag, std::random_access_iterator_tag, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, int root){
		return gather_n_randomaccess(bm, first, std::distance(first, last), d_first, std::distance(d_first, d_last), root);
	}
	template<class BlockingMode, class RandomAccessIt1, class RandomAccessIt2, class Size>
	auto gather_n_randomaccess(BlockingMode bm, RandomAccessIt1 first, Size count, RandomAccessIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_contiguity(detail::is_contiguous<RandomAccessIt1>{}, detail::is_contiguous<RandomAccessIt2>{}, first, count, d_first, d_count, root);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_contiguity(BlockingMode bm, std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_contiguous_builtin(bm,
			detail::is_basic<typename std::iterator_traits<ContiguousIt1>::value_type>{}, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt2>::value_type>{}, 
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_contiguity(BlockingMode bm, std::false_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_noncontiguous_contiguous_strided(bm, 
			detail::is_strided<ContiguousIt1>{}, detail::is_strided<ContiguousIt2>{},
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_noncontiguous_contiguous_strided(BlockingMode bm, std::true_type, std::false_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(bm, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt1>::value_type>{}, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt2>::value_type>{}, 
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(blocking_mode,std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		mpi3::type datatype1_strided = mpi3::type(datatype1{}).vector(count, 1, first.stride());
		int s = MPI_Gather(detail::data(first.base()), 1, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all gather");
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(nonblocking_mode,std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		mpi3::type datatype1_strided = mpi3::type(datatype1{}).vector(count, 1, first.stride());
		request ret;
		int s = MPI_Igather(detail::data(first.base()), 1, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot Igather");
		return ret;
	}
*/
/*	template<class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_contiguous_builtin(std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		int s = MPI_Gather(detail::data(first), count, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
	}*/
private:
	template<class It1, class It2>
	auto scatter_builtinQ(std::true_type, It1 first, It1 last, It2 d_first, int root){
		return scatter_builtin(
			detail::iterator_category<It1>{}, 
			detail::iterator_category<It2>{}, first, last, d_first, root
		);
	}
	template<class It1, class It2, 
		class sendtype = detail::basic_datatype<typename std::iterator_traits<It1>::value_type>,
		class recvtype = detail::basic_datatype<typename std::iterator_traits<It2>::value_type>
	>
	auto scatter_builtin(detail::contiguous_iterator_tag, detail::contiguous_iterator_tag, It1 first, It1 last, It2 d_first, int root){
		int s = MPI_Scatter( detail::data(first), std::distance(first, last), 
			sendtype::value,
			detail::data(d_first), std::distance(first, last),
			recvtype::value,
			root,
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class Iterator1, class Iterator2>
	auto scatter_builtinQ(std::false_type, Iterator1 first, Iterator2 last, Iterator1 d_first, int root)
//	{ TODO implement }
	;
public:
	template<class T>
	auto broadcast_value(T& t, int root = 0){return broadcast_n(std::addressof(t), 1, root);}
	template<class It, typename Size>
	auto ibroadcast_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int root
	){
		request r;
		int s = MPI_Ibcast(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			root, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot ibroadcast"};
		return r;
	}
	template<class It, typename Size>
	void broadcast_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int root
	){
		auto e = static_cast<enum error>(
			MPI_Bcast(
				detail::data(first), count, 
				detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			root, impl_)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot broadcast"};
	}
	template<class It>
	auto ibroadcast(
		It first, It last, 
			detail::random_access_iterator_tag,
		int root
	){
		return ibroadcast_n(
			first, std::distance(first, last), root
		);
	}
	template<class It>
	void broadcast(
		It first, It last, 
			detail::random_access_iterator_tag,
		int root
	){broadcast_n(first, std::distance(first, last), root);}
public:
	template<class It, typename Size>
	auto ibroadcast_n(It first, Size count, int root = 0){
		return ibroadcast_n(
			first,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			root
		);
	}
	template<class It, typename Size>
	void broadcast_n(It first, Size count, int root = 0){
		return broadcast_n(
			first,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count, 
			root
		);
	}
	template<class It>
	auto ibroadcast(It first, It last, int root = 0){
		return ibroadcast(
			first, last,
				detail::iterator_category_t<It>{},
			root
		);
	}
	template<class It>
	void broadcast(It first, It last, int root = 0){
		return broadcast(
			first, last,
				detail::iterator_category_t<It>{},
			root
		);
	}
	template<class T, class Op = std::plus<> >
	void reduce_value(T const& t, T& ret, Op op = {}, int root = 0){
		reduce_n(std::addressof(t), 1, std::addressof(ret), op, root); 
	}
	template<class T, class Op = std::plus<> >
	auto reduce_value(T const& t, Op op = {}, int root = 0){
		T ret = T(0);
		reduce_value(t, ret, op, root); // if(rank() == root) return optional<T>(ret);
		return ret;
	}
	template<class It1, class Size, class It2, class Op, class PredefinedOp>
	It2 reduce_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Op,
		PredefinedOp,
		int root
	){
		static_assert(std::is_same<typename std::iterator_traits<It1>::value_type, typename std::iterator_traits<It2>::value_type>{}, "!");
		int s = MPI_Reduce(
			detail::data(first), detail::data(d_first), count, 
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, PredefinedOp{}, 
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce n");
		return rank()==root?d_first + count:d_first;
	}
			
	template<class It1, class Size, class It2, class Op>
	It2 reduce_n(It1 first, Size count, It2 d_first, Op op, int root = 0){
		return reduce_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			op, 
				predefined_operation<Op>{},
			root
		);			
	}
//protected: // these functions have been made public by popular demand, they might be deprecated in the future by a more composable idiom
public:
	template<class T, class Op = std::plus<> >
	void all_reduce_value(T const& t, T& ret, Op op = {}){
		auto e = all_reduce_n(std::addressof(t), 1, std::addressof(ret), op); (void)e;
		assert( e = std::addressof(ret) + 1 );
	}
	template<class T, class Op = std::plus<>, typename = decltype(T{T(0)})>
	auto all_reduce_value(T const& t, Op op = {}){
		T ret = T(0);
		all_reduce_value(t, ret, op); // if(rank() == root) return optional<T>(ret);
		return ret;
	}
	template<class Op = std::logical_and<>>
	bool all_reduce_value(bool t, Op op={}){
		int ret; all_reduce_value(int{t}, ret, op); return ret;
	}
	template<class T>
	auto max(T const& t){
		auto ret = std::numeric_limits<T>::lowest(); 
		all_reduce_value(t, ret, mpi3::max<>{}); return ret;
	}
	template<class T>
	auto min(T const& t){
		auto ret = std::numeric_limits<T>::max();
		all_reduce_value(t, ret, mpi3::min<>{}); return ret;
	}

public:
	template<
		class It1, class Size, class It2, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class V2 = typename std::iterator_traits<It2>::value_type,
		class P1 = decltype(detail::data(It1{})), class P2 = decltype(detail::data(It2{})),
		typename = std::enable_if_t<std::is_same<V1, V2>{}>,
		class PredefinedOp = predefined_operation<Op>,
		typename = decltype(predefined_operation<Op>::value)
	>
	It2 all_reduce_n(It1 first, Size count, It2 d_first, Op /*op*/ = {}){
		using detail::data;
		int s = MPI_Allreduce(data(first), detail::data(d_first), count, detail::basic_datatype<V1>{}, PredefinedOp{}/*op*/, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce n");
	}
	template<class It1, class Size, class It2, class Op = std::plus<>, class PredefinedOp>
	auto all_reduce_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Op,
			PredefinedOp
	){
		static_assert(std::is_same<typename std::iterator_traits<It1>::value_type, typename std::iterator_traits<It2>::value_type>{}, "!");
		using detail::data;
		int s = MPI_Allreduce(detail::data(first), detail::data(d_first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{}, PredefinedOp{}/*op*/, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot reduce n"};
		return d_first + count;
	}
	template<class It1, class Size, class It2, class Op = std::plus<>>
	It2 all_reduce_n(It1 first, Size count, It2 d_first, Op op = {}){
		return all_reduce_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			op, 
				predefined_operation<Op>{}
		);
	}
	template<typename It1, typename It2, class Op = std::plus<>>
	It2 all_reduce(It1 first, It1 last, It2 d_first, Op op = {}){
		return all_reduce_n(first, std::distance(first, last), d_first, op);
	}
	template<class T> static auto data_(T&& t){
		using detail::data;
		return data(std::forward<T>(t));
	}
	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(data_(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto all_reduce_in_place_n(It1 first, Size count, Op /*op*/){ // TODO check why op is not used 
		int s = MPI_Allreduce(MPI_IN_PLACE, data_(first), count, detail::basic_datatype<V1>{}, PredefinedOp{}, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all reduce n");
	}
	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>,
		typename = decltype(predefined_operation<Op>::value)
	>
	void all_reduce_n(It1 first, Size count, Op op){
		return all_reduce_in_place_n(first, count, op);
	}
	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(data_(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_in_place_n(It1 first, Size count, Op /*op*/, int root = 0){
		int s = MPI_SUCCESS;
		if (rank() == root) {
			s = MPI_Reduce(MPI_IN_PLACE, data_(first), count, detail::basic_datatype<V1>{}, PredefinedOp{}, root, impl_);
		} else {
			s = MPI_Reduce(data_(first), NULL, count, detail::basic_datatype<V1>{}, PredefinedOp{}, root, impl_);
		}
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot reduce n"};
	}
	template<
		class It1, class Size, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_n(It1 first, Size count, Op op = {}, int root = 0){
		if(rank() == root) return reduce_in_place_n(first, count, op, root);
		return reduce_n(first, 0, nullptr, op, root);
	}
	template<
		class It1, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_in_place(It1 first, It1 last, Op op, int root = 0){
		return reduce_in_place_n(first, std::distance(first, last), op, root);
	}
	template<
		class It1, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce(It1 first, It1 last, Op op = {}, int root = 0){
		return reduce_n(first, std::distance(first, last), op, root);
	}
	
private:

	template<class ReducePolicy, class It1, class Size, class It2, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n(ReducePolicy rp, It1 first, Size count, It2 d_first, Op op, int root = 0){
		return reduce_n_dispatch(rp, 
			std::integral_constant<bool, detail::is_basic<V1>{} and detail::is_basic<V2>{} and detail::is_contiguous<It1>{} and detail::is_contiguous<It2>{}>{}, 
			first, count, d_first, op, root
		);
	}
	template<class ReducePolicy, class RandomAccessIt1, class It2, class Op>
	auto reduce_category(ReducePolicy rp, std::random_access_iterator_tag, RandomAccessIt1 first, RandomAccessIt1 last, It2 d_first, Op op, int root){
		return reduce_n(rp, first, std::distance(first, last), d_first, op, root);
	}
	template<class It1, class Size, class It2, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(all_reduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int /*root*/ = 0){
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce");
	}
	template<class It1, class Size, class It2, class Op,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(reduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int root = 0){
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_,
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2, class Op,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(ireduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int root = 0){
		request ret;
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_,
			root, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
		return ret;
	}
public:

#if 0
	template<class InputIt1, class Size, class InputIt2, 
		class V1 = typename std::iterator_traits<InputIt1>::value_type, 
		class V2 = typename std::iterator_traits<InputIt2>::value_type 
	>
	void scatter_n(InputIt1 first, Size count, InputIt2 d_first, int root = 0){
		return scatter_n_dispatch( 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<InputIt1>::value and detail::is_contiguous<InputIt2>::value
			>{}, first, count, d_first, root
		);
	}
	template<class InputIt1, class Size, class InputIt2, 
		class V1 = typename std::iterator_traits<InputIt1>::value_type, 
		class V2 = typename std::iterator_traits<InputIt2>::value_type 
	>
	void scatter_from_n(InputIt2 d_first, Size count, InputIt1 first, int root = 0){
		return scatter_from_n_dispatch( 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<InputIt1>::value and detail::is_contiguous<InputIt2>::value
			>{}, d_first, count, first, root
		);
	}
	template<class InputIt1, class InputIt2>
	void scatter(InputIt1 first, InputIt1 last, InputIt2 d_first, int root = 0){
		return scatter_category(typename std::iterator_traits<InputIt1>::iterator_category{}, first, last, d_first, root);
	}
	template<class It1, class OIt2>
	auto scatter(It1 first, It1 last, It2 d_first, detail::output_iterator_tag, int root){
		vector<typename detail::iterator_traits<It2>::value_type> buff(std::distance(first, last)/size());
		auto e = scatter(first, last, typename detail::iterator_traits<It1>::iterator_category{}, buff.data()); assert( e == buff.data() + buff.size() );
		return move(begin(buff), end(buff), d_first);
	}
	template<class It1, class T>
	auto scatter(It1 first, It1 last, detail::input_iterator_tag, T const* d_first, int root){
		vector<typename detail::iterator_traits<It1>::value_type> buff(first, last);
		scatter_n(buff.data(), buff.size(),  
		return scatter(first, last, typename detail::iterator_traits<It1>::iterator_category{}, d_first);
	}
	template<class In1, class Size, class It2>
	auto scatter_n(In1 first, Size n, detail::input_iterator_tag, It2 d_first){
		vector<typename detail::iterator_traits<In1>::value_type> buff; buff.reserve(n);
		std::copy_n(first, n, begin(buff));
		scatter_n(buff.data(), n, d_first, detail::iterator_traits<It2>::iterator_category{});
	}
	template<class In1, class Size, class It2>
	auto scatter_n(In1 first, Size n, detail::input_iterator_tag, It2 d_first){
		vector<typename detail::iterator_traits<In1>::value_type> buff; buff.reserve(n);
		std::copy_n(first, n, begin(buff));
		scatter_n(buff.data(), n, d_first, detail::iterator_traits<It2>::iterator_category{});
	}
	template<class It1, class Size, class It2>
	auto scatter_n(It1 first, Size n, It2 d_first){
		assert( n % size() == 0 );
		return scatter_n(first, n, typename detail::iterator_traits<It1>::iterator_category, d_first);
	}
	template<class In1, class It2>
	auto scatter(In1 first, In1 last, detail::input_iterator_tag, It2 d_first, int root){
		vector<typename detail::iterator_traits<In1>::value_type> buff(first, last);
		return scatter_n(buff.data(), buff.size(), d_first, root);
	}
	template<class It1, class It2>
	auto scatter(It1 first, It1 last, It2 d_first, int root = 0){
		return scatter(first, last, detail::iterator_traits<It1>::iterator_category{}, root);
	}
	template<class It1, class It2>
	auto scatter(It1 first, It1 last, It2 d_first, int root = 0){
		return scatter(first, last, d_first, detail::iterator_traits<It2>::iterator_category{}, root);
	}
	template<class InputIt1, class InputIt2>
	void scatter_from(InputIt1 d_first, InputIt1 d_last, InputIt2 first, int root = 0){
		return scatter_from_category(typename std::iterator_traits<InputIt1>::iterator_category{}, d_first, d_last, first, root);
	}
	template<class InputIt1, class T>
	void scatter_value(InputIt1 first, T& t, int root = 0){
		return scatter_n(first, size(), std::addressof(t), root);
	}
	template<class InputIt1>
	typename std::iterator_traits<InputIt1>::value_type 
	scatter_value(InputIt1 first, int root = 0){
		typename std::iterator_traits<InputIt1>::value_type ret;
		scatter_value(first, ret, root);
		return ret;
	}
	template<class T>
	T scatter_value(std::vector<T> const& v, int root = 0){
		if(rank() == root) assert( static_cast<std::ptrdiff_t>(v.size()) == size() );
		return scatter_value(v.begin(), root);
	}
private:
	template<class InputIt1, class InputIt2>
	void scatter_category(std::random_access_iterator_tag, InputIt1 first, InputIt1 last, InputIt2 d_first, int root){
		return scatter_n(first, std::distance(first, last), d_first, root);
	}
	template<class InputIt1, class InputIt2>
	void scatter_from_category(std::random_access_iterator_tag, InputIt1 d_first, InputIt1 d_last, InputIt2 first, int root){
		return scatter_from_n(d_first, std::distance(d_first, d_last), first, root);
	}
	template<class ContIt1, class Size, class ContIt2,
		class V1 = typename std::iterator_traits<ContIt1>::value_type, 
		class V2 = typename std::iterator_traits<ContIt2>::value_type 
	>
	void scatter_n_dispatch(std::true_type, ContIt1 first, Size count, ContIt2 d_first, int root){
		if(count % size() != 0) throw std::runtime_error("not matching size for scatter, elements count is " + std::to_string(count) + ", comm size is " + std::to_string(size()));
		int s = MPI_Scatter(
			detail::data(first), count/size(), detail::basic_datatype<V1>{},
			detail::data(d_first), count/size(), detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class ContIt1, class Size, class ContIt2,
		class V1 = typename std::iterator_traits<ContIt1>::value_type, 
		class V2 = typename std::iterator_traits<ContIt2>::value_type 
	>
	void scatter_from_n_dispatch(std::true_type, ContIt2 d_first, Size count, ContIt1 first,  int root){
		if(count % size() != 0) throw std::runtime_error("not matching size for scatter");
		int s = MPI_Scatter(
			detail::data(first), count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
#endif
//////////////////////////////
	template<class CIt1, class Size, class CIt2>
	auto scatter_n(
		CIt1 first, Size n, detail::contiguous_iterator_tag, detail::basic_tag,
		CIt2 d_first,       detail::contiguous_iterator_tag, detail::basic_tag,
		int root
	){
		auto e = static_cast<enum error>( MPI_Scatter(
			detail::data(first  ), n/size(), detail::basic_datatype<typename std::iterator_traits<CIt1>::value_type>{},
			detail::data(d_first), n/size(), detail::basic_datatype<typename std::iterator_traits<CIt2>::value_type>{},
			root, impl_
		) );
		if(e != mpi3::error::success) throw std::system_error{e, "cannot scatter"};
		if(rank()==root) std::advance(d_first, n);
		return d_first;
	}
	template<class In1, class Size, class It2, class Any2, class Any3>
	It2 scatter_n(
		In1 first, Size n, detail::input_iterator_tag, detail::basic_tag,
		It2 d_first,       Any2                      , Any3,
		int root
	){
		vector<typename std::iterator_traits<In1>::value_type> buff; buff.reserve(n);
		using std::copy_n;
		copy_n(first, n, std::back_inserter(buff));
	//	auto e = 
		scatter_n(buff.begin(), n, d_first, root);
		std::advance(d_first, n/size());
		return d_first;
	}
	template<class It1, class Size, class It2, class V1 = typename std::iterator_traits<It1>::value_type, class V2 = typename std::iterator_traits<It2>::value_type>
	auto scatter_n(It1 first, Size n, It2 d_first, int root = 0){
		return scatter_n(
			first, n, detail::iterator_category_t<It1>{}, detail::value_category_t<V1>{},
			d_first,  detail::iterator_category_t<It2>{}, detail::value_category_t<V2>{},
			root
		);
	}
	template<class It1, class Size, class It2>
	auto scatter_n_from(It2 d_first, Size n, It1 first, int root = 0){
		return scatter_n(first, n*size(), d_first, root);
	//	std::advance(d_first, n);
	//	return d_first;
	}
	template<class RA1, class It2>
	auto scatter(RA1 first, RA1 last, detail::random_access_iterator_tag, It2 d_first, int root){
		return scatter_n(first, std::distance(first, last), d_first, root);
	}
	template<class It1, class It2>
	auto scatter(It1 first, It1 last, It2 d_first, int root = 0){
		return scatter(first, last, detail::iterator_category_t<It1>{}, d_first, root);
	}
	template<class RA2, class It1>
	auto scatter_from(RA2 first, RA2 last, detail::random_access_iterator_tag, It1 d_first, int root){
		return scatter_n_from(first, std::distance(first, last), d_first, root);
	}
	template<class V, class It1>
	auto scatter_value_from(V& v, It1 first, int root = 0){
		return scatter_n_from(std::addressof(v), 1, first, root); 
	}
	template<class It1, class V = typename std::iterator_traits<It1>::value_type>
	V scatter(It1 first, It1 last, int root = 0){
		if(rank()==root) assert(std::distance(first, last) == size());
		V v;
	//	auto e = 
		scatter_value_from(v, first, root);
	//	std::cerr << std::distance(e, last) << std::endl;
	//	assert(e == last);
		return v;
	}
	template<class It, class V = typename std::iterator_traits<It>::value_type>
	V scatter(It first, int root = 0){
		V v;  
		scatter_value_from(v, first, root); 
		return v;
	}
	template<class Container, class V = typename std::iterator_traits<typename Container::iterator>::value_type>
	V scatter(Container c, int root = 0, void* = 0){
		assert( (int)c.size() == (rank()==root?size():0) );
		using std::begin;
		return scatter(begin(c), root);
	}

	template<class It1, class It2>
	auto scatter_from(It2 d_first, It2 d_last, It1 first, int root = 0){
		return scatter_from(d_first, d_last, detail::iterator_category_t<It2>{}, first, root);
	}
////////////////////////////////////////////////////////////////////////////////
	template<class CIt1, class Size, class CIt2>
	void scatterv_n(
		CIt1 first, int* counts, int* displs,
		detail::contiguous_iterator_tag, detail::basic_tag,
		CIt2 d_first, Size n, 
		detail::contiguous_iterator_tag, detail::basic_tag,
		int root = 0
	){
		MPI_(Scatterv)( 
			detail::data(first), counts, displs,
			detail::basic_datatype<typename std::iterator_traits<CIt1>::value_type>{},
			detail::data(d_first), n,
			detail::basic_datatype<typename std::iterator_traits<CIt2>::value_type>{},
			root, impl_
		);
	}
	template<class CItCIt1, class CItN, class CIt2>
	auto scatterv_n(
		CItCIt1 citcit1, CItN ns, detail::contiguous_iterator_tag, detail::basic_tag, 
		CIt2 it2,                 detail::contiguous_iterator_tag, detail::basic_tag,
		int root
	){
		std::vector<int> counts(ns, ns + size());
		std::vector<int> displs(size());
		for(int i = 0; i != size(); ++i) displs[i+1] = detail::data(citcit1[i+1]) - detail::data(citcit1[i]);//adjacent_difference doesn't work here becuase of type incompatibility
	//	std::adjacent_difference(citcit1, citcit1 + size(), begin(displs) + 1, [](auto& a, auto& b){return detail::data(a) - detail::data(b);});
		int n = scatter(counts.begin(), counts.end());
		scatterv_n(
			detail::data(*citcit1), counts.data(), displs.data(), detail::contiguous_iterator_tag{}, detail::basic_tag{},
			it2, n                                              , detail::contiguous_iterator_tag{}, detail::basic_tag{},
			root
		);
		return it2;
	}
	template<class CItIt1, class CItN, class It2>
	auto scatterv_n(CItIt1 citit1, CItN ns, It2 it2, int root = 0){
		scatterv_n(
			citit1, ns, detail::iterator_category_t<typename std::iterator_traits<CItIt1>::value_type>{}, detail::basic_tag{},
			it2       , detail::iterator_category_t<It2>{},                                               detail::basic_tag{},
			root
		);
	}
	template<class Container, class It2, typename = typename Container::iterator>
	auto scatterv(Container const& c, It2 it2, int root = 0){
		assert( (int)c.size() == ((rank()==root)?size():0) );
		if(rank()==root){
			std::cerr<< "in scatterv " << detail::data(c[1].begin()) - detail::data(c[0].begin()) << " " << detail::data(c[2].begin()) - detail::data(c[0].begin()) << std::endl;
		}
		
		using std::begin; using std::end; 
		std::vector<int> displs(c.size());
		for(int i = 0; i!=(int)displs.size(); ++i){
			std::ptrdiff_t d = detail::data(c[i].begin()) - detail::data(c[0].begin());
			assert( d <= static_cast<std::ptrdiff_t>(std::numeric_limits<int>::max()) );
			displs[i] = d;
		}
		if(rank()==root){
			std::cerr<< "in scatterv 2 " << displs[0] << " " <<  displs[1] << " " <<  displs[2] << std::endl;
		}
		std::vector<int> counts(c.size());
		std::transform(counts.begin(), counts.end(), begin(c), counts.begin(), [](auto&, auto& b){return std::distance(begin(b), end(b));}); 
		int n = scatter(counts);
		scatterv_n(
			detail::data(begin(*begin(c))), counts.data(), displs.data(), detail::contiguous_iterator_tag{}, detail::basic_tag{}, 
			detail::data(it2), n                                        , detail::contiguous_iterator_tag{}, detail::basic_tag{}, 
			root
		);
		return it2 + n;
	}
	template<class T> class id{using type = T;};
	
	template<class MultiIt, class Size, class MultiIt2, class=std::enable_if_t<(MultiIt::dimensionality>=1)> >
	auto gather_n(MultiIt first, Size count, MultiIt2 d_first, int root=0)
	->decltype( MPI_Gather (mpi3::base(first), count, mpi3::type{first}, mpi3::base(d_first), count, mpi3::type{d_first}, root, impl_), d_first+count){
		return MPI_(Gather)(mpi3::base(first), count, mpi3::type{first}, mpi3::base(d_first), count, mpi3::type{d_first}, root, impl_), d_first+count;}

	template<class It1, typename Size1, class It2, typename Size2>
	It2 gather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count, 
		int root
	){
		MPI_(Gather)(
				detail::data(first), count, 
				detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
				detail::data(d_first), d_count,
				detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
				root, impl_
		);
		return d_first + ((rank()==root)?d_count*size():0);
	}
	template<class It1, class Size1, class It2, class Size2>
	auto gather_n(It1 first, Size1 count, It2 d_first, Size2 d_count, int root){
		return gather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count,
			root
		);
	}
	template<class It1, class Size1, class It2, class=std::enable_if_t<(not has_dimensionality<It1>{})>>
	auto gather_n(It1 first, Size1 count, It2 d_first, int root = 0){
		return gather_n(first, count, d_first, count, root);
	}
	
	template<class It2, class Size, class It1>
	auto scatterv_n_from(It2 d_first, Size n, It1 first, int root = 0){
		std::vector<int> counts(size());//rank()==root?size():0);
		int nn = n;
		gather_n(std::addressof(nn), 1, counts.data(), root);
		std::vector<int> displs(size());//rank()==root?size():0);
		partial_sum(counts.begin(), counts.end(), displs.begin()+1);
		scatterv_n(
			first, counts.data(), displs.data(), 
			detail::iterator_category_t<It1>{}, detail::value_category_t<typename std::iterator_traits<It1>::value_type>{}, 
			d_first, n, 
			detail::iterator_category_t<It2>{}, detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}, 
			root
		);
		std::advance(first, rank()==root?displs.back() + counts.back():0);
		return first;
	}
	template<class RA2, class It1>
	auto scatterv_from(RA2 d_first, RA2 d_last, detail::random_access_iterator_tag, It1 first, int root){
		return scatterv_n_from(d_first, std::distance(d_first, d_last), first, root);
	}
	template<class It1, class It2>
	auto scatterv_from(It2 d_first, It2 d_last, It1 first, int root = 0){
		return scatterv_from(d_first, d_last, detail::iterator_category_t<It2>{}, first, root);
	}
public:
	template<class T, class It> 
	void all_gather_value(T const& t, It first){all_gather_n(std::addressof(t), 1, first);}
	template<class T> std::vector<T> 
	all_gather_value(T const& t){
		std::vector<T> ret(size());
		all_gather_value(t, ret.begin());
		return ret;
	}

	template<typename T, typename It> 
	It gather_value(T const& t, It first, int root){
		return gather_n(std::addressof(t), 1, first, root);
	}
	template<class T> 
	std::vector<T> gather_value(T const& t, int root = 0){
		std::vector<T> ret((rank() == root)?size():0);
		gather_value(t, ret.begin(), root);
		return ret;
	}
	protected:
	template<class It, typename Size> 
	void advance(It& it, Size count){std::advance(it, count);}
	template<class It, typename Size> 
	void advance(It& it, Size s, int r){std::advance(it, rank()==r?s:0);}
	template<
		class GatherMode,
		typename It1, typename Size, typename It2, 
		typename = std::enable_if_t<
			std::is_same<
				typename std::iterator_traits<It1>::value_type, 
				typename std::iterator_traits<It2>::value_type
			>{}
		>, class... Root
	>
	auto a_gather_n(
		GatherMode gm,
		It1 first, 
			detail::contiguous_iterator_tag, 
			detail::memcopyable_tag,
		Size count, 
		It2 d_first, 
			detail::contiguous_iterator_tag,
			detail::memcopyable_tag,
		Root... root
	){
		int s = gm(
			detail::data(first), count*sizeof(*first), MPI_BYTE, 
			detail::data(d_first), count*sizeof(*d_first), MPI_BYTE, 
			root..., impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
		advance(d_first, count*size(), root...);
	//	std::advance(d_first, count);
		return d_first;
	}
	public:
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(
		It1 first  , Size count  , 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, Size d_count,
			detail::contiguous_iterator_tag,
			detail::basic_tag
	){
		auto e = static_cast<enum error>(
			MPI_Allgather(
				detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
				detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
				impl_
			)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot allgather, all_gather works if all sending sizes are equal"};
		return d_first + d_count*size();
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(It1 first, Size count, It2 d_first, Size d_count){
		return all_gather_n(
			first, count,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first, d_count,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}
		);
	}
	template<class It1, typename It2>
	It2 all_gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag, 
		It2 d_first
	){return all_gather_n(first, std::distance(first, last), d_first);}
	template<class It1, typename It2>
	It2 all_gather(It1 first, It1 last, It2 d_first){
		return all_gather(first, last, detail::iterator_category_t<It1>{}, d_first);
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(
		It1 first, Size count, 
		It2 d_first, 
			detail::output_iterator_tag,
		CountsIt counts, DisplsIt displs
	){
		auto s = std::accumulate(counts, counts + size(), typename std::iterator_traits<CountsIt>::value_type{0});
		std::vector<typename std::iterator_traits<It1>::value_type> buff(s);
		auto e = all_gatherv_n(first, count, buff.data(), counts, displs);
		assert( e == buff.data() + buff.size() );
		using std::move;
		return move(buff.begin(), buff.end(), d_first);
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(
		It1 first, Size count, 
		It2 d_first, 
			detail::forward_iterator_tag,
		CountsIt counts, DisplsIt displs){
		return all_gatherv_n(
			first, count, 
				detail::iterator_category_t<It1>{}, 
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			counts, displs); 
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(It1 first, Size count, It2 d_first, CountsIt counts, DisplsIt displs){
		return all_gatherv_n(
			first, count, 
			d_first,
				detail::iterator_category_t<It2>{},
			counts, displs
		);
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(
		It1 first, Size count, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		CountsIt counts, DisplsIt displs
	){
		int s = MPI_Allgatherv(
			detail::data(first), count,
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), detail::data(counts), detail::data(displs),
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all_gatherv");
		return d_first + detail::data(displs)[size()-1] + detail::data(counts)[size()-1];
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(
		It1 first  , Size count  , 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, Size d_count,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(size());
		std::vector<int> displs(size());
		all_gather_n(&posize, 1, sizes.begin(), 1);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		all_gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data());
		package_iarchive pia(pi);
		d_count*=size();
		while(d_count--) pia >> *(d_first++);
		return d_first;
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(It1 first, Size count, It2 d_first){
		return all_gather_n(first, count, d_first, count);
	}
	template<typename V, typename It2>
	auto all_gather_v(V const& v, It2 d_first){
		return all_gather_n(std::addressof(v), 1, d_first);
	}
	template<class Vector, typename V>
	auto all_gather_as(V const& v){
		Vector ret(size());
		all_gather_v(v, ret.begin());
		return ret;
	}
	template<typename It1, typename Size, typename It2>
	auto all_gatherv_n(It1 first, Size count, It2 d_first){
		std::vector<int> counts(size());
		std::vector<int> displs(size()+1);
		int c = count;
		all_gather_n(&c, 1, counts.begin());
		partial_sum(counts.begin(), counts.end(), displs.begin()+1);
		return all_gatherv_n(first, count, d_first, counts.begin(), displs.begin());
	}
	template<typename It1, typename Size, typename It2>
	auto gatherv_n(It1 first, Size count, It2 d_first){
		std::vector<int> counts(size());
		std::vector<int> displs(size()+1);
		int c = count;
		all_gather_n(&c, 1, counts.begin());
		partial_sum(counts.begin(), counts.end(), displs.begin()+1);
		return gatherv_n(first, count, d_first, counts.begin(), displs.begin());
	}		
	template<class It1, typename Size1, class It2, class Size2>
	auto all_gather_n(
		It1 first,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size1 count, 
		It2 d_first,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size2 d_count
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(size());
		std::vector<int> displs(size());
		all_gather_n(&posize, 1, sizes.begin(), 1);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		all_gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data());
		package_iarchive pia(pi);
		d_count*=size();
		while(d_count--) pia >> *(d_first++);
		return d_first;
	}
	template<typename It1, typename It2>
	auto all_gatherv(It1 first, It1 last, detail::random_access_iterator_tag, It2 d_first){
		return all_gatherv_n(first, std::distance(first, last), d_first);
	}
	template<typename It1, typename It2>
	auto all_gatherv(It1 first, It1 last, It2 d_first){
		return all_gatherv(
			first, last,
				detail::iterator_category_t<It1>{},
			d_first
		);
	}
public:
	template<class It1, typename Size1, class It2, class Itc, class Itd>
	auto gatherv_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Itc counts,
			detail::contiguous_iterator_tag,
		Itd displs,
			detail::contiguous_iterator_tag,
		int root
	){
		MPI_Gatherv(
			detail::data(first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), detail::data(counts), detail::data(displs), detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			root, impl_
		);
	}
	template<class It1, typename Size1, class It2, class Itc, class Itd>
	auto gatherv_n(
		It1 first, Size1 count,
		It2 d_first, Itc counts, Itd displs = 0,
		int root = 0
	){
		return gatherv_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count,
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			counts, 
				detail::iterator_category_t<Itc>{},
			displs, 
				detail::iterator_category_t<Itd>{},
			root
		);
	}
	template<class It1, typename Size1, class It2, typename Size2>
	auto igather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count, 
		int root
	){
		request ret;
		int s = MPI_Igather(
			detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, 
			root, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot Igather"};
		return ret;
	}
	template<class It1, typename Size1, class It2, typename Size2>
	auto iall_gather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count
	){
		request ret;
		int s = MPI_Iallgather(
			detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, 
			impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot Iallgather");
		return ret;
	}
	template<class It1, typename Size1, class It2, class Size2>
	auto gather_n(
		It1 first,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		Size1 count, 
		It2 d_first,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		Size2 d_count,
		int root
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(rank()==root?size():0);
		gather_n(&posize, 1, sizes.begin(), 1, root);
		std::vector<int> displs(sizes.size()+1);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data(), root);
		if(rank() == root){
			package_iarchive pia(pi);
			d_count*=size();
			while(d_count--) pia >> *(d_first++);
		}
		return d_first;
	}
	template<class It1, class Size1, class It2, class Size2>
	auto igather_n(It1 first, Size1 count, It2 d_first, Size2 d_count, int root){
		return igather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count,
			root
		);
	}
	template<class It1, class Size1, class It2, class Size2>
	auto iall_gather_n(It1 first, Size1 count, It2 d_first, Size2 d_count){
		return iall_gather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count
		);
	}
	template<class It1, class Size1, class It2>
	auto igather_n(It1 first, Size1 count, It2 d_first, int root = 0){
		return igather_n(first, count, d_first, count, root);
	}
	template<class It1, class Size1, class It2>
	auto iall_gather_n(It1 first, Size1 count, It2 d_first){
		return iall_gather_n(first, count, d_first, count);
	}
	template<typename It1, typename It2>
	auto gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, int root
	){
		return gather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto gather(
		It1 first, It1 last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, int root
	){
		return gather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto igather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, int root
	){
		return igather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto iall_gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first
	){
		return iall_gather_n(first, std::distance(first, last), d_first);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::input_iterator_tag,
			detail::basic_tag,
		It2 d_first, int root
	){
		mpi3::vector<typename std::iterator_traits<It1>::value_type> buffer(first, last);
		return gather_n(buffer.data(), buffer.size(), d_first, root);
	}
	template<typename It1, typename It2>
	auto gather(It1 first, It1 last, It2 d_first, int root = 0){
		return gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
			root
		);
	}
	template<typename It1, typename It2>
	auto igather(It1 first, It1 last, It2 d_first, int root){
		return igather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
			root
		);
	}
	template<typename It1, typename It2>
	auto iall_gather(It1 first, It1 last, It2 d_first){
		return iall_gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, It2 d_last,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int root
	){
		return gather_n(
			first, std::distance(first, last), 
			d_first, std::distance(d_last, d_last), 
			root
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, It2 d_last,
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int root
	){
		return gather_n(
			first, std::distance(first, last), 
			d_first, std::distance(d_last, d_last), 
			root
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		It2 d_first, It2 d_last,
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int root
	){
		mpi3::vector<typename std::iterator_traits<It1>::value_type> v(first, last);
		return gather_n(
			v.data(), v.size(), 
			d_first, std::distance(d_first, d_last), 
			root
		);
	}
	public:
	template<typename It1, typename It2>
	auto gather(It1 first, It1 last, It2 d_first, It2 d_last, int root){
		return gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first, d_last,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			root
		);
	}
//	template<class It1, class Size, class It2>
//	auto igather_n(It1 first, Size count, It2 d_first, int root = 0){return gather_n(igather_mode{}, first, count, d_first, root);}
/*
	template<class It1, class It2>
	//	[[nodiscard]] 
	auto igather(It1 first, It1 last, It2 d_first, int root = 0){
		return igather_category(
			typename std::iterator_traits<It1>::iterator_category{},
			first, last, d_first, root
		);
	}
	template<class It1, class It2>
	auto igather_category(std::random_access_iterator_tag, It1 first, It1 last, It2 d_first, int root = 0){
		return igather_n(first, std::distance(first, last), d_first, root);
	}*/
private:
/*	template<class GatherPolicy, class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n(GatherPolicy gp, It1 first, Size count, It2 d_first, int root = 0){
		return gather_n_dispatch(gp, 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<It1>::value and detail::is_contiguous<It2>::value
			>{}, first, count, d_first, root);
	}*/
/*	template<class GatherPolicy, class It1, class It2>
	auto gather(GatherPolicy gp, It1 first, It1 last, It2 d_first, int root){
		return gather_category(gp, typename std::iterator_traits<It1>::iterator_category{}, first, last, d_first, root);
	}*/
/*	template<class GatherPolicy, class RandomAccessIt1, class It2>
	auto gather_category(GatherPolicy gp, std::random_access_iterator_tag, RandomAccessIt1 first, RandomAccessIt1 last, It2 d_first, int root){
		return gather_n(gp, first, std::distance(first, last), d_first, root);
	}*/
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(gather_mode g, std::true_type, It1 first, Size count, It2 d_first, int root = 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(gather_mode g, std::false_type, It1 first, Size count, It2 d_first, int root = 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(igather_mode g, std::true_type, It1 first, Size count, It2 d_first, int root = 0){
		request r;
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
		return r;
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(all_gather_mode g, std::true_type, It1 first, Size count, It2 d_first, int /*root*/= 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}

//	template<class IContiguousIterator, class OContiguousIterator>
//	void gather_n(IContiguousIterator first, std::ptrdiff_t send_count, OContiguousIterator result, int root = 0){
//		gather_n(first, send_count, result, send_count, root);
//	}
//	template<class IContiguousIterator, class OContiguousIterator>
//	void scatter(IContiguousIterator first, IContiguousIterator last, OContiguousIterator result, int root = 0){
//		scatter_n(first, std::distance(first, last), result, root);
//	}
//	template<class IContiguousIterator, class OContiguousIterator>
//	void gather (IContiguousIterator first, IContiguousIterator last, OContiguousIterator result, int root = 0){
//		gather_n(first, std::distance(first, last), result, root);
//	}
public:
	std::string get_name() const{
		int len;
		char comm_name[MPI_MAX_OBJECT_NAME];
		int status = MPI_Comm_get_name(impl_, comm_name, &len);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot set name");
		return std::string(comm_name, len);
	}
	void set_name(std::string const& s){
		int status = MPI_Comm_set_name(impl_, s.c_str());
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot get name");
	}
	std::string name() const{
		return get_name();
	}
	void name(std::string const& s){set_name(s);}

//	communicator parent() const{
//		communicator ret;
//		MPI_Comm_get_parent(&ret.impl_);
//		return ret;
//	}
	mpi3::communicator& parent(){
		static_assert(sizeof(MPI_Comm) == sizeof(mpi3::communicator), "!");
		static_assert(std::is_same<decltype(impl_), MPI_Comm>{}, "!");
		MPI_Comm* p{}; MPI_Comm_get_parent(p); assert(p);
		return reinterpret_cast<mpi3::communicator&>(*p);
	}
	communicator spawn(std::string const& argv0, int np) const{
		communicator intercomm;
		MPI_Comm_spawn(argv0.data(), MPI_ARGV_NULL, np, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm.impl_, MPI_ERRCODES_IGNORE );
		return intercomm;
	}
	communicator intercommunicator_create(int local_leader, communicator const& peer, int remote_leader, int tag = 0) const{
		communicator ret;
		int s = MPI_Intercomm_create(impl_, local_leader, peer.impl_, remote_leader, tag, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot create intercommunicator");
		return ret;
	}
	communicator create(int local_leader, communicator const& peer, int remote_leader, int tag = 0) const{
		return intercommunicator_create(local_leader, peer, remote_leader, tag);
	}
	communicator create(group const& g) const;
	communicator create_group(group const& g, int tag) const;
	FILE* fopen(const char* filename, int amode = MPI_MODE_RDWR | MPI_MODE_CREATE);

inline std::string const& name(communicator::topology const& t){
	static std::map<communicator::topology, std::string> const names = {
		{communicator::topology::undefined, "undefined"}, 
		{communicator::topology::graph, "graph"},
		{communicator::topology::cartesian, "cartesian"}};
	return names.find(t)->second;
}

//template<class T>
//friend auto operator,(communicator& comm, T const& t){
//	std::vector<T> ret(comm.size());
//	comm.all_gather_n(std::addressof(t), 1, first, root); 
//}

template<class T>
friend T operator+=(communicator& comm, T const& t){
	return comm.all_reduce_value(t, std::plus<>{});
}

template<class T>
friend T operator&=(communicator& comm, T const& t){
	return comm.all_reduce_value(t, std::logical_and<>{});
}

template<class T>
friend communicator& operator<<(communicator& comm, T const& t){
	comm.send_value(t);
	return comm;
}


};

inline void barrier(communicator const& self){self.barrier();}

inline communicator::communicator(group const& g, int tag){
	MPI_(Comm_create_group)(MPI_COMM_WORLD, &const_cast<group&>(g), tag, &impl_);
}

inline communicator::communicator(group const& g){
	MPI_(Comm_create)(MPI_COMM_WORLD, &const_cast<group&>(g), &impl_);
}
// https://www.open-mpi.org/doc/v3.0/man3/MPI_Comm_create_group.3.php
// MPI_Comm_create_group is similar to MPI_Comm_create; however, MPI_Comm_create must be called by all processes in the group of comm, whereas MPI_Comm_create_group must be called by all processes in group, which is a subgroup of the group of comm. In addition, MPI_Comm_create_group requires that comm is an intracommunicator. MPI_Comm_create_group returns a new intracommunicator, newcomm, for which the group argument defines the communication group. No cached information propagates from comm to newcomm. 

inline communicator::communicator(communicator const& o, group const& g){
	MPI_(Comm_create)(o.impl_, &const_cast<group&>(g), &impl_);
}

inline communicator::operator group() const{
	group ret; MPI_(Comm_group)(impl_, &ret.impl_); return ret;
}
//inline group::group(communicator const& c){MPI_(Comm_group)(&const_cast<communicator&>(c), &impl_);}

inline communicator communicator::create(group const& g) const{
	communicator ret;
	int s = MPI_Comm_create(impl_, &const_cast<group&>(g), &ret.impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot crate group"};
	return ret;
}

inline communicator communicator::create_group(class group const& g, int tag = 0) const{
	communicator ret;
	MPI_(Comm_create_group)(impl_, &const_cast<group&>(g), tag, &ret.impl_);
	return ret;
}

template<class T>
inline void communicator::deallocate_shared(pointer<T>){
//	MPI_Free_mem(p.base_ptr(rank()));
}

template<class T>
inline void communicator::deallocate(pointer<T>&, MPI_Aint){
//	p.pimpl_->fence();
//	MPI_Free_mem(p.local_ptr());
//	MPI_Win_free(&p.pimpl_->impl_);
//	delete p.pimpl_;
//	p.pimpl_ == nullptr;
}

#if 0
template<class T>
inline window<T> communicator::make_window(mpi3::size_t size){
	mpi3::info inf;
	void* ptr;
	window<T> ret;
	int s = MPI_Win_allocate(size*sizeof(T), sizeof(T), inf.impl_, this->impl_, &ptr, &ret.impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot window_allocate");
	return ret;
}
#endif

struct strided_range{
	int first;
	int last;
	int stride = 1;
	int front() const{return first;}
	int back() const{return last - 1;}
	int size() const{return (last - first) / stride;}
};

#if 0
struct group{
	MPI_Group impl_;
//	static group empty(){return group();}
	group() : impl_{MPI_GROUP_EMPTY}{}
	group(communicator const& comm){MPI_Comm_group(&comm, &impl_);}
//	template<class ContiguousIntIterator>
//	group(group const& other, ContiguousIntIterator ranks_begin, std::size_t n){
//		MPI_Group_incl(other.impl_, n, ranks_begin, &impl_);
//	}
	group(group const& other, std::vector<int> ranks){// : group(other, ranks.data(), ranks.size()){
		MPI_Group_incl(other.impl_, ranks.size(), ranks.data(), &impl_);
	}
	~group(){MPI_Group_free(&impl_);}
	int rank() const{
		int rank = -1; 
		int s = MPI_Group_rank(impl_, &rank);
		if(s != MPI_SUCCESS) throw std::runtime_error("rank not available");
		return rank;
	}
//	bool is_null() const{return MPI_GROUP_NULL == impl_;}
//	operator bool() const{return not is_null();}
	int size() const{
		int size = -1;
		MPI_Group_size(impl_, &size); 
		return size;
	}
	bool empty() const{
		return size() == 0;
	}

	mpi3::equality compare(group const& other) const{
		int result;
		MPI_Group_compare(impl_, other.impl_, &result);
		return static_cast<boost::mpi3::equality>(result);
	}
	friend mpi3::equality compare(group const& self, group const& other){return self.compare(other);}

	bool operator==(group const& other) const{return compare(other) == boost::mpi3::identical;}
	bool operator!=(group const& other) const{return not operator==(other);} 

	std::vector<int> translate_ranks(std::vector<int> const& ranks1, group const& other) const{
		std::vector<int> ranks2(ranks1.size());
		MPI_Group_translate_ranks(impl_, ranks1.size(), ranks1.data(), other.impl_, ranks2.data());
		return ranks2;
	}
	std::vector<int> translate_ranks(group const& other) const{
		std::vector<int> ranks1(other.size());
		for(int i = 0; i != (int)ranks1.size(); ++i) ranks1[i] = i;
	//	std::iota(ranks1.begin(), ranks1.end(), 0);
		return translate_ranks(ranks1, other);
	}
public:
	template<class Iterator, class Size>
	group include_n(Iterator it, Size count) const{
		return include_n_dispatch(typename std::iterator_traits<Iterator>::value_type{}, detail::iterator_category<Iterator>{}, it, count); 
	}
	template<class Iterator>
	group include(Iterator first, Iterator last) const{return include_n(first, std::distance(first, last));}
	template<class Range>
	group include(Range const& r) const{using std::begin; using std::end; return include(begin(r), end(r));}
private:
	template<class Iterator, class Size>
	group include_n_dispatch(int, detail::contiguous_iterator_tag, Iterator first, Size count) const{
		group ret;
		MPI_Group_incl(impl_, count, detail::data(first), &ret.impl_);
		return ret;
	}
	template<class Iterator, class Size, class Value, class = decltype(int(Value{}))>
	group include_n_dispatch(Value, std::input_iterator_tag, Iterator first, Size count) const{
		std::vector<int> v(count); std::copy_n(first, count, v.begin());
		return include_n(v.data(), v.size()); 
	}
public:
	group include(std::vector<int> const& ranks){return group(*this, ranks);}
	template<class ContiguousIntIterator, class Size>
	group exclude_n(ContiguousIntIterator it, Size count) const{
		group ret;
		MPI_Group_excl(impl_, count, boost::mpi3::detail::data(it), &ret.impl_);
		return ret;
	}
	template<class ContiguousIntIterator>
	group exclude(ContiguousIntIterator first, ContiguousIntIterator last) const{
		return exclude_n(first, std::distance(first, last));
	}
	template<class ContiguousRangeIterator, class Size>
	group range_exclude_n(ContiguousRangeIterator it, Size n) const{
		group ret;
		using int3 = int[3];
		auto ptr = const_cast<int3*>(reinterpret_cast<int3 const*>(boost::mpi3::detail::data(it)));
		int status = MPI_Group_range_excl(impl_, n, ptr, &ret.impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot exclude");
		return ret;
	}
	template<class ContiguousRangeIterator, class Size>
	group range_include_n(ContiguousRangeIterator it, Size n) const{
		group ret;
		using int3 = int[3];
		auto ptr = const_cast<int3*>(reinterpret_cast<int3 const*>(boost::mpi3::detail::data(it)));
		int status = MPI_Group_range_incl(impl_, n, ptr, &ret.impl_);
		if(status  != MPI_SUCCESS) throw std::runtime_error("cannot include");
		return ret;
	}

	template<class ContiguousRangeIterator>
	group range_exclude(ContiguousRangeIterator first, ContiguousRangeIterator last) const{
		return range_exclude_n(first, std::distance(first, last));
	}
	group range_exclude(std::vector<boost::mpi3::strided_range> const& ranges){
		return range_exclude(ranges.begin(), ranges.end());
	}
	template<class ContiguousRangeIterator>
	group range_include(ContiguousRangeIterator first, ContiguousRangeIterator last) const{
		return range_include_n(first, std::distance(first, last));
	}
	group range_include(std::vector<boost::mpi3::strided_range> const& ranges){
		return range_include(ranges.begin(), ranges.end());
	}
	group set_union(group const& other) const{
		group ret;
		MPI_Group_union(impl_, other.impl_, &ret.impl_); 
		return ret;
	}
	friend group set_union(group const& self, group const& other){
		return self.set_union(other);
	}
	group set_intersection(group const& other) const{
		group ret;
		MPI_Group_intersection(impl_, other.impl_, &ret.impl_);
		return ret;
	}
	group intersection(group const& other) const{return set_intersection(other);}
	friend group set_intersection(group const& self, group const& other){
		return self.set_intersection(other);
	}
	friend group intersection(group const& self, group const& other){
		return self.intersection(other);
	}
//	template<class Iterator, class Size>
//	group include_n(Iterator it, Size n) const{
//		group ret;
//		MPI_Group_incl(impl_, n, std::addressof(*it), &ret.impl_);
//		return ret;
//	}
#if 0
	template<class ContiguousIntIterator>
	group include(ContiguousIntIterator it1, ContiguousIntIterator it2) const{
		return include_n(it1, std::distance(it1, it2));
	}
	template<class Range>
	group include(Range const& r) const{
		using std::begin;
		using std::end;
		return include(begin(r), end(r));
	}
#endif
};
#endif

template<class Range>
auto operator/(Range const& r, communicator& self)
	->decltype(self.scatter(begin(r), end(r))){
		return self.scatter(begin(r), end(r));}

//inline communicator::communicator(communicator const& other, struct group const& g, int tag) : communicator(){
//	MPI_Comm_create_group(other.impl_, g.impl_, tag, &impl_);
//}

//package communicator::make_package(int n){return package(*this, n);}

/*
inline void communicator::send_value(package const& p, int dest, int tag){
	int s = MPI_Send(p.data(), p.buffer_.size(), MPI_PACKED, dest, tag, impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannotsend_value package");
}*/


/*
template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type>
void communicator::send_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag){
	package p(*this);
	detail::package_oarchive poa(p);
	while( n ){
		poa << *first;
		++first;
		--n;
	}
	p.send(dest, tag);
//	return first
}*/

/*
template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type>
void communicator::receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag){
//	assert(0);
	package p(*this);
	p.receive(dest, tag);
	detail::package_iarchive pia(p);
	while(n){
		pia >> *first;
		++first;
		--n;
	}
}*/

/*
template<class ContiguousIt, typename Size>
void communicator::broadcast_n_contiguous_builtinQ(std::false_type, ContiguousIt first, Size count, int root){
	package p(*this);
	if(rank() == root){
		detail::package_oarchive poa(p);
		while(count){
			poa << *first;
			++first;
			--count;
		}
	}
	p.broadcast(root);
	if(rank() != root){
		detail::package_iarchive pia(p);
		while(count){
			pia >> *first;
			++first;
			--count;
		}
	}
}
*/

}}

//BOOST_SERIALIZATION_REGISTER_ARCHIVE(boost::mpi3::package_oarchive)
//BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(boost::mpi3::detail::package_oarchive)
//BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(boost::mpi3::detail::package_iarchive)

#if not __INCLUDE_LEVEL__

#include "../mpi3/main.hpp"
#include "../mpi3/version.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

class V{
	mpi3::communicator comm_;
	public:
	V(mpi3::communicator const& c) : comm_(c){}
	V(mpi3::communicator&& c) : comm_(std::move(c)){}
};

int mpi3::main(int, char*[], mpi3::communicator world){
	std::cout << mpi3::undefined << std::endl;

	static_assert(std::is_nothrow_constructible<mpi3::communicator>::value, "MyType should be noexcept MoveConstructible");

//	auto worldcopy1 = world;
//	auto worldcopy2 = std::move(worldcopy1);
//	V v(worldcopy);
//	V v2(std::move(v));

	if(world.rank() == 0) cout << "MPI version " <<  mpi3::version() << '\n';
//	if(world.rank() == 0) cout << "Topology: " << name(world.topo()) << '\n';

	cout << "MPI_ERR_COMM = " << MPI_ERR_COMM << '\n';

	mpi3::communicator comm;
	assert(!comm);
//	cout << comm.rank() << '\n';

	mpi3::communicator comm2 = world;
	assert(comm2);
	assert(comm2.size() == world.size());
	assert(comm2 == world);
	assert(&comm2 != &world);

	mpi3::communicator comm3 = world;//.duplicate();
	assert(comm3);
	assert(comm3 == world);
	assert(&comm3 != &world);
	comm = comm2;
	assert(&comm != &comm2);

//	world2 = world;

	return 0;
#if 0
//	boost::mpi3::communicator newcomm = world;
	{
		int color = world.rank()/3;
		communicator row_comm;
		row_comm = world.split(color);
		world.barrier();
		std::cout << std::to_string(world.rank()) + " " + std::to_string(row_comm.rank()) + "\n";// << std::endl;
		world.barrier();
	}
	{
		communicator row_comm = world/3;
		world.barrier();
		std::cout << std::to_string(world.rank()) + " " + std::to_string(row_comm.rank()) + "\n";// << std::endl;
		world.barrier();
	}

	world.barrier();
	if(world.rank() == 0) cout << "prime communicator" << '\n';
	world.barrier();

	{
	//	group world_group(world);
	//	const int ranks[4] = {2, 3, 5, 7};
	//	group prime = world_group.include(ranks, ranks + 4);
	//	communicator prime_comm(world, prime);
		auto prime_comm = world.subcomm({2,3,5,7});
		cout << world.rank() << " -> " << prime_comm.rank() << "/" << prime_comm.size() << '\n';
#if 0
		if(communicator::null != prime_comm){
			cout << world.rank() << " -> " << prime_comm.rank() << "/" << prime_comm.size() << '\n';
		}else{
			cout << world.rank() << " not in prime comm\n";
		}
#endif
	}

	world.barrier();
	if(world.rank() == 0) cout << "prime communicator" << '\n';
	world.barrier();

	if(0){
		auto prime = world.subcomm({2,3,5,7});
		if(prime.is_empty()){
	//	if (communicator::null != prime){
			cout << world.rank() << " -> " << prime.rank() << "/" << prime.size() << '\n';
		}else{
			cout << world.rank() << " not in prime comm\n";
		}
	}
#endif
}

#endif
#endif


