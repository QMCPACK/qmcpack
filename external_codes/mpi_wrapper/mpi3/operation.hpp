#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -D_TEST_BOOST_MPI3_OPERATION $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_OPERATION_HPP
#define BOOST_MPI3_OPERATION_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include "../mpi3/detail/datatype.hpp"
#include "../mpi3/handle.hpp"

#include<utility> // forward

namespace boost{
namespace mpi3{

struct commutative_operation;
struct builtin_operation;

struct operation : detail::nondefault_handle<operation, MPI_Op, MPI_Op_free>{ // user_operation, operator_ . operator is a C++ keyword
	using base = detail::nondefault_handle<operation, MPI_Op, MPI_Op_free>;
	using detail::nondefault_handle<operation, MPI_Op, MPI_Op_free>::nondefault_handle;
	public:
	template<class F, typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation>{}> >
	operation(F&& f, bool commutative) : operation(detail::uninitialized{}){
		MPI_Op_create(
			&f,
		//	reinterpret_cast<void (*)(void*, void*, int*, int*)>(&f),
			commutative, 
			&impl_
		);
	}

	operation() = delete;
	operation(operation const&) = delete;
	operation& operator=(operation const&) = delete;
	~operation() = default;

#if 0
	enum struct code : MPI_Op{
		maximum = MPI_MAX, minimum = MPI_MIN, 
		sum = MPI_SUM, product = MPI_PROD, 
		logical_and = MPI_LAND, bitwise_and = MPI_BAND, 
		logical_or  = MPI_LOR,   bitwise_or = MPI_BOR,
		logical_xor = MPI_LXOR, bitwise_xor = MPI_BXOR,
		max_value_location = MPI_MAXLOC,
		min_value_location = MPI_MINLOC
	};
#endif
//	operation(operation::code c) : base((MPI_Op)c){}

//	static operation const sum;//(operation::code::sum);
//	static operation const product;
//	static operation const maximum;
//	static operation const minimum;

};

//operation const sum    (operation::code::sum);
//operation const product(operation::code::product);
//operation const maximum(operation::code::maximum);
//operation const minimum(operation::code::minimum);

template<class T = void>
using plus = std::plus<T>;
template<class T = void>
using minus = std::minus<T>;
template<class T = void>
using multiplies = std::multiplies<T>;

template<class T = void> struct min{
	T const& operator()(T const& t1, T const& t2) const{return std::min(t1, t2);}
};
template<> struct min<void>{
	template<class T1, class T2> decltype(auto) operator()(T1&& t1, T2&& t2) const{return std::min(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<class T = void> struct max{
	T const& operator()(T const& t1, T const& t2) const{return std::max(t1, t2);}
};
template<> struct max<void>{
	template<class T1, class T2> decltype(auto) operator()(T1&& t1, T2&& t2) const{return std::max(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<class Op> struct predefined_operation;

#define BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(CppoP, MpinamE, NamE) \
template<> struct predefined_operation<CppoP>{ \
/*	constexpr*/ operator MPI_Op() const{return MpinamE;} \
/*	static constexpr MPI_Op value = MpinamE;*/ \
}; \
using NamE = predefined_operation<CppoP>

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::plus<>       , MPI_SUM , sum        );
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::multiplies<> , MPI_PROD, product    );
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::logical_and<>, MPI_LAND, logical_and);

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(max<>, MPI_MAX, maximum);
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(min<>, MPI_MIN, minimum);

#undef BOOST_MPI3_DECLARE_PREDEFINED_OPERATION

struct commutative_operation : operation{
	template<class F,  typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation>{}> >
	commutative_operation(F&& f) : operation(std::forward<F>(f), true){}
};

struct non_commutative_operation : operation{
	template<class F,  typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation>{}>>
	non_commutative_operation(F&& f) : operation(std::forward<F>(f), false){}
};

}}

#ifdef _TEST_BOOST_MPI3_OPERATION

#include "../mpi3/main.hpp"
#include "../mpi3/error_handler.hpp"

void addem_int(int const* invec, int *inoutvec, int *len, int* f){
	for(int i=0; i<*len; i++) inoutvec[i] += invec[i];
}

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	int correct_result = world.size()*(world.size()-1)/2;

	int data = world.rank();
	{
		int result = -1;
		world.reduce_n(&data, 1, &result, std::plus<>{}, 0);
		if(world.root()) assert( result == correct_result ); else assert( result == -1 );
		world.broadcast_n(&result, 1, 0);
		assert(result == correct_result);
	}
	{
		int result = -1;
		world.all_reduce_n(&data, 1, &result, std::plus<>{});
		assert(result == correct_result);
	}
	{
	//	int result = world.all_reduce_value<std::plus<>>(data);
	//	assert( result == correct_result );
	}
	{
	//	int result = world.all_reduce_value(world.rank(), std::plus<>{});
	//	assert( result == correct_result );
	}

	return 0;
}

#endif
#endif

