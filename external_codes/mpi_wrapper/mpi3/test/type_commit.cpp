// #define BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/type.hpp"

#include "../../mpi3/detail/tuple_offset.hpp"

#include<tuple>

namespace mpi3 = boost::mpi3;

//template<class T>
//auto gen_type(T const&);

//template<> auto gen_type<int >(int  const&) {return type{MPI_INT };}
//template<> auto gen_type<char>(char const&) {return type{MPI_CHAR};}

//template<class... Ts, std::size_t... Is>
//auto gen_type_aux<std::tuple<T...>>(std::tuple<Ts...> const& t, std::index_sequence<Is...>) {
//	std::vector<int> blocklen = { (il.size(), 1);
//	std::vector<MPI_Aint> disp;
//	return type::struct_({get_type(std::get<Is>(t), ...});
//}

//template<class... Ts>
//auto gen_type<std::tuple<T...>>(std::tuple<Ts...> const& t) {
//	return gen_type_aux(t, std::index_sequence_for<Ts...>{});
//}

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator /*world*/) try {

	{
		using Tuple = std::tuple<int, double, int, char, float, long double>;
		Tuple tup;
		auto offset4 = mpi3::detail::element_offset<4, Tuple>();
		assert( reinterpret_cast<char*>(&tup) + offset4 == reinterpret_cast<char*>(&std::get<4>(tup)) );  // NOLINT(cert-dcl03-c,hicpp-static-assert,misc-static-assert,cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-bounds-pointer-arithmetic) for some compiler this is not a constexpr
	}

	mpi3::type t = mpi3::int_[100]; // mpi3::type::int_.contiguous(100);


#if 0
	t.commit_as<std::array<int, 100>>();
	t.commit_as<int[100]>();
//	std::array<int, 100> buffer;
	int buffer[100];

	if(world.rank() == 0) world.send_n(&buffer, 1, 1, 123);
	else if(world.rank() == 1) world.receive_n(&buffer, 1, 0, 123);
#endif

	return 0;
} catch(...) {return 1;}

