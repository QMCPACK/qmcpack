/* -*- indent-tabs-mode: t -*- */

#ifndef BOOST_MPI3_INFO_HPP
#define BOOST_MPI3_INFO_HPP

#include "../mpi3/handle.hpp"

#include <algorithm>  // for std::for_each
#include <iostream>
#include <tuple>  // for std::tie

namespace boost {
namespace mpi3 {

struct info :
	detail::regular_handle<
		info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free
	>
{
	using base = detail::regular_handle<info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free>;
	using detail::regular_handle<info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free>::call;

	info() = default;
	info(info const&) = default;
	info(info     &&) = delete;  // TODO(correaa) consider relation with default constructor

	info& operator=(info const&) = default;
	info& operator=(info     &&) = delete;  // TODO(correaa) consider relation with default constructor

	~info() = default;

	// cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3, initialize_list ctor must be implicit
	info(std::initializer_list<std::pair<std::string, std::string>> il) {
		std::for_each(il.begin(), il.end(), [this](auto const& e) {set(e.first, e.second);});
	//  for(auto const& e : il) {set(e.first, e.second);}
	}

//  void                        delete_(std::string const& key) {call<&MPI_Info_delete>(key);}
	std::pair<std::string, int> get(std::string const& key, int valuelen) const{return base::call<&MPI_Info_get>(key, valuelen);}
	int                         get_nkeys() const{return call<&MPI_Info_get_nkeys>();}
	std::string                 get_nthkey(int n) const{return call<&MPI_Info_get_nthkey>(n);}
	std::pair<int, int>         get_valuelen(std::string const& key) const{return call<&MPI_Info_get_valuelen>(key);}
	void                        set(std::string const& key, std::string const& value){call<&MPI_Info_set>(key, value);}

	void insert(std::string const& key, std::string const& value){return set(key, value);}
	void erase(std::string const& key) {call<&MPI_Info_delete>(key);}
	int size() const{return get_nkeys();}
	std::string operator[](std::string const& key) const{
		int valuelen = MPI_MAX_INFO_VAL;
		int flag;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		std::tie(valuelen, flag) = get_valuelen(key);
		if(flag == 0) {throw std::logic_error{"key `" + key + "' not found in info set"};}
		return get(key, valuelen).first;
	};
	std::pair<std::string, std::string> operator[](int n) const{
		auto key = call<&MPI_Info_get_nthkey>(n);
		return {key, operator[](key)};
	}
	friend std::ostream& operator<<(std::ostream& os, info const& self) {
		for(int i = 0; i != self.get_nkeys(); ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
			auto p = self[i];
			os<< p.first <<" : "<< p.second <<'\n';
		}
		return os;
	}
};

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_INFO

//#include "../mpi3/main.hpp"
//#include<iostream>

//using std::cout;
//using std::endl;

//int boost::mpi3::main(int, char*[], boost::mpi3::communicator world){
//	if(world.rank() == 0){
//		boost::mpi3::info nfo;
//		nfo.set("file", "runfile.txt");
//		nfo.set("soft", "host");
//		cout << nfo.get_nkeys() << '\n';
//		cout << nfo << '\n';
//		nfo.delete_("soft");
//		cout << nfo << '\n';
//		assert( nfo["file"] == "runfile.txt" );
//		boost::mpi3::info nfo2 = nfo;
//		cout << nfo2 << '\n';
//	}
//	return 0;
//}

//#endif
#endif
