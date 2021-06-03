#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 `#-Wfatal-errors` -D_TEST_BOOST_MPI3_INFO $0x.cpp -o $0x.x && time mpirun -np 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_INFO_HPP
#define BOOST_MPI3_INFO_HPP

#include "../mpi3/handle.hpp"

#include <tuple> // tie
#include<iostream>

namespace boost{
namespace mpi3{

struct info : 
	detail::regular_handle<
		info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free
	>
{
	using base = detail::regular_handle<info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free>;
	using detail::regular_handle<info, MPI_Info, MPI_Info_create, MPI_Info_dup, MPI_Info_free>::call;

	info() = default;
	info(info const& other) = default;
	info& operator=(info const&) = default;
	~info() = default;

	info(std::initializer_list<std::pair<std::string, std::string>> il) : info(){
		for(auto& e : il) set(e.first, e.second);
	}

	void                        delete_(std::string const& key){call<&MPI_Info_delete>(key);}
	std::pair<std::string, int> get(std::string const& key, int valuelen) const{return base::call<&MPI_Info_get>(key, valuelen);}
	int                         get_nkeys() const{return call<&MPI_Info_get_nkeys>();}
	std::string                 get_nthkey(int n) const{return call<&MPI_Info_get_nthkey>(n);}
	std::pair<int, int>         get_valuelen(std::string const& key) const{return call<&MPI_Info_get_valuelen>(key);}
	void                        set(std::string const& key, std::string const& value){call<&MPI_Info_set>(key, value);}

	void insert(std::string const& key, std::string const& value){return set(key, value);}
	void erase(std::string const& key){delete_(key);}
	int size() const{return get_nkeys();}
	std::string operator[](std::string const& key) const{
		int valuelen = MPI_MAX_INFO_VAL;
		int flag;
		std::tie(valuelen, flag) = get_valuelen(key);
		if(flag == false) throw std::logic_error("key `" + key + "' not found in info set");
		return get(key, valuelen).first;
	};
	std::pair<std::string, std::string> operator[](int n) const{
		auto key = call<&MPI_Info_get_nthkey>(n);
		return {key, operator[](key)};
	}
	friend std::ostream& operator<<(std::ostream& os, info const& self){
		for(int i = 0; i != self.get_nkeys(); ++i){
			auto p = self[i];
			os << p.first << " : " << p.second << '\n';
		}
		return os;
	}
};

}}

#ifdef _TEST_BOOST_MPI3_INFO

#include "../mpi3/main.hpp"
#include<iostream>

using std::cout;
using std::endl;

int boost::mpi3::main(int, char*[], boost::mpi3::communicator world){
	if(world.rank() == 0){
		boost::mpi3::info nfo;
		nfo.set("file", "runfile.txt");
		nfo.set("soft", "host");
		cout << nfo.get_nkeys() << '\n';
		cout << nfo << '\n';
		nfo.delete_("soft");
		cout << nfo << '\n';
		assert( nfo["file"] == "runfile.txt" );
		boost::mpi3::info nfo2 = nfo;
		cout << nfo2 << '\n';
	}
	return 0;
}

#endif
#endif

