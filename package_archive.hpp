#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo '#include"'$0'"'>$0.cpp)&&mpic++ -std=c++14 -O3 -Wall -Wextra -Wfatal-errors -D_TEST_MPI3_PACKAGE_ARCHIVE $0.cpp -o$0x -lboost_serialization&&mpirun -n 2 $0x && rm $0x $0.cpp; exit
#endif

#ifndef MPI3_PACKAGE_ARCHIVE_HPP
#define MPI3_PACKAGE_ARCHIVE_HPP

#include"../mpi3/detail/package.hpp"

#include <boost/archive/detail/common_iarchive.hpp>
#include <boost/archive/detail/common_oarchive.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/item_version_type.hpp>

#include <boost/version.hpp>

#include <sstream>

namespace boost{
namespace mpi3{

namespace detail{

class basic_package_iprimitive{
protected:
	package& p_;
public:
	// we provide an optimized save for all basic (and fundamental) types
	// this note is taken from boost serialization:
	// typedef serialization::is_bitwise_serializable<mpl::_1> 
	// use_array_optimization;
	// workaround without using mpl lambdas
	struct use_array_optimization {
		template<class T>  
		struct apply : 
			public mpl::bool_<mpi3::detail::is_basic<T>::value>
		{};  
	};
	template<class T>
#if(BOOST_VERSION < 106100)
	void load_array(boost::serialization::array<T>& t, unsigned int = 0){ // for boost pre 1.63
#else
	void load_array(boost::serialization::array_wrapper<T>& t, unsigned int = 0){
#endif
		p_.unpack_n(t.address(), t.count()); 
	}
    template<class T>
    void load(T& t){p_ >> t;}
	basic_package_iprimitive(mpi3::detail::package& p) : p_(p){}
};

class basic_package_oprimitive{
protected:
	package& p_;
public:
	struct use_array_optimization {
		template <class T>
		struct apply : public boost::serialization::is_bitwise_serializable<T>{};  
	};
	template<class T>
	void save(const T& t, unsigned int = 0){p_.pack_n(std::addressof(t), 1);}//p_ << t;}
	template<class T>
#if(BOOST_VERSION < 106100)
	void save_array(boost::serialization::array<T> const& t, unsigned int = 0){
#else
	void save_array(boost::serialization::array_wrapper<T> const& t, unsigned int = 0){
#endif
		p_.pack_n(t.address(), t.count()); 
	}
#if 0
	void save(const boost::archive::object_id_type&){}
	void save(const boost::archive::object_reference_type&){}
	void save(const boost::archive::class_id_type&){}
	void save(const boost::archive::class_id_optional_type&){}
	basic_memory_oprimitive(size_t& os) : os_(os){}
#endif
	basic_package_oprimitive(mpi3::detail::package& p) : p_(p){}
};

template<class Archive>
class basic_package_iarchive : public boost::archive::detail::common_iarchive<Archive>{
	friend class boost::archive::detail::interface_iarchive<Archive>;
	typedef boost::archive::detail::common_iarchive<Archive> detail_common_iarchive;
	template<class T>
	void load_override(T& t, /*BOOST_PFTO*/ int = 0){
#if(BOOST_VERSION < 105900)
		this->detail_common_iarchive::load_override(t, 0);
#else
		this->detail_common_iarchive::load_override(t);//, 0);
#endif
	}
	protected:
	basic_package_iarchive(unsigned int flags) : boost::archive::detail::common_iarchive<Archive>(flags){}
};

template<class Archive>
class basic_package_oarchive : public boost::archive::detail::common_oarchive<Archive>{
	friend class boost::archive::detail::interface_oarchive<Archive>;
	typedef boost::archive::detail::common_oarchive<Archive> detail_common_oarchive;
protected:
	template<class T>
	void save_override(T& t, /*BOOST_PFTO*/ int = 0){
#if(BOOST_VERSION < 105900)
	  this->detail_common_oarchive::save_override(t, 0);//, 0);
#else
	  this->detail_common_oarchive::save_override(t);
#endif
	}
#if 0
	void save_override(const object_id_type&, int){/* this->This()->newline(); this->detail_common_oarchive::save_override(t, 0);*/}
	void save_override(const class_id_optional_type&, int){}
	void save_override(const class_name_type&, int){/*  const std::string s(t); * this->This() << s;*/}
#endif
	protected:
	basic_package_oarchive(unsigned int flags) : 
		boost::archive::detail::common_oarchive<Archive>(flags)
	{}
};

template<class Archive>
class package_iarchive_impl : 
	public basic_package_iprimitive, 
	public basic_package_iarchive<Archive>
{
	public:
	template<class T>
	void load(T& t){basic_package_iprimitive::load(t);}
// empty functions follow, so that metadata is communicated
	void load(boost::archive::version_type&){}
//	void save(const boost::serialization::item_version_type&){/*save(static_cast<const unsigned int>(t));*/}
	void load(boost::archive::tracking_type&){/*save(static_cast<const unsigned int>(t));*/}
	void load(boost::archive::object_id_type&){}
	void load(boost::archive::object_reference_type&){}
	void load(boost::archive::class_id_type&){}
	void load(boost::archive::class_id_optional_type&){}
	void load(boost::archive::class_id_reference_type&){}
	void load(boost::archive::class_name_type&){}

	void load(boost::serialization::collection_size_type& t){
		unsigned int x = 0;
		load(x);
		t = serialization::collection_size_type(x);
	}
	void load(boost::serialization::item_version_type&){}

	void load(char* s){
		assert(0);
		const std::size_t len = std::ostream::traits_type::length(s);
		*this->This() << len;
		p_.pack_n(s, len);
	}
	void load(wchar_t * ws){
		const std::size_t l = std::wcslen(ws);
		*this->This() << l;
		assert(0);
	}
	void load(std::string &s){
		std::size_t size; //  *this->This() >> size;
		p_.unpack_n(&size, 1);
		s.resize(size);
		p_.unpack_n(const_cast<char*>(s.c_str()), size);
	}
	void load(std::wstring &ws){
    	const std::size_t size = ws.size();
		*this->This() << size;
	//	++tokens_; //	this->This()->newtoken();
	//	os_ += ws.size()*sizeof(wchar_t);//	os << s;
		assert(0);
	}
	public:
	package_iarchive_impl(mpi3::detail::package& p, unsigned int flags) : // size_t& os, size_t& tokens, unsigned int flags) :
		basic_package_iprimitive(p),
		basic_package_iarchive<Archive>(flags)
	{}
};

template<class Archive>
class package_oarchive_impl : public basic_package_oprimitive, public basic_package_oarchive<Archive>{
public:
	template<class T>
	void save(const T& t){basic_package_oprimitive::save(t);}
#if(BOOST_VERSION < 106100)
	void save(boost::serialization::array<double>&){
#else
	void save(boost::serialization::array_wrapper<double>&){
#endif
		assert(0);
	}
    void save(const boost::archive::version_type&){}
//	void save(const boost::serialization::item_version_type&){/*save(static_cast<const unsigned int>(t));*/}
    void save(const boost::archive::tracking_type&){/*save(static_cast<const unsigned int>(t));*/}
	void save(const boost::archive::object_id_type&){}
	void save(const boost::archive::object_reference_type&){}
	void save(const boost::archive::class_id_type&){}
	void save(const boost::archive::class_id_optional_type&){}
	void save(const boost::archive::class_id_reference_type&){}
	void save(const boost::archive::class_name_type&){}

	void save(const boost::serialization::collection_size_type& t){
		save(static_cast<      unsigned int>(t));
//		save(static_cast<const unsigned int>(t));
	}
	void save(const boost::serialization::item_version_type&){}

	// string types (like char*, string, etc) have special handling
	// types that need special handling
	void save(const char * s){
		assert(0);
		const std::size_t len = std::ostream::traits_type::length(s);
		*this->This() << len;
	//	++tokens_;//	this->This()->newtoken();
	//	os_ += len*sizeof(char);//	os << s;
		p_.pack_n(s, len);
	}
	void save(const wchar_t * ws){
		const std::size_t l = std::wcslen(ws);
		*this->This() << l;
		assert(0);
	//	++tokens_; // this->This()->newtoken();
	//	os_ += l*sizeof(wchar_t);//	os.write((const char *)ws, l * sizeof(wchar_t)/sizeof(char));
	}
	void save(std::string const& s){
    	const std::size_t size = s.size();
	//	*this->This() << size;
		p_.pack_n(&size, 1);
	//	std::cout << " packed size = " << size << '\n';
	//	++tokens_; // this->This()->newtoken();
	//	os_ += s.size()*sizeof(char);//	os << s;
		p_.pack_n(s.c_str(), size);
	}
	void save(const std::wstring &ws){
    	const std::size_t size = ws.size();
		*this->This() << size;
	//	++tokens_; //	this->This()->newtoken();
	//	os_ += ws.size()*sizeof(wchar_t);//	os << s;
		assert(0);
	}
//	using package_oarchive_impl<package_oarchive>::save_override;
#if 1
	// Save all supported datatypes directly
	template<class T>
#if(BOOST_VERSION < 106100)
	void save(boost::serialization::array<T> const& t, unsigned int){
#else
	void save(boost::serialization::array_wrapper<T> const& t, unsigned int){
#endif
		assert(0);
		save_override(t, boost::mpl::bool_<true>{});//std::true_type{});
	}
#endif
	public:
    package_oarchive_impl(mpi3::detail::package& p, unsigned int flags) : // size_t& os, size_t& tokens, unsigned int flags) :
		basic_package_oprimitive(p),
		basic_package_oarchive<Archive>(flags)
	{}
};

} // boost::mpi3::detail

struct package_iarchive : public detail::package_iarchive_impl<package_iarchive>{
    package_iarchive(mpi3::detail::package& p, unsigned int flags = 0) 
    : package_iarchive_impl<package_iarchive>(p, flags){}
};

struct package_oarchive : public detail::package_oarchive_impl<package_oarchive>{
	package_oarchive(mpi3::detail::package& p, unsigned int flags = 0) : 
		package_oarchive_impl<package_oarchive>(p, flags)
	{}
	using package_oarchive_impl<package_oarchive>::operator&;
#if(BOOST_VERSION < 106100)
	package_oarchive& operator & (boost::serialization::array<double>&){
#else
	package_oarchive& operator & (boost::serialization::array_wrapper<double>&){
#endif
		assert(0);
		return *this;
	}
};

}}

// maybe needed for optimization to take effect?
// BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(boost::archive::package_oarchive)

#ifdef _TEST_MPI3_PACKAGE_ARCHIVE

#include "../mpi3/main.hpp"
#include "../mpi3/process.hpp"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

namespace mpi3 = boost::mpi3;
using std::cout; 

int mpi3::main(int, char*[], mpi3::communicator world){
	assert(world.size() > 1);
	switch(world.rank()){
		case 0: {
			mpi3::detail::package p(world);
			mpi3::package_oarchive poa(p);
			std::string s("hello");
			int 
				i = 12, 
				j = 13
			;
			std::vector<double> v(20, 5.);
			std::map<int, int> m = {{1,2},{2,4},{3,4}};
			poa 
				<< s 
				<< i 
				<< j
				<< v
				<< 5
				<< m
			;
			p.send(1);
		} break;
		case 1: {
			mpi3::detail::package p(world);
			mpi3::package_iarchive pia(p);
			p.receive(0);
			std::string s;
			int 
				i,
				j
			;
			std::vector<double> v;
			int c;
			std::map<int, int> m;
			pia 
				>> s
				>> i
				>> j
				>> v
				>> c
				>> m
			;
			assert( s == "hello" );
			assert( i == 12 );
			assert( j == 13 );
			assert(v.size() == 20);
			assert(c == 5);
			assert( m[3] == 4 );
		}
	}
	return 0;
}
#endif
#endif

