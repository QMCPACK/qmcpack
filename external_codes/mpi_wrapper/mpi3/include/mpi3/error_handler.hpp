// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2017-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_ERROR_HANDLER_HPP
#define BOOST_MPI3_ERROR_HANDLER_HPP

#include<mpi.h>

#include "../mpi3/communicator.hpp"

namespace boost {
namespace mpi3 {

class fatal {
//	explicit operator MPI_Errhandler() const {return MPI_ERRORS_ARE_FATAL;}
};

class code {
//	explicit operator MPI_Errhandler() const {return MPI_ERRORS_RETURN;}
};

struct error_handler {
	MPI_Errhandler impl_ = MPI_ERRORS_ARE_FATAL;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa)

	explicit constexpr error_handler(MPI_Errhandler impl) noexcept : impl_{impl} {}

	error_handler() = default;
	error_handler(error_handler const&) = delete;
	friend mpi3::communicator;

	error_handler& operator=(error_handler const&) = delete;
	error_handler& operator=(error_handler     &&) = delete;

 private:
	error_handler(error_handler&&) = default;  // this is necessary in C++14 (to return from function)

 public:
//	error_handler(void(*fn)(MPI_Comm*, int* err, ...)){
//		MPI_Comm_create_errhandler(fn, &impl_);
//	}
//	void operator()(communicator& comm, int error) const{comm.call_error_handler(error);}
	~error_handler() {
		if(impl_ != MPI_ERRORS_ARE_FATAL and impl_ != MPI_ERRORS_RETURN) {
			MPI_Errhandler_free(&impl_);
		}
	}
//	static error_handler const exception;
	static void exception(MPI_Comm* /*comm*/, int const* err) {//, ...){
		int len;  // NOLINT(cppcoreguidelines-init-variables,-warnings-as-errors) delayed init
		std::array<char, MPI_MAX_ERROR_STRING> estring{};
		MPI_Error_string(*err, estring.data(), &len);
		std::string const w(estring.data(), static_cast<std::string::size_type>(len));
		throw std::runtime_error{"error code"+ std::to_string(*err) +" "+ w};
//		throw boost::mpi3::exception("error code " + std::to_string(*err) + " from comm " + std::to_string(*comm) + ": " + w);
//		throw std::runtime_error("error code " + std::to_string(*err) + " from comm " + std::to_string(*comm) + ": " + w);
	}

	static error_handler const fatal;
	static error_handler const code;
};

error_handler const error_handler::fatal{MPI_ERRORS_ARE_FATAL};  // NOLINT(misc-definitions-in-headers,fuchsia-statically-constructed-objects) TODO(correaa)
error_handler const error_handler::code{MPI_ERRORS_RETURN};  // NOLINT(misc-definitions-in-headers,fuchsia-statically-constructed-objects) TODO(correaa)

inline void communicator::set_error_handler(error_handler const& eh) {
	MPI_(Comm_set_errhandler)(impl_, eh.impl_);
}

inline error_handler communicator::get_error_handler() const {
	error_handler ret;
	MPI_(Comm_get_errhandler)(impl_, &ret.impl_);
	return ret;
}

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_ERROR_HANDLER

//#include "../mpi3/main.hpp"

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) {

////	error_counter_t ec;
////	mpi3::error_handler<&ehh> ehherr;
////	mpi3::error_handler newerr(eh);
////	world.set_error_handler(ec);
////	world.set_error_handler(newerr);
////	world.set_error_handler(
////	world.error(error::other);
////	cout << ec.count() << '\n';
////	newerr(world, MPI_ERR_OTHER);

////	auto f = world.get_error_handler();

//	return 0;
//}

//#endif
#endif
