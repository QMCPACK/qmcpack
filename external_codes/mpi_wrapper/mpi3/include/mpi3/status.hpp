/* -*- indent-tabs-mode: t -*- */
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_STATUS_HPP
#define BOOST_MPI3_STATUS_HPP

//#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "../mpi3/detail/datatype.hpp"

#include<stdexcept>

namespace boost {
namespace mpi3 {

struct [[nodiscard]] status {
	MPI_Status impl_;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa)

	status() = default;
	status(status const&) = default; //  TODO(correaa) check
	status(status     &&) = default; //  TODO(correaa) check

	status& operator=(status const&) = default;
	status& operator=(status     &&) = default;

	~status() noexcept = default; //  TODO(correaa) use MPI_Status_free
	//{
	//	if(impl_ != MPI_STATUS_NULL) 
	//	MPI_Status_free(&impl_);
	//}

	template<class T>  // = char>
	int count() const {  // entries received of datatype T
		int ret = -1;
		MPI_Get_count(&impl_, datatype<T>{}(), &ret);  // can't use MPI_(Get_count)
		return ret;
	}

	template<class T>
	int elements() const {
		int ret = -1;
		int const s = MPI_Get_elements(&impl_, datatype<T>{}(), &ret);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {throw std::runtime_error{"cannot elements"};}
		return ret;
	}
	template<class T>
	void set_elements(int count) {
		int const s = MPI_Status_set_elements(&impl_, datatype<T>{}(), count);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {throw std::runtime_error{"cannot set elements"};}
	}

	int source() const {return impl_.MPI_SOURCE;}
	void set_source(int s) {impl_.MPI_SOURCE = s;}

	int tag() const {return impl_.MPI_TAG;}
	void set_tag(int t) {impl_.MPI_TAG = t;}

	void set_cancelled(bool flag = true) {
		int const s = MPI_Status_set_cancelled(&impl_, flag?1:0);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {throw std::runtime_error{"cannot set cancelled"};}
	}
//	bool cancelled() const {
//		int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
//		int s = MPI_Test_cancelled(&impl_, &ret);
//		if(s != MPI_SUCCESS) {throw std::runtime_error{"cannot test cancelled"};}
//		return ret != 0;
//	}
//	constexpr static auto const ignore = MPI_STATUS_IGNORE;
};

}  // end namespace mpi3
}  // end namespace boost
#endif
