// Copyright 2019-2025 Alfredo A. Correa

#ifndef BOOST_MPI3_DETAIL_CALL_HPP
#define BOOST_MPI3_DETAIL_CALL_HPP

#include <mpi3/error.hpp>
#include <mpi3/status.hpp>

#include <mpi3/config/NODISCARD.hpp>

#include <mpi3/detail/mpi_impl.h>

#include<string>

#ifdef __PRETTY_FUNCTION__
#define BOOST_MPI3_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#define BOOST_MPI3_PRETTY_FUNCTION "some_MPI_function"  // NOLINT(cppcoreguidelines-macro-usage)
#endif

namespace boost {
namespace mpi3 {
namespace detail {

template<int(*F)(int*)>
int call() {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
	auto const e = static_cast<enum error>((*F)(&ret));
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};} 
	return ret;
}

template<int(*F)(char*, int*)>
std::string call() {
	int len = -1;
	std::string name(MPI_MAX_PROCESSOR_NAME, '\0');  // std::array<char, MPI_MAX_PROCESSOR_NAME> name{};
	auto const e = static_cast<enum error>((*F)(name.data(), &len));
	assert(len >= 0);
	name.resize(static_cast<std::size_t>(len));
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};}
	return name;
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()...)))* = nullptr>
void call(Args... args) {
	// NOLINTNEXTLINE(bugprone-multi-level-implicit-pointer-conversion)
	auto const e = static_cast<enum error>((*F)(args...));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};}
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<MPI_Status*>())))* = nullptr>
mpi3::status call(Args... args) {
	mpi3::status ret;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret.impl_));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};}
	return ret;
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<int*>(), std::declval<int*>())))* = nullptr>
[[nodiscard]] auto call(Args... args) {
	std::pair<int, int> ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret.first, &ret.second));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};}
	// cppcheck-suppress uninitvar
	return ret;
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<int*>())))* = nullptr>
[[nodiscard]] int call(Args... args) {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{BOOST_MPI3_PRETTY_FUNCTION}};}
	// cppcheck-suppress uninitvar
	return ret;
}

#define MPI3_CALL(F) detail::call<decltype(F), F>  // NOLINT(cppcoreguidelines-macro-usage)
#define MPI_(F) MPI3_CALL(MPI_##F)  // NOLINT(cppcoreguidelines-macro-usage): name concatenation

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

#undef BOOST_MPI3_PRETTY_FUNCTION

#endif
