// © Alfredo A. Correa 2019-2021
#ifndef BOOST_MPI3_DETAIL_CALL_HPP
#define BOOST_MPI3_DETAIL_CALL_HPP

#include "../error.hpp"
#include "../status.hpp"

#include "../config/NODISCARD.hpp"

#include<mpi.h>  // MPI_MAX_PROCESSOR_NAME

#include<string>

namespace boost {
namespace mpi3 {
namespace detail {

template<int(*F)(int*)>
int call() {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
	auto const e = static_cast<enum error>((*F)(&ret));
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}	
	return ret;
}

template<int(*F)(char*, int*)>
std::string call() {
	int len = -1;
	std::array<char, MPI_MAX_PROCESSOR_NAME> name{};
	auto const e = static_cast<enum error>((*F)(name.data(), &len));
	assert(len >= 0);
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}
	return {name.data(), static_cast<std::size_t>(len)};
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()...)))* = nullptr>
void call(Args... args) {
	auto const e = static_cast<enum error>((*F)(args...));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<MPI_Status*>())))* = nullptr>
mpi3::status call(Args... args) {
	mpi3::status ret;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret.impl_));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}
	return ret;
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<int*>(), std::declval<int*>())))* = nullptr>
[[nodiscard]] auto call(Args... args) {
	std::pair<int, int> ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret.first, &ret.second));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}
	// cppcheck-suppress uninitvar
	return ret;
}

template<class FT, FT* F, class... Args, decltype(static_cast<enum error>((*F)(std::declval<Args>()..., std::declval<int*>())))* = nullptr>
[[nodiscard]] int call(Args... args) {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	auto const e = static_cast<enum error>((*F)(args..., &ret));  // NOLINT(clang-analyzer-optin.mpi.MPI-Checker) // non-blocking calls have wait in request destructor
	if(e != mpi3::error::success) {throw std::system_error{e, "cannot call function " + std::string{__PRETTY_FUNCTION__}};}
	// cppcheck-suppress uninitvar
	return ret;
}

#define MPI3_CALL(F) detail::call<decltype(F), F>  // NOLINT(cppcoreguidelines-macro-usage)
#define MPI_(F) MPI3_CALL(MPI_##F)  // NOLINT(cppcoreguidelines-macro-usage): name concatenation

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

#endif
