/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// archive_exception.cpp:

// (C) Copyright 2009 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
// NOLINTBEGIN(altera-unroll-loops,cppcoreguidelines-pro-type-member-init,hicpp-member-init,misc-const-correctness) external code
#if (defined _MSC_VER) && (_MSC_VER == 1200)
#  pragma warning (disable : 4786) // too long name, harmless warning
#endif

#include <cstring>
#include <exception>
#include <string>

#define BOOST_ARCHIVE_SOURCE
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/config.hpp>
#include <boost/version.hpp>

namespace boost {
namespace archive {

#if BOOST_VERSION < 105900
#define BOOST_NOEXCEPT
#endif

BOOST_ARCHIVE_DECL
unsigned int
archive_exception::append(unsigned int l, const char * a){
    while(l < (sizeof(m_buffer) - 1)){
        char c = *a++;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) external code
        if('\0' == c) {
            break; }
        m_buffer[l++] = c;  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index) external code
    }
    m_buffer[l] = '\0';  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index) external code
    return l;
}

BOOST_ARCHIVE_DECL
archive_exception::archive_exception(  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) external code
    exception_code c, 
    const char * e1,
    const char * e2
) BOOST_NOEXCEPT :
    code(c)
{
    unsigned int length = 0;
    switch(code){
    case no_exception:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "uninitialized exception");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case unregistered_class:
        length = append(length, "unregistered class");
        if(nullptr != e1){
            length = append(length, " - ");
            // cppcheck-suppress [unreadVariable] ; external code
            length = append(length, e1);  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        }    
        break;
    case invalid_signature:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "invalid signature");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case unsupported_version:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "unsupported version");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case pointer_conflict:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "pointer conflict");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case incompatible_native_format:
        length = append(length, "incompatible native format");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        if(nullptr != e1){
            length = append(length, " - ");
            // cppcheck-suppress [unreadVariable] ; external code
            length = append(length, e1);  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        }
        break;
    case array_size_too_short:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "array size too short");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case input_stream_error:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "input stream error");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case invalid_class_name:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "class name too long");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case unregistered_cast:
        length = append(length, "unregistered void cast ");
        length = append(length, (nullptr != e1) ? e1 : "?");
        length = append(length, "<-");
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, (nullptr != e2) ? e2 : "?");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case unsupported_class_version:
        length = append(length, "class version ");
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, (nullptr != e1) ? e1 : "<unknown class>");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case other_exception:
        // if get here - it indicates a derived exception 
        // was sliced by passing by value in catch
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "unknown derived exception");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    case multiple_code_instantiation:
        length = append(length, "code instantiated in more than one module");
        if(nullptr != e1){
            length = append(length, " - ");
            // cppcheck-suppress [unreadVariable] ; external code
            length = append(length, e1);  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        }
        break;
    case output_stream_error:
        // cppcheck-suppress [unreadVariable] ; external code
        length = append(length, "output stream error");  // NOLINT(clang-analyzer-deadcode.DeadStores) external code
        break;
    default:
        BOOST_ASSERT(false);  // NOLINT(cert-dcl03-c,hicpp-static-assert,misc-static-assert) external code
        length = append(length, "programming error");  // cppcheck-suppress [unreadVariable] ; external code
        break;
    }
}

#if BOOST_VERSION > 105900
BOOST_ARCHIVE_DECL
archive_exception::archive_exception(archive_exception const & oth) BOOST_NOEXCEPT :  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) external code
	std::exception(oth),
	code(oth.code)
{
	std::memcpy(m_buffer,oth.m_buffer,sizeof m_buffer);
}
#endif

BOOST_ARCHIVE_DECL
archive_exception::~archive_exception() BOOST_NOEXCEPT_OR_NOTHROW {}  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) external code

BOOST_ARCHIVE_DECL const char *
archive_exception::what() const BOOST_NOEXCEPT_OR_NOTHROW {
    return m_buffer;
}

BOOST_ARCHIVE_DECL
archive_exception::archive_exception() BOOST_NOEXCEPT :
    code(no_exception)
{}

} // end namespace archive
} // end namespace boost
// NOLINTEND(altera-unroll-loops,cppcoreguidelines-pro-type-member-init,hicpp-member-init,misc-const-correctness) external code
