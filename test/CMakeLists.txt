# Copyright 2018-2026 Alfredo A. Correa
# Copyright 2024 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.16) # 3.16 for precompiled_headers

if (DEFINED BOOST_SUPERPROJECT_VERSION AND NOT DEFINED BOOST_MULTI_STANDALONE AND BUILD_TESTING)

	include(BoostTestJamfile OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

	if(HAVE_BOOST_TEST)

		boost_test_jamfile(FILE Jamfile LINK_LIBRARIES Boost::multi Boost::multi_array Boost::iterator)

	endif()

else()

	enable_testing()

	if(MSVC AND CMAKE_CUDA_COMPILER)
		set(CMAKE_CUDA_RUNTIME_LIBRARY "Shared")
	endif()

	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	# option(CMAKE_COMPILE_WARNING_AS_ERROR "Compile warnings as errors" ON)

	find_package(Boost CONFIG)
# find_package(Boost COMPONENTS headers)  # should be last


#libomp is keg-only, which means it was not symlinked into /opt/homebrew,
#because it can override GCC headers and result in broken builds.
#
#For compilers to find libomp you may need to set:
#  export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
#  export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

	if(NOT Boost_FOUND)
		message(WARNING "preliminary: Cannot find Boost, Multi library will have a very minimal test. If you want to test the library install Boost.Test, for example please run:\n    sudo apt install libboost-test-dev\n    sudo dnf install boost-devel")
	endif()

	if(ENABLE_CUDA)
		enable_language(CUDA)
		set(CMAKE_CUDA_EXTENSIONS OFF)
		if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		  set(CMAKE_CUDA_ARCHITECTURES native)
		endif()
		foreach(_thrust_dir ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} /usr/include)
			if(EXISTS "${_thrust_dir}/thrust/version.h")
				file(STRINGS "${_thrust_dir}/thrust/version.h" _thrust_ver_line REGEX "#define THRUST_VERSION ")
				string(REGEX MATCH "[0-9]+" _thrust_ver "${_thrust_ver_line}")
				math(EXPR _thrust_major "${_thrust_ver} / 100000")
				math(EXPR _thrust_minor "(${_thrust_ver} / 100) % 1000")
				math(EXPR _thrust_patch "${_thrust_ver} % 100")
				message(STATUS "Thrust ${_thrust_major}.${_thrust_minor}.${_thrust_patch} (${_thrust_dir})")
				break()
			endif()
		endforeach()
	endif()

	find_package(TBB)

	enable_testing()
	include(CTest)

	if (CMAKE_CUDA_HOST_COMPILER_ID STREQUAL "MSVC")
		message(STATUS "NVCC's host compiler is MSVC.")
		# Add specific actions or configurations for MSVC host compiler
	else()
		message(STATUS "NVCC's host compiler is not MSVC (ID: ${CMAKE_CUDA_HOST_COMPILER_ID}).")
		# Add specific actions or configurations for other host compilers
	endif()

	set(WARNS
	$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
		#remove Xcompile for Windows?
		#-Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wdouble-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
		--Werror=cross-execution-space-call
		--Werror=ext-lambda-captures-this
		-Xcudafe=--display_error_number
		#--use_fast_math --extra-device-vectorization --restrict
		#-Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
		# $<$<VERSION_LESS:$<CUDA_COMPILER_VERSION>,12>:
		# 	-Xcudafe=--diag_suppress=20011 # for std::string in test/allocator.cpp, more?
		# 	-Xcudafe=--diag_suppress=20014
		# >
		#-Xcudafe=--diag_suppress=1427  # for member_array_cast.cpp in msvc
	>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>,$<NOT:$<CUDA_COMPILER_ID:Clang>>>: # -Werror
		$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-ftrivial-auto-var-init=pattern>
	# -WNSObject-attribute # (gcc 12, not in 11)
	# -Wabi=13 -Wabi-tag # (maybe important when linking with very old libraries)
	-Waddress
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Waddress-of-packed-member>
	-Wno-aggregate-return # (disallow return classes or structs, seems a C-compatibility warning)
	-Waggressive-loop-optimizations
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wanalyzer-out-of-bounds>
	-Waligned-new=all
	# -Walloc-size-larger-than=<bytes> # (gcc 12, not in 11)
	-Walloc-zero
	-Walloca  # -Walloca-larger-than=<number>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Warith-conversion>
	-Warray-bounds # -Warray-bounds=<0,2>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Warray-compare>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Warray-parameter>  #=<0,2> # (gcc 11, not in gcc 10)
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wattribute-alias>  #=<0,2> # (gcc 12, not in 11, not in 7)
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,9>:-Wattribute-warning>
	-Wattributes
	-Wbool-compare -Wbool-operation
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wbidi-chars -Wbidi-chars=any>
	-Wbuiltin-declaration-mismatch -Wbuiltin-macro-redefined
	#-Wc++-compat
	-Wc++0x-compat -Wc++11-compat -Wc++14-compat -Wc++17-compat
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wc++17-extensions>
	-Wc++1z-compat
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wc++20-compat>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wc++20-extensions -Wc++23-extensions>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wc++2a-compat>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,9>:-Wcannot-profile>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8>:-Wcast-align=strict -Wcast-function-type>  # (gcc 8, not in gcc 7)
	-Wcast-qual
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8>:-Wcatch-value>  #=<0, 3>
	-Wchar-subscripts
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wclass-conversion>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wclass-memaccess>
	-Wclobbered
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wcomma-subscript>
	-Wcomment  # (same as -Wcomments)
	-Wconditionally-supported
	-Wconversion -Wconversion-null
	-Wcoverage-mismatch -Wcpp
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wctad-maybe-unsupported>
	-Wctor-dtor-privacy
	-Wdangling-else
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wdangling-pointer>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,13>:-Wdangling-reference>
	-Wdate-time
	-Wdelete-incomplete -Wdelete-non-virtual-dtor
	-Wdeprecated
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wdeprecated-copy -Wdeprecated-copy-dtor>  # (gcc 11, not in gcc 8)
	-Wdeprecated-declarations
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wdeprecated-enum-enum-conversion -Wdeprecated-enum-float-conversion>
	-Wdisabled-optimization
	-Wdiv-by-zero -Wdouble-promotion
	-Wduplicated-branches -Wduplicated-cond
	# -Weffc++ # (doesn't allow some advanced techniques, such as CRTP)
	-Wempty-body -Wendif-labels
	-Wenum-compare
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wenum-conversion>  # (gcc 11, not in 10)
	-Wexpansion-to-defined
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wexceptions>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8>:-Wextra-semi>
	-Wfloat-conversion -Wfloat-equal
	-Wformat=2
	-Wformat-contains-nul
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10>:-Wformat-diag>
	-Wformat-extra-args -Wformat-nonliteral
	-Wformat-overflow=1
	-Wformat-security -Wformat-signedness -Wformat-truncation -Wformat-y2k -Wformat-zero-length
	-Wframe-address  # -Wframe-larger-than=<byte-size>
	-Wfree-nonheap-object -Whsa
	# -Wif-not-aligned # (gcc 8, not in 7)
	-Wignored-attributes
	-Wignored-qualifiers
	-Wimplicit-fallthrough #=3  # -Wimplicit-fallthrough=<0,5>
	# -Winaccessible-base # (gcc 12, not in 11)
	# -Winfinite-recursion # (gcc 12, not in 9)
	-Winherited-variadic-ctor
	#-Winit-list-lifetime # (gcc 12, not in 11)
	-Winit-self
	-Wno-inline  # not all inline attempts can be successful
	-Wint-in-bool-context -Wint-to-pointer-cast
	# -Winterference-size # (gcc 12, not in 11)
	# -Winvalid-imported-macros # (gcc 11, not in gcc 10)
	-Winvalid-memory-model -Winvalid-offsetof -Winvalid-pch
	# -Wlarger-than=<byte-size>  # (disallow large objects types? in executable)
	-Wliteral-suffix
	-Wlogical-not-parentheses -Wlogical-op
	# -Wlong-long # (C++98 warning)
	-Wlto-type-mismatch -Wmain -Wmaybe-uninitialized
	-Wmemset-elt-size -Wmemset-transposed-args
	-Wmisleading-indentation
	# -Wmismatched-dealloc -Wmismatched-new-delete # (gcc 11, not in gcc 10)
	# -Wmismatched-tags # (gcc 11, not in 9)
	# -Wmissing-attributes # (gcc 8, not in 8)
	-Wmissing-braces -Wmissing-declarations -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn
	#-Wmissing-profile # (gcc 11, not in gcc 8)
	# -Wmissing-requires -Wmissing-template-keyword # (gcc 12, not in 11)
	-Wmultichar
	-Wno-multiple-inheritance # (disallows composition by inheritance)
	# -Wmultistatement-macros # (gcc 8, not in 7)
	-Wno-namespaces # (disallows use of namespaces, seems a C-tool)
	-Wnarrowing
	# -Wno-alloc-size-larger-than=<bytes> -Wframe-larger-than=<bytes> -Wno-larger-than<bytes> -Wstack-usage=<bytes> # (gcc 112, not in 11)
	-Wnoexcept -Wnoexcept-type
	-Wnon-template-friend -Wnon-virtual-dtor
	-Wnonnull -Wnonnull-compare
	-Wnormalized  #=nfc -Wnormalized=[none|id|nfc|nfkc]
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14>:-Wnrvo>  # (not in gcc 13)
	-Wnull-dereference
	-Wodr
	-Wold-style-cast  # TODO(correaa) fix this
	# -Wopenacc-parallelism # (gcc 12, not in 11)
	-Wopenmp-simd -Woverflow
	-Woverlength-strings -Woverloaded-virtual
	-Wpacked -Wpacked-bitfield-compat
	# -Wpacked-not-aligned # (gcc 8, not 7)
	# -Wpadded # (disallows structs that need padding for alignment)
	-Wparentheses
	#-Wpessimizing-move # (gcc 11, not in gcc 8)
	-Wplacement-new  #=1  -Wplacement-new=<0,2>
	-Wpmf-conversions
	-Wpointer-arith -Wpointer-compare
	-Wpragmas
	#-Wprio-ctor-dtor # (gcc 11, not in gcc 8)
	-Wpsabi
	# -Wrange-loop-construct # (gcc 11, not in gcc 10)
	-Wredundant-decls
	#-Wredundant-move # (gcc 11, not in gcc 8)
	# -Wredundant-tags # (gcc 11, not in gcc 9)
	-Wregister
	-Wreorder
	-Wreturn-local-addr -Wreturn-type
	-Wrestrict
	-Wscalar-storage-order -Wsequence-point
	-Wshadow -Wshadow-compatible-local -Wshadow-local -Wshadow=compatible-local -Wshadow=local
	-Wshift-count-negative -Wshift-count-overflow -Wshift-negative-value -Wshift-overflow  #=1 -Wshift-overflow=<0,2>
	-Wsign-compare -Wsign-conversion -Wsign-promo
	-Wsized-deallocation
	-Wsizeof-array-argument
	# -Wsizeof-array-div # (gcc 11, not in 10)
	-Wsizeof-pointer-memaccess
	-Wstack-protector  # -Wstack-usage=<byte-size>
	-Wstrict-aliasing #=3 # -Wstrict-aliasing=<0,3>
	-Wstrict-null-sentinel  
	# -Wstrict-overflow=<0,5>
	-Wstrict-overflow=1
	# -Wstrict-overflow=2  # generates false positives when comparing pointers
	# -Wstrict-overflow=3  # incompatible with index_range operator-
	# -Wstrict-overflow=4  # not compatible with stl heap
	# -Wstrict-overflow=5  # not compatible with Boost.Test
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wstring-compare>  # (not in gcc 9)
	-Wstringop-overflow  #=2  -Wstringop-overflow=<0,4>
	# -Wstringop-overread # (gcc 11, not in 10)
	# -Wstringop-truncation # (gcc 8, not in 7)
	-Wsubobject-linkage
	# -Wsuggest-attribute=cold # (gcc 12, not in 11, not in 7)
	-Wsuggest-attribute=const -Wsuggest-attribute=format
	# -Wsuggest-attribute=malloc # (gcc 8, not in 7)
	-Wsuggest-attribute=noreturn # -Wsuggest-attribute=pure # (false positive in gcc 8.5.0 RedHat)
	-Wsuggest-final-methods -Wsuggest-final-types
	-Wsuggest-override
	-Wswitch -Wswitch-bool
	-Wno-switch-default
	-Wswitch-enum
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10>:-Wswitch-outside-range>
	-Wswitch-unreachable
	-Wsync-nand -Wsynth
	-Wno-system-headers  # (expects system headers to be warning-compliant which they are not)
	-Wtautological-compare
	-Wno-templates # (disallows templates, C-tool)
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wterminate>  # (gcc 12, not in 11)
	-Wtrampolines -Wtrigraphs
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wtrivial-auto-var-init>  # (not in gcc 11)
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Wtsan>  # (gcc 11, not in 10)
	-Wtype-limits -Wundef -Wuninitialized
	-Wno-unknown-pragmas  # (see above) -Wunknown-pragmas (other compilers need their own pragmas for their warnings)
	-Wunreachable-code -Wunsafe-loop-optimizations
	-Wunused -Wunused-but-set-parameter -Wunused-but-set-variable
	-Wunused-const-variable  #=2 TODO(correaa) add [[maybe_unused]] to niebloids
	-Wunused-function -Wunused-label -Wunused-local-typedefs -Wunused-macros -Wunused-parameter -Wunused-result -Wunused-value -Wunused-variable
	# -Wuse-after-free  # =<0,3> # (gcc 12, not in 11)
	-Wuseless-cast
	-Wvarargs -Wvariadic-macros -Wvector-operation-performance
	# -Wvexing-parse # (gcc 11, not in gcc 10)
	-Wvirtual-inheritance -Wvirtual-move-assign
	-Wvla # -Wvla-larger-than=<number> # (gcc 12, not in 11)
	# -Wvla-parameter # (gcc 11, not in gcc 10)
	# -Wvolatile # (gcc 11, not in gcc 9)
	-Wvolatile-register-var
	-Wwrite-strings
	-Wzero-as-null-pointer-constant
	# -Wzero-length-bounds # (gcc 12, not in 11)
	>
	$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CUDA_COMPILER_ID:Clang>>:
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-ftrivial-auto-var-init=pattern>
	-Wall -Wextra
	-Wno-error=\#warnings  # -Wno-\#warnings disable this working for a bug in rocm ~5.6
	-Wno-c++98-compat -Wno-c++98-compat-bind-to-temporary-copy -Wno-c++98-compat-pedantic
	-Wno-poison-system-directories
	#$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wno-pre-c++17-compat>
	# -Wno-exit-time-destructors  # Boost.LightweightTest uses static objects with destructor
	-Wno-large-by-value-copy  # subarray references can be large
	#-Wno-padded  # most classes, e.g. layouts, need padding
	-Wstrict-aliasing -Wstrict-aliasing=0 -Wstrict-aliasing=1 -Wstrict-aliasing=2
	-Wstrict-overflow -Wstrict-overflow=0 -Wstrict-overflow=1 -Wstrict-overflow=2 -Wstrict-overflow=3 -Wstrict-overflow=4 -Wstrict-overflow=5
	-Wno-undefined-reinterpret-cast  # needed to cast to c-array
	-Wunknown-pragmas  # other compilers need their own pragmas
	# -Rpass -Rpass-analysis -Rpass-missed
	# -Rsanitize-address # -Rsearch-path-usage
	>
	$<$<CXX_COMPILER_ID:Intel>: # -Werror
	-Wall -Wextra # also IntelLLVM, XL (ibm), XLClang (ibm)
	-diag-disable=remark
	-diag-error:3846
	-diag-disable=1011  # disables warning missing return at the end of non-void function
	-diag-disable=2196  # disables "error #2196: routine is both "inline" and "noinline"" in icpc 2021.5.0
	-wd161
	-Wabi
	-Warray-bounds
	-Wcast-qual
	-Wchar-subscripts
	-Wcomment
	-Wdeprecated
	-Wenum-compare
	-Wextra-tokens
	-Wformat -Wformat=2 -Wformat-security
	-Wic-pointer
	-Wnarrowing
	-Wreturn-type
	-Wnon-virtual-dtor
	-Wnonnull
	-Wmaybe-uninitialized
	-Wmain
	-Wmissing-declarations -Wmissing-prototypes
	-Wmultichar
	-Woverloaded-virtual
	-Woverflow
	-Wp64
	-Wparentheses
	-Wpointer-arith
	-Wpointer-sign
	-Wreorder
	-Wreturn-type
	-Wsequence-point
	-Wshadow
	-Wsign-compare
	-Wshorten-64-to-32
	-Wmissing-prototypes
	-Wstrict-aliasing
	-Wno-unknown-pragmas
	-Wstrict-prototypes
	-Wtrigraphs
	-Wtype-limits
	-Wuninitialized
	-Wunused -Wunused-but-set-variable -Wunused-function -Wunused-parameter -Wunused-variable
	-Wwrite-strings
	>
	$<$<CXX_COMPILER_ID:NVHPC>:
		$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,23>:
			#-Wfatal-errors
			-fmax-errors=2
		>
		-Wall -Wextra  # -Wpedantic
		#-Wcast-qual -Wformat=2 -Wshadow -Wsign-compare -Wundef -Wuninitialized -Wunused
	>
	$<$<CXX_COMPILER_ID:MSVC>:
		$<$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:
			/W4  # /volatile:iso /WX /Zc:wchar_t /Zc:forScope /Zc:inline
			/wd4514  # unreferenced inline function has been removed
			/wd4866  # compiler may not enforce left-to-right evaluation order for call to 'boost::multi::subarray<int,1,int *,boost::multi::layout_t<1,__int64> >::operator[]' [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\test\array_ref.cpp.x.vcxproj]
			/wd4868  # compiler may not enforce left-to-right evaluation order in braced initializer list
			/wd5045  # Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\include\boost\multi\adaptors\thrust\omp\test\omp.cpp.x.vcxproj]
			/wd5031  # likely mismatch, popping warning state pushed in different file
		>
		$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
			-Xcompiler=/W4  # /volatile:iso /WX /Zc:wchar_t /Zc:forScope /Zc:inline /W4
			-Xcompiler=/wd4464  # relative include path contains '..'
			-Xcompiler=/wd4514  # unreferenced inline function has been removed
			-Xcompiler=/wd4668  # warning C4668: '__NV_NO_HOST_COMPILER_CHECK' is not defined as a preprocessor macro, replacing with '0' for '#if/#elif'
			-Xcompiler=/wd4866  # compiler may not enforce left-to-right evaluation order for call
			-Xcompiler=/wd4868  # compiler may not enforce left-to-right evaluation order in braced initializer list
			-Xcompiler=/wd5045  # Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\include\boost\multi\adaptors\thrust\omp\test\omp.cpp.x.vcxproj]
		>
	>
	)

	file(
		GLOB TEST_SRCS
		RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
		*.cpp
	)

	if(NOT Boost_FOUND)
		message(WARNING "post: Cannot find Boost, Multi library will have a very minimal test. If you want to test the library install Boost.Core (for Boost.LightweightTest), for example please run:\n    sudo apt install libboost-dev\n    sudo dnf install boost-devel")

		add_executable(main main.cpp)

		target_link_libraries(main PUBLIC multi)
		target_compile_features   (main PUBLIC cxx_std_17)
		target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/include)

		if(NOT DEFINED ENABLE_CIRCLE)
			target_compile_options(main PRIVATE ${WARNS})
		else()
			target_compile_options    (main PRIVATE -Wall)  # -Werror
		endif()
		add_test(NAME main COMMAND $<TARGET_FILE:main>)

	else()

		foreach(TEST_FILE ${TEST_SRCS})
			set(TEST_EXE "${TEST_FILE}.x")
			add_executable(${TEST_EXE} ${TEST_FILE})

			if(ENABLE_CUDA)
				set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
			endif()

			target_link_libraries(${TEST_EXE} PRIVATE multi)
			target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

			if(TBB_FOUND)
				target_link_libraries    (${TEST_EXE} PRIVATE TBB::tbb)
				target_compile_definitions(${TEST_EXE} PUBLIC TBB_FOUND=1)
			endif()

			# set_property(TARGET ${TEST_EXE} PROPERTY COMPILE_WARNING_AS_ERROR ON)

			if(NOT DEFINED ENABLE_CIRCLE)
				target_compile_options(${TEST_EXE} PRIVATE ${WARNS})
			else()
				target_compile_options(${TEST_EXE} PRIVATE -Wall)  # -Werror
			endif()

			set_target_properties(${TEST_EXE} PROPERTIES CXX_EXTENSIONS OFF)

			add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
		endforeach()
	endif()
endif()
