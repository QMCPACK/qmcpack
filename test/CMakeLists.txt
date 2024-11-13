# Copyright 2018-2024 Alfredo A. Correa
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

	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_CXX_EXTENSIONS OFF)

	find_package(Boost)  #  COMPONENTS boost unit_test_framework)  # tests require Boost.Test `sudo apt install libboost-test-dev`

	if(ENABLE_CUDA)
		enable_language(CUDA)
		set(CMAKE_CUDA_EXTENSIONS OFF)
		if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		  set(CMAKE_CUDA_ARCHITECTURES native)
		endif()
	endif()

	find_package(TBB)

	enable_testing()
	include(CTest)

	set(WARNS
	$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
		-Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wdouble-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
		#--Werror=cross-execution-space-call
		--Werror ext-lambda-captures-this
		--use_fast_math
		--extra-device-vectorization
		--restrict
		-Xcudafe=--display_error_number -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
	>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>,$<NOT:$<CUDA_COMPILER_ID:Clang>>>: -Werror
		$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-ftrivial-auto-var-init=pattern>
	# -WNSObject-attribute # (gcc 12, not in 11)
	# -Wabi=13 -Wabi-tag # (maybe important when linking with very old libraries)
	-Waddress
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Waddress-of-packed-member>
	-Wno-aggregate-return # (disallow return classes or structs, seems a C-compatibility warning)
	-Waggressive-loop-optimizations
	-Waligned-new=all # (gcc 12, not in 11)
	# -Walloc-size-larger-than=<bytes> # (gcc 12, not in 11)
	-Walloc-zero
	-Walloca  # -Walloca-larger-than=<number>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>:-Warith-conversion>
	-Warray-bounds # -Warray-bounds=<0,2>
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Warray-compare>
	# -Warray-parameter #=<0,2> # (gcc 11, not in gcc 10)
	# -Wattribute-alias #=<0,2> # (gcc 12, not in 11, not in 7)
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
	#-Wcannot-profile # (gcc 9, not in gcc 8)
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8>:-Wcast-align=strict -Wcast-function-type> # (gcc 8, not in gcc 7)
	-Wcast-qual
	#-Wcatch-value  #=<0, 3> # (gcc 8, not in gcc 7)
	-Wchar-subscripts
	#-Wclass-conversion # (gcc 11, not in 8)
	#-Wclass-memaccess # (gcc 12, not in 11, not in 7)
	-Wclobbered
	# -Wcomma-subscript # (gcc 12, not in 11)
	-Wcomment  # (same as -Wcomments)
	-Wconditionally-supported
	-Wconversion -Wconversion-null
	-Wcoverage-mismatch -Wcpp
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wctad-maybe-unsupported>
	-Wctor-dtor-privacy
	-Wdangling-else
	# -Wdangling-pointer # (gcc 12, not in 11)
	# -Wdangling-reference # (gcc 13, not in 12)
	-Wdate-time
	-Wdelete-incomplete -Wdelete-non-virtual-dtor
	-Wdeprecated
	#-Wdeprecated-copy -Wdeprecated-copy-dtor # (gcc 11, not in gcc 8)
	-Wdeprecated-declarations
	# -Wdeprecated-enum-enum-conversion -Wdeprecated-enum-float-conversion # (gcc 11, not in gcc 10)
	-Wdisabled-optimization
	-Wdiv-by-zero -Wdouble-promotion
	-Wduplicated-branches -Wduplicated-cond
	# -Weffc++ # (doesn't allow some advanced techniques, such as CRTP)
	-Wempty-body -Wendif-labels
	-Wenum-compare
	# -Wenum-conversion # (gcc 11, not in 10)
	-Wexpansion-to-defined
	# -Wexceptions  # (gcc 11, not in 10)
	# -Wextra-semi # (gcc 8, not in 7)
	-Wfloat-conversion -Wfloat-equal
	-Wformat=2
	-Wformat-contains-nul # (gcc 12, not in 11)
	# -Wformat-diag # (gcc 10, not in 9)
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
	# -Wnrvo # (not in gcc 12)
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
	-Wreorder # (gcc 12, not in 11)
	-Wreturn-local-addr -Wreturn-type
	-Wrestrict -Wreorder
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
	-Wstrict-null-sentinel  #=1  -Wstrict-overflow=<0,5>
	-Wstrict-overflow=2
	# -Wstrict-overflow=3  # incompatible with index_range operator-
	# -Wstrict-overflow=4  # not compatible with stl heap
	# -Wstrict-overflow=5  # not compatible with Boost.Test
	# -Wstring-compare # (gcc 11, not in 9)
	-Wstringop-overflow  #=2  -Wstringop-overflow=<0,4>
	# -Wstringop-overread # (gcc 11, not in 10)
	# -Wstringop-truncation # (gcc 8, not in 7)
	-Wsubobject-linkage
	# -Wsuggest-attribute=cold # (gcc 12, not in 11, not in 7)
	-Wsuggest-attribute=const -Wsuggest-attribute=format
	# -Wsuggest-attribute=malloc # (gcc 8, not in 7)
	-Wsuggest-attribute=noreturn # -Wsuggest-attribute=pure # (false positive in gcc 8.5.0 RedHat)
	-Wsuggest-final-methods -Wsuggest-final-types
	-Wsuggest-override # (gcc 12, not in gcc 11)
	-Wswitch -Wswitch-bool
	-Wno-switch-default
	-Wswitch-enum
	# -Wswitch-outside-range # (gcc 11, not in gcc 9)
	-Wswitch-unreachable
	-Wsync-nand -Wsynth
	-Wno-system-headers  # (expects system headers to be warning-compliant which they are not)
	-Wtautological-compare
	-Wno-templates # (disallows templates, C-tool)
	-Wterminate # (gcc 12, not in 11)
	-Wtrampolines -Wtrigraphs
	# -Wtrivial-auto-var-init # (gcc 12, not in 11)
	# -Wtsan # (gcc 11, not in 10)
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
	$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CUDA_COMPILER_ID:Clang>>:-Werror -Wall -Wpedantic -Wextra
	-Wno-error=\#warnings  # -Wno-\#warnings disable this working for a bug in rocm ~5.6
	-WCFString-literal
	-WCL4
	-WIndependentClass-attribute
	-WNSObject-attribute
	-Wabi
	-Wabsolute-value
	-Wabstract-final-class -Wabstract-vbase-init
	-Waddress -Waddress-of-packed-member -Waddress-of-temporary
	-Waggregate-return
	# -Waix-compat (clang 13, not in clang 11)
	# -Walign-mismatch (clang 13, not in clang 11)
	# -Walloca (clang 11, not in clang 9)
	-Walloca-with-align-alignof
	# -Walways-inline-coroutine (clang 15, not in clang 13)
	-Wambiguous-delete -Wambiguous-ellipsis -Wambiguous-macro -Wambiguous-member-template
	# -Wambiguous-reversed-operator (clang 11, not in clang 9)
	-Wanalyzer-incompatible-plugin
	# -Wanon-enum-enum-conversion (clang 11, not in clang 9)
	-Wanonymous-pack-parens
	-Warc -Warc-bridge-casts-disallowed-in-nonarc -Warc-maybe-repeated-use-of-weak -Warc-non-pod-memaccess -Warc-performSelector-leaks -Warc-repeated-use-of-weak -Warc-retain-cycles -Warc-unsafe-retained-assign
	# -Wargument-outside-range (clang 15, not in clang 13)
	# -Wargument-undefined-behaviour (clang 15, not in clang 13)
	-Warray-bounds -Warray-bounds-pointer-arithmetic
	-Wasm -Wasm-operand-widths
	-Wassign-enum
	-Wassume
	-Wat-protocol
	-Watimport-in-framework-header
	# -Watomic-access (clang 15, not in clang 13)
	-Watomic-alignment
	# -Watomic-implicit-seq-cst (not in clang 7)
	-Watomic-memory-ordering -Watomic-properties -Watomic-property-with-user-defined-accessor
	-Wattribute-packed-for-bitfield
	# -Wattribute-warning  (clang 15, not in clang 13)
	-Wattributes
	-Wauto-disable-vptr-sanitizer -Wauto-import -Wauto-storage-class -Wauto-var-id
	-Wavailability
	# -Wavr-rtlib-linking-quirks (not in clang 7)
	-Wbackend-plugin
	-Wbackslash-newline-escape
	-Wbad-function-cast
	-Wbinary-literal
	-Wbind-to-temporary-copy
	-Wbinding-in-condition
	# -Wbit-int-extension (clang 15, not in clang 13)
	-Wbitfield-constant-conversion -Wbitfield-enum-conversion -Wbitfield-width
	# -Wbitwise-conditional-parentheses (clang 11, not in clang 9)
	# -Wbitwise-instead-of-logical (clang 15, not in clang 13)
	-Wbitwise-op-parentheses
	-Wblock-capture-autoreleasing
	-Wbool-conversion -Wbool-conversions
	# -Wbool-operation (clang 11, not in clang 9)
	-Wbraced-scalar-init
	# -Wbranch-protection (clang 15, not in clang 13)
	-Wbridge-cast
	# -Wbuiltin-assume-aligned-alignment (clang 11, not in clang 9)
	-Wbuiltin-macro-redefined -Wbuiltin-memcpy-chk-size -Wbuiltin-requires-header
	# -Wc++-compat
	# -Wc++0x-compat -Wc++0x-extensions -Wc++0x-narrowing
	# -Wc++11-compat -Wc++11-compat-deprecated-writable-strings -Wc++11-compat-pedantic -Wc++11-compat-reserved-user-defined-literal -Wc++11-extensions -Wc++11-extra-semi -Wc++11-inline-namespace -Wc++11-long-long -Wc++11-narrowing
	# -Wc++14-attribute-extensions (clang 15, not in clang 13)
	-Wc++14-binary-literal
	# -Wc++14-compat -Wc++14-compat-pedantic -Wc++14-extensions
	# -Wc++17-attribute-extensions (clang 15, not in clang 13)
	# -Wc++17-compat -Wc++17-compat-mangling -Wc++17-compat-pedantic -Wc++17-extensions
	# -Wc++1y-extensions
	# -Wc++1z-compat -Wc++1z-compat-mangling -Wc++1z-extensions
	# -Wc++20-attribute-extensions (clang 15, not in clang 13)
	# -Wc++20-compat -Wc++20-compat-pedantic
	# -Wc++20-designator -Wc++20-extensions
	# -Wc++2a-compat -Wc++2a-compat-pedantic
	# -Wc++2a-extensions
	# -Wc++2b-extensions
	# -Wc++98-c++11-c++14-c++17-compat -Wc++98-c++11-c++14-c++17-compat-pedantic -Wc++98-c++11-c++14-compat -Wc++98-c++11-c++14-compat-pedantic
	# -Wc++98-c++11-compat -Wc++98-c++11-compat-binary-literal -Wc++98-c++11-compat-pedantic
	# -Wc++98-compat
	-Wno-c++98-compat-bind-to-temporary-copy
	#-Wc++98-compat-extra-semi -Wc++98-compat-local-type-template-args -Wc++98-compat-pedantic -Wc++98-compat-unnamed-type-template-args
	-Wno-c++98-compat-unnamed-type-template-args
	# -Wc11-extensions
	# -Wc2x-extensions
	# -Wc99-compat
	# -Wc99-designator -Wc99-extensions
	-Wcast-align -Wcast-qual
	# -Wcall-to-pure-virtual-from-ctor-dtor (not in clang 7)
	# -Wcalled-once-parameter (clang 13, not in clang 11)
	-Wcast-align -Wcast-calling-convention
	# -Wcast-function-type (clang 13, not in clang 11)
	-Wcast-of-sel-type -Wcast-qual -Wcast-qual-unrelated
	-Wchar-align -Wchar-subscripts
	-Wclang-cl-pch
	# -Wclass-conversion (clang 11, not in clang 9)
	-Wclass-varargs
	# -Wcmse-union-leak (clang 11, not in clang 9)
	-Wcomma
	-Wcomment -Wcomments
	-Wcompare-distinct-pointer-types
	# -Wcompletion-handler (clang 13, not in clang 11)
	-Wcomplex-component-init
	# -Wcompound-token-split -Wcompound-token-split-by-macro -Wcompound-token-split-by-space (clang 13, not in clang 11)
	-Wconversion
	#-Wconcepts-ts-compat (not working in clang++ 16)
	-Wconditional-type-mismatch -Wconditional-uninitialized
	-Wconfig-macros
	-Wconstant-conversion
	# -Wconstant-evaluated (clang 11, not in clang 9)
	-Wconstant-logical-operand
	-Wconstexpr-not-const
	-Wconsumed
	-Wconversion -Wconversion-null
	-Wcoroutine -Wcoroutine-missing-unhandled-exception
	-Wcovered-switch-default
	-Wcpp
	-Wcstring-format-directive
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,18>:-Wctad-maybe-unsupported>
	-Wctor-dtor-privacy
	# -Wctu # (not in clang 7)
	-Wcuda-compat
	-Wcustom-atomic-properties
	# -Wcxx-attribute-extension (clang 13, not in clang 11)
	-Wdangling -Wdangling-else -Wdangling-field
	# -Wdangling-gsl (clang 13, not in clang 9)
	-Wdangling-initializer-list
	# -Wdarwin-sdk-settings (not in clang 7)
	-Wdate-time
	-Wdealloc-in-category
	-Wdebug-compression-unavailable
	-Wdeclaration-after-statement
	# -Wdefaulted-function-deleted (not in clang 7)
	-Wdelegating-ctor-cycles
	# -Wdelete-abstract-non-virtual-dtor (not in clang 7)
	-Wdelete-incomplete
	# -Wdelete-non-abstract-non-virtual-dtor (not in clang 7)
	# -Wdelete-non-virtual-dtor (not in clang 7)
	# -Wdelimited-escape-sequence-extension (clang 15, not in clang 13)
	-Wdeprecated
	# -Wdeprecated-altivec-src-compat (clang 13, not in clang 11)
	# -Wdeprecated-anon-enum-enum-conversion -Wdeprecated-array-compare (clang 11, not in clang 9)
	-Wdeprecated-attributes
	# -Wdeprecated-comma-subscript (clang 11, not in clang 9)
	# -Wdeprecated-copy -Wdeprecated-copy-dtor (clang 11, not in clang 9)
	# -Wdeprecated-copy-with-dtor (clang 13, not in clang 11)
	# -Wdeprecated-copy-with-user-provided-copy (clang 13, not in clang 11)
	# -Wdeprecated-copy-with-user-provided-dtor (clang 13, not in clang 11)
	# -Wdeprecated-coroutine (clang 15, not in clang 13)
	-Wdeprecated-declarations -Wdeprecated-dynamic-exception-spec
	# -Wdeprecated-enum-compare -Wdeprecated-enum-compare-conditional -Wdeprecated-enum-enum-conversion (clang 11, not in clang 9)
	# -Wdeprecated-enum-float-conversion (clang 11, not in clang 9)
	# -Wdeprecated-experimental-coroutine (clang 15, not in clang 13)
	-Wdeprecated-implementations -Wdeprecated-increment-bool
	# -Wdeprecated-non-prototype (clang 15, not in clang 13)
	-Wdeprecated-objc-isa-usage -Wdeprecated-objc-pointer-introspection -Wdeprecated-objc-pointer-introspection-performSelector
	# -Wdeprecated-pragma (clang 15, not in clang 13)
	-Wdeprecated-register -Wdeprecated-this-capture
	# -Wdeprecated-type (clang 15, not in clang 13)
	# -Wdeprecated-volatile (clang 11, not in clang 9)
	-Wdeprecated-writable-strings
	-Wdirect-ivar-access
	-Wdisabled-macro-expansion
	-Wdisabled-optimization
	-Wdiscard-qual
	-Wdistributed-object-modifiers
	-Wdiv-by-zero
	-Wdivision-by-zero
	-Wdll-attribute-on-redeclaration -Wdllexport-explicit-instantiation-decl -Wdllimport-static-field-def
	-Wdocumentation -Wdocumentation-deprecated-sync -Wdocumentation-html -Wdocumentation-pedantic -Wdocumentation-unknown-command
	-Wdollar-in-identifier-extension
	-Wdouble-promotion
	# -Wdtor-name -Wdtor-typedef (clang 11, not in clang 9)
	-Wduplicate-decl-specifier -Wduplicate-enum -Wduplicate-method-arg -Wduplicate-method-match -Wduplicate-protocol
	-Wdynamic-class-memaccess -Wdynamic-exception-spec
	# -Weffc++
	# -Welaborated-enum-base -Welaborated-enum-class (clang 11, not in clang 9)
	-Wembedded-directive
	-Wempty-body -Wempty-decomposition
	# -Wempty-init-stmt (not in clang 7)
	-Wempty-translation-unit
	-Wencode-type
	-Wendif-labels
	-Wenum-compare
	# -Wenum-compare-conditional (clang 11, not in clang 9)
	-Wenum-compare-switch -Wenum-conversion
	# -Wenum-enum-conversion -Wenum-float-conversion (clang 11, not in clang 9)
	-Wenum-too-large
	-Wexceptions
	# -Wexcess-initializers (clang 11, not in clang 9)
	-Wexit-time-destructors
	-Wexpansion-to-defined
	-Wexplicit-initialize-call -Wexplicit-ownership-type
	# -Wexport-unnamed (not in clang 7)
	# -Wexport-using-directive (not in clang 7)
	-Wextern-c-compat -Wextern-initializer
	#-Wextra
	-Wextra-qualification
	# -Wextra-semi (clang 8, not in 7)
	# -Wextra-semi-stmt (not in clang 7)
	-Wextra-tokens
	# -Wfinal-dtor-non-final-class (clang 11, not in clang 9)
	# -Wfinal-macro (clang 15, not in clang 13)
	# -Wfixed-enum-extension (not in clang 7)
	# -Wfixed-point-overflow (clang 11, not in clang 9)
	-Wflag-enum
	-Wflexible-array-extensions
	-Wfloat-conversion -Wfloat-equal -Wfloat-overflow-conversion -Wfloat-zero-conversion
	-Wfor-loop-analysis
	-Wformat -Wformat-extra-args
	# -Wformat-insufficient-args (clang 13, not in clang 11)
	-Wformat-invalid-specifier -Wformat-non-iso -Wformat-nonliteral -Wformat-pedantic -Wformat-security
	# -Wformat-type-confusion (clang 11, not in clang 9)
	-Wformat-y2k -Wformat-zero-length -Wformat=2
	# -Wfortify-source (not in clang 7)
	-Wfour-char-constants
	# -Wframe-address (clang 11, not in clang 9)
	# -Wframe-larger-than # -Wframe-larger-than=
	-Wframework-include-private-from-public
	# -Wfree-nonheap-object (clang 13, not in clang 11)
	-Wfunction-def-in-objc-container -Wfunction-multiversion
	# -Wfuse-ld-path (clang 13, not in clang 11)
	# -Wfuture-attribute-extensions (clang 15, not in clang 13)
	# -Wfuture-compat
	# -Wgcc-compat
	-Wglobal-constructors
	# -Wglobal-isel (clang 11, not in clang 9)
	-Wgnu -Wgnu-alignof-expression -Wgnu-anonymous-struct -Wgnu-array-member-paren-init -Wgnu-auto-type -Wgnu-binary-literal -Wgnu-case-range -Wgnu-complex-integer -Wgnu-compound-literal-initializer -Wgnu-conditional-omitted-operand -Wgnu-designator -Wgnu-empty-initializer -Wgnu-empty-struct -Wgnu-flexible-array-initializer -Wgnu-flexible-array-union-member -Wgnu-folding-constant -Wgnu-imaginary-constant -Wgnu-include-next
	# -Wgnu-inline-cpp-without-extern (clang 11, not in clang 9)
	-Wgnu-label-as-value
	# -Wgnu-null-pointer-arithmetic
	# -Wgnu-pointer-arith (clang 15, not in clang 13)
	-Wgnu-redeclared-enum -Wgnu-statement-expression -Wgnu-static-float-init -Wgnu-string-literal-operator-template -Wgnu-union-cast -Wgnu-variable-sized-type-not-at-end -Wgnu-zero-line-directive -Wgnu-zero-variadic-macro-arguments
	# -Wgpu-maybe-wrong-side (clang 15, not in clang 13)
	-Wheader-guard -Wheader-hygiene
	# -Whip-only (clang 11, not in clang 9)
	-Widiomatic-parentheses
	-Wignored-attributes
	# -Wignored-availability-without-sdk-settings (clang 13, not in clang 11)
	-Wignored-optimization-argument
	#-Wignored-pragma-intrinsic
	#-Wignored-pragma-optimize
	#-Wignored-pragmas (other compilers need pragmas)
	-Wignored-qualifiers  # turns out that returning const types is not really ignored by the compiler and prevents certain misleading syntaxes
	# -Wignored-reference-qualifiers (clang 15, not in clang 13)
	-Wimplicit -Wimplicit-atomic-properties
	# -Wimplicit-const-int-float-conversion (clang 11, not in clang 9)
	-Wimplicit-conversion-floating-point-to-bool -Wimplicit-exception-spec-mismatch -Wimplicit-fallthrough -Wimplicit-fallthrough-per-function
	# -Wimplicit-fixed-point-conversion -Wimplicit-float-conversion (not in clang 7)
	-Wimplicit-function-declaration -Wimplicit-int
	# -Wimplicit-int-conversion (not in clang 7)
	# -Wimplicit-int-float-conversion (clang 11, not in clang 9)
	-Wimplicit-retain-self -Wimplicitly-unsigned-literal
	-Wimport -Wimport-preprocessor-directive-pedantic
	-Winaccessible-base
	-Winclude-next-absolute-path -Winclude-next-outside-header
	-Wincompatible-exception-spec -Wincompatible-function-pointer-types -Wincompatible-library-redeclaration -Wincompatible-ms-struct -Wincompatible-pointer-types -Wincompatible-pointer-types-discards-qualifiers -Wincompatible-property-type -Wincompatible-sysroot -Wincomplete-framework-module-declaration -Wincomplete-implementation -Wincomplete-module
	# -Wincomplete-setjmp-declaration (not in clang 7)
	-Wincomplete-umbrella
	-Winconsistent-dllimport -Winconsistent-missing-destructor-override -Winconsistent-missing-override
	-Wincrement-bool
	-Winfinite-recursion
	-Winit-self
	-Winitializer-overrides
	-Winjected-class-name
	-Winline -Winline-asm
	# -Winline-namespace-reopened-noninline (clang 11, not in clang 9)
	-Winline-new-delete
	-Winstantiation-after-specialization
	-Wint-conversion -Wint-conversions
	# -Wint-in-bool-context (clang 11, not in clang 9)
	-Wint-to-pointer-cast -Wint-to-void-pointer-cast
	-Winteger-overflow
	# -Winterrupt-service-routine (clang 13, not in clang 11)
	-Winvalid-command-line-argument -Winvalid-constexpr -Winvalid-iboutlet -Winvalid-initializer-from-system-header -Winvalid-ios-deployment-target
	# -Winvalid-no-builtin-names (clang 11, not in clang 9)
	-Winvalid-noreturn -Winvalid-offsetof -Winvalid-or-nonexistent-directory -Winvalid-partial-specialization -Winvalid-pch -Winvalid-pp-token -Winvalid-source-encoding -Winvalid-token-paste
	-Wjump-seh-finally
	-Wkeyword-compat -Wkeyword-macro
	-Wknr-promoted-parameter
	-Wlanguage-extension-token
	-Wno-large-by-value-copy  # subarray references can be large
	-Wliblto
	# -Wlinker-warnings (clang 15, not in clang 13)
	-Wliteral-conversion -Wliteral-range
	# -Wlocal-type-template-args
	-Wlogical-not-parentheses -Wlogical-op-parentheses
	#-Wlong-long (C++98 warning)
	-Wloop-analysis
	-Wmacro-redefined
	-Wmain -Wmain-return-type
	-Wmalformed-warning-check
	-Wmany-braces-around-scalar-init
	# -Wmax-tokens (clang 11, not in clang 9)
	-Wmax-unsigned-zero
	-Wmemset-transposed-args -Wmemsize-comparison
	-Wmethod-signatures
	-Wmicrosoft
	# -Wmicrosoft-abstract (clang 13, not in clang 11)
	-Wmicrosoft-anon-tag -Wmicrosoft-cast -Wmicrosoft-charize -Wmicrosoft-comment-paste -Wmicrosoft-const-init -Wmicrosoft-cpp-macro -Wmicrosoft-default-arg-redefinition
	# -Wmicrosoft-drectve-section (not in clang 7)
	-Wmicrosoft-end-of-file -Wmicrosoft-enum-forward-reference -Wmicrosoft-enum-value -Wmicrosoft-exception-spec -Wmicrosoft-exists -Wmicrosoft-explicit-constructor-call -Wmicrosoft-extra-qualification -Wmicrosoft-fixed-enum -Wmicrosoft-flexible-array -Wmicrosoft-goto -Wmicrosoft-inaccessible-base -Wmicrosoft-include  -Wmicrosoft-mutable-reference -Wmicrosoft-pure-definition -Wmicrosoft-redeclare-static -Wmicrosoft-sealed
	# -Wmicrosoft-static-assert (clang 13, not in clang 11)
	-Wmicrosoft-template
	# -Wmicrosoft-template-shadow (clang 11, not in clang 9)
	-Wmicrosoft-union-member-reference -Wmicrosoft-unqualified-friend -Wmicrosoft-using-decl -Wmicrosoft-void-pseudo-dtor
	# -Wmisexpect (clang 15, not in clang 13)
	# -Wmisleading-indentation (clang 11, not in clang 9)
	-Wmismatched-new-delete -Wmismatched-parameter-types -Wmismatched-return-types -Wmismatched-tags
	-Wmissing-braces
	# -Wmissing-constinit (clang 11, not in clang 9)
	-Wmissing-declarations -Wmissing-exception-spec -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-method-return-type -Wmissing-noescape -Wmissing-noreturn -Wmissing-prototype-for-cc -Wmissing-prototypes -Wmissing-selector-name -Wmissing-sysroot -Wmissing-variable-declarations
	# -Wmisspelled-assumption (clang 13, not in clang 11)
	# -Rmodule-build
	-Wmodule-conflict -Wmodule-file-config-mismatch -Wmodule-file-extension
	# -Rmodule-import (not in clang 7)
	-Wmodule-import-in-extern-c
	# -Rmodule-lock
	-Wmodules-ambiguous-internal-linkage -Wmodules-import-nested-redundant
	-Wmost
	-Wmove
	-Wmsvc-include -Wmsvc-not-found
	-Wmultichar
	-Wmultiple-move-vbase
	-Wnon-virtual-dtor
	-Wnarrowing
	-Wnested-anon-types -Wnested-externs
	-Wnew-returns-null
	-Wnewline-eof
	# -Wnoderef (not in clang 7)
	-Wnoexcept-type
	# -Wnon-c-typedef-for-linkage (clang 11, not in clang 9)
	-Wnon-gcc
	-Wnon-literal-null-conversion
	-Wnon-modular-include-in-framework-module -Wnon-modular-include-in-module
	-Wnon-pod-varargs
	# -Wnon-power-of-two-alignment (clang 11, not in clang 9)
	-Wnon-virtual-dtor
	-Wnonnull
	#-Wnonportable-cfstrings
	-Wnonportable-include-path
	-Wnonportable-system-include-path
	-Wnonportable-vector-initialization
	-Wnontrivial-memaccess
	-Wnsconsumed-mismatch
	-Wnsreturns-mismatch
	-Wnull-arithmetic -Wnull-character -Wnull-conversion -Wnull-dereference -Wnull-pointer-arithmetic
	# -Wnull-pointer-subtraction (clang 13, not in clang 11)
	-Wnullability -Wnullability-completeness -Wnullability-completeness-on-arrays -Wnullability-declspec -Wnullability-extension -Wnullability-inferred-on-nested-type
	-Wnullable-to-nonnull-conversion
	-Wodr
	-Wold-style-cast
	-Wold-style-definition
	-Wopencl-unsupported-rgba
	# -Wopenmp (clang 11, not in clang 9)
	# -Wopenmp-51-extensions (clang 13, not in clang 11)
	-Wopenmp-clauses -Wopenmp-loop-form
	# -Wopenmp-mapping (clang 11, not in clang 9)
	-Wopenmp-target
	-Woption-ignored
	-Wordered-compare-function-pointers
	-Wout-of-line-declaration -Wout-of-scope-function
	-Wover-aligned
	-Woverflow
	-Woverlength-strings
	-Woverloaded-shift-op-parentheses -Woverloaded-virtual
	# -Woverride-init (not in clang 7)
	-Woverride-module -Woverriding-method-mismatch
	# -Woverriding-option (not in clang 14)
	# -Woverriding-t-option (not in clang 18)
	-Wpacked
	# -Wpadded (ask for aligment)
	-Wparentheses
	-Wparentheses-equality
	-Wpartial-availability
	# -Rpass -Rpass-analysis
	-Wpass-failed
	# -Rpass-missed
	-Wpch-date-time
	-Wpedantic -Wpedantic-core-features
	# -Wpedantic-macros (clang 15, not in clang 13)
	-Wpessimizing-move
	-Wpointer-arith -Wpointer-bool-conversion
	# -Wpointer-compare (clang 11, not in clang 9)
	# -Wpointer-integer-compare (not in clang 7)
	-Wpointer-sign
	# -Wpointer-to-enum-cast (clang 11, not in clang 9)
	-Wpointer-to-int-cast -Wpointer-type-mismatch
	# -Wpoison-system-directories (clang 11, not in clang 9)
	# -Wpotentially-direct-selector (clang 11, not in clang 9)
	-Wpotentially-evaluated-expression
	# -Wpragma-clang-attribute -Wpragma-once-outside-header -Wpragma-pack -Wpragma-pack-suspicious-include -Wpragma-system-header-outside-header -Wpragmas
	# -Wpre-c++14-compat -Wpre-c++14-compat-pedantic
	# -Wpre-c++17-compat -Wpre-c++17-compat-pedantic
	# -Wpre-c++20-compat -Wpre-c++20-compat-pedantic -Wpre-c++2b-compat -Wpre-c++2b-compat-pedantic (clang 13, not in clang 11)
	# -Wpre-c2x-compat -Wpre-c2x-compat-pedantic (clang 13, not in clang 11)
	# -Wpre-openmp-51-compat (clang 13, not in clang 11)
	-Wpredefined-identifier-outside-function
	-Wprivate-extern -Wprivate-header -Wprivate-module
	-Wprofile-instr-missing -Wprofile-instr-out-of-date -Wprofile-instr-unprofiled
	-Wproperty-access-dot-syntax -Wproperty-attribute-mismatch
	-Wprotocol -Wprotocol-property-synthesis-ambiguity
	# -Wpsabi (clang 11, not in clang 9)
	-Wqualified-void-return-type
	-Wquoted-include-in-framework-header
	# -Wrange-loop-analysis (this check is overzealous in clang 9)
	-Wrange-loop-analysis
	# -Wrange-loop-bind-reference -Wrange-loop-construct (clang 11, not in clang 9)
	-Wreadonly-iboutlet-property
	-Wreceiver-expr -Wreceiver-forward-class
	-Wredeclared-class-member
	# -Wredundant-consteval-if (clang 15, not in clang 13)
	-Wredundant-decls -Wredundant-move -Wredundant-parens
	-Wregister
	-Wreinterpret-base-class
	-Rremark-backend-plugin
	-Wreorder
	# -Wreorder-ctor -Wreorder-init-list (clang 11, not in clang 9)
	-Wrequires-super-attribute
	-Wreserved-id-macro
	# -Wreserved-identifier -Wreserved-macro-identifier (clang 13, not in clang 11)
	-Wreserved-user-defined-literal
	# -Wrestrict-expansion (clang 15, not in clang 13)
	-Wretained-language-linkage
	-Wreturn-stack-address -Wreturn-std-move -Wreturn-type -Wreturn-type-c-linkage
	# -Wrewrite-not-bool (clang 11, not in clang 9)
	# -Rround-trip-cc1-args (clang 13, not in clang 11)
	# -Wrtti (clang 13, not in clang 11)
	-Wshadow
	-Wold-style-cast
	-Woverloaded-virtual
	-Rsanitize-address
	# -Rsearch-path-usage
	-Wsection
	-Wselector -Wselector-type-mismatch
	-Wself-assign -Wself-assign-field -Wself-assign-overloaded
	-Wself-move
	-Wsemicolon-before-method-body
	-Wsentinel
	-Wsequence-point
	-Wserialized-diagnostics
	-Wshadow -Wshadow-all -Wshadow-field -Wshadow-field-in-constructor -Wshadow-field-in-constructor-modified -Wshadow-ivar -Wshadow-uncaptured-local
	-Wshift-count-negative -Wshift-count-overflow -Wshift-negative-value -Wshift-op-parentheses -Wshift-overflow -Wshift-sign-overflow
	-Wshorten-64-to-32
	-Wsign-compare -Wsign-conversion -Wsign-promo
	-Wsigned-enum-bitfield
	# -Wsigned-unsigned-wchar (clang 11, not in clang 9)
	-Wsizeof-array-argument -Wsizeof-array-decay
	# -Wsizeof-array-div (clang 11, not in clang 9)
	# -Wsizeof-pointer-div (not in clang 7)
	-Wsizeof-pointer-memaccess
	-Wslash-u-filename
	# -Wslh-asm-goto (clang 11, not in clang 9)
	-Wsometimes-uninitialized
	# -Wsource-mgr (clang 13, not in clang 11)
	-Wsource-uses-openmp
	-Wspir-compat
	# -Wspirv-compat (clang 15, not in clang 13)
	# -Wstack-exhausted (clang 11, not in clang 9)
	-Wstack-protector
	-Wstatic-float-init -Wstatic-in-inline -Wstatic-inline-explicit-instantiation -Wstatic-local-in-inline -Wstatic-self-init
	-Wstdlibcxx-not-found
	-Wstrict-aliasing -Wstrict-aliasing=0 -Wstrict-aliasing=1 -Wstrict-aliasing=2
	# -Wstrict-overflow 
	# -Wstrict-overflow=0
	# -Wstrict-overflow=1
	# -Wstrict-overflow=2
	# -Wstrict-overflow=3
	# -Wstrict-overflow=4
	-Wstrict-overflow=5
	# -Wstrict-potentially-direct-selector (clang 11, not in clang 9)
	-Wstrict-prototypes -Wstrict-selector-match
	-Wstring-compare
	# -Wstring-concatenation (clang 13, not in clang 11)
	-Wstring-conversion -Wstring-plus-char -Wstring-plus-int
	-Wstrlcpy-strlcat-size
	-Wstrncat-size
	# -Wsuggest-destructor-override (clang 11, not in clang 9)
	# -Wsuggest-override (clang 11, not in clang 9)
	-Wsuper-class-method-mismatch
	-Wsuspicious-bzero -Wsuspicious-memaccess
	# -Wswift-name-attribute (clang 13, not in clang 11)
	-Wswitch -Wswitch-bool
	-Wno-switch-default
	-Wswitch-enum
	-Wsync-fetch-and-nand-semantics-changed
	-Wsynth
	# -Wtarget-clones-mixed-specifiers (clang 15, not in clang 13)
	# -Wtautological-bitwise-compare (clang 11, not in clang 9)
	-Wtautological-compare -Wtautological-constant-compare -Wtautological-constant-in-range-compare -Wtautological-constant-out-of-range-compare
	# -Wtautological-objc-bool-compare (not in clang 7)
	-Wtautological-overlap-compare -Wtautological-pointer-compare -Wtautological-type-limit-compare -Wtautological-undefined-compare
	# -Wtautological-unsigned-char-zero-compare -Wtautological-unsigned-enum-zero-compare
	# -Wtautological-unsigned-zero-compare -Wtautological-value-range-compare (clang 13, not in clang 11)
	# -Wtcb-enforcement (clang 13, not in clang 11)
	-Wtentative-definition-incomplete-type
	-Wthread-safety -Wthread-safety-analysis -Wthread-safety-attributes -Wthread-safety-beta -Wthread-safety-negative -Wthread-safety-precise -Wthread-safety-reference -Wthread-safety-verbose
	-Wtrigraphs
	-Wtype-limits -Wtype-safety
	-Wtypedef-redefinition
	-Wtypename-missing
	-Wunable-to-open-stats-file
	# -Wunaligned-access (clang 15, not in clang 13)
	# -Wunaligned-qualifier-implicit-cast (clang 15, not in clang 13)
	-Wunavailable-declarations
	-Wundeclared-selector
	-Wundef
	# -Wundef-prefix (clang 11, not in clang 9)
	-Wundefined-bool-conversion -Wundefined-func-template -Wundefined-inline -Wundefined-internal -Wundefined-internal-type
	-Wno-undefined-reinterpret-cast  # needed to cast to c-array
	-Wundefined-var-template
	# -Wunderaligned-exception-object (not in clang 7)
	-Wunevaluated-expression
	-Wunguarded-availability -Wunguarded-availability-new
	-Wunicode -Wunicode-homoglyph -Wunicode-whitespace
	# -Wunicode-zero-width (not in clang 7)
	-Wuninitialized
	# -Wuninitialized-const-reference (clang 11, not in clang 9)
	-Wunknown-argument
	# -Wunknown-assumption (clang 13, not in clang 11)
	-Wunknown-attributes
	# -Wunknown-cuda-version (clang 11, not in clang 9)
	-Wunknown-escape-sequence
	# -Wunknown-pragmas # (other compilers need their own pragmas)
	-Wno-unknown-pragmas # (other compilers need their own pragmas)
	# -Wunknown-sanitizers
	# -Wunknown-warning-option
	-Wunnamed-type-template-args
	-Wunneeded-internal-declaration -Wunneeded-member-function
	# -Wunqualified-std-cast-call (clang 15, not in clang 13)
	-Wunreachable-code -Wunreachable-code-aggressive -Wunreachable-code-break
	# -Wunreachable-code-fallthrough (clang 15, not in clang 13)
	-Wunreachable-code-loop-increment -Wunreachable-code-return
	$<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,15>:-Wunsafe-buffer-usage>
	-Wunsequenced
	# -Wunsupported-abi (clang 15, not in clang 13)
	-Wunsupported-abs -Wunsupported-availability-guard -Wunsupported-cb -Wunsupported-dll-base-class-template
	# -Wunsupported-floating-point-opt (clang 11, not in clang 9)
	-Wunsupported-friend -Wunsupported-gpopt -Wunsupported-nan -Wunsupported-target-opt -Wunsupported-visibility
	-Wunusable-partial-specialization
	-Wunused -Wunused-argument
	# -Wunused-but-set-parameter -Wunused-but-set-variable (clang 13, not in clang 11)
	-Wunused-command-line-argument -Wunused-comparison -Wunused-const-variable -Wunused-exception-parameter -Wunused-function -Wunused-getter-return-value -Wunused-label -Wunused-lambda-capture -Wunused-local-typedef -Wunused-local-typedefs -Wunused-macros -Wunused-member-function -Wunused-parameter -Wunused-private-field -Wunused-property-ivar -Wunused-result -Wunused-template -Wunused-value -Wunused-variable -Wunused-volatile-lvalue
	-Wused-but-marked-unused
	-Wuser-defined-literals -Wuser-defined-warnings
	-Wvarargs
	-Wvariadic-macros
	-Wvec-elem-size
	-Wvector-conversion -Wvector-conversions
	-Wvexing-parse
	-Wvisibility
	-Wvla -Wvla-extension
	# -Wvoid-pointer-to-enum-cast -Wvoid-pointer-to-int-cast (clang 11, not in clang 9)
	-Wvoid-ptr-dereference
	-Wvolatile-register-var
	# -Wwasm-exception-spec (clang 11, not in clang 9)
	-Wweak-template-vtables -Wweak-vtables
	-Wwritable-strings
	-Wwrite-strings
	# -Wxor-used-as-pow (clang 11, not in clang 9)
	-Wzero-as-null-pointer-constant -Wzero-length-array
	>
	$<$<CXX_COMPILER_ID:Intel>: -Werror -Wall -Wextra # also IntelLLVM, XL (ibm), XLClang (ibm)
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
	# $<$<OR:$<CXX_COMPILER_ID:PGI>,$<CXX_COMPILER_ID:NVHPC>>:
	$<$<CXX_COMPILER_ID:NVHPC>:
	-Werror -Wall
	-Wcast-qual
	-Wformat=2
	-Wshadow
	-Wuninitialized
	--diag_suppress unrecognized_gcc_pragma
	>
	$<$<CXX_COMPILER_ID:MSVC>:
	/GS  # Enables buffer security check (can help detect stack-based buffer overflows)
	/Wall /EHsc  # /volatile:iso /WX /Zc:wchar_t /Zc:forScope /Zc:inline /permissive- /W4
	/wd4371  # layout of class may have changed from a previous version of the compiler due to better packing of member
	/wd4514  # unreferenced inline function has been removed
	/wd4623  # default constructor was implicitly defined as deleted #TODO(correaa) investigate
	/wd4625  # copy constructor was implicitly defined as deleted
	/wd4626  # assignment operator was implicitly defined as deleted #TODO(correaa) investigate
	/wd4710  # function not inlined
	/wd4711  # function selected for automatic inline expansion
	/wd4820  # bytes padding added after base class
	/wd4866  # compiler may not enforce left-to-right evaluation order for call
	/wd4848  # support for attribute [[msvc::no_unique_address]] in C++17 and earlier is a vendor extension
	/wd4868  # compiler may not enforce left-to-right evaluation order in braced initializer list
	/wd5026  # move constructor was implicitly defined as deleted
	/wd5027  # move assignment operator was implicitly defined as deleted
	>
	)

	file(
		GLOB TEST_SRCS
		RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
		*.cpp
	)

	if(NOT Boost_FOUND)
		message(WARNING "Cannot find Boost, Multi library will have a very minimal test. If you want to test the library install Boost.Test, for example please run:\n    sudo apt install libboost-test-dev\n    sudo dnf install boost-devel")

		add_executable(main main.cpp)

		# target_link_libraries(main PRIVATE multi)
		target_compile_features   (main PUBLIC cxx_std_17)
		target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/include)

		# target_compile_definitions(main PRIVATE "$<$<CUDA_COMPILER_ID:NVIDIA>:BOOST_PP_VARIADICS>")

		if(NOT DEFINED ENABLE_CIRCLE)
			target_compile_options(main PRIVATE ${WARNS})
		else()
			target_compile_options    (main PRIVATE -Werror -Wall)
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
			# target_link_libraries     (${TEST_EXE} PRIVATE multi Boost::boost   )
			if(TBB_FOUND)
				target_link_libraries     (${TEST_EXE} PRIVATE TBB::tbb)
				target_compile_definitions(${TEST_EXE} PUBLIC TBB_FOUND=1)
			endif()

			if(NOT DEFINED ENABLE_CIRCLE)
				target_compile_options(${TEST_EXE} PRIVATE ${WARNS})
			else()
				target_compile_options    (${TEST_EXE} PRIVATE -Werror -Wall)
			endif()
			add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
		endforeach()
        
		# find_package(BLAS REQUIRED)
		# target_link_libraries(redux.cpp.x PRIVATE BLAS::BLAS)

	endif()

endif()
