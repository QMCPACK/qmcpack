set(CMAKE_CXX_COMPILER "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "2021.1.0.20201208")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "14")
set(CMAKE_CXX_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17;cxx_std_20")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "cxx_std_20")
set(CMAKE_CXX23_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "GNU")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_CXX_SIMULATE_VERSION "10.3.0")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_CXX_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;mpp;CPP;ixx;cppm)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/hdf5-1.10.6-6wt7jwobpoa3iwvoyabbtoilhpsguglf/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-vtune-2022.0.0-iggnhkd5gvt6es6xuj7kvvsjndwspbxy/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/ninja-1.10.2-4xa5uhs7kbc7egby2vhkts2yl3lvovfx/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-lfs-2.11.0-ey7djahsqy2b52fkw2murkwak4gag2jp/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-2.35.2-tf3qm22shnekjfqgkzk435xf7tjvuj3g/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/include/intel64;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/include/icc;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/include/c++/10.3.0;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/include/c++/10.3.0/x86_64-pc-linux-gnu;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/include/c++/10.3.0/backward;/usr/local/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/lib/gcc/x86_64-pc-linux-gnu/10.3.0/include;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/lib/gcc/x86_64-pc-linux-gnu/10.3.0/include-fixed;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/include;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpi;imf;svml;irng;stdc++;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/hwloc-2.7.1-ighhiqrlckk5ptrzojvkkyesnxkjkuyx/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/libevent-2.1.12-ooh2skj64cuq3awhabsjtz63gdpdxyhd/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/pmix-4.1.2-oxvur7ycvpgiyl73us2ywa3d6sw3cgib/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/hdf5-1.10.6-6wt7jwobpoa3iwvoyabbtoilhpsguglf/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/lib/intel64_lin;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-vtune-2022.0.0-iggnhkd5gvt6es6xuj7kvvsjndwspbxy/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/ninja-1.10.2-4xa5uhs7kbc7egby2vhkts2yl3lvovfx/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-lfs-2.11.0-ey7djahsqy2b52fkw2murkwak4gag2jp/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-2.35.2-tf3qm22shnekjfqgkzk435xf7tjvuj3g/lib;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/lib/gcc/x86_64-pc-linux-gnu/10.3.0;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/lib64;/lib64;/usr/lib64;/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/lib;/lib;/usr/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
