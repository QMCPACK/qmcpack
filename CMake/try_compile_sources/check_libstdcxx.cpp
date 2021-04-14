// Test the compatibility between the compiler and the libstdc++ from GNU
#include <cstdio>

int main(int argc, char **argv)
{
// Unfortunately this check doesn't work for compilers <=v7.0 because _GLIBCXX_RELEASE appeared in the GCC 7.1 release.
// It is kept here as an example for the future.
#if ( defined(__INTEL_COMPILER) && (  _GLIBCXX_RELEASE < 7 ) )
#error You are using an Intel compiler. They obtain libstdc++ from a GNU compiler installation. For Intel compilers, you must use a GNU version >= 7. Found version <7.
#endif
// libstdc++ from GCC 8 is bad for Intel 19 in both C++14 and C++17
#if ( ( __INTEL_COMPILER == 1900 ) && (  _GLIBCXX_RELEASE > 7 ) )
#error You are using the Intel compiler v19 which obtains libstdc++ from a GNU compiler installation. You must use GNU version 7 with this Intel compiler. Found version >7. Alternatively (preferred route), use a more recent Intel compiler.
#endif
#if ( ( __INTEL_COMPILER == 1910 ) && (  _GLIBCXX_RELEASE > 9 ) )
#error You are using the Intel compiler v19.1 ("20") which obtains libstdc++ from a GNU compiler installation. Due to incompatibilities, you must use a GNU version <= 9 with this Intel compiler version. Found version >9.
#endif
    return 0;
}
