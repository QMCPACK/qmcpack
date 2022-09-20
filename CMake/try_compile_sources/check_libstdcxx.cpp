// Test the compatibility between the compiler and the libstdc++ from GNU
#include <cstdio>

int main(int argc, char **argv)
{
#if ( ( __INTEL_COMPILER == 1910 ) && ( _GLIBCXX_RELEASE < 9 || _GLIBCXX_RELEASE > 9 ) )
#error You are using the Intel compiler v19.1 ("20") which obtains libstdc++ from a GCC installation. Due to incompatibilities, you must use a GCC version 9 with this Intel compiler version. Check with 'icpc -v'.
#endif
    return 0;
}
