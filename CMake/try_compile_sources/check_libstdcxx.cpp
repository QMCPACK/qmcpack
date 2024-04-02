// Test the compatibility between the compiler and the libstdc++ from GNU
#include <cstdio>

int main(int argc, char **argv)
{
#if ( defined(_GLIBCXX_RELEASE) && ( _GLIBCXX_RELEASE < 9 ) )
#error Detected libstdc++. libstdc++ from a GCC version lower than 9 is not supported.
#endif
    return 0;
}
