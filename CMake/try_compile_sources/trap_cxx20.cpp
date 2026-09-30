// Trap C++20 or newer

int main(int argc, char **argv)
{
#if __cplusplus >= 202002L
#error Compiling in C++20 or newer.
#endif
    return 0;
}
