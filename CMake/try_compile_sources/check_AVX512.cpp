// Check if AVX512 is activated by the compiler
int main(int argc, char **argv)
{
#if !defined(__AVX512F__)
#error "AVX512 not found"
#endif
  return 0;
}
