// Test the OpenMP taskloop construct
int main(int argc, char **argv) {
  const int size = 10;
  int a[size];
#pragma omp taskloop default(shared)
  for (int i = 0; i < 10; i++) {
    a[i] = i * 2;
  }

  for (int i = 0; i < 10; i++)
    a[i]++;

  return 0;
}
