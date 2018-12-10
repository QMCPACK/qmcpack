

// Code to create coverage example
/* Steps to compare a single file
   1. g++ -O0 --coverage sample.cpp
      This should create sample.cpp.gcno (and a.out)
   2. Run a.out
      This should create sample.cpp.gcda
   3. Run gcov -b sample.cpp.gcda
      This should create sample.cpp.gcov
   4. Run python read-gcov.py -a test sample.cpp.gcov to ensure that read/write
      of gcov files works correctly.

   Steps to create a unit-test-like file and compare
   1. Make directories base/, unit/, and diff/
   2. Move previous sample.cpp.gcov file to base/ directory
   3. Remove sample.cpp.gcda
   4. g++ -O0 -DUNIT --coverage sample.cpp -o sample_unit
   5. Run sample_unit
   6. Run gcov -b sample.cpp.gcda
   7. Move sample.cpp to unit/ directory
   8. Use the 'diff' action to generate text-based comparisons and summaries:
        python compare_gcov.py --action diff --base-dir base --unit-dir unit

   9. Use the 'compare' action to create comparison gcov files:
        python compare_gcov.py --action compare --base-dir base --unit-dir unit --output-dir diff
       Now there should be sample.cpp.gcov in the diff/ directory.
       The code in func5 (and the call site) should be listed as uncovered.
*/


// Should be listed as unexecutable code
#if 0
int unused_func()
{
  int i = 1;
  return i;
}
#endif

int func1()
{
  int j = 0;
  int k = 0;
  for (int i = 0; i < 10; i++) {
    if (i < 5) { // branch coverage
      j += 4;
    }
    k += i;
  }
  return k;
}

// Should be uncovered
void func2()
{
  int i = 1;
  int j = i+2;
}

int func3(bool should_throw)
{
  if (should_throw) {
    throw  "An exception";
  }
  return 3;
}

void func4(int i)
{
  int j = 0;
  try {
    j = func3(false);
  } catch (...) {
    j = 3;
  }

}

// covered in base, uncovered in unit
void func5(int i)
{
  int j = i+1;
}


int main()
{
  int k = func1();
  // Should not be executed, so func2 is uncovered
  if (k < 0) {
    func2();
  }
  func4(0);
#ifndef UNIT
  func5(1);
#endif
  return 0;
}
