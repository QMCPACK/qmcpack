#include <iostream>
#include <complex>
#include "Blitz.h"

int main() {
  typedef Array<complex<double>,1> zVec;
  zVec a(4);
  a(0) = complex<double>(1,2);
  a(1) = complex<double>(3,4);
  a(2) = complex<double>(5,6);
  a(3) = complex<double>(7,8);

  double ** values= new double*[2];
	
  values[0] = reinterpret_cast<double*>(&a(0));
  values[1] = reinterpret_cast<double*>(&a(2));

  cout << "The first value is " << (values[0])[0] << endl;
  cout << "The second value is " << (values[1])[0] << endl;
}
