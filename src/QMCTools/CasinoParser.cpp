#include "QMCTools/CasinoParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

void CasinoParser::parse(const std::string& fname) {

  std::ifstream fin(fname.c_str());

  std::string spin_unrestricted;
  search(fin,"Spin");
  getValue(fin,spin_unrestricted);
  if(spin_unrestricted.find("false")<spin_unrestricted.size()) {
    SpinRestricted=true;
  } else {
    SpinRestricted=false;
  }

  search(fin,"electrons");
  getValue(fin,NumberOfEls);

  search(fin, "GEOMETRY");
  getGeometry(fin);

  search(fin, "BASIS");
  getGaussianCenters(fin);

  //search(fin, "MULTIDETERMINANT");
  search(fin, "EIGENVECTOR");
  int nline = (SizeOfBasisSet*SizeOfBasisSet)/4;
  //skip the line
  fin.getline( dbuffer, sizeof ( dbuffer ));
  EigVecU="\n";
  while(nline) {
    fin.getline( dbuffer, sizeof ( dbuffer ));
    EigVecU.append(dbuffer);
    EigVecU.append("\n");
    nline--;
  }
}

//std::copy(Qv.begin(), Qv.end(),ostream_iterator<double>(cout, ","));
void CasinoParser::getGeometry(std::istream& is) {
  //Number of atoms
  getNumberOfAtoms(is);
  //Atomic positions (au):
  getAtomicPositions(is);
  //Atomic numbers for each atom:
  getAtomicNumbers(is);
  //Valence charges for each atom:
  getValenceCharges(is);
}

void CasinoParser::getNumberOfAtoms(std::istream& is) {
  search(is,"Number");
  getValue(is,NumberOfAtoms);
  R.resize(NumberOfAtoms);
  GroupID.resize(NumberOfAtoms);
  Qv.resize(NumberOfAtoms);
  gBound.resize(NumberOfAtoms+1);
}
void CasinoParser::getAtomicPositions(std::istream& is) {
  search(is,"Atomic");
  getValues(is,R.begin(),R.end());
}
void CasinoParser::getAtomicNumbers(std::istream& is) {
  search(is,"Atomic");
  getValues(is,GroupID.begin(),GroupID.end());
}
void CasinoParser::getValenceCharges(std::istream& is) {
  search(is,"Valence");
  getValues(is,Qv.begin(),Qv.end());
}

void CasinoParser::getGaussianCenters(std::istream& is) {
  int n=0;
  search(is, "shells");
  getValue(is,n);
  gShell.resize(n); gNumber.resize(n);

  search(is, "AO");
  getValue(is,SizeOfBasisSet);

  search(is, "Gaussian");
  getValue(is,n);
  gExp.resize(n); gC0.resize(n); gC1.resize(n);

  search(is, "Code");
  getValues(is,gShell.begin(), gShell.end());

  search(is, "Number");
  getValues(is,gNumber.begin(), gNumber.end());

  search(is, "Sequence");
  getValues(is,gBound.begin(), gBound.end());
  for(int i=0; i<gBound.size(); i++) gBound[i]-= 1;

  search(is, "Exponents");
  getValues(is,gExp.begin(), gExp.end());

  search(is, "Correctly");
  getValues(is,gC0.begin(), gC0.end());

  search(is, "2nd");
  getValues(is,gC1.begin(), gC1.end());
}
