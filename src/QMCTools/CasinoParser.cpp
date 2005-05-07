#include "QMCTools/CasinoParser.h"
#include "Utilities/OhmmsInfo.h"
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>

using namespace std;

CasinoParser::CasinoParser() {
  basisName = "casino-G2";
  Normalized = "yes";
}

CasinoParser::CasinoParser(int argc, char** argv): 
  QMCGaussianParserBase(argc,argv) {
  basisName = "casino-G2";
  Normalized = "yes";
}

void CasinoParser::parse(const std::string& fname) {

  std::ifstream fin(fname.c_str());

  //Grep the first word of the first line to assign a file name
  fin.getline(dbuffer,sizeof(dbuffer));
  std::istringstream a(dbuffer);
  a>>Title;

  std::string spin_unrestricted;
  LOGMSG("Looking for Spin ")
  search(fin,"Spin");
  getValue(fin,spin_unrestricted);
  if(spin_unrestricted.find("false")<spin_unrestricted.size()) {
    SpinRestricted=true;
  } else {
    SpinRestricted=false;
  }

  LOGMSG("Looking for electrons ")
  search(fin,"electrons");
  getValue(fin,NumberOfEls);

  LOGMSG("Looking for GEOMETRY ")
  search(fin, "GEOMETRY");
  getGeometry(fin);

  GroupName.resize(GroupID.size());
  for(int i=0; i<GroupID.size(); i++) GroupName[i]=IonName[GroupID[i]];

  LOGMSG("Looking for BASIS ")
  search(fin, "BASIS");
  getGaussianCenters(fin);

  //search(fin, "MULTIDETERMINANT");

  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
  vector<value_type> etemp;

  LOGMSG("Looking for EIGENVECTOR ")
  search(fin, "EIGENVECTOR");
  int nstates=SizeOfBasisSet;
  if(SpinRestricted) {
    EigVec.resize(SizeOfBasisSet*SizeOfBasisSet);
    etemp.resize(SizeOfBasisSet);
  } else {
    EigVec.resize(2*SizeOfBasisSet*SizeOfBasisSet);
    etemp.resize(2*SizeOfBasisSet);
  }
  getValues(fin,EigVec.begin(), EigVec.end());

  std::string aline;
  getline(fin,aline,'\n');
  getline(fin,aline,'\n');
  if(aline.find("EIGENVALUES")<aline.size()) {
    LOGMSG("Looking for EIGENVALUES ")
    skiplines(fin,2);
    getValues(fin,etemp.begin(), etemp.end());
    std::copy(etemp.begin(), etemp.begin()+SizeOfBasisSet, EigVal_alpha.begin());
    if(SpinRestricted) {
      EigVal_beta = EigVal_alpha;
    }else {
      std::copy(etemp.begin()+SizeOfBasisSet, etemp.end(), EigVal_beta.begin());
    }
  } else {
    WARNMSG("Missing EIGENVALUES. Modify determinant/@orbitals")
  }

  //Need to make correction to the eigen vectors
  int tot=0;
  for(int i=0; i<nstates; i++) {
    for(int j=0; j<SizeOfBasisSet; j++, tot++) EigVec[tot]*=BasisCorrection[j];
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
  search(is,"positions");
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
  LOGMSG("Checking shells")
  search(is, "shells");
  getValue(is,n);
  gShell.resize(n); gNumber.resize(n);

  LOGMSG("Checking SO")
  search(is, "AO");
  getValue(is,SizeOfBasisSet);

  LOGMSG("Checking Gaussian")
  search(is, "Gaussian");
  getValue(is,n);
  gExp.resize(n); gC0.resize(n); gC1.resize(n);

  LOGMSG("Checking shell types")
  search(is, "Code");
  getValues(is,gShell.begin(), gShell.end());

  LOGMSG("Checking the number of gaussians per shell")
  search(is, "Number");
  getValues(is,gNumber.begin(), gNumber.end());

  LOGMSG("Checking the bound of shells")
  search(is, "Sequence");
  getValues(is,gBound.begin(), gBound.end());
  for(int i=0; i<gBound.size(); i++) gBound[i]-= 1;

  LOGMSG("Checking the gaussian exponents")
  search(is, "Exponents");
  getValues(is,gExp.begin(), gExp.end());

  LOGMSG("Checking the gaussian contractions")
  search(is, "contraction");
  getValues(is,gC0.begin(), gC0.end());

  LOGMSG("Checking the gaussian contractions (sp)")
  search(is, "2nd");
  getValues(is,gC1.begin(), gC1.end());

  makeCorrections();
}

double CasinoParser::contractionCorrection(int shell_id, double alpha) {
  const double pi = 4.0*atan(1.0);
  const double sqrtpi=sqrt(pi);

  double fac;

  switch(shell_id) {
    case(1): //s, 1/(2 sqrt(pi))
      fac = 2.0*sqrtpi; break; 
    case(2): // sp
      fac = 2.0*sqrtpi; break;
    case(3): // p
      fac = sqrt(4.0/3.0)*sqrtpi; break;
    case(4): // d
      fac = sqrt(16.0/15.0)*sqrtpi; break;
    case(5): // f
      fac = 1.0e0;
      //fac *= pow(2.0e0*alpha,2)/sqrt(pi); break;
    default: // others, return 1 for now
      fac = 1.0e0; break;
  }

  return fac;
}

/** make corrections to gC0 and gC1 and tabulate the corrections to the eigen vectors
 */
void CasinoParser::makeCorrections() {

  int n=gC0.size();
  for(int i=0; i<n; i++) {
    gC0[i] *= contractionCorrection(gShell[i],gExp[i]);
    if(gShell[i] == 2) gC1[i] *= contractionCorrection(3,gExp[i]);
  }

  //s, sp and p do not need any correction
  BasisCorrection.resize(SizeOfBasisSet,1.0);
  std::vector<int> offset(10,0);
  offset[1]=1; //s
  offset[2]=4; //sp
  offset[3]=3; //p
  offset[4]=5; //d
  offset[5]=7; //f

  const double m40=sqrt(3.0e0);
  const double m50=sqrt(15.0e0/8.0e0);
  const double m51=1.5e0*sqrt(5.0e0);
  const double m52=15.0e0/sqrt(2.0e0);
  const double m53=15.0e0*sqrt(3.0e0);

  int basisCount=0;
  for(int i=0; i<n; i++) {
    if(gShell[i] == 4) {//d-orbital corrections
      BasisCorrection[basisCount++]=m40;//m=0
      BasisCorrection[basisCount++]=0.5;//m=1
      BasisCorrection[basisCount++]=0.5;//m=-1
      BasisCorrection[basisCount++]=1.0;//m=2
      BasisCorrection[basisCount++]=0.5;//m=-2
    } else if(gShell[i] == 5) {//f-orbital corrections
      BasisCorrection[basisCount++]=m50;//m=0
      BasisCorrection[basisCount++]=m51;//m=1
      BasisCorrection[basisCount++]=m51;//m=-1
      BasisCorrection[basisCount++]=m52;//m=2
      BasisCorrection[basisCount++]=m52;//m=-2
      BasisCorrection[basisCount++]=m53;//m=3
      BasisCorrection[basisCount++]=m53;//m=-3
    } else {//increase the count
      basisCount += offset[gShell[i]];
    }
  }
}
