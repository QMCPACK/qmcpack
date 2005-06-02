#include "QMCTools/GaussianFCHKParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

GaussianFCHKParser::GaussianFCHKParser() {
  basisName = "Gaussian-G2";
  Normalized = "no";
}

GaussianFCHKParser::GaussianFCHKParser(int argc, char** argv): 
  QMCGaussianParserBase(argc,argv) {
  basisName = "Gaussian-G2";
  Normalized = "no";
}

void GaussianFCHKParser::parse(const std::string& fname) {

  std::ifstream fin(fname.c_str());

  getwords(currentWords,fin); //1
  Title = currentWords[0];

  getwords(currentWords,fin);//2  SP RHF Gen
  //if(currentWords[1]=="ROHF" || currentWords[1]=="UHF") {
  if(currentWords[1]=="UHF") {
    SpinRestricted=false;
  } else {
    SpinRestricted=true;
  }

  getwords(currentWords,fin);//3  Number of atoms
  NumberOfAtoms = atoi(currentWords.back().c_str());

  getwords(currentWords,fin); //4 Charge
  getwords(currentWords,fin); //5 Multiplicity

  getwords(currentWords,fin); //6 Number of electrons
  NumberOfEls=atoi(currentWords.back().c_str());

  getwords(currentWords,fin); //7 Number of alpha electrons
  int nup=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //8 Number of beta electrons 
  int ndown=atoi(currentWords.back().c_str());

  //NumberOfAlpha=nup-ndown;
  NumberOfAlpha=nup;
  NumberOfBeta=ndown;

  getwords(currentWords,fin); //9 Number of basis functions
  SizeOfBasisSet=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //10 Number of independant functions 
   ///
  int ng;
  bool notfound = true;
  while(notfound) {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++){
      if("contracted" == currentWords[i]){
	ng=atoi(currentWords.back().c_str());
	notfound = false;
      }
    }
  }
  //getwords(currentWords,fin); //11 Number of contracted shells
  //int ng=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //12 Number of contracted shells
  getwords(currentWords,fin); //13 Number of contracted shells
  getwords(currentWords,fin); //14 Number of contracted shells
  int nx=atoi(currentWords.back().c_str()); //the number of exponents

  //allocate everything here
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  Qv.resize(NumberOfAtoms);

  gBound.resize(NumberOfAtoms+1);
  gShell.resize(ng); 
  gNumber.resize(ng);
  gExp.resize(nx); 
  gC0.resize(nx); 
  gC1.resize(nx);


  search(fin, "Atomic numbers");//search for Atomic numbers
  getGeometry(fin);

  std::cout << "Number of gaussians " << ng << std::endl;
  std::cout << "Number of primitives " << nx << std::endl;
  std::cout << "Number of atoms " << NumberOfAtoms << std::endl;

  search(fin, "Shell types");
  getGaussianCenters(fin);

  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
  search(fin, "Alpha"); //search "Alpha Orbital Energies"
  getValues(fin,EigVal_alpha.begin(), EigVal_alpha.end());

  if(SpinRestricted) {
    EigVec.resize(SizeOfBasisSet*SizeOfBasisSet);
    EigVal_beta=EigVal_alpha;
  }
  else {
    EigVec.resize(2*SizeOfBasisSet*SizeOfBasisSet);
  }

  search(fin, "Alpha MO");
  getValues(fin,EigVec.begin(), EigVec.end());
}

void GaussianFCHKParser::getGeometry(std::istream& is) {
  //atomic numbers
  vector<int> q(NumberOfAtoms);
  getValues(is,q.begin(),q.end());
  search(is,"Nuclear");
  getValues(is,Qv.begin(),Qv.end());
  search(is,"coordinates");
  getValues(is,IonSystem.R.begin(),IonSystem.R.end());

  for(int i=0; i<NumberOfAtoms; i++) {
    GroupName[i]=IonName[q[i]];
    IonSystem.GroupID[i]=IonSystem.getSpeciesSet().addSpecies(GroupName[i]);
  }
  //for(int i=0; i<GroupID.size(); i++) GroupName[i]=IonName[GroupID[i]];
}

void GaussianFCHKParser::getGaussianCenters(std::istream& is) {

  //map between Gaussian to Casino Shell notation
  std::map<int,int> gsMap;
  gsMap[0] =1; //s
  gsMap[-1]=2; //sp
  gsMap[1] =3; //p
  gsMap[-2]=4; //d

  vector<int> n(gShell.size()), dn(NumberOfAtoms,0);

  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++) gShell[i]=gsMap[n[i]];

  search(is, "Number");
  getValues(is,gNumber.begin(), gNumber.end());

  search(is, "Shell");
  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++) dn[n[i]-1]+=1;

  gBound[0]=0;
  for(int i=0; i<NumberOfAtoms; i++) {
    gBound[i+1]=gBound[i]+dn[i];
  }

  search(is, "Primitive");
  getValues(is,gExp.begin(), gExp.end());

  search(is, "Contraction");
  getValues(is,gC0.begin(), gC0.end());

  search(is, "P(S=P)");
  getValues(is,gC1.begin(), gC1.end());
}

