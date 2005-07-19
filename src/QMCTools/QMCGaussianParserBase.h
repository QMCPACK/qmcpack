#ifndef OHMMS_TOOLS_EXTERNAL_GAUSSIANPARSERBASE_H
#define OHMMS_TOOLS_EXTERNAL_GAUSSIANPARSERBASE_H
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include "OhmmsData/OhmmsElementBase.h"
#include "Utilities/SimpleParser.h"
#include "Particle/ParticleSet.h"

using namespace ohmmsqmc;

struct OhmmsAsciiParser {

  static const int bufferSize=200;
  char dbuffer[bufferSize];
  vector<string> currentWords;

  inline void skiplines(std::istream& is, int n) {
    while(n>0) {
      is.getline( dbuffer, bufferSize);--n;
    }
  }
  template<class T>
  inline void getValue(std::istream& is, T& aval) {
    is.getline( dbuffer,bufferSize);
    std::istringstream a(dbuffer); a>> aval;
  }

  template<class IT>
  inline void getValues(std::istream& is, IT first, IT last) {
 
    while(first != last) {
      is.getline( dbuffer,bufferSize);
      std::istringstream a(dbuffer);
      while(first != last && a >> *first){first++;}
    }

  }

  int search(std::istream& is, const std::string& keyword) {
    bool notfound = true;
    while(notfound) {
      std::string aline;
      getline(is,aline,'\n');
      if(! is){
	cout << "KEYWORD " << keyword << " : NOT FOUND. " << endl;
	abort();
      }
      if(aline.find(keyword) < aline.size()) {
        notfound = false;
      } 
    }
    return 1;
  }
};

struct QMCGaussianParserBase {

  typedef double value_type;
  typedef ParticleSet::SingleParticlePos_t SingleParticlePos_t;

  bool BohrUnit;
  bool SpinRestricted;
  int IonChargeIndex;
  int ValenceChargeIndex;
  int AtomicNumberIndex;
  int NumberOfAtoms;
  int NumberOfEls;
  int SpinMultiplicity;
  int NumberOfAlpha, NumberOfBeta;
  int SizeOfBasisSet;
  std::string Title;
  std::string basisType;
  std::string basisName;
  std::string Normalized;
  std::string CurrentCenter;

  ParticleSet IonSystem;

  std::vector<string> GroupName;
  std::vector<int> gShell, gNumber, gBound;
  std::vector<int> Occ_alpha, Occ_beta;
  std::vector<value_type> Qv;
  std::vector<value_type> gExp, gC0, gC1;
  std::vector<value_type> EigVal_alpha, EigVal_beta;
  std::vector<value_type> EigVec;
  //std::vector<GaussianCombo<value_type> > gExp, gC0, gC1;
  //std::string EigVecU, EigVecD;
  xmlNodePtr gridPtr;

  QMCGaussianParserBase();
  QMCGaussianParserBase(int argc, char** argv);

  void setOccupationNumbers();

  void createGridNode(int argc, char** argv);

  xmlNodePtr createElectronSet();
  xmlNodePtr createIonSet();
  xmlNodePtr createBasisSet();
  xmlNodePtr createCenter(int iat, int _off);
  void createShell(int n, int ig, int off_, xmlNodePtr abasis);
  xmlNodePtr createDeterminantSet();

  void map2GridFunctors(xmlNodePtr cur);
  virtual void parse(const std::string& fname) = 0;

  //static std::vector<std::string> IonName;
  static std::map<int,std::string> IonName;
  static std::vector<std::string> gShellType;
  static std::vector<int>         gShellID;

  static void init();
};
#endif
