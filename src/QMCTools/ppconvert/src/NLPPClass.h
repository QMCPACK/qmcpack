#ifndef PSEUDO_CLASS_H
#define PSEUDO_CLASS_H

#include <map>
#include <string>
#include <vector>
#include "CubicSpline.h"
#include "XMLWriterClass2.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "common/PotentialBase.h"

struct State
{
  int n, l;
  double occ;
  State(int n_, int l_, double occ_) : n(n_), l(l_), occ(occ_)
  { }
  State(string s);
};

class Config
{
  vector<State> States;
  map<char,int> ChannelRevMap;
public:
  inline const State& operator[](int i)
  { return States[i]; }

  size_t size() { return States.size(); }

  Config() {}
  Config (string s);
};



using namespace std;

class ChannelPotentialClass
{
public:
  int l, n_principal;
  // Vl is stored in hartrees
  CubSpline Vl, Vlr;
  CubSpline ul;
  bool HasProjector;
  double Cutoff, Occupation, Eigenvalue;
  double FindCutoff();
  void WriteChannelLog (XMLWriterClass &writer, bool writeVl);
  void WriteChannelLinear (XMLWriterClass &writer, 
			   double dr, double rmax, bool writeVl);
  ChannelPotentialClass() :
    n_principal(0)
  {
    HasProjector= false;
    Occupation = 0.0;
    Eigenvalue = 0.0;
  }
};

typedef enum { XC_LDA, XC_GGA, XC_HF, XC_DF, XC_NONE} XCType;

class PseudoClass
  : Potential
{
private:
  map<string,double> UnitToHartreeMap;
  map<string,double> UnitToBohrMap;
  map<XCType,string> XCMap;
  map<string,XCType> XCRevMap;
  map<int,string> ChannelMap;
  map<string,int> ChannelRevMap;
  map<int, string> ZToSymbolMap;
  map<int, double> ZToMassMap;
  map<string, int> SymbolToZMap;
  void SetupMaps();

  vector<ChannelPotentialClass> ChannelPotentials;
  CubSpline RhoAtom;
  int AtomicNumber;
  double PseudoCharge, TotalEnergy;
  string EnergyUnit, LengthUnit;
  int LocalChannel;
  // The grid is stored in bohr
  SimpleGrid PotentialGrid;
  bool Relativistic;
  XCType XC;
  void Write4Block(FILE *fout, vector<double> &data, int indent=2);
  bool GetNextState (string &state, int &n, int &l, double &occ);
  //
  const double grid_delta;
public:
  bool WriteLogGrid;

  int GetNumChannels();
  bool HaveProjectors();
  bool ReadBFD_PP (string fileName);
  bool ReadCASINO_PP (string fileName);
  bool ReadCASINO_WF (string fileName, int l);
  bool ReadFHI_PP (string fileName);
  bool ReadGAMESS_PP (string fileName);
  void WriteXML (string fileName);
  void WriteABINIT (string fileName="");
  void WriteUPF (string fileName);
  void WriteASCII();
  void WriteFHI(string filename);
  void WriteFPMD(string filename);
  void WriteCASINO (string filename);
  void WriteHDF5 (string filename);
  void CalcProjector (string refstate, int lchannel);
  double V     (double r);   double V     (int l, double r);
  double dVdr  (double r);   double dVdr  (int l, double r);
  double d2Vdr2(double r);   double d2Vdr2(int l, double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);

  PseudoClass() : XC(XC_NONE), Relativistic(false),
		  grid_delta(0.001), WriteLogGrid(false)
  {
    SetupMaps();
  }
};

#endif
