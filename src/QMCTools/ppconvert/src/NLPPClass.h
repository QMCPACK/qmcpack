//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
  State( std::string s);
};

class Config
{
  std::vector<State> States;
  std::map<char,int> ChannelRevMap;
public:
  inline const State& operator[](int i)
  { return States[i]; }

  size_t size() { return States.size(); }

  Config() {}
  Config ( std::string s);
};




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
  std::map<std::string,double> UnitToHartreeMap;
  std::map<std::string,double> UnitToBohrMap;
  std::map<XCType,std::string> XCMap;
  std::map<std::string,XCType> XCRevMap;
  std::map<int,std::string> ChannelMap;
  std::map<std::string,int> ChannelRevMap;
  std::map<int, string> ZToSymbolMap;
  std::map<int, double> ZToMassMap;
  std::map<std::string, int> SymbolToZMap;
  void SetupMaps();

  std::vector<ChannelPotentialClass> ChannelPotentials;
  CubSpline RhoAtom;
  int AtomicNumber;
  double PseudoCharge, TotalEnergy;
  std::string EnergyUnit, LengthUnit;
  int LocalChannel;
  // The grid is stored in bohr
  SimpleGrid PotentialGrid;
  bool Relativistic;
  XCType XC;
  void Write4Block(FILE *fout, std::vector<double> &data, int indent=2);
  bool GetNextState ( std::string &state, int &n, int &l, double &occ);
  //
  const double grid_delta;
public:
  bool WriteLogGrid;

  int GetNumChannels();
  bool HaveProjectors();
  bool ReadBFD_PP ( std::string fileName);
  bool ReadCASINO_PP ( std::string fileName);
  bool ReadCASINO_WF ( std::string fileName, int l);
  bool ReadFHI_PP ( std::string fileName);
  bool ReadGAMESS_PP ( std::string fileName);
  bool ReadUPF_PP ( std::string fileName);
  void WriteXML ( std::string fileName);
  void WriteABINIT ( std::string fileName="");
  void WriteUPF ( std::string fileName);
  void WriteASCII();
  void WriteFHI( std::string filename);
  void WriteFPMD( std::string filename);
  void WriteCASINO ( std::string filename);
  void WriteHDF5 ( std::string filename);
  void CalcProjector ( std::string refstate, int lchannel);
  double V     (double r);   double V     (int l, double r);
  double dVdr  (double r);   double dVdr  (int l, double r);
  double d2Vdr2(double r);   double d2Vdr2(int l, double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);

  void SetLocalChannel (int local)
  {
    LocalChannel = local;
  }

  PseudoClass() : XC(XC_NONE), Relativistic(false), LocalChannel(-1),
		  grid_delta(0.001), WriteLogGrid(false)
  {
    SetupMaps();
  }
};

#endif
