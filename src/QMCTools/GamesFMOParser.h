//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



//#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
//#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"
class Monomer;
class IDmer;

class GamesFMOParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  GamesFMOParser();

  GamesFMOParser(int argc, char** argv);


  std::streampos pivot_begin;
  std::string psi_tag, ion_tag; 
  std::vector<std::string> tags;
  std::string *IDMonomer; 

  bool Mono; 
  std::string MOtype;


  //int nCartMO;
  int readtype;
  int TotNumAtom;
  int NumMonomer;
  int NumDimer;
  int NumTrimer;
  int FMOMethod;
  int MonomerID;

  Monomer * Tot_Monomer;
  
  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getAllMonomer(std::istream& is,int Index);

  void getMO(std::istream& is, std::string temp_tag);

  void getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal, std::string temp_tag);

  void getGaussianCenters(std::istream& is);

};



class Monomer: public QMCGaussianParserBase
{
private:
     int test;
public:
      Monomer();


  int q,AtomIndex;
  std::string tags;
  double X, Y, Z,ESPq;

  void print_Geometry();
  void parse(const std::string& fname){ };


};

class IDmer
{

public:
  IDmer();
  std::string MyId;
  int IndexI;
  int IndexJ;
  int IndexK;
  
};
