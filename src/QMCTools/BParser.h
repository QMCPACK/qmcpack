//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file BParser.h
 * @brief Declaration of BParser
 */
#ifndef QMCPLUSPLUS_TOOLS_BPARSER_H
#define QMCPLUSPLUS_TOOLS_BPARSER_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class BMakeFuncBase;

struct AGPLambda
{
  int I;
  int J;
  double X;
  AGPLambda(int i, int j, double x): I(i), J(j), X(x) {}
  AGPLambda(std::vector<std::string>& w)
  {
    I=atoi(w[0].c_str());
    J=atoi(w[1].c_str());
    X=atof(w[2].c_str());
  }
  xmlNodePtr createNode();
};

class BParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  int DetShells;
  int J3Shells;
  int J2Index;
  int DetSize;
  int J3Size;
  int DetNonZero;
  int J3NonZero;


  /** Size of the basis set per atom for the determinant */
  std::vector<int> detBasisPerAtom;
  /** Size of the basis set per atom for the three-body jastrow */
  std::vector<int> j3BasisPerAtom;
  /** Basis set per atom for the determinant*/
  std::map<int,std::vector<BMakeFuncBase*>*> detBasisSet;
  /** Basis set per atom for the three-body jastrow*/
  std::map<int,std::vector<BMakeFuncBase*>*> j3BasisSet;
  /** Occupation mask for the expanded basis set for the determinant */
  std::vector<int> detOcc;
  /** Occupation mask for the expanded basis set for the three-body jastrow */
  std::vector<int> j3Occ;

  /** non-zero paired Lambda elements for the determinant **/
  std::vector<AGPLambda> detPairedLambda;
  /** non-zero un-paired Lambda elements for the determinant **/
  std::vector<AGPLambda> detUnPairedLambda;
  /** non-zero Lambda elements for the three-body jastrow **/
  std::vector<AGPLambda> j3Lambda;

  ///default constructor
  BParser();

  ///another constructor
  BParser(int argc, char** argv);

  ///overwrite the virtual function
  void parse(const std::string& fname);
  void dump(const std::string& psi_tag, const std::string& ion_tag);

  void getGeometry(std::istream& is);
  void getBasisSetForDet(std::istream& is);
  void getBasisSetForJ3(std::istream& is);
  void getOccupationForDet(std::istream& is);
  void getOccupationForJ3(std::istream& is);
  void getLambdaForDet(std::istream& is);
  void getLambdaForJ3(std::istream& is);

  xmlNodePtr createDeterminantSet();
  xmlNodePtr createJ3();
  xmlNodePtr createBasisSet(std::map<int,std::vector<BMakeFuncBase*>*>& bset,
                            std::vector<int>& basisPerAtom, std::vector<int>& occ,
                            bool jastrow);
};
#endif
