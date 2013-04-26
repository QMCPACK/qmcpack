//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  AGPLambda(vector<string>& w)
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
  vector<int> detBasisPerAtom;
  /** Size of the basis set per atom for the three-body jastrow */
  vector<int> j3BasisPerAtom;
  /** Basis set per atom for the determinant*/
  map<int,vector<BMakeFuncBase*>*> detBasisSet;
  /** Basis set per atom for the three-body jastrow*/
  map<int,vector<BMakeFuncBase*>*> j3BasisSet;
  /** Occupation mask for the expanded basis set for the determinant */
  vector<int> detOcc;
  /** Occupation mask for the expanded basis set for the three-body jastrow */
  vector<int> j3Occ;

  /** non-zero paired Lambda elements for the determinant **/
  vector<AGPLambda> detPairedLambda;
  /** non-zero un-paired Lambda elements for the determinant **/
  vector<AGPLambda> detUnPairedLambda;
  /** non-zero Lambda elements for the three-body jastrow **/
  vector<AGPLambda> j3Lambda;

  ///default constructor
  BParser();

  ///another constructor
  BParser(int argc, char** argv);

  ///overwrite the virtual function
  void parse(const std::string& fname);
  void dump(const string& psi_tag, const string& ion_tag);

  void getGeometry(std::istream& is);
  void getBasisSetForDet(std::istream& is);
  void getBasisSetForJ3(std::istream& is);
  void getOccupationForDet(std::istream& is);
  void getOccupationForJ3(std::istream& is);
  void getLambdaForDet(std::istream& is);
  void getLambdaForJ3(std::istream& is);

  xmlNodePtr createDeterminantSet();
  xmlNodePtr createJ3();
  xmlNodePtr createBasisSet(map<int,vector<BMakeFuncBase*>*>& bset,
                            vector<int>& basisPerAtom, vector<int>& occ,
                            bool jastrow);
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
