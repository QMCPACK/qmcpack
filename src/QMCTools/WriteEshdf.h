/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////


#ifndef WRITE_ESHDF_H
#define WRITE_ESHDF_H
#include "hdf5.h"
#include <string>
#include <vector>
#include <cmath>

class XmlNode;
class FftContainer;

class EshdfFile 
{
private:
  hid_t file_;
  herr_t error_;

  int wrapped(int i, int size) const;
  void writeApplication(const std::string& appName, int major, int minor, int sub);
  void writeCreator();
  void writeFormat();
  void writeVersion();

  void readInEigFcn(const XmlNode& nd, FftContainer& cont);
  void handleSpinGroup(const XmlNode* nd, hid_t groupLoc, double& nocc, FftContainer& cont);
  double getOccupation(const XmlNode* nd) const;

  EshdfFile(const EshdfFile& f); // no copy constructor
  EshdfFile& operator=(const EshdfFile& f); // operator= not allowed
public:
  EshdfFile(const std::string& hdfFileName);
  ~EshdfFile();

  void writeBoilerPlate(const std::string& appName, const XmlNode& qboxSample);
  void writeSupercell(const XmlNode& qboxSample);
  void writeAtoms(const XmlNode& qboxSample);
  void writeElectrons(const XmlNode& qboxSample);
};

#endif
